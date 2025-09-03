#!/bin/bash

# ORF Coordinate-Based Annotation Script
# Usage: ./orf_coordinate_annotate.sh INPUT_FILE SAMPLE_NAME OUTPUT_FILE

# Parse command line arguments
INPUT_FILE=$1
SAMPLE_NAME=$2
OUTPUT_FILE=$3

if [ $# -ne 3 ]; then
    echo "Usage: $0 INPUT_FILE SAMPLE_NAME OUTPUT_FILE"
    echo "Example: $0 /path/to/genome.dna B006 /path/to/output.gbk"
    echo ""
    echo "This script finds ORF coordinates using BLAST and adds them directly to annotations"
    exit 1
fi

# Get script directory and source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Define required files
ORF_DNA_FILE="${ORF_DATABASE_DIR}/${SAMPLE_NAME}_ORFs.fasta"
ORF_PROTEIN_FILE="${PROTEIN_DATABASE_DIR}/${SAMPLE_NAME}_protein.fasta"

# Check required files
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

if [ ! -f "$ORF_DNA_FILE" ]; then
    echo "Error: ORF DNA file not found: $ORF_DNA_FILE"
    echo "Expected: ${SAMPLE_NAME}_ORFs.fasta in $ORF_DATABASE_DIR"
    exit 1
fi

if [ ! -f "$ORF_PROTEIN_FILE" ]; then
    echo "Error: ORF protein file not found: $ORF_PROTEIN_FILE"
    echo "Expected: ${SAMPLE_NAME}_protein.fasta in $PROTEIN_DATABASE_DIR"
    exit 1
fi

echo "Processing: $SAMPLE_NAME"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "ORF DNA: $ORF_DNA_FILE"
echo "ORF Proteins: $ORF_PROTEIN_FILE"

# Create temporary working directory
WORK_DIR=$(mktemp -d)
echo "Working in: $WORK_DIR"

# Convert input file to FASTA format for BLAST
GENOME_FASTA="${WORK_DIR}/${SAMPLE_NAME}_genome.fasta"

echo "Converting input file to FASTA..."
$PYTHON3 -c "
from Bio import SeqIO
import sys

try:
    input_file = '$INPUT_FILE'
    output_file = '$GENOME_FASTA'
    
    # Read input file
    if input_file.lower().endswith('.dna'):
        try:
            record = SeqIO.read(input_file, 'snapgene')
        except:
            record = SeqIO.read(input_file, 'genbank')
    elif input_file.lower().endswith(('.gbk', '.gb', '.genbank')):
        record = SeqIO.read(input_file, 'genbank')
    else:
        record = SeqIO.read(input_file, 'fasta')
    
    # Write as FASTA
    SeqIO.write(record, output_file, 'fasta')
    print(f'Converted to FASTA: {len(record.seq)} bp')
    
except Exception as e:
    print(f'Error converting file: {e}')
    sys.exit(1)
"

if [ $? -ne 0 ]; then
    echo "Error: Failed to convert input file"
    rm -rf "$WORK_DIR"
    exit 1
fi

# Activate conda environment for BLAST
source $CONDA_SETUP_PATH
conda activate $CONDA_ENV_ITR  # Uses itr_analysis environment which has BLAST

# Create BLAST database from genome
echo "Creating BLAST database..."
DB_NAME="${WORK_DIR}/genome_db"
$MAKEBLASTDB -in "$GENOME_FASTA" -dbtype nucl -out "$DB_NAME" -title "${SAMPLE_NAME}_genome"

if [ $? -ne 0 ]; then
    echo "Error: Failed to create BLAST database"
    rm -rf "$WORK_DIR"
    exit 1
fi

# Run BLAST to find ORF coordinates
BLAST_RESULTS="${WORK_DIR}/orf_blast_results.txt"
echo "Running BLAST to find ORF coordinates..."

$BLASTN -query "$ORF_DNA_FILE" \
        -db "$DB_NAME" \
        -out "$BLAST_RESULTS" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -perc_identity 98 \
        -qcov_hsp_perc 90 \
        -word_size 11 \
        -max_target_seqs 10 \
        -max_hsps 10 \
        -dust no \
        -soft_masking false

if [ $? -ne 0 ]; then
    echo "Error: BLAST search failed"
    rm -rf "$WORK_DIR"
    exit 1
fi

# Check if we got results
if [ ! -s "$BLAST_RESULTS" ]; then
    echo "Warning: No BLAST matches found for ORFs"
    echo "Check that your ORF sequences match the genome sequence"
    rm -rf "$WORK_DIR"
    exit 1
fi

echo "BLAST search completed. Found matches:"
wc -l < "$BLAST_RESULTS"

# Create annotations using Python script
echo "Creating ORF annotations..."
$PYTHON3 "${SCRIPT_DIR}/create_orf_annotations.py" \
    "$INPUT_FILE" \
    "$BLAST_RESULTS" \
    "$ORF_PROTEIN_FILE" \
    "$OUTPUT_FILE" \
    --sample-name "$SAMPLE_NAME" \
    --min-identity 95

if [ $? -ne 0 ]; then
    echo "Error: Failed to create annotations"
    rm -rf "$WORK_DIR"
    exit 1
fi

# Show summary
echo ""
echo "=== ORF ANNOTATION SUMMARY ==="
if [ -f "$OUTPUT_FILE" ]; then
    echo "Successfully created: $OUTPUT_FILE"
    
    # Count features added
    FEATURE_COUNT=$(grep -c "/gene=" "$OUTPUT_FILE" 2>/dev/null || echo "0")
    echo "ORF features added: $FEATURE_COUNT"
    
    # Show some examples
    echo ""
    echo "Sample ORF annotations:"
    grep "/gene=" "$OUTPUT_FILE" | head -5 | sed 's/.*\/gene="/  - /' | sed 's/".*//'
else
    echo "Error: Output file was not created"
    rm -rf "$WORK_DIR"
    exit 1
fi

# Cleanup
rm -rf "$WORK_DIR"
echo ""
echo "ORF coordinate annotation completed successfully!"