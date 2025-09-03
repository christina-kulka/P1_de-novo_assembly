#!/bin/bash

# Prokka ORF-Only Annotation Script
# Usage: ./prokka_orf_only.sh INPUT_FILE SAMPLE_NAME OUTPUT_FILE

# Parse command line arguments
INPUT_FILE=$1
SAMPLE_NAME=$2
OUTPUT_FILE=$3

if [ $# -ne 3 ]; then
    echo "Usage: $0 INPUT_FILE SAMPLE_NAME OUTPUT_FILE"
    echo "Example: $0 /path/to/genome.dna B006 /path/to/output.gbk"
    exit 1
fi

# Get script directory and source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
ORF_PROTEIN_FILE="${PROTEIN_DATABASE_DIR}/${SAMPLE_NAME}_protein.fasta"
PROKKA_PROTEIN_FILE="${PROTEIN_DATABASE_DIR}/${SAMPLE_NAME}_protein_prokka.faa"

# Check required files
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
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

# Create temporary working directory
TEMP_DIR=$(mktemp -d)
echo "Working in: $TEMP_DIR"

# Convert input file to FASTA format for Prokka
INPUT_FASTA="${TEMP_DIR}/${SAMPLE_NAME}_input.fasta"

# Extract sequence from input file (works for .dna, .gbk, .fasta)
$PYTHON3 -c "
from Bio import SeqIO
import sys

try:
    # Try different formats
    input_file = '$INPUT_FILE'
    output_file = '$INPUT_FASTA'
    
    # Determine format and read
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
    print(f'Converted {input_file} to FASTA format')
    print(f'Sequence length: {len(record.seq)} bp')
    
except Exception as e:
    print(f'Error converting file: {e}')
    sys.exit(1)
"

if [ $? -ne 0 ]; then
    echo "Error: Failed to convert input file to FASTA"
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Check if pre-formatted protein file exists, if not create it
if [ ! -f "$PROKKA_PROTEIN_FILE" ]; then
    echo "Creating pre-formatted protein file..."
    if [ -f "${SCRIPT_DIR}/4a_format_proteins.py" ]; then
        $PYTHON3 "${SCRIPT_DIR}/4a_format_proteins.py" "$ORF_PROTEIN_FILE" "$PROKKA_PROTEIN_FILE"
    else
        echo "Error: 4a_format_proteins.py not found in $SCRIPT_DIR"
        rm -rf "$TEMP_DIR"
        exit 1
    fi
fi

echo "Using protein file: $PROKKA_PROTEIN_FILE"

# Activate conda environment
source $CONDA_SETUP_PATH
conda activate prokka

# Change to temp directory
cd "$TEMP_DIR"

# Run Prokka with ORF proteins only
echo "Running Prokka with ORF database only..."
prokka \
    --kingdom Viruses \
    --genus Parapoxvirus \
    --strain "$SAMPLE_NAME" \
    --locustag "$SAMPLE_NAME" \
    --proteins "$PROKKA_PROTEIN_FILE" \
    --prefix "$SAMPLE_NAME" \
    --outdir . \
    --force \
    --cpus $THREADS \
    --evalue 1e-6 \
    --coverage 50 \
    --norrna \
    --notrna \
    --rfam \
    "$INPUT_FASTA"

if [ $? -ne 0 ]; then
    echo "Error: Prokka failed"
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "Prokka completed successfully!"

# Restore ORF names using existing script
if [ -f "${SCRIPT_DIR}/4b_restore_orf_names.py" ]; then
    echo "Restoring ORF names..."
    # $PYTHON3 "${SCRIPT_DIR}/4b_restore_orf_names.py" "$SAMPLE_NAME" "$ORF_PROTEIN_FILE"
    # TODO change this back when using the pipeline!
    $PYTHON3 "${SCRIPT_DIR}/4b_restore_orf_names.py" "$SAMPLE_NAME" "$ORF_PROTEIN_FILE" "$TEMP_DIR" # for proposal annotation

    # Use the ORF-named file if it was created
    if [ -f "${SAMPLE_NAME}.orf_named.gbk" ]; then
        FINAL_GBK="${SAMPLE_NAME}.orf_named.gbk"
    else
        FINAL_GBK="${SAMPLE_NAME}.gbk"
    fi
else
    echo "Warning: 4b_restore_orf_names.py not found, using Prokka output as-is"
    FINAL_GBK="${SAMPLE_NAME}.gbk"
fi

cp "${SAMPLE_NAME}".* /tmp/ 2>/dev/null
if [ -f "${SAMPLE_NAME}.orf_named.gbk" ]; then
    cp "${SAMPLE_NAME}.orf_named.gbk" /tmp/
fi

# Copy result to final location
if [ -f "$FINAL_GBK" ]; then
    cp "$FINAL_GBK" "$OUTPUT_FILE"
    echo "Output saved to: $OUTPUT_FILE"
    
    # Show summary
    echo ""
    echo "=== SUMMARY ==="
    if [ -f "${SAMPLE_NAME}.txt" ]; then
        cat "${SAMPLE_NAME}.txt"
    fi
    
    # Show some ORF annotations found
    echo ""
    echo "ORF annotations added:"
    grep "product=" "$FINAL_GBK" | head -10 | sed 's/.*product="/  - /' | sed 's/".*//'
    
else
    echo "Error: Expected output file not found"
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Cleanup
rm -rf "$TEMP_DIR"
echo "Annotation completed!"