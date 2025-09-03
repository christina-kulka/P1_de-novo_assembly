#!/bin/bash

# Viral Verification Script (Step 3c)
# Usage: ./03b_viralverify.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# Check if sample name is provided
SAMPLE=$1
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Check if Miniasm assembly exists
MINIASM_DIR="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"
if [ ! -f "$MINIASM_DIR/polished_assembly.fasta" ]; then
    echo "Error: Miniasm assembly not found at $MINIASM_DIR/polished_assembly.fasta"
    echo "Run 03_assembly_miniasm.sh first!"
    exit 1
fi

# Check if clean reads file exists
CLEAN_READS="${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq"
if [ ! -f "$CLEAN_READS" ]; then
    echo "Error: Clean reads file not found: $CLEAN_READS"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input Miniasm assembly: $MINIASM_DIR"
echo "Input clean reads: $CLEAN_READS"
echo "Output directory: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification"

# Activate conda environment
source $CONDA_SETUP_PATH
conda activate $CONDA_ENV_QC

# Check if Pfam database exists
if [ ! -f "$PFAM_DB" ]; then
    echo "Error: Pfam database not found at $PFAM_DB"
    echo "Download with: wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz"
    exit 1
fi

echo "Running viral verification with viralVerify..."

# Run viralVerify to confirm viral identity
$VIRALVERIFY -f "$MINIASM_DIR/polished_assembly.fasta" \
           --hmm "$PFAM_DB" \
           -o "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification"

echo "Viral verification completed!"

# Parse and display results
VERIFICATION_DIR="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification"
VERIFICATION_FILE="${VERIFICATION_DIR}/result_table.csv"

if [ -f "$VERIFICATION_FILE" ]; then
    echo "=== VIRAL VERIFICATION RESULTS ==="
    # Display the results table
    cat "$VERIFICATION_FILE"
    
    # Check if assembly is classified as viral
    PREDICTION=$(tail -n 1 "$VERIFICATION_FILE" | cut -d',' -f2)
    if [ "$PREDICTION" = "Virus" ]; then
        echo ""
        echo "✅ CONFIRMED: Assembly is classified as VIRAL"
        echo "✅ Your ORFV genome assembly is verified!"
        
        # Show assembly statistics
        echo ""
        echo "=== FINAL ASSEMBLY STATISTICS ==="
        grep "^>" "$MINIASM_DIR/polished_assembly.fasta" | wc -l | xargs echo "Number of contigs:"
        GENOME_LENGTH=$(grep -v "^>" "$MINIASM_DIR/polished_assembly.fasta" | wc -c)
        echo "Genome length: $((GENOME_LENGTH-1)) bp"
        echo "Final verified ORFV genome: $MINIASM_DIR/polished_assembly.fasta"
    else
        echo ""
        echo "⚠️  WARNING: Assembly not classified as viral - check results"
        echo "Prediction: $PREDICTION"
    fi
else
    echo "Warning: Verification results file not found"
    echo "Check output in: $VERIFICATION_DIR"
fi

echo ""
echo "Verification results saved in: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification/"
echo "Original assembly file: $MINIASM_DIR/polished_assembly.fasta"