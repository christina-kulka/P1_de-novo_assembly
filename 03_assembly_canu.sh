#!/bin/bash

# Canu Assembly Script (Step 3c)
# Usage: ./03_assembly_canu.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Check if sample name is provided
SAMPLE=$1
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Create output directory for this sample
mkdir -p "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}"

# Check if processed reads file exists
CLEAN_READS="${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_super_clean_final.fastq"

if [ ! -f "$CLEAN_READS" ]; then
    echo "Error: Clean reads file not found: $CLEAN_READS"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input clean reads: $CLEAN_READS"
echo "Output directory: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out"

# Activate conda environment and run Canu
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_CANU

echo "Running Canu assembly..."
canu \
    -p "$SAMPLE" \
    -d "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out" \
    genomeSize=$GENOME_SIZE \
    -nanopore "$CLEAN_READS" \
    maxThreads=$THREADS \
    maxMemory=$CANU_MEMORY \
    correctedErrorRate=$CORRECTED_ERROR_RATE \
    minReadLength=$MIN_READ_LENGTH \
    minOverlapLength=$MIN_OVERLAP \
    corMaxEvidenceErate=0.20 \
    #utgRepeats=$UTG_REPEATS \
    "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"

# Check if assembly succeeded
if [ -f "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.contigs.fasta" ]; then
    echo "Canu assembly completed successfully!"
    grep "^>" "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.contigs.fasta" | wc -l | xargs echo "Number of contigs:"
    echo "Assembly file: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.contigs.fasta"
    
    # Also check for unitigs (uncorrected assembly)
    if [ -f "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.unitigs.fasta" ]; then
        echo "Unitigs file: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.unitigs.fasta"
    fi
else
    echo "Error: Canu assembly failed - no contigs.fasta found"
    exit 1
fi