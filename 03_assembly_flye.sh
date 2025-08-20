#!/bin/bash

# Flye Assembly Script (Step 3a)
# Usage: ./03_assembly_flye.sh SAMPLE_NAME

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

# Check if clean reads file exists
CLEAN_READS="${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq"

if [ ! -f "$CLEAN_READS" ]; then
    echo "Error: Clean reads file not found: $CLEAN_READS"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input clean reads: $CLEAN_READS"
echo "Output directory: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out"

# Activate conda environment and run Flye
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_FLYE

echo "Running Flye assembly..."
flye --nano-raw "$CLEAN_READS" \
     --out-dir "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out" \
     --threads $THREADS \
     --genome-size $GENOME_SIZE \
     --min-overlap $MIN_OVERLAP \
     --iterations 3

# Check if assembly succeeded
if [ -f "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta" ]; then
    echo "Flye assembly completed successfully!"
    grep "^>" "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta" | wc -l | xargs echo "Number of contigs:"
    echo "Assembly file: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta"
else
    echo "Error: Flye assembly failed - no assembly.fasta found"
    exit 1
fi