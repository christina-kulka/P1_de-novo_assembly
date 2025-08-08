#!/bin/bash

# Initial Assembly Script for Viral Genome (Step 3)
# Usage: ./03_assembly.sh SAMPLE_NAME

# Set parameters
SAMPLE=$1
BASE_DIR="data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
PREPROC_OUTPUT_DIR="${BASE_DIR}/02_preprocessing"
ASSEMBLY_OUTPUT_DIR="${BASE_DIR}/03_assembly"

# Check if sample name is provided
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Create output directory for this sample
mkdir -p "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate flye

# Check if clean reads file exists
CLEAN_READS="${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq"
if [ ! -f "$CLEAN_READS" ]; then
    echo "Error: Clean reads file not found: $CLEAN_READS"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input clean reads: $CLEAN_READS"
echo "Output directory: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}"

# Run Flye assembly (following paper's parameters)
echo "Running Flye assembly..."
flye --nano-raw "$CLEAN_READS" \
     --out-dir "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}" \
     --threads 12 \
     --meta \
     --keep-haplotypes \
     --min-overlap 5000 \
     --iterations 3

echo "Assembly completed for $SAMPLE"
echo "Results saved in: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}"

# Show assembly statistics
if [ -f "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/assembly.fasta" ]; then
    echo "Assembly statistics:"
    grep "^>" "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/assembly.fasta" | wc -l | xargs echo "Number of contigs:"
    echo "Assembly file: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/assembly.fasta"
else
    echo "Warning: assembly.fasta not found - check for errors"
fi