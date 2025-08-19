#!/bin/bash

# Flye Assembly Script (Step 3a)
# Usage: ./03a_flye_assembly.sh SAMPLE_NAME

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
conda activate flye

echo "Running Flye assembly..."
flye --nano-raw "$CLEAN_READS" \
     --out-dir "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out" \
     --threads 20 \
     --genome-size 140k \
     --min-overlap 1000 \
     --iterations 1

# Check if assembly succeeded
if [ -f "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta" ]; then
    echo "Flye assembly completed successfully!"
    grep "^>" "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta" | wc -l | xargs echo "Number of contigs:"
    echo "Assembly file: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/flye_out/assembly.fasta"
else
    echo "Error: Flye assembly failed - no assembly.fasta found"
    exit 1
fi