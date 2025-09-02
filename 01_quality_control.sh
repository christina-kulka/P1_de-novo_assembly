#!/bin/bash

# Quality Control Script for Nanopore Data (Step 1)
# Usage: ./01_quality_control.sh SAMPLE_NAME

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
mkdir -p "${QC_OUTPUT_DIR}/${SAMPLE}"

# Activate conda environment
source $CONDA_SETUP_PATH
conda activate $CONDA_ENV_QC

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "NanoPlot location: $(which NanoPlot)"
echo "Python location: $(which python)"
echo "PATH: $PATH"

# Find the raw fastq file for this sample (contains all reads)
SAMPLE_DIR="${RAW_DATA_DIR}/${SAMPLE}_results/${SAMPLE}_results/Sequencing/"
echo "Processing sample: $SAMPLE"
echo "Looking for raw fastq file in: $SAMPLE_DIR"

# Find the raw fastq file specifically (handles .gz files)
RAW_FASTQ=$(find "$SAMPLE_DIR" -name "*raw*.fastq*" -type f | head -1)

if [ -z "$RAW_FASTQ" ]; then
    echo "Error: No raw fastq file found for sample $SAMPLE"
    exit 1
fi

echo "Found raw fastq file: $RAW_FASTQ"

# Run NanoPlot on the raw fastq file
NanoPlot --fastq "$RAW_FASTQ" \
         --outdir "${QC_OUTPUT_DIR}/${SAMPLE}" \
         --prefix "${SAMPLE}_raw" \
         --plots dot --plots kde

echo "Quality control completed for $SAMPLE using raw reads"
echo "Results saved in: ${QC_OUTPUT_DIR}/${SAMPLE}"