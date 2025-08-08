#!/bin/bash

# Data Cleaning/Preprocessing Script for Nanopore Data
# Usage: ./02_preprocessing.sh SAMPLE_NAME

# Set parameters
SAMPLE=$1
BASE_DIR="data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"  
QC_OUTPUT_DIR="${BASE_DIR}/01_quality_control"
PREPROC_OUTPUT_DIR="${BASE_DIR}/02_preprocessing"

# Tool paths
PORECHOP="/home/ubuntu/miniconda3/bin/porechop"
NANOFILT="/home/ubuntu/.local/bin/NanoFilt"

# Check if sample name is provided
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Create output directory for this sample
mkdir -p "${PREPROC_OUTPUT_DIR}/${SAMPLE}"

# Find the raw fastq file from the original data
RAW_DATA_DIR="${BASE_DIR}/00_raw_data_microsynth"
SAMPLE_DIR="${RAW_DATA_DIR}/${SAMPLE}_results"
RAW_FASTQ=$(find "$SAMPLE_DIR" -name "*raw*.fastq*" -type f | head -1)

if [ -z "$RAW_FASTQ" ]; then
    echo "Error: No raw fastq file found for sample $SAMPLE"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input file: $RAW_FASTQ"
echo "Output directory: ${PREPROC_OUTPUT_DIR}/${SAMPLE}"

# Step 1: Adapter trimming with Porechop
echo "Step 1: Running Porechop for adapter trimming..."
$PORECHOP -i "$RAW_FASTQ" \
          -o "${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_trimmed.fastq" \
          --threads 4

# Step 2: Quality and length filtering with NanoFilt
echo "Step 2: Running NanoFilt for quality/length filtering..."
$NANOFILT --length 1000 \
          --quality 8 \
          --headcrop 10 \
          --tailcrop 10 \
          < "${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_trimmed.fastq" \
          > "${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq"

echo "Preprocessing completed for $SAMPLE"
echo "Final clean reads: ${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq"

# Optional: Get stats on clean reads
echo "Clean reads statistics:"
echo "Number of reads: $(grep -c '^@' ${PREPROC_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_clean.fastq)"