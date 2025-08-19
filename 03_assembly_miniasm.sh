#!/bin/bash

# Miniasm Assembly Script
# Usage: ./03_assembly_miniasm.sh SAMPLE_NAME

# Set parameters
SAMPLE=$1
BASE_DIR="$HOME/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
RAW_DATA_DIR="${BASE_DIR}/00_raw_data_microsynth"
ASSEMBLY_OUTPUT_DIR="${BASE_DIR}/03_assembly"

# Check if sample name is provided
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Create output directory
mkdir -p "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"

# Find raw reads file (using the same approach as your working script)
SAMPLE_DIR="${RAW_DATA_DIR}/${SAMPLE}_results"
RAW_READS=$(find "$SAMPLE_DIR" -name "*raw*.fastq*" -type f | head -1)

if [ -z "$RAW_READS" ]; then
    echo "Error: No raw fastq file found for sample $SAMPLE"
    echo "Searched in: $SAMPLE_DIR"
    exit 1
fi

echo "Found raw reads file: $RAW_READS"

echo "Processing sample: $SAMPLE"
echo "Input raw reads: $RAW_READS"
echo "Output directory: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"

# Activate miniasm environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate miniasm

cd "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"

# Step 1: Generate overlaps with minimap2
echo "Step 1: Finding overlaps with minimap2..."
minimap2 -x ava-ont -t 20 -K 1G "$RAW_READS" "$RAW_READS" > overlaps.paf

# Step 2: Run miniasm assembly
echo "Step 2: Running miniasm assembly..."
miniasm -f "$RAW_READS" overlaps.paf > assembly.gfa

# Step 3: Convert GFA to FASTA
echo "Step 3: Converting to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' assembly.gfa > assembly.fasta

# Step 4: Polish with minipolish
echo "Step 4: Polishing with minipolish..."
minipolish --threads 20 "$RAW_READS" assembly.gfa > polished.gfa
awk '/^S/{print ">"$2"\n"$3}' polished.gfa > polished_assembly.fasta

echo "Miniasm assembly completed!"

# Show assembly statistics
if [ -f "polished_assembly.fasta" ]; then
    echo "=== ASSEMBLY STATISTICS ==="
    grep "^>" polished_assembly.fasta | wc -l | xargs echo "Number of contigs:"
    echo "Final assembly: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
fi