#!/bin/bash

# Miniasm Assembly Script (Step 3b)
# Usage: ./03_assembly_miniasm.sh SAMPLE_NAME

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
source $CONDA_SETUP_PATH
conda activate $CONDA_ENV_MINIASM

cd "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"

# Step 1: Generate overlaps with minimap2
echo "Step 1: Finding overlaps with minimap2..."
$MINIMAP2 -x ava-ont -t $THREADS -K 1G "$RAW_READS" "$RAW_READS" > overlaps.paf

# Step 2: Run miniasm assembly
echo "Step 2: Running miniasm assembly..."
miniasm -f "$RAW_READS" overlaps.paf > assembly.gfa

# Step 3: Convert GFA to FASTA
echo "Step 3: Converting to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' assembly.gfa > assembly.fasta

# Step 4: Polish with minipolish
echo "Step 4: Polishing with minipolish..."
$MINIPOLISH --threads $THREADS "$RAW_READS" assembly.gfa > polished.gfa
awk '/^S/{print ">"$2"\n"$3}' polished.gfa > polished_assembly.fasta

echo "Miniasm assembly completed!"

# Show assembly statistics
if [ -f "polished_assembly.fasta" ]; then
    echo "=== ASSEMBLY STATISTICS ==="
    grep "^>" polished_assembly.fasta | wc -l | xargs echo "Number of contigs:"
    echo "Final assembly: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
fi