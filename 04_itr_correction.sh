#!/bin/bash

# ITR Correction Script (Step 4 - Final)
# Usage: ./04_itr_correction.sh SAMPLE_NAME

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

mkdir -p "${ITR_OUTPUT_DIR}/${SAMPLE}"

# Find polished assembly and raw reads
MINIASM_DIR="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"
POLISHED_ASSEMBLY=$(find "$MINIASM_DIR" -name "*polished*.fasta" -type f | head -1)
RAW_READS=$(find "${RAW_DATA_DIR}/${SAMPLE}_results" -name "*raw*.fastq.gz" -type f | head -1)

if [ -z "$POLISHED_ASSEMBLY" ]; then
    echo "Error: No polished assembly found in $MINIASM_DIR"
    exit 1
fi

if [ -z "$RAW_READS" ]; then
    echo "Error: No raw reads found for sample $SAMPLE"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Input polished assembly: $POLISHED_ASSEMBLY"
echo "Input raw reads: $RAW_READS"
echo "Output directory: ${ITR_OUTPUT_DIR}/${SAMPLE}"

# Activate environment (TandemTools needs Python)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_MINIASM

cd "${ITR_OUTPUT_DIR}/${SAMPLE}"

echo "Running TandemTools ITR correction..."

# Run TandemTools with proper parameters
python "${TANDEMTOOLS}/tandemquast.py" \
    --only-polish \
    --nano "$RAW_READS" \
    -t $THREADS \
    -o tandem_output \
    "$POLISHED_ASSEMBLY"

echo "ITR correction completed!"

# Copy the final result to a clear filename
if [ -f "tandem_output/polished_sequences.fasta" ]; then
    cp tandem_output/polished_sequences.fasta final_orfv_genome.fasta
    echo "✅ Final ORFV genome with corrected ITRs: ${ITR_OUTPUT_DIR}/${SAMPLE}/final_orfv_genome.fasta"
else
    echo "⚠️  Check TandemTools output in: ${ITR_OUTPUT_DIR}/${SAMPLE}/tandem_output/"
fi