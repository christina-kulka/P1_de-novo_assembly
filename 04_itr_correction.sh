#!/bin/bash

# ITR Correction Script (Step 4 - Final)
# Usage: ./04_itr_correction.sh SAMPLE_NAME

SAMPLE=$1
BASE_DIR="$HOME/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
ASSEMBLY_OUTPUT_DIR="${BASE_DIR}/03_assembly"
ITR_OUTPUT_DIR="${BASE_DIR}/04_final_results"

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

mkdir -p "${ITR_OUTPUT_DIR}/${SAMPLE}"

# Find polished assembly and raw reads
MINIASM_DIR="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"
POLISHED_ASSEMBLY=$(find "$MINIASM_DIR" -name "*polished*.fasta" -type f | head -1)
RAW_READS=$(find "${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results" -name "*raw*.fastq.gz" -type f | head -1)

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
conda activate miniasm

cd "${ITR_OUTPUT_DIR}/${SAMPLE}"

echo "Running TandemTools ITR correction..."

# Run TandemTools with proper parameters
python /home/ubuntu/TandemTools/tandemquast.py \
    --only-polish \
    --nano "$RAW_READS" \
    -t 20 \
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