#!/bin/bash

# Prokka Annotation Script (Step 5)
# Usage: ./05_annotation_prokka.sh SAMPLE_NAME

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

# Create output directory
ANNOTATION_OUTPUT_DIR="${BASE_DIR}/05_annotation"
mkdir -p "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}"

# Find assembly file (looking in miniasm output)
ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found for sample $SAMPLE"
    echo "Expected: $ASSEMBLY_FILE"
    exit 1
fi

echo "Found assembly file: $ASSEMBLY_FILE"
echo "Processing sample: $SAMPLE"
echo "Output directory: ${ANNOTATION_OUTPUT_DIR}/${SAMPLE}"

# Activate conda environment (assuming prokka is in your miniasm env or create new one)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_PROKKA

cd "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}"

# Find the corresponding ORF reference file for this sample
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
ORF_PROTEIN_FILE="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_ORFs_proteins.faa"

# Check if ORF reference exists and convert to proteins
if [ ! -f "$ORF_REFERENCE_FILE" ]; then
    echo "Warning: ORF reference not found at $ORF_REFERENCE_FILE"
    echo "Running Prokka with default viral database..."
    PROTEIN_FLAG=""
else
    echo "Using custom ORF reference: $ORF_REFERENCE_FILE"
    echo "Converting nucleotide ORFs to proteins..."
    
    # Convert nucleotide ORFs to protein sequences
    transeq -sequence "$ORF_REFERENCE_FILE" -outseq "$ORF_PROTEIN_FILE" -frame 1 -trim
    
    if [ -f "$ORF_PROTEIN_FILE" ]; then
        echo "Protein file created: $ORF_PROTEIN_FILE"
        PROTEIN_FLAG="--proteins $ORF_PROTEIN_FILE"
    else
        echo "Warning: Failed to create protein file. Using default database..."
        PROTEIN_FLAG=""
    fi
fi

# Run Prokka annotation
echo "Running Prokka annotation..."
prokka \
    --kingdom Viruses \
    --genus Parapoxvirus \
    --strain "$SAMPLE" \
    --locustag "$SAMPLE" \
    --prefix "$SAMPLE" \
    --outdir . \
    --force \
    --cpus $THREADS \
    $PROTEIN_FLAG \
    "$ASSEMBLY_FILE"

echo "Prokka annotation completed!"

# Show annotation summary
if [ -f "${SAMPLE}.txt" ]; then
    echo "=== ANNOTATION SUMMARY ==="
    cat "${SAMPLE}.txt"
    echo "Key output files:"
    echo "  - ${SAMPLE}.gff: Gene annotations"
    echo "  - ${SAMPLE}.faa: Protein sequences"
    echo "  - ${SAMPLE}.gbk: GenBank format"
fi