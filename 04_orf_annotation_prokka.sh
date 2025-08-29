#!/bin/bash

# Prokka Annotation Script (Step 5)
# Usage: ./05_annotation_prokka.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Parse command line arguments
SAMPLE=$1
ASSEMBLY_TYPE=${2:-"miniasm_viral"}  # Default to miniasm_viral

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    echo "Example: $0 B006"
    echo "Example: $0 B006 canu"
    echo "Assembly types: miniasm_viral, miniasm_raw, canu, canu_ultra, canu_super, microsynth"
    exit 1
fi

# Define assembly file paths based on type
case $ASSEMBLY_TYPE in
    "miniasm_viral")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta"
        ;;
    "miniasm_raw")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
        ;;
    "canu")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_ultra_output/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra_trimmed")
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed.fasta"
        ;;
    "canu_super")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "microsynth")
        ASSEMBLY_FILE="${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"
        ;;
    *)
        echo "Error: Unknown assembly type '$ASSEMBLY_TYPE'"
        echo "Valid types: miniasm_viral, miniasm_raw, canu, microsynth"
        exit 1
        ;;
esac

# Create output directory with assembly type subfolder
mkdir -p "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found for sample $SAMPLE, type $ASSEMBLY_TYPE"
    echo "Expected: $ASSEMBLY_FILE"
    exit 1
fi

echo "Found assembly file: $ASSEMBLY_FILE"
echo "Processing sample: $SAMPLE"
echo "Assembly type: $ASSEMBLY_TYPE"
echo "Output directory: ${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

# Activate conda environment (assuming prokka is in your miniasm env or create new one)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_PROKKA

cd "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

# Find the corresponding ORF reference file for this sample
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
ORF_PROTEIN_FILE="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}/${SAMPLE}_ORFs_proteins.faa"

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