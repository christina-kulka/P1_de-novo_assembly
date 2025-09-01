#!/bin/bash

# Prokka Annotation Script (Step 5) - Clean Version
# Usage: ./05_annotation_prokka.sh SAMPLE_NAME [ASSEMBLY_TYPE]

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
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_150.fasta" #TODO fix this
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

ORF_PROTEIN_FILE="/home/ubuntu/data-volume/001_Raw_Data/Databases/Proteins/${SAMPLE}_protein.fasta"

# Create output directory
mkdir -p "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

# Validate input files
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found for sample $SAMPLE, type $ASSEMBLY_TYPE"
    echo "Expected: $ASSEMBLY_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "Assembly type: $ASSEMBLY_TYPE"
echo "Assembly file: $ASSEMBLY_FILE"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_PROKKA

cd "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

# Run Prokka with or without custom proteins
if [ ! -f "$ORF_PROTEIN_FILE" ]; then
    echo "ORF reference not found. Running Prokka with default viral database..."
    
    prokka \
        --kingdom Viruses \
        --genus Parapoxvirus \
        --strain "$SAMPLE" \
        --locustag "$SAMPLE" \
        --prefix "$SAMPLE" \
        --outdir . \
        --force \
        --cpus $THREADS \
        "$ASSEMBLY_FILE"
else
    echo "Using custom ORF reference: $ORF_PROTEIN_FILE"
    
    # Format protein file using external script
    PROKKA_PROTEIN_FILE="${SAMPLE}_prokka_proteins.faa"

    if [ ! -f "${SCRIPT_DIR}/4a_format_proteins.py" ]; then
        echo "Error: 4a_format_proteins.py not found in ${SCRIPT_DIR}"
        exit 1
    fi
    
    python3 "${SCRIPT_DIR}/4a_format_proteins.py" "$ORF_PROTEIN_FILE" "$PROKKA_PROTEIN_FILE"
    
    # Run Prokka with custom proteins
    prokka \
        --kingdom Viruses \
        --genus Parapoxvirus \
        --strain "$SAMPLE" \
        --locustag "$SAMPLE" \
        --proteins "$PROKKA_PROTEIN_FILE" \
        --prefix "$SAMPLE" \
        --outdir . \
        --force \
        --cpus $THREADS \
        --evalue 1e-6 \
        --coverage 50 \
        --notrna \
        "$ASSEMBLY_FILE"
fi

echo "Prokka annotation completed!"

# Restore ORF names using external script
if [ -f "$ORF_PROTEIN_FILE" ] && [ -f "${SCRIPT_DIR}/4b_restore_orf_names.py" ]; then
    echo "Restoring ORF names..."
    python3 "${SCRIPT_DIR}/4b_restore_orf_names.py" "$SAMPLE" "$ORF_PROTEIN_FILE"
elif [ -f "$ORF_PROTEIN_FILE" ]; then
    echo "Warning: 4b_restore_orf_names.py not found in ${SCRIPT_DIR}"
fi

# Show results
echo ""
echo "=== ANNOTATION SUMMARY ==="
if [ -f "${SAMPLE}.txt" ]; then
    cat "${SAMPLE}.txt"
fi

echo ""
echo "=== OUTPUT FILES ==="
echo "Main files:"
echo "  - ${SAMPLE}.gff (Gene annotations)"
echo "  - ${SAMPLE}.gbk (GenBank format)"
echo "  - ${SAMPLE}.faa (Protein sequences)"
echo "  - ${SAMPLE}.fna (Nucleotide sequences)"

if [ -f "${SAMPLE}.orf_named.gff" ]; then
    echo ""
    echo "ORF-named files:"
    echo "  - ${SAMPLE}.orf_named.gff (GFF with ORF names)"
    echo "  - ${SAMPLE}.orf_named.gbk (GenBank with ORF names)"
    echo "  - ${SAMPLE}.orf_named.faa (Proteins with ORF names)"
    echo "  - ${SAMPLE}_orf_mapping.txt (Mapping report)"
fi

echo ""
echo "Annotation pipeline completed!"