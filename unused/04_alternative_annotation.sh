#!/bin/bash

# Alternative Annotation Script - Creates DNA files with ORF annotations
# Usage: ./05_annotation_alternative.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SAMPLE=$1
ASSEMBLY_TYPE=${2:-"miniasm_viral"}

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    exit 1
fi

# Define assembly file paths
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
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_1000.fasta"
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

mkdir -p "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
cd "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE with assembly type: $ASSEMBLY_TYPE"

# Locate ORF files
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
ORF_PROTEIN_FILE="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}/${SAMPLE}_ORFs_proteins.faa"

# Ensure ORF protein file exists
if [ ! -f "$ORF_PROTEIN_FILE" ]; then
    if [ -f "$ORF_REFERENCE_FILE" ]; then
        echo "Translating ORF reference to protein file..."
        transeq -sequence "$ORF_REFERENCE_FILE" -outseq "$ORF_PROTEIN_FILE" -frame 1 -trim
    else
        echo "Error: ORF reference not found at $ORF_REFERENCE_FILE"
        exit 1
    fi
fi

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_PROKKA

# Run Prokka annotation with ORF protein file
echo "Running Prokka annotation with $ORF_PROTEIN_FILE..."
prokka \
    --kingdom Viruses \
    --genus Parapoxvirus \
    --strain "$SAMPLE" \
    --locustag "$SAMPLE" \
    --prefix "$SAMPLE" \
    --proteins "$ORF_PROTEIN_FILE" \
    --outdir . \
    --force \
    --cpus $THREADS \
    "$ASSEMBLY_FILE"

echo "Prokka annotation completed!"

# === Post-process GFF and GBK to ensure ORF IDs are used ===
if [ -f "${SAMPLE}.gff" ]; then
    echo "Replacing locus_tag/gene with ORF IDs..."

    awk -v OFS="\t" '
    $3=="CDS"{
      if($9 ~ /product=ORF_[0-9]+/){
        match($9,/product=(ORF_[0-9]+)/,m)
        id=m[1]
        gsub(/locus_tag=[^;]+/,"locus_tag="id,$9)
        if($9 ~ /gene=/){ gsub(/gene=[^;]+/,"gene="id,$9) } else { $9=$9";gene="id }
      }
    }1' ${SAMPLE}.gff > ${SAMPLE}.tmp.gff

    agat_sp_gff2gbk.pl --gff ${SAMPLE}.tmp.gff --fasta ${SAMPLE}.fna -o ${SAMPLE}.tmp.gbk

    mv ${SAMPLE}.tmp.gff ${SAMPLE}.gff
    mv ${SAMPLE}.tmp.gbk ${SAMPLE}.gbk
    echo "GBK fi
