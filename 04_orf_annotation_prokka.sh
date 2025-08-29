#!/bin/bash

# Prokka Annotation Script (Step 5) - Enhanced ORF Naming
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
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_1000.fasta" #TODO fix this
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

# === Enhanced ORF ID mapping using BLAST + awk ===
if [ -f "${SAMPLE}.gff" ] && [ -f "$ORF_PROTEIN_FILE" ] && command -v blastp &> /dev/null; then
    echo "Mapping Prokka genes to reference ORFs..."
    
    # Create BLAST database
    makeblastdb -in "$ORF_PROTEIN_FILE" -dbtype prot -out orf_db -logfile /dev/null
    
    # BLAST Prokka proteins against ORF database
    blastp -query "${SAMPLE}.faa" -db orf_db \
        -outfmt "6 qseqid sseqid pident evalue" \
        -max_target_seqs 1 -evalue 1e-3 -num_threads $THREADS > blast_results.tmp
    
    # Create mapping file (only high-confidence matches: >70% identity)
    awk '$3 >= 70 && $4 <= 1e-3 {print $1"\t"$2"\t"$3}' blast_results.tmp > gene_to_orf_map.txt
    
    # Count mappings
    MAPPED_GENES=$(wc -l < gene_to_orf_map.txt)
    echo "Found $MAPPED_GENES high-confidence ORF mappings"
    
    if [ $MAPPED_GENES -gt 0 ]; then
        # Backup original GFF
        cp "${SAMPLE}.gff" "${SAMPLE}.gff.backup"
        
        # Update GFF with ORF names using awk
        awk -F'\t' '
        BEGIN {
            # Read the mapping file
            while ((getline line < "gene_to_orf_map.txt") > 0) {
                split(line, fields, "\t")
                map[fields[1]] = fields[2]
                identity[fields[1]] = fields[3]
            }
            close("gene_to_orf_map.txt")
        }
        /^#/ || NF != 9 { print; next }  # Print headers and non-feature lines as-is
        $3 == "CDS" {
            # Extract locus_tag from attributes
            if (match($9, /locus_tag=([^;]+)/, locus)) {
                locus_tag = locus[1]
                if (locus_tag in map) {
                    orf_name = map[locus_tag]
                    # Replace locus_tag and gene with ORF name
                    gsub(/locus_tag=[^;]+/, "locus_tag=" orf_name, $9)
                    gsub(/gene=[^;]+/, "gene=" orf_name, $9)
                    # Add gene if it doesnt exist
                    if ($9 !~ /gene=/) {
                        $9 = $9 ";gene=" orf_name
                    }
                    # Update product if its generic
                    if ($9 ~ /product=hypothetical protein/) {
                        gsub(/product=hypothetical protein/, "product=" orf_name " (" identity[locus_tag] "% id)", $9)
                    }
                }
            }
        }
        { print }
        ' "${SAMPLE}.gff.backup" > "${SAMPLE}.gff"
        
        # Update protein sequences (FAA file) with ORF names
        echo "Updating protein sequences with ORF names..."
        cp "${SAMPLE}.faa" "${SAMPLE}.faa.backup"
        
        # Update FAA headers using the mapping
        awk '
        BEGIN {
            # Read the mapping file
            while ((getline line < "gene_to_orf_map.txt") > 0) {
                split(line, fields, "\t")
                map[fields[1]] = fields[2]
            }
            close("gene_to_orf_map.txt")
        }
        /^>/ {
            # Extract gene ID from header
            if (match($0, />([^ ]+)/, gene)) {
                gene_id = gene[1]
                if (gene_id in map) {
                    # Replace gene ID with ORF name in header
                    gsub(gene_id, map[gene_id], $0)
                }
            }
        }
        { print }
        ' "${SAMPLE}.faa.backup" > "${SAMPLE}.faa"
        
        # Update GenBank file if AGAT is available
        if command -v agat_sp_gff2gbk.pl &> /dev/null; then
            echo "Regenerating GenBank file with ORF names..."
            agat_sp_gff2gbk.pl --gff "${SAMPLE}.gff" --fasta "${SAMPLE}.fna" -o "${SAMPLE}.gbk" 2>/dev/null
            
            # Verify GBK file was created successfully
            if [ -f "${SAMPLE}.gbk" ]; then
                echo "GenBank file successfully updated with ORF names"
            else
                echo "Warning: GenBank file generation failed"
            fi
        else
            echo "Warning: AGAT not available - GenBank file not updated"
            echo "Install with: conda install -c bioconda agat"
        fi
        
        # Create mapping report
        echo "=== ORF MAPPING REPORT ===" > "${SAMPLE}_orf_mapping.txt"
        echo "Sample: $SAMPLE" >> "${SAMPLE}_orf_mapping.txt"
        echo "Genes mapped to ORFs: $MAPPED_GENES" >> "${SAMPLE}_orf_mapping.txt"
        echo "" >> "${SAMPLE}_orf_mapping.txt"
        echo "Prokka_Gene\tORF_Name\tIdentity%" >> "${SAMPLE}_orf_mapping.txt"
        cat gene_to_orf_map.txt >> "${SAMPLE}_orf_mapping.txt"
        
        echo "Successfully updated annotations with ORF names!"
        echo "Mapping report: ${SAMPLE}_orf_mapping.txt"
        
        # Clean up
        rm -f blast_results.tmp gene_to_orf_map.txt orf_db.*
    else
        echo "No high-confidence ORF mappings found. Keeping original Prokka annotations."
        rm -f blast_results.tmp gene_to_orf_map.txt orf_db.*
    fi
    
elif [ ! -f "$ORF_PROTEIN_FILE" ]; then
    echo "Skipping ORF mapping: No ORF protein file found"
elif ! command -v blastp &> /dev/null; then
    echo "Skipping ORF mapping: BLAST not available"
else
    echo "Skipping ORF mapping: Required files not found"
fi

# Show annotation summary
echo ""
echo "=== ANNOTATION SUMMARY ==="
if [ -f "${SAMPLE}.txt" ]; then
    cat "${SAMPLE}.txt"
fi

echo ""
echo "=== OUTPUT FILES ==="
echo "  - ${SAMPLE}.gff: Gene annotations (with ORF names)"
echo "  - ${SAMPLE}.gbk: GenBank file (with ORF names)"  
echo "  - ${SAMPLE}.faa: Protein sequences (with ORF names)"
echo "  - ${SAMPLE}.fna: Nucleotide sequences"
if [ -f "${SAMPLE}_orf_mapping.txt" ]; then
    echo "  - ${SAMPLE}_orf_mapping.txt: ORF mapping report"
fi

echo ""
echo "Annotation pipeline completed!"