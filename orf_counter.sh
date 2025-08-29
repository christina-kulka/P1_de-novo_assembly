#!/bin/bash

# Simple ORF Counter Script
# Usage: ./count_orfs.sh SAMPLE_NAME [ASSEMBLY_TYPE]

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Parse command line arguments
SAMPLE=$1
ASSEMBLY_TYPE=${2:-"canu_super"}

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    echo "Example: $0 B006"
    echo "Example: $0 B006 canu_super"
    exit 1
fi

# Define assembly file paths using config
case $ASSEMBLY_TYPE in
    "canu")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_ultra_output/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra_trimmed")
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed.fasta"
        ;;
    "microsynth")
        ASSEMBLY_FILE="${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"
        ;;
    *)
        echo "Error: Unknown assembly type '$ASSEMBLY_TYPE'"
        exit 1
        ;;
esac

# Check if files exist
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

# Use ORF reference from config
ORF_REF_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
if [ ! -f "$ORF_REF_FILE" ]; then
    echo "Error: ORF reference file not found: $ORF_REF_FILE"
    exit 1
fi

# Create output directory using config
OUTPUT_DIR="${ITR_ANALYSIS_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
mkdir -p "$OUTPUT_DIR"

echo "=== ORF COUNTER ==="
echo "Sample: $SAMPLE"
echo "Assembly: $ASSEMBLY_TYPE"
echo "Assembly file: $ASSEMBLY_FILE"
echo "ORF reference: $ORF_REF_FILE"
echo ""

cd "$OUTPUT_DIR"

# Get assembly stats
ASSEMBLY_LENGTH=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)
TOTAL_ORFS_IN_REF=$(grep -c "^>" "$ORF_REF_FILE")

echo "Assembly length: $ASSEMBLY_LENGTH bp"
echo "Total ORFs in reference: $TOTAL_ORFS_IN_REF"
echo ""

# Create BLAST database
echo "Creating BLAST database..."
$MAKEBLASTDB -in "$ASSEMBLY_FILE" -dbtype nucl -out assembly_db -logfile makeblastdb.log

# BLAST ORFs against assembly
echo "Searching ORFs in assembly..."
$BLASTN -query "$ORF_REF_FILE" \
        -db assembly_db \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -num_threads $THREADS \
        -out orf_blast_results.txt

# Create output file
OUTPUT_FILE="orf_counts_${SAMPLE}_${ASSEMBLY_TYPE}.txt"

echo "# ORF Count Analysis" > "$OUTPUT_FILE"
echo "# Sample: $SAMPLE" >> "$OUTPUT_FILE"
echo "# Assembly: $ASSEMBLY_TYPE" >> "$OUTPUT_FILE"
echo "# Assembly length: $ASSEMBLY_LENGTH bp" >> "$OUTPUT_FILE"
echo "# Total ORFs in reference: $TOTAL_ORFS_IN_REF" >> "$OUTPUT_FILE"
echo "# Generated: $(date)" >> "$OUTPUT_FILE"
echo "#" >> "$OUTPUT_FILE"

# Filter high-quality hits and count occurrences
echo "# High-quality hits (≥95% identity, ≥90% coverage)" >> "$OUTPUT_FILE"
echo "# Format: ORF_NAME COUNT POSITIONS" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

awk '$3 >= 95 && $13 >= 90 {
    start = ($9 < $10) ? $9 : $10
    end = ($9 < $10) ? $10 : $9
    orientation = ($9 < $10) ? "+" : "-"
    base_orf = $1
    gsub(/_[FR]$/, "", base_orf)
    print base_orf "\t" start "-" end "(" orientation ")"
}' orf_blast_results.txt | sort -u | cut -f1 | uniq -c | sort -nr > temp_counts.txt

# Process and format results
echo "Processing results..."
found_orfs=0
total_occurrences=0
duplicated_orfs=0

while read count positions orf_name; do
    found_orfs=$((found_orfs + 1))
    total_occurrences=$((total_occurrences + count))
    
    if [ $count -gt 1 ]; then
        duplicated_orfs=$((duplicated_orfs + 1))
        status="[DUPLICATED]"
    else
        status=""
    fi
    
    echo "$orf_name $count $positions $status" >> "$OUTPUT_FILE"
done < temp_counts.txt

# Add summary statistics
echo "" >> "$OUTPUT_FILE"
echo "# SUMMARY STATISTICS" >> "$OUTPUT_FILE"
echo "# Unique ORFs found: $found_orfs" >> "$OUTPUT_FILE"
echo "# Total ORF occurrences: $total_occurrences" >> "$OUTPUT_FILE"
echo "# ORFs with multiple copies: $duplicated_orfs" >> "$OUTPUT_FILE"

missing_orfs=$((TOTAL_ORFS_IN_REF - found_orfs))
echo "# ORFs not found: $missing_orfs" >> "$OUTPUT_FILE"

if [ $total_occurrences -gt $found_orfs ]; then
    concatemer_ratio=$(echo "scale=2; $total_occurrences / $found_orfs" | bc -l)
    echo "# Average copies per ORF: $concatemer_ratio" >> "$OUTPUT_FILE"
    if (( $(echo "$concatemer_ratio > 1.5" | bc -l) )); then
        echo "# WARNING: High duplication suggests concatemers" >> "$OUTPUT_FILE"
    fi
fi

# Create simple summary for console
echo "=== RESULTS ==="
echo "Unique ORFs found: $found_orfs/$TOTAL_ORFS_IN_REF"
echo "Total occurrences: $total_occurrences"
echo "Duplicated ORFs: $duplicated_orfs"

if [ $duplicated_orfs -gt 0 ]; then
    echo ""
    echo "Most duplicated ORFs:"
    head -n5 temp_counts.txt | while read count positions orf_name; do
        if [ $count -gt 1 ]; then
            echo "  $orf_name: $count copies"
        fi
    done
fi

echo ""
echo "Results saved to: $OUTPUT_FILE"

# Clean up
rm -f temp_counts.txt assembly_db.n*

echo "=== ANALYSIS COMPLETE ==="