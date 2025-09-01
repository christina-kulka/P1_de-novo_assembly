#!/bin/bash

# ITR Analysis Script (Step 6)
# Usage: ./06_itr_analysis.sh SAMPLE_NAME [ASSEMBLY_TYPE]

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Parse command line arguments
SAMPLE=$1
ASSEMBLY_TYPE=${2:-"canu_super"}  # Default to canu_super

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    echo "Example: $0 B006"
    echo "Example: $0 B006 canu"
    exit 1
fi

# Define assembly file paths based on type
    #"miniasm_viral")
    #    ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta"
    #    ;;
    #"miniasm_raw")
    #    ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
    #    ;;
    #"canu_super")
    #    ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.longest_contig.fasta"
    #    ;;
case $ASSEMBLY_TYPE in
    "canu")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_ultra_output/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra_trimmed")
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_150.fasta"
        ;;
    "microsynth")
        ASSEMBLY_FILE="${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"
        ;;
    *)
        echo "Error: Unknown assembly type '$ASSEMBLY_TYPE'"
        exit 1
        ;;
esac

# Check if required files exist
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
PROKKA_GFF="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}/${SAMPLE}.gff"
PROKKA_FAA="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}/${SAMPLE}.faa"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

if [ ! -f "$ORF_REFERENCE_FILE" ]; then
    echo "Error: ORF reference file not found: $ORF_REFERENCE_FILE"
    exit 1
fi

if [ ! -f "$PROKKA_GFF" ]; then
    echo "Error: Prokka GFF file not found: $PROKKA_GFF"
    echo "Please run prokka annotation first"
    exit 1
fi

# Create output directory
OUTPUT_DIR="${ITR_ANALYSIS_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
mkdir -p "$OUTPUT_DIR"

echo "=== ITR ANALYSIS FOR $SAMPLE ($ASSEMBLY_TYPE) ==="
echo "Assembly file: $ASSEMBLY_FILE"
echo "ORF reference: $ORF_REFERENCE_FILE"
echo "Output directory: $OUTPUT_DIR"

cd "$OUTPUT_DIR"

# Get assembly length
echo "Getting assembly statistics..."
ASSEMBLY_LENGTH=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)
CONTIG_COUNT=$(grep -c "^>" "$ASSEMBLY_FILE")

echo "Assembly length: $ASSEMBLY_LENGTH bp"
echo "Number of contigs: $CONTIG_COUNT"

# Extract terminal regions for hairpin search
echo "Extracting terminal regions..."
head -n1 "$ASSEMBLY_FILE" > header.tmp
grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' > sequence.tmp

# Left terminal region
head -c $HAIRPIN_SEARCH_REGION sequence.tmp > left_terminal.fasta
echo ">Left_terminal_${HAIRPIN_SEARCH_REGION}bp" > left_terminal_formatted.fasta
echo $(cat left_terminal.fasta) >> left_terminal_formatted.fasta

# Right terminal region  
tail -c $HAIRPIN_SEARCH_REGION sequence.tmp > right_terminal.fasta
echo ">Right_terminal_${HAIRPIN_SEARCH_REGION}bp" > right_terminal_formatted.fasta
echo $(cat right_terminal.fasta) >> right_terminal_formatted.fasta

# Search for hairpin structures using RNAfold
echo "Searching for hairpin structures with RNAfold..."
echo "=== HAIRPIN ANALYSIS ===" > hairpin_report.txt
echo "Search parameters:" >> hairpin_report.txt
echo "  Search region: $HAIRPIN_SEARCH_REGION bp from each end" >> hairpin_report.txt
echo "  Tool: Vienna RNAfold" >> hairpin_report.txt
echo "" >> hairpin_report.txt

# Search for telomere resolution sequence (CRS) patterns
echo "Searching for telomere resolution sequences (CRS)..." >> hairpin_report.txt

# Known CRS consensus: 5'-T/A-T6-N8-TAAAT-3'
# More flexible pattern for detection
crs_patterns=(
    "TTTTTT.{8}TAAAT"    # T6-N8-TAAAT
    "ATTTTT.{8}TAAAT"    # A-T5-N8-TAAAT  
    "TTTTTTT.{7}TAAAT"   # T7-N7-TAAAT (slightly different)
)

full_sequence=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n')


# Search for similar sequences to known hairpin with fuzzy matching
echo "Searching for sequences similar to known hairpin..." >> hairpin_report.txt
known_hairpin="GGAGGCCGTCCTCCCTCCAAAACTTTTCGTAAAATCTCTTCGGAGGCCGTCCTCCCTCCAAAACTTTTCGTAAAATCT"

# First try exact match
position=$(echo "$full_sequence" | grep -b -o "$known_hairpin" | cut -d: -f1)

if [ ! -z "$position" ]; then
    end_position=$((position + ${#known_hairpin}))
    echo "Known hairpin FOUND (exact match) at position: $((position + 1))-$end_position" >> hairpin_report.txt
else
    echo "Exact hairpin not found - searching for similar structures..." >> hairpin_report.txt
    
    # Method 1: Search for shorter characteristic motifs (more flexible)
    # Break the hairpin into meaningful chunks
    motifs=(
        "GGAGGCCGTCCTCCCTCC"      # First part of repeat
        "AAAACTTTTCGTAAAATCT"     # Second part of repeat
        "CTTCGGAGGCCGTCCTCC"      # Bridge region
        "GGAGGCCGTCCTCCCTCCAAAACTTTTCGTAAAATCT"  # Full repeat unit
    )
    
    echo "  Searching for characteristic motifs:" >> hairpin_report.txt
    
    motif_positions=()
    for motif in "${motifs[@]}"; do
        positions=$(echo "$full_sequence" | grep -aob "$motif" | cut -d: -f1)
        if [ ! -z "$positions" ]; then
            echo "    Motif '$motif' found at:" >> hairpin_report.txt
            for pos in $positions; do
                echo "      Position $((pos+1))" >> hairpin_report.txt
                motif_positions+=($pos)
            done
        fi
    done
    
    # Method 2: Look for sequences with similar length and structure
    echo "  Searching for similar-length AT-rich regions with hairpin potential:" >> hairpin_report.txt
    
    # Search for regions of similar length (±20bp) that are AT-rich
    target_length=${#known_hairpin}
    min_length=$((target_length - 20))
    max_length=$((target_length + 20))
    
    # Use sliding window to find AT-rich regions of similar length
    candidates_found=0
    for start in $(seq 1 50 $((${#full_sequence} - max_length))); do
        for length in $(seq $min_length 10 $max_length); do
            end=$((start + length - 1))
            if [ $end -gt ${#full_sequence} ]; then
                break
            fi
            
            candidate_seq=$(echo "$full_sequence" | cut -c$start-$end)
            
            # Check AT content (should be similar to known hairpin)
            at_count=$(echo "$candidate_seq" | grep -o '[AT]' | wc -l)
            at_content=$((at_count * 100 / length))
            
            # Calculate AT content of known hairpin for comparison
            known_at_count=$(echo "$known_hairpin" | grep -o '[AT]' | wc -l)
            known_at_content=$((known_at_count * 100 / ${#known_hairpin}))
            
            # If AT content is similar (±15%), check if it can form hairpins
            at_diff=$((at_content > known_at_content ? at_content - known_at_content : known_at_content - at_content))
            
            if [ $at_diff -le 15 ]; then
                # Test if this sequence can form a hairpin structure
                result=$(echo "$candidate_seq" | RNAfold --noPS 2>/dev/null)
                structure=$(echo "$result" | tail -n1 | awk '{print $1}')
                energy=$(echo "$result" | tail -n1 | awk '{print $2}' | tr -d '()')
                stem_count=$(echo "$structure" | grep -o '(' | wc -l)
                
                # If it forms a decent hairpin structure
                if [ $stem_count -ge 8 ] && (( $(echo "$energy < -10" | bc -l) )); then
                    candidates_found=$((candidates_found + 1))
                    
                    echo "    Candidate similar hairpin at position $start-$end:" >> hairpin_report.txt
                    echo "      Length: ${length}bp (target: ${target_length}bp)" >> hairpin_report.txt
                    echo "      AT content: ${at_content}% (target: ${known_at_content}%)" >> hairpin_report.txt
                    echo "      Structure: $structure" >> hairpin_report.txt
                    echo "      Energy: $energy kcal/mol" >> hairpin_report.txt
                    echo "      Stem pairs: $stem_count" >> hairpin_report.txt
                    
                    # Calculate rough sequence similarity (simple approach)
                    # Count matching positions in overlapping regions
                    if [ $length -ge $target_length ]; then
                        compare_seq=$(echo "$candidate_seq" | cut -c1-$target_length)
                    else
                        compare_seq=$candidate_seq
                    fi
                    
                    # Simple character-by-character comparison (first 50 chars for speed)
                    compare_len=$((${#compare_seq} < 50 ? ${#compare_seq} : 50))
                    known_compare=$(echo "$known_hairpin" | cut -c1-$compare_len)
                    candidate_compare=$(echo "$compare_seq" | cut -c1-$compare_len)
                    
                    matches=0
                    for i in $(seq 1 $compare_len); do
                        known_char=$(echo "$known_compare" | cut -c$i)
                        candidate_char=$(echo "$candidate_compare" | cut -c$i)
                        if [ "$known_char" = "$candidate_char" ]; then
                            matches=$((matches + 1))
                        fi
                    done
                    
                    similarity=$((matches * 100 / compare_len))
                    echo "      Sequence similarity (first ${compare_len}bp): ${similarity}%" >> hairpin_report.txt
                    echo "" >> hairpin_report.txt
                    
                    # Stop after finding 10 candidates to avoid too much output
                    if [ $candidates_found -ge 10 ]; then
                        echo "    (Stopping after 10 candidates to avoid excessive output)" >> hairpin_report.txt
                        break 2
                    fi
                fi
            fi
        done
    done
    
    if [ $candidates_found -eq 0 ]; then
        echo "    No structurally similar hairpin candidates found" >> hairpin_report.txt
    fi
fi

# Predict structure for reference hairpin
echo "Structure prediction for reference hairpin:" >> hairpin_report.txt
echo "$known_hairpin" | RNAfold --noPS >> hairpin_report.txt 2>/dev/null
echo "" >> hairpin_report.txt

# Predict structure for any found patterns
if [ ! -z "$position" ] || [ ! -z "$positions" ]; then
    echo "Structure prediction for reference hairpin:" >> hairpin_report.txt
    echo "$known_hairpin" | RNAfold --noPS >> hairpin_report.txt 2>/dev/null
fi
echo "" >> hairpin_report.txt



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
# Clean up temporary files
rm -f high_quality_hits.txt orf_counts_detailed.txt

# Generate summary report
echo "=== SUMMARY REPORT ===" > summary_report.txt
echo "Sample: $SAMPLE" >> summary_report.txt
echo "Assembly type: $ASSEMBLY_TYPE" >> summary_report.txt
echo "Assembly length: $ASSEMBLY_LENGTH bp" >> summary_report.txt
echo "Number of contigs: $CONTIG_COUNT" >> summary_report.txt
echo "Unique ORFs found: $unique_orfs" >> summary_report.txt
echo "" >> summary_report.txt

# Count duplicated ORFs
duplicated_count=$(awk '$1 > 1' orf_counts.txt | wc -l)
echo "ORFs appearing multiple times: $duplicated_count" >> summary_report.txt

if [ $duplicated_count -gt 0 ]; then
    echo "WARNING: Multiple ORF copies detected - possible concatemers" >> summary_report.txt
    echo "See orf_analysis.txt for details" >> summary_report.txt
fi

echo "" >> summary_report.txt
echo "Hairpin analysis completed - see hairpin_report.txt" >> summary_report.txt

# Clean up temporary files
rm -f header.tmp sequence.tmp left_terminal.fasta right_terminal.fasta
rm -f assembly_db.n* left_hairpins.txt right_hairpins.txt

echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "Output files created in: $OUTPUT_DIR"
echo "  - summary_report.txt: Overview of findings"
echo "  - orf_analysis.txt: Detailed ORF analysis"
echo "  - hairpin_report.txt: Hairpin structure analysis"
echo "  - orf_hits.txt: Raw BLAST results"
echo "  - orf_counts.txt: ORF occurrence counts"
echo ""
echo "Summary:"
cat "$OUTPUT_DIR/summary_report.txt"