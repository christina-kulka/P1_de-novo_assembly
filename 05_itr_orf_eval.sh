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

# Function to find hairpins around CRS positions
find_hairpins_around_crs() {
    local file=$1
    local offset=$2  # Offset to adjust coordinates back to assembly
    local output_file=$3
    
    # Extract sequence without header
    seq=$(grep -v "^>" $file | tr -d '\n')
    
    window_size=100
    step_size=20
    seq_len=${#seq}
    
    for start in $(seq 1 $step_size $((seq_len - window_size + 1))); do
        end=$((start + window_size - 1))
        if [ $end -gt $seq_len ]; then
            end=$seq_len
        fi
        
        subseq=${seq:$((start-1)):$window_size}
        
        # Calculate AT content
        at_count=$(echo "$subseq" | grep -o '[AT]' | wc -l)
        at_content=$((at_count * 100 / ${#subseq}))
        
        # Run RNAfold
        result=$(echo "$subseq" | RNAfold --noPS 2>/dev/null)
        structure=$(echo "$result" | tail -n1 | awk '{print $1}')
        energy=$(echo "$result" | tail -n1 | awk '{print $2}' | tr -d '()')
        
        # Count hairpin stems
        stem_count=$(echo "$structure" | grep -o '(' | wc -l)
        
        # Lower threshold for CRS-associated hairpins
        if [ $stem_count -ge 6 ] && (( $(echo "$energy < -4" | bc -l) )); then
            # Adjust coordinates back to assembly position
            assembly_start=$((offset + start))
            assembly_end=$((offset + end))
            
            echo "    Hairpin near CRS at position $assembly_start-$assembly_end:" >> $output_file
            echo "      Structure: $structure" >> $output_file
            echo "      Energy: $energy kcal/mol" >> $output_file
            echo "      Stem pairs: $stem_count" >> $output_file
            echo "      AT content: ${at_content}%" >> $output_file
            echo "" >> $output_file
        fi
    done
}

# Function to find hairpins with RNAfold - IMPROVED VERSION
find_hairpins_rnafold() {
    local file=$1
    local region=$2
    local output_file=$3
    
    echo "Analyzing $region region for parapoxvirus-like hairpins..." >> $output_file
    
    # Extract sequence without header
    seq=$(grep -v "^>" $file | tr -d '\n')
    
    # Use sliding window optimized for parapoxvirus hairpins
    window_size=150  # Reduced since hairpins are ~100-120bp
    step_size=25     # Smaller steps for better coverage
    seq_len=${#seq}
    
    # Debug info
    echo "  Sequence length: $seq_len bp" >> $output_file
    
    candidates_found=0
    at_filtered=0
    
    for start in $(seq 1 $step_size $((seq_len - window_size + 1))); do
        end=$((start + window_size - 1))
        if [ $end -gt $seq_len ]; then
            end=$seq_len
        fi
        
        subseq=${seq:$((start-1)):$window_size}
        
        # Calculate AT content (should be high for parapoxvirus hairpins)
        at_count=$(echo "$subseq" | grep -o '[AT]' | wc -l)
        at_content=$((at_count * 100 / ${#subseq}))
        
        # Lower AT threshold and report skipped sequences
        if [ $at_content -lt 50 ]; then
            at_filtered=$((at_filtered + 1))
            continue
        fi
        
        # Run RNAfold
        result=$(echo "$subseq" | RNAfold --noPS 2>/dev/null)
        structure=$(echo "$result" | tail -n1 | awk '{print $1}')
        energy=$(echo "$result" | tail -n1 | awk '{print $2}' | tr -d '()')
        
        # Count hairpin stems
        stem_count=$(echo "$structure" | grep -o '(' | wc -l)
        
        # More lenient criteria for detection
        # - 8+ base pairs (lowered from 10)
        # - AT-rich (>50%, lowered from 60%)  
        # - Any reasonable folding energy
        if [ $stem_count -ge 8 ] && (( $(echo "$energy < -5" | bc -l) )); then
            candidates_found=$((candidates_found + 1))
            
            # Mark high-confidence vs low-confidence
            confidence="LOW"
            if [ $stem_count -ge 10 ] && [ $stem_count -le 30 ] && [ $at_content -ge 60 ] && (( $(echo "$energy < -8" | bc -l) )); then
                confidence="HIGH"
            fi
            
            echo "  Potential hairpin at position $start-$end ($confidence confidence):" >> $output_file
            echo "    Structure: $structure" >> $output_file
            echo "    Energy: $energy kcal/mol" >> $output_file
            echo "    Stem pairs: $stem_count" >> $output_file
            echo "    AT content: ${at_content}%" >> $output_file
            echo "" >> $output_file
        fi
    done
    
    # Summary
    echo "  Summary: $candidates_found candidates found, $at_filtered regions filtered by AT content" >> $output_file
    echo "" >> $output_file
}

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

# Also search around CRS positions for hairpins
crs_positions=()

for pattern in "${crs_patterns[@]}"; do
    positions=$(echo "$full_sequence" | grep -aob -E "$pattern" | cut -d: -f1)
    if [ ! -z "$positions" ]; then
        echo "CRS-like pattern '$pattern' found at positions:" >> hairpin_report.txt
        for pos in $positions; do
            # Store CRS positions for later hairpin analysis
            crs_positions+=($pos)
            
            # Extract the matching sequence
            match_seq=$(echo "$full_sequence" | cut -c$((pos+1))-$((pos+20)))
            echo "  Position $((pos+1)): $match_seq" >> hairpin_report.txt
            
            # Check if near assembly ends (typical for terminal hairpins)
            distance_from_start=$((pos+1))
            distance_from_end=$((ASSEMBLY_LENGTH - pos - 20))
            echo "    (${distance_from_start}bp from start, ${distance_from_end}bp from end)" >> hairpin_report.txt
        done
        echo "" >> hairpin_report.txt
    fi
done

# Search for hairpins around CRS positions
if [ ${#crs_positions[@]} -gt 0 ]; then
    echo "Analyzing regions around CRS positions for hairpins..." >> hairpin_report.txt
    
    for crs_pos in "${crs_positions[@]}"; do
        # Extract 500bp region around CRS
        region_start=$((crs_pos > 250 ? crs_pos - 250 : 0))
        region_end=$((crs_pos + 250))
        if [ $region_end -gt $ASSEMBLY_LENGTH ]; then
            region_end=$ASSEMBLY_LENGTH
        fi
        
        region_seq=$(echo "$full_sequence" | cut -c$((region_start+1))-$region_end)
        
        # Create temp file for this region
        echo ">CRS_region_${crs_pos}" > crs_region.tmp
        echo "$region_seq" >> crs_region.tmp
        
        echo "  Analyzing 500bp region around CRS at position $((crs_pos+1)):" >> hairpin_report.txt
        
        # Use the same hairpin finding function but adjust coordinates
        find_hairpins_around_crs "crs_region.tmp" "$region_start" "hairpin_report.txt"
        
        rm -f crs_region.tmp
    done
fi

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


# Filter hairpins near ORF 134
echo "Filtering hairpins near ORF 134..." >> hairpin_report.txt

# Find ORF 134 position from BLAST results
orf134_position=""
if grep -qi "134" orf_hits.txt; then
    orf134_line=$(grep -i "134" orf_hits.txt | head -n1)
    orf134_start=$(echo "$orf134_line" | awk '{print $9}')
    orf134_end=$(echo "$orf134_line" | awk '{print $10}')
    orf134_position="$orf134_start-$orf134_end"
    echo "ORF 134 found at position: $orf134_position" >> hairpin_report.txt
else
    echo "ORF 134 not found - checking terminal ORFs instead" >> hairpin_report.txt
    # Get ORFs in terminal regions (first/last 10kb)
    terminal_threshold=10000
    awk -v thresh=$terminal_threshold -v len=$ASSEMBLY_LENGTH '
        $3 >= 95 && $13 >= 90 && ($9 < thresh || $9 > (len - thresh)) {
            print $1, $9, $10
        }' orf_hits.txt > terminal_orfs.txt
    
    if [ -s terminal_orfs.txt ]; then
        echo "Terminal ORFs found:" >> hairpin_report.txt
        cat terminal_orfs.txt >> hairpin_report.txt
    fi
fi

# Define proximity threshold (e.g., within 2kb of ORF)

echo "" >> hairpin_report.txt
echo "Hairpins within ${PROXIMITY_THRESHOLD}bp of ORF 134 or terminal ORFs:" >> hairpin_report.txt

# Parse hairpin positions and check proximity
grep -n "Potential hairpin at position" hairpin_report.txt | while read line; do
    # Extract hairpin coordinates
    hairpin_coords=$(echo "$line" | grep -o '[0-9]*-[0-9]*')
    hairpin_start=$(echo "$hairpin_coords" | cut -d'-' -f1)
    hairpin_end=$(echo "$hairpin_coords" | cut -d'-' -f2)
    
    # Check if near ORF 134
    if [ ! -z "$orf134_position" ]; then
        orf_start=$(echo "$orf134_position" | cut -d'-' -f1)
        orf_end=$(echo "$orf134_position" | cut -d'-' -f2)
        
        # Calculate distance
        distance_to_start=$((hairpin_start > orf_start ? hairpin_start - orf_end : orf_start - hairpin_end))
        
        if [ $distance_to_start -lt $PROXIMITY_THRESHOLD ]; then
            echo "  CLOSE TO ORF 134: $line" >> hairpin_report.txt
        fi
    fi
    
    # Also check proximity to assembly ends (potential cutting sites)
    distance_to_left_end=$hairpin_start
    distance_to_right_end=$((ASSEMBLY_LENGTH - hairpin_end))
    
    if [ $distance_to_left_end -lt 5000 ] || [ $distance_to_right_end -lt 5000 ]; then
        echo "  NEAR ASSEMBLY END: $line (distances: ${distance_to_left_end}bp from start, ${distance_to_right_end}bp from end)" >> hairpin_report.txt
    fi
done




# ORF analysis
echo "Analyzing ORF content..."
echo "=== ORF ANALYSIS ===" > orf_analysis.txt

# Create BLAST database from assembly
makeblastdb -in "$ASSEMBLY_FILE" -dbtype nucl -out assembly_db

# BLAST ORFs against assembly
blastn -query "$ORF_REFERENCE_FILE" \
       -db assembly_db \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
       -out orf_hits.txt

# Count ORF occurrences
echo "ORF occurrence analysis:" >> orf_analysis.txt
echo "Total ORFs in reference: $(grep -c "^>" $ORF_REFERENCE_FILE)" >> orf_analysis.txt
echo "" >> orf_analysis.txt

# Filter high-quality hits first (>95% identity, >90% coverage)
awk '$3 >= 95 && $13 >= 90 {print $0}' orf_hits.txt > high_quality_hits.txt

# Create detailed ORF occurrence report with positions and orientations
echo "ORFs found multiple times (potential concatemers):" >> orf_analysis.txt
awk '$3 >= 95 && $13 >= 90 {
    orientation = ($9 < $10) ? "forward" : "reverse"
    start = ($9 < $10) ? $9 : $10
    end = ($9 < $10) ? $10 : $9
    print $1 "\t" orientation "\t" start "-" end "\t" $3 "%" 
}' orf_hits.txt | sort | uniq -c | sort -nr > orf_counts_detailed.txt

# Report duplicated ORFs with full details
while read count orf orientation positions identity; do
    if [ $count -gt 1 ]; then
        echo "  $orf: $count occurrences" >> orf_analysis.txt
        
        # Get all occurrences of this ORF
        echo "    Positions and orientations:" >> orf_analysis.txt
        awk -v orf="$orf" '$1 == orf && $3 >= 95 && $13 >= 90 {
            orientation = ($9 < $10) ? "forward" : "reverse"
            start = ($9 < $10) ? $9 : $10
            end = ($9 < $10) ? $10 : $9
            printf "      %s: %d-%d (%s, %.1f%% identity, %.1f%% coverage)\n", orientation, start, end, orientation, $3, $13
        }' orf_hits.txt >> orf_analysis.txt
        echo "" >> orf_analysis.txt
    fi
done < orf_counts_detailed.txt

# Count total unique ORFs found (not occurrences, but unique ORF types)
unique_orfs=$(awk '$3 >= 95 && $13 >= 90 {print $1}' orf_hits.txt | sort -u | wc -l)
total_occurrences=$(awk '$3 >= 95 && $13 >= 90' orf_hits.txt | wc -l)

echo "Unique ORFs found in assembly: $unique_orfs" >> orf_analysis.txt
echo "Total ORF occurrences: $total_occurrences" >> orf_analysis.txt

# Calculate concatemer evidence
duplicated_count=$(awk '$1 > 1' orf_counts_detailed.txt | wc -l)
if [ $duplicated_count -gt 0 ]; then
    echo "ORFs with multiple copies: $duplicated_count" >> orf_analysis.txt
    concatemer_ratio=$(echo "scale=2; $total_occurrences / $unique_orfs" | bc -l)
    echo "Average copies per ORF: $concatemer_ratio" >> orf_analysis.txt
    
    if (( $(echo "$concatemer_ratio > 1.5" | bc -l) )); then
        echo "WARNING: High concatemer ratio suggests genome concatenation" >> orf_analysis.txt
    fi
fi

echo "" >> orf_analysis.txt

# Look for ORF 134 specifically with all occurrences
echo "ORF 134 analysis:" >> orf_analysis.txt
orf134_hits=$(grep -i "134" orf_hits.txt | awk '$3 >= 95 && $13 >= 90')

if [ ! -z "$orf134_hits" ]; then
    echo "ORF 134 found:" >> orf_analysis.txt
    echo "$orf134_hits" | while read line; do
        orf_name=$(echo "$line" | awk '{print $1}')
        start_pos=$(echo "$line" | awk '{print $9}')
        end_pos=$(echo "$line" | awk '{print $10}')
        identity=$(echo "$line" | awk '{print $3}')
        coverage=$(echo "$line" | awk '{print $13}')
        
        if [ $start_pos -lt $end_pos ]; then
            orientation="forward"
            position="${start_pos}-${end_pos}"
        else
            orientation="reverse"
            position="${end_pos}-${start_pos}"
        fi
        
        echo "  $orf_name: $position ($orientation, ${identity}% identity, ${coverage}% coverage)" >> orf_analysis.txt
        
        # Check if near assembly ends
        distance_from_start=$([ $start_pos -lt $end_pos ] && echo $start_pos || echo $end_pos)
        distance_from_end=$((ASSEMBLY_LENGTH - $([ $start_pos -gt $end_pos ] && echo $start_pos || echo $end_pos)))
        
        echo "    Distance from assembly start: ${distance_from_start}bp" >> orf_analysis.txt
        echo "    Distance from assembly end: ${distance_from_end}bp" >> orf_analysis.txt
        echo "" >> orf_analysis.txt
    done
else
    echo "  ORF 134 not found with high confidence (>95% identity, >90% coverage)" >> orf_analysis.txt
    echo "  Checking for terminal ORFs..." >> orf_analysis.txt
    
    # Get ORFs in terminal regions (first/last 10kb)
    terminal_threshold=10000
    echo "  ORFs in terminal regions (first/last ${terminal_threshold}bp):" >> orf_analysis.txt
    
    awk -v thresh=$terminal_threshold -v len=$ASSEMBLY_LENGTH '$3 >= 95 && $13 >= 90 {
        start = ($9 < $10) ? $9 : $10
        end = ($9 < $10) ? $10 : $9
        orientation = ($9 < $10) ? "forward" : "reverse"
        
        if (start < thresh || end > (len - thresh)) {
            printf "    %s: %d-%d (%s, %.1f%% identity, %.1f%% coverage)\n", $1, start, end, orientation, $3, $13
        }
    }' orf_hits.txt >> orf_analysis.txt
fi

# Analysis of ORF distribution along genome
echo "" >> orf_analysis.txt
echo "ORF distribution analysis:" >> orf_analysis.txt

# Divide genome into 10 segments and count ORFs in each
segment_size=$((ASSEMBLY_LENGTH / 10))
echo "Genome divided into 10 segments of ~${segment_size}bp each:" >> orf_analysis.txt

for i in $(seq 0 9); do
    segment_start=$((i * segment_size))
    segment_end=$(((i + 1) * segment_size))
    
    if [ $i -eq 9 ]; then
        segment_end=$ASSEMBLY_LENGTH  # Ensure last segment goes to end
    fi
    
    # Count ORFs in this segment
    orf_count=$(awk -v start=$segment_start -v end=$segment_end '$3 >= 95 && $13 >= 90 {
        hit_start = ($9 < $10) ? $9 : $10
        hit_end = ($9 < $10) ? $10 : $9
        if (hit_start >= start && hit_end <= end) print $0
    }' orf_hits.txt | wc -l)
    
    echo "  Segment $((i+1)) (${segment_start}-${segment_end}): $orf_count ORFs" >> orf_analysis.txt
done

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