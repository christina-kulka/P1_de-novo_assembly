#!/bin/bash

# CRS-Based Genome Trimmer Script
# Usage: ./trim_genome_by_crs.sh SAMPLE_NAME [ASSEMBLY_TYPE] [MIN_DISTANCE]
# Trims genome assemblies based on CRS positions found in hairpin reports

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Parse command line arguments
SAMPLE=$1
ASSEMBLY_TYPE=${2:-"canu_super"}
MIN_DISTANCE=${3:-$MIN_CRS_DISTANCE}  # Minimum distance between CRS positions to consider them at different ends

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE] [MIN_DISTANCE]"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME   - Sample identifier (e.g., B006)"
    echo "  ASSEMBLY_TYPE - Assembly type (default: canu_super)" 
    echo "  MIN_DISTANCE  - Minimum distance between CRS to consider different ends (default: ${MIN_CRS_DISTANCE}bp)"
    echo ""
    echo "Examples:"
    echo "  $0 B006"
    echo "  $0 B006 canu_super 15000"
    exit 1
fi

# Define file paths
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
    "canu_super")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.longest_contig.fasta"
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
ITR_DIR="${ITR_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
HAIRPIN_REPORT="${ITR_DIR}/hairpin_report.txt"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

if [ ! -f "$HAIRPIN_REPORT" ]; then
    echo "Error: Hairpin report not found: $HAIRPIN_REPORT"
    echo "Please run ITR analysis first to generate CRS positions"
    exit 1
fi

# Get assembly length
ASSEMBLY_LENGTH=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)

echo "=== CRS-BASED GENOME TRIMMER ==="
echo "Sample: $SAMPLE"
echo "Assembly type: $ASSEMBLY_TYPE"
echo "Assembly file: $ASSEMBLY_FILE"
echo "Assembly length: $ASSEMBLY_LENGTH bp"
echo "Minimum distance between CRS: $MIN_DISTANCE bp"
echo ""

# Extract CRS positions from hairpin report and deduplicate
echo "Extracting CRS positions..."
crs_positions=$(grep -A1 "Position [0-9]*:" "$HAIRPIN_REPORT" | \
                grep "Position" | \
                grep -o "Position [0-9]*" | \
                cut -d' ' -f2 | \
                sort -n | \
                uniq)

if [ -z "$crs_positions" ]; then
    echo "Error: No CRS positions found in hairpin report"
    echo "Please check that ITR analysis found CRS sequences"
    exit 1
fi

echo "Raw CRS positions found:"
for pos in $crs_positions; do
    distance_from_start=$pos
    distance_from_end=$((ASSEMBLY_LENGTH - pos))
    echo "  Position $pos (${distance_from_start}bp from start, ${distance_from_end}bp from end)"
done
echo ""

# Group CRS positions by proximity to assembly ends
left_crs=()
right_crs=()
middle_crs=()

for pos in $crs_positions; do
    distance_from_start=$pos
    distance_from_end=$((ASSEMBLY_LENGTH - pos))
    
    # Consider positions in first quarter as "left", last quarter as "right"
    quarter_length=$((ASSEMBLY_LENGTH / 4))
    
    if [ $distance_from_start -lt $quarter_length ]; then
        left_crs+=($pos)
    elif [ $distance_from_end -lt $quarter_length ]; then
        right_crs+=($pos)
    else
        middle_crs+=($pos)
    fi
done

echo "CRS position analysis:"
echo "  Left end CRS: ${left_crs[*]}"
echo "  Right end CRS: ${right_crs[*]}"
echo "  Middle CRS: ${middle_crs[*]}"
echo ""

# Determine cutting strategy
cut_start=""
cut_end=""

# Select leftmost CRS for left cutting (if multiple at left end)
if [ ${#left_crs[@]} -gt 0 ]; then
    cut_start=${left_crs[0]}  # Take leftmost (smallest position)
    echo "Will cut from start to position $cut_start (removing left concatemer)"
fi

# Select rightmost CRS for right cutting (if multiple at right end) 
if [ ${#right_crs[@]} -gt 0 ]; then
    # Sort right CRS positions and take rightmost (largest position)
    IFS=$'\n' right_sorted=($(sort -nr <<<"${right_crs[*]}"))
    cut_end=${right_sorted[0]}
    echo "Will cut from position $cut_end to end (removing right concatemer)"
fi

# Handle middle CRS positions
if [ ${#middle_crs[@]} -gt 0 ]; then
    echo "Warning: Found CRS in middle of genome at positions: ${middle_crs[*]}"
    echo "This might indicate complex concatemer structure"
    
    # If no end cuts planned, use middle CRS positions
    if [ -z "$cut_start" ] && [ -z "$cut_end" ]; then
        if [ ${#middle_crs[@]} -eq 1 ]; then
            cut_end=${middle_crs[0]}
            echo "Using middle CRS at $cut_end as right cut point"
        elif [ ${#middle_crs[@]} -eq 2 ]; then
            # Check if they're far apart enough to be different ends
            pos1=${middle_crs[0]}
            pos2=${middle_crs[1]}
            distance=$((pos2 - pos1))
            
            if [ $distance -gt $MIN_DISTANCE ]; then
                cut_start=$pos1
                cut_end=$pos2
                echo "Using middle CRS positions as cut points: $cut_start (left) and $cut_end (right)"
            else
                echo "Middle CRS positions too close ($distance bp) - using rightmost"
                cut_end=$pos2
            fi
        fi
    fi
fi

# Validate cutting plan
if [ -z "$cut_start" ] && [ -z "$cut_end" ]; then
    echo "Error: No valid cutting points identified"
    exit 1
fi

echo ""
echo "CUTTING PLAN:"
if [ ! -z "$cut_start" ]; then
    echo "  Remove bases 1 to $cut_start (${cut_start}bp from left end)"
fi
if [ ! -z "$cut_end" ]; then
    echo "  Remove bases $cut_end to $ASSEMBLY_LENGTH ($((ASSEMBLY_LENGTH - cut_end))bp from right end)"
fi

# Calculate final genome length
crs_length=$CRS_LENGTH
final_start=$((cut_start - crs_length + 1))
if [ -z "$cut_start" ]; then
    final_start=1
fi

final_end=$((cut_end + crs_length - 1))
if [ -z "$cut_end" ]; then
    final_end=$ASSEMBLY_LENGTH
fi

final_length=$((final_end - final_start + 1))
echo "  Final genome length: ${final_length}bp (was ${ASSEMBLY_LENGTH}bp)"
echo ""

# Create output directory
OUTPUT_DIR="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
mkdir -p "$OUTPUT_DIR"

# Extract the trimmed sequence
echo "Extracting trimmed genome..."
full_sequence=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n')

# Extract trimmed portion
trimmed_sequence=$(echo "$full_sequence" | cut -c${final_start}-${final_end})

# Create output file
output_file="${OUTPUT_DIR}/${SAMPLE}_${ASSEMBLY_TYPE}_trimmed_${crs_length}.fasta"
header=$(head -n1 "$ASSEMBLY_FILE")

echo "$header [trimmed ${final_start}-${final_end}]" > "$output_file"
echo "$trimmed_sequence" >> "$output_file"

# Create summary report
summary_file="${OUTPUT_DIR}/trimming_summary_${crs_length}.txt"
echo "# Genome Trimming Summary" > "$summary_file"
echo "Sample: $SAMPLE" >> "$summary_file"
echo "Assembly type: $ASSEMBLY_TYPE" >> "$summary_file"
echo "Original length: $ASSEMBLY_LENGTH bp" >> "$summary_file"
echo "Trimmed length: $final_length bp" >> "$summary_file"
echo "Bases removed: $((ASSEMBLY_LENGTH - final_length)) bp" >> "$summary_file"
echo "" >> "$summary_file"
echo "CRS positions used:" >> "$summary_file"
for pos in $crs_positions; do
    echo "  $pos" >> "$summary_file"
done
echo "" >> "$summary_file"
echo "Cutting coordinates:" >> "$summary_file"
echo "  Start: $final_start" >> "$summary_file"
echo "  End: $final_end" >> "$summary_file"
if [ ! -z "$cut_start" ]; then
    echo "  Left cut at: $cut_start" >> "$summary_file"
fi
if [ ! -z "$cut_end" ]; then
    echo "  Right cut at: $cut_end" >> "$summary_file"
fi
echo "Generated: $(date)" >> "$summary_file"

echo "=== TRIMMING COMPLETE ==="
echo "Output files:"
echo "  Trimmed genome: $output_file"
echo "  Summary report: $summary_file"
echo ""
echo "Verification:"
echo "  Original: $ASSEMBLY_LENGTH bp"
echo "  Trimmed: $final_length bp" 
echo "  Removed: $((ASSEMBLY_LENGTH - final_length)) bp"