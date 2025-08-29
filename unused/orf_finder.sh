#!/bin/bash

# ORF Position Finder Script
# Usage: ./find_orf_positions.sh SAMPLE_NAME ORF_NUMBER DIRECTION [ASSEMBLY_TYPE]
# Example: ./find_orf_positions.sh B006 134 leftmost
# Example: ./find_orf_positions.sh B006 134 rightmost canu_super

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Parse command line arguments
SAMPLE=$1
ORF_NUMBER=$2
DIRECTION=$3
ASSEMBLY_TYPE=${4:-"canu_super"}  # Default to canu_super

if [ -z "$SAMPLE" ] || [ -z "$ORF_NUMBER" ] || [ -z "$DIRECTION" ]; then
    echo "Usage: $0 SAMPLE_NAME ORF_NUMBER DIRECTION [ASSEMBLY_TYPE]"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME   - Sample identifier (e.g., B006)"
    echo "  ORF_NUMBER    - ORF number to search for (e.g., 134)"
    echo "  DIRECTION     - leftmost or rightmost"
    echo "  ASSEMBLY_TYPE - Assembly type (default: canu_super)"
    echo ""
    echo "Examples:"
    echo "  $0 B006 134 leftmost"
    echo "  $0 B006 134 rightmost canu"
    echo "  $0 B007 130 leftmost miniasm_viral"
    exit 1
fi

# Validate direction parameter
if [ "$DIRECTION" != "leftmost" ] && [ "$DIRECTION" != "rightmost" ]; then
    echo "Error: DIRECTION must be 'leftmost' or 'rightmost'"
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

# Check if required files exist
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

if [ ! -f "$ORF_REFERENCE_FILE" ]; then
    echo "Error: ORF reference file not found: $ORF_REFERENCE_FILE"
    exit 1
fi

# Create output directory
OUTPUT_DIR="${ITR_ANALYSIS_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}/orf_finder"
mkdir -p "$OUTPUT_DIR"

echo "=== ORF POSITION FINDER ==="
echo "Sample: $SAMPLE"
echo "ORF: $ORF_NUMBER"
echo "Direction: $DIRECTION"
echo "Assembly type: $ASSEMBLY_TYPE"
echo "Assembly file: $ASSEMBLY_FILE"
echo ""

cd "$OUTPUT_DIR"

# Get assembly length for context
ASSEMBLY_LENGTH=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)

# Create BLAST database from assembly (if not exists)
if [ ! -f "assembly_db.nin" ]; then
    echo "Creating BLAST database..."
    makeblastdb -in "$ASSEMBLY_FILE" -dbtype nucl -out assembly_db -logfile makeblastdb.log
fi

# BLAST ORFs against assembly
echo "Searching for ORF $ORF_NUMBER..."
blastn -query "$ORF_REFERENCE_FILE" \
       -db assembly_db \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
       -out temp_orf_hits.txt

# Filter for the specific ORF with high quality hits
grep -i "$ORF_NUMBER" temp_orf_hits.txt | awk '$3 >= 95 && $13 >= 90' > orf_${ORF_NUMBER}_hits.txt

# Check if ORF was found
if [ ! -s "orf_${ORF_NUMBER}_hits.txt" ]; then
    echo "Error: ORF $ORF_NUMBER not found with high confidence (>95% identity, >90% coverage)"
    echo "Checking for lower confidence matches..."
    grep -i "$ORF_NUMBER" temp_orf_hits.txt | awk '$3 >= 80 && $13 >= 70' > orf_${ORF_NUMBER}_low_conf.txt
    
    if [ -s "orf_${ORF_NUMBER}_low_conf.txt" ]; then
        echo "Found lower confidence matches:"
        cat orf_${ORF_NUMBER}_low_conf.txt
    else
        echo "No matches found for ORF $ORF_NUMBER"
    fi
    exit 1
fi

# Process the hits to find leftmost or rightmost occurrence
echo "Processing ORF $ORF_NUMBER hits..."

OUTPUT_FILE="orf_${ORF_NUMBER}_${DIRECTION}_position.txt"

# Create output file with header
echo "# ORF Position Finder Results" > "$OUTPUT_FILE"
echo "# Sample: $SAMPLE" >> "$OUTPUT_FILE"
echo "# ORF: $ORF_NUMBER" >> "$OUTPUT_FILE"
echo "# Direction: $DIRECTION" >> "$OUTPUT_FILE"
echo "# Assembly: $ASSEMBLY_TYPE" >> "$OUTPUT_FILE"
echo "# Assembly length: $ASSEMBLY_LENGTH bp" >> "$OUTPUT_FILE"
echo "# Generated: $(date)" >> "$OUTPUT_FILE"
echo "#" >> "$OUTPUT_FILE"
echo "# Format: ORF_NAME ORIENTATION START END IDENTITY COVERAGE DISTANCE_FROM_START DISTANCE_FROM_END" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Find leftmost or rightmost occurrence
if [ "$DIRECTION" = "leftmost" ]; then
    # Sort by start position (leftmost first)
    selected_hit=$(awk '{
        start = ($9 < $10) ? $9 : $10
        end = ($9 < $10) ? $10 : $9
        orientation = ($9 < $10) ? "forward" : "reverse"
        print $1 "\t" orientation "\t" start "\t" end "\t" $3 "\t" $13 "\t" start
    }' orf_${ORF_NUMBER}_hits.txt | sort -k7,7n | head -n1)
    
    echo "Found leftmost ORF $ORF_NUMBER:"
else
    # Sort by start position (rightmost first) 
    selected_hit=$(awk '{
        start = ($9 < $10) ? $9 : $10
        end = ($9 < $10) ? $10 : $9
        orientation = ($9 < $10) ? "forward" : "reverse"
        print $1 "\t" orientation "\t" start "\t" end "\t" $3 "\t" $13 "\t" start
    }' orf_${ORF_NUMBER}_hits.txt | sort -k7,7nr | head -n1)
    
    echo "Found rightmost ORF $ORF_NUMBER:"
fi

if [ ! -z "$selected_hit" ]; then
    # Parse the selected hit
    orf_name=$(echo "$selected_hit" | cut -f1)
    orientation=$(echo "$selected_hit" | cut -f2)
    start_pos=$(echo "$selected_hit" | cut -f3)
    end_pos=$(echo "$selected_hit" | cut -f4)
    identity=$(echo "$selected_hit" | cut -f5)
    coverage=$(echo "$selected_hit" | cut -f6)
    
    # Calculate distances from assembly ends
    distance_from_start=$start_pos
    distance_from_end=$((ASSEMBLY_LENGTH - end_pos))
    
    # Write to output file
    echo "$orf_name $orientation $start_pos $end_pos $identity $coverage $distance_from_start $distance_from_end" >> "$OUTPUT_FILE"
    
    # Display results
    echo "  ORF: $orf_name"
    echo "  Orientation: $orientation"
    echo "  Position: $start_pos-$end_pos"
    echo "  Identity: ${identity}%"
    echo "  Coverage: ${coverage}%"
    echo "  Distance from assembly start: ${distance_from_start}bp"
    echo "  Distance from assembly end: ${distance_from_end}bp"
    echo ""
    echo "Results saved to: $OUTPUT_FILE"
    
    # Also create a simple position file for easy parsing by other scripts
    echo "$start_pos $end_pos" > "orf_${ORF_NUMBER}_${DIRECTION}_coords.txt"
    
    # Show all occurrences for context
    echo ""
    echo "All ORF $ORF_NUMBER occurrences found:"
    echo "Position Range    Orientation  Identity  Coverage"
    echo "----------------  -----------  --------  --------"
    
    awk '{
        start = ($9 < $10) ? $9 : $10
        end = ($9 < $10) ? $10 : $9
        orientation = ($9 < $10) ? "forward" : "reverse"
        printf "%-8d-%-7d  %-11s  %6.1f%%    %6.1f%%\n", start, end, orientation, $3, $13
    }' orf_${ORF_NUMBER}_hits.txt | sort -k1,1n
    
else
    echo "Error: Could not process ORF hits"
    exit 1
fi

# Clean up temporary files
rm -f temp_orf_hits.txt

echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "Output files:"
echo "  - $OUTPUT_FILE (detailed results)"
echo "  - orf_${ORF_NUMBER}_${DIRECTION}_coords.txt (coordinates only)"
echo "  - orf_${ORF_NUMBER}_hits.txt (all high-quality hits)"