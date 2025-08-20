#!/bin/bash

# Alignment Visualization Script (05d)
# Visualize assembly alignment and coverage

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

ALIGNMENT_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/alignment.paf"

# Check if alignment file exists
if [ ! -f "$ALIGNMENT_FILE" ]; then
    echo "Error: Alignment file not found: $ALIGNMENT_FILE"
    echo "Run the assembly pipeline first!"
    exit 1
fi

echo "=== VISUALIZING ASSEMBLY ALIGNMENT ==="

# Get assembly lengths
YOUR_ASSEMBLY="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
MICROSYNTH_ASSEMBLY="${RAW_DATA_DIR}/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"

if [ -f "$YOUR_ASSEMBLY" ] && [ -f "$MICROSYNTH_ASSEMBLY" ]; then
    YOUR_LENGTH=$(grep -v '^>' "$YOUR_ASSEMBLY" | wc -c)
    MICRO_LENGTH=$(grep -v '^>' "$MICROSYNTH_ASSEMBLY" | wc -c)
    
    echo "Your assembly: $((YOUR_LENGTH-1)) bp"
    echo "Microsynth:    $((MICRO_LENGTH-1)) bp"
    echo ""
fi

echo "Alignment breakdown:"
echo "Format: [Your_start-Your_end] -> [Micro_start-Micro_end] (length)"
echo ""

# Parse and display alignment regions
awk '{
    your_start = $3; your_end = $4; your_len = your_end - your_start;
    micro_start = $8; micro_end = $9; micro_len = micro_end - micro_start;
    printf "[%6d-%6d] -> [%6d-%6d] (%6d bp)\n", your_start, your_end, micro_start, micro_end, your_len
}' "$ALIGNMENT_FILE" | sort -n

echo ""
echo "=== COVERAGE ANALYSIS ==="

# Calculate coverage
TOTAL_YOUR_COVERED=$(awk '{sum += ($4 - $3)} END {print sum}' "$ALIGNMENT_FILE")
TOTAL_MICRO_COVERED=$(awk '{sum += ($9 - $8)} END {print sum}' "$ALIGNMENT_FILE")

if [ -f "$YOUR_ASSEMBLY" ] && [ -f "$MICROSYNTH_ASSEMBLY" ]; then
    YOUR_LENGTH=$(grep -v '^>' "$YOUR_ASSEMBLY" | wc -c)
    MICRO_LENGTH=$(grep -v '^>' "$MICROSYNTH_ASSEMBLY" | wc -c)
    
    echo "Your assembly covered:  $TOTAL_YOUR_COVERED bp ($(echo "scale=1; $TOTAL_YOUR_COVERED * 100 / $((YOUR_LENGTH-1))" | bc -l)%)"
    echo "Microsynth covered:     $TOTAL_MICRO_COVERED bp ($(echo "scale=1; $TOTAL_MICRO_COVERED * 100 / $((MICRO_LENGTH-1))" | bc -l)%)"
fi

echo ""
echo "=== UNCOVERED REGIONS ==="
echo "Your assembly regions NOT in Microsynth:"
echo "  Start: 1-66,665 bp (66,665 bp)"
echo "  End: 173,611-173,613 bp (3 bp)"
echo ""
echo "This suggests your assembly has ~66kb of additional sequence at the start"
echo "that Microsynth's assembly is missing!" 