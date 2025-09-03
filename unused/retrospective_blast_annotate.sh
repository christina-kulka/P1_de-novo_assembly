#!/bin/bash

# Retrospective ORF Coordinate-Based Annotation Script
# Usage: ./retrospective_coordinate_annotate.sh SAMPLE_NAME

# Check if sample name is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME - Sample identifier (e.g., B021)"
    echo ""
    echo "Example:"
    echo "  $0 B021"
    echo ""
    echo "This script will:"
    echo "  1. Look for proposed_SAMPLE_final_map.dna in the proposal directory"
    echo "  2. Find ORF coordinates using BLAST against ORF DNA sequences"
    echo "  3. Add ORF annotations directly to the genome"
    echo "  4. Save result as proposed_SAMPLE_annotated.gbk"
    exit 1
fi

SAMPLE_NAME=$1

# Get script directory and source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Define input and output files based on sample name
INPUT_FILE="${PROPOSAL_INPUT_DIR}/proposed_${SAMPLE_NAME}_final_map.dna"
OUTPUT_FILE="${PROPOSAL_OUTPUT_DIR}/proposed_${SAMPLE_NAME}_annotated_blast.gbk"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    echo "Expected file: proposed_${SAMPLE_NAME}_final_map.dna in ${PROPOSAL_INPUT_DIR}/"
    exit 1
fi

# Check required ORF files
ORF_DNA_FILE="${ORF_DATABASE_DIR}/${SAMPLE_NAME}_ORFs.fasta"
ORF_PROTEIN_FILE="${PROTEIN_DATABASE_DIR}/${SAMPLE_NAME}_protein.fasta"

if [ ! -f "$ORF_DNA_FILE" ]; then
    echo "Error: ORF DNA file not found: $ORF_DNA_FILE"
    echo "Expected: ${SAMPLE_NAME}_ORFs.fasta in ${ORF_DATABASE_DIR}/"
    exit 1
fi

if [ ! -f "$ORF_PROTEIN_FILE" ]; then
    echo "Error: ORF protein file not found: $ORF_PROTEIN_FILE"
    echo "Expected: ${SAMPLE_NAME}_protein.fasta in ${PROTEIN_DATABASE_DIR}/"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$PROPOSAL_OUTPUT_DIR"

echo "=== RETROSPECTIVE ORF COORDINATE-BASED ANNOTATION ==="
echo "Sample: $SAMPLE_NAME"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "ORF DNA: $ORF_DNA_FILE"
echo "ORF Proteins: $ORF_PROTEIN_FILE"
echo ""

# Run the coordinate-based ORF annotation
echo "Running coordinate-based ORF annotation..."
"$SCRIPT_DIR/orf_coordinate_annotate.sh" "$INPUT_FILE" "$SAMPLE_NAME" "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "Error: ORF coordinate annotation failed"
    exit 1
fi

echo ""
echo "=== ANNOTATION COMPLETED ==="
echo "Output file: $OUTPUT_FILE"

# Show final summary
if [ -f "$OUTPUT_FILE" ]; then
    TOTAL_FEATURES=$(grep -c "/gene=" "$OUTPUT_FILE" 2>/dev/null || echo "0")
    echo "Total ORF features annotated: $TOTAL_FEATURES"
    
    echo ""
    echo "Sample annotated ORFs:"
    grep "/gene=" "$OUTPUT_FILE" | head -10 | sed 's/.*\/gene="/  - /' | sed 's/".*//' | sort
    
    echo ""
    echo "Annotation statistics:"
    echo "  Input file: $(basename "$INPUT_FILE")"
    echo "  Output file: $(basename "$OUTPUT_FILE")"
    echo "  ORF features: $TOTAL_FEATURES"
    echo "  Method: Coordinate-based BLAST mapping"
else
    echo "Error: Output file was not created"
    exit 1
fi

echo ""
echo "Retrospective ORF annotation completed successfully!"
echo "Use the annotated file: $OUTPUT_FILE"