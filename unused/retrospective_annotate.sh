#!/bin/bash

# Retrospective Annotation Script
# Usage: ./retrospective_annotate.sh SAMPLE_NAME

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
    echo "  2. Run ORF annotation pipeline"
    echo "  3. Save result as proposed_SAMPLE_annotated.gbk"
    exit 1
fi

SAMPLE_NAME=$1

# Get script directory and source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Define input and output files based on sample name
INPUT_FILE="${PROPOSAL_INPUT_DIR}/proposed_${SAMPLE_NAME}_final_map.dna"
OUTPUT_FILE="${PROPOSAL_OUTPUT_DIR}/proposed_${SAMPLE_NAME}_annotated.gbk"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    echo "Expected file: proposed_${SAMPLE_NAME}_final_map.dna in ${PROPOSAL_INPUT_DIR}/"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$PROPOSAL_OUTPUT_DIR"

echo "=== RETROSPECTIVE ANNOTATION PIPELINE ==="
echo "Sample: $SAMPLE_NAME"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Step 1: Run Prokka to get ORF annotations (saves to temp directory)
echo "Step 1: Running Prokka for ORF annotation..."
"$SCRIPT_DIR/prokka_annotation_proposal.sh" "$INPUT_FILE" "$SAMPLE_NAME" "${TEMP_DIR}/${SAMPLE_NAME}_orfs.gbk"

if [ $? -ne 0 ]; then
    echo "Error: Prokka annotation failed"
    exit 1
fi

# Step 2: Restore ORF names in the Prokka output
echo "Step 2: Restoring ORF names..."
$PYTHON3 "$SCRIPT_DIR/4b_restore_orf_names.py" "$SAMPLE_NAME" "${PROTEIN_DATABASE_DIR}/${SAMPLE_NAME}_protein.fasta" "$TEMP_DIR"

if [ $? -ne 0 ]; then
    echo "Error: ORF name restoration failed"
    exit 1
fi

# Step 3: Merge with original file - ONLY true ORFs, no generic hypothetical proteins
echo "Step 3: Merging annotations..."
$PYTHON3 "$SCRIPT_DIR/merge_annotations.py" "$INPUT_FILE" "${TEMP_DIR}/${SAMPLE_NAME}.orf_named.gbk" "$OUTPUT_FILE" --allow-overlaps

if [ $? -ne 0 ]; then
    echo "Error: Annotation merging failed"
    exit 1
fi

echo ""
echo "=== RETROSPECTIVE ANNOTATION COMPLETED ==="
echo "Output file: $OUTPUT_FILE"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -f "${TEMP_DIR}/${SAMPLE_NAME}"*

echo "Annotation pipeline completed successfully!"
