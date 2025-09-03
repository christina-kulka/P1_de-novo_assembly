#!/bin/bash

# Improved Retrospective ORF Annotation Script
# Usage: ./improved_retrospective_annotate.sh SAMPLE_NAME

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
    echo "  2. Find ORF locations using direct sequence matching (no BLAST coordinate issues)"
    echo "  3. Add ORF annotations directly to the genome with proper translation"
    echo "  4. Save result as proposed_SAMPLE_annotated.gbk"
    exit 1
fi

SAMPLE_NAME=$1

# Get script directory and source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Define input and output files based on sample name
# Check for both .dna and .fasta/.fa files
INPUT_FILE_DNA="${PROPOSAL_INPUT_DIR}/proposed_${SAMPLE_NAME}_final_map.dna"
INPUT_FILE_FASTA="${PROPOSAL_INPUT_DIR}/proposed_${SAMPLE_NAME}_final_map.fasta"
INPUT_FILE_FA="${PROPOSAL_INPUT_DIR}/proposed_${SAMPLE_NAME}_final_map.fa"

# Find which input file exists
if [ -f "$INPUT_FILE_DNA" ]; then
    INPUT_FILE="$INPUT_FILE_DNA"
    echo "Found .dna file: $INPUT_FILE"
elif [ -f "$INPUT_FILE_FASTA" ]; then
    INPUT_FILE="$INPUT_FILE_FASTA"
    echo "Found .fasta file: $INPUT_FILE"
elif [ -f "$INPUT_FILE_FA" ]; then
    INPUT_FILE="$INPUT_FILE_FA"
    echo "Found .fa file: $INPUT_FILE"
else
    echo "Error: No input file found. Looked for:"
    echo "  $INPUT_FILE_DNA"
    echo "  $INPUT_FILE_FASTA"
    echo "  $INPUT_FILE_FA"
    echo "Expected file: proposed_${SAMPLE_NAME}_final_map.[dna|fasta|fa] in ${PROPOSAL_INPUT_DIR}/"
    exit 1
fi

OUTPUT_FILE="${PROPOSAL_OUTPUT_DIR}/proposed_${SAMPLE_NAME}_annotated_improved.gbk"

# Check required ORF files
ORF_DNA_FILE="${ORF_DATABASE_DIR}/${SAMPLE_NAME}_ORFs.fasta"

if [ ! -f "$ORF_DNA_FILE" ]; then
    echo "Error: ORF DNA file not found: $ORF_DNA_FILE"
    echo "Expected: ${SAMPLE_NAME}_ORFs.fasta in ${ORF_DATABASE_DIR}/"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$PROPOSAL_OUTPUT_DIR"

echo "=== IMPROVED RETROSPECTIVE ORF ANNOTATION ==="
echo "Sample: $SAMPLE_NAME"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "ORF DNA: $ORF_DNA_FILE"
echo ""

# Run the improved ORF annotation (no BLAST, direct sequence matching)
echo "Running improved ORF annotation (direct sequence matching)..."
$PYTHON3 "${SCRIPT_DIR}/improved_orf_annotator.py" \
    "$INPUT_FILE" \
    "$ORF_DNA_FILE" \
    "$OUTPUT_FILE" \
    --sample-name "$SAMPLE_NAME" \
    --min-identity 95.0

if [ $? -ne 0 ]; then
    echo "Error: ORF annotation failed"
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
    echo "  Method: Direct sequence matching (no BLAST)"
    
    # Check for any translation issues
    echo ""
    echo "Quality check - looking for potential issues:"
    
    # Check for very short translations (might indicate frame issues)
    SHORT_TRANSLATIONS=$(grep -A1 "/translation=" "$OUTPUT_FILE" | grep -v "^--$" | grep -v "/translation=" | grep -E "^[A-Z*]{1,10}$" | wc -l)
    if [ "$SHORT_TRANSLATIONS" -gt 0 ]; then
        echo "  Warning: Found $SHORT_TRANSLATIONS very short protein translations"
    else
        echo "  ✓ All protein translations appear reasonable length"
    fi
    
    # Check for excessive stop codons in translations
    STOP_CODONS=$(grep -A1 "/translation=" "$OUTPUT_FILE" | grep -v "^--$" | grep -v "/translation=" | grep -o "\*" | wc -l)
    TOTAL_TRANSLATIONS=$(grep -c "/translation=" "$OUTPUT_FILE")
    if [ "$TOTAL_TRANSLATIONS" -gt 0 ]; then
        AVG_STOPS=$(echo "scale=2; $STOP_CODONS / $TOTAL_TRANSLATIONS" | bc 2>/dev/null || echo "0")
        if (( $(echo "$AVG_STOPS > 2" | bc -l 2>/dev/null) )); then
            echo "  Warning: Average $AVG_STOPS stop codons per translation (may indicate frame issues)"
        else
            echo "  ✓ Stop codon frequency appears normal"
        fi
    fi
    
else
    echo "Error: Output file was not created"
    exit 1
fi

echo ""
echo "Improved retrospective ORF annotation completed successfully!"
echo "Use the annotated file: $OUTPUT_FILE"
