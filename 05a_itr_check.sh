#!/bin/bash

# ITR Check Script (05a)
# Check for mirror-image sequences in ORFV genome

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

ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"

# Check if assembly file exists
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    echo "Run the assembly pipeline first!"
    exit 1
fi

cd "${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out"

echo "Checking ITRs (mirror-image sequences) in ORFV genome..."

# Extract DNA sequence
grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' > genome_sequence.txt

GENOME_LENGTH=$(wc -c < genome_sequence.txt)
echo "Genome length: $GENOME_LENGTH bp"

# Check ITRs
for LENGTH in 100 500 1000 2000 5000; do
    echo "Checking ITRs of $LENGTH bp..."
    
    # Extract start and end
    START_SEQ=$(head -c $LENGTH genome_sequence.txt)
    END_SEQ=$(tail -c $LENGTH genome_sequence.txt)
    
    # Reverse the end sequence (mirror it)
    END_REVERSED=$(echo "$END_SEQ" | rev)
    
    # Compare start with reversed end
    if [ "$START_SEQ" = "$END_REVERSED" ]; then
        echo "✅ PERFECT ITR: $LENGTH bp mirror-image repeats found!"
    else
        echo "❌ NO ITR: $LENGTH bp sequences are not mirror images"
        
        # Show first 30 characters for comparison
        echo "   Start:     $(echo "$START_SEQ" | head -c 30)..."
        echo "   End(rev):  $(echo "$END_REVERSED" | head -c 30)..."
    fi
done

# Clean up
rm -f genome_sequence.txt
echo "ITR analysis completed!" 