#!/bin/bash

ASSEMBLY_FILE="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out/polished_assembly.fasta"

cd /home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out

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