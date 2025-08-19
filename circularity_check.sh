#!/bin/bash

# Circularity Checker Script
ASSEMBLY="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out/polished_assembly.fasta"

echo "=== CHECKING FOR CIRCULAR GENOME ==="

# Extract genome sequence
grep -v '^>' "$ASSEMBLY" | tr -d '\n' > genome_seq.txt
GENOME_LENGTH=$(wc -c < genome_seq.txt)
echo "Genome length: $GENOME_LENGTH bp"

# Check if start of genome matches end of genome (circular overlap)
echo "Checking for circular overlaps..."

for OVERLAP in 1000 2000 5000 10000; do
    if [ $OVERLAP -lt $((GENOME_LENGTH / 2)) ]; then
        echo "Testing $OVERLAP bp overlap..."
        
        # Get start and end sequences
        head -c $OVERLAP genome_seq.txt > start_seq.txt
        tail -c $OVERLAP genome_seq.txt > end_seq.txt
        
        # Check if start matches end (indicates circular genome)
        if cmp -s start_seq.txt end_seq.txt; then
            echo "🔄 CIRCULAR GENOME DETECTED: $OVERLAP bp overlap found!"
            echo "   Start: $(head -c 50 start_seq.txt)..."
            echo "   End:   $(head -c 50 end_seq.txt)..."
        else
            echo "❌ No circular overlap at $OVERLAP bp"
        fi
    fi
done

# Check if genome can be trimmed to remove overlap
echo ""
echo "If circular, the genome should be trimmed to remove redundant sequence"

# Clean up
rm -f genome_seq.txt start_seq.txt end_seq.txt
echo "Circularity check completed!"