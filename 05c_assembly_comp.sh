#!/bin/bash

# Assembly Comparison Script (05c)
# Compare your assembly vs Microsynth assembly

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

YOUR_ASSEMBLY="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
MICROSYNTH_ASSEMBLY="${RAW_DATA_DIR}/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"

# Check if assembly files exist
if [ ! -f "$YOUR_ASSEMBLY" ]; then
    echo "Error: Your assembly file not found: $YOUR_ASSEMBLY"
    echo "Run the assembly pipeline first!"
    exit 1
fi

if [ ! -f "$MICROSYNTH_ASSEMBLY" ]; then
    echo "Error: Microsynth assembly file not found: $MICROSYNTH_ASSEMBLY"
    exit 1
fi

echo "=== ASSEMBLY COMPARISON ==="

# Basic stats
echo "Your assembly:"
echo "  Length: $(grep -v '^>' $YOUR_ASSEMBLY | wc -c) bp"
echo "  Contigs: $(grep -c '^>' $YOUR_ASSEMBLY)"

echo "Microsynth assembly:"
echo "  Length: $(grep -v '^>' $MICROSYNTH_ASSEMBLY | wc -c) bp" 
echo "  Contigs: $(grep -c '^>' $MICROSYNTH_ASSEMBLY)"

# Show contig headers
echo "Your contigs:"
grep '^>' $YOUR_ASSEMBLY

echo "Microsynth contigs:"
grep '^>' $MICROSYNTH_ASSEMBLY

# Compare ITRs between assemblies
echo "=== ITR COMPARISON ==="

# Check your assembly ITRs
echo "Your assembly ITRs:"
grep -v '^>' $YOUR_ASSEMBLY | tr -d '\n' > your_genome.txt
head -c 1000 your_genome.txt > your_start.txt
tail -c 1000 your_genome.txt > your_end.txt

# Check Microsynth ITRs  
echo "Microsynth assembly ITRs:"
grep -v '^>' $MICROSYNTH_ASSEMBLY | tr -d '\n' > micro_genome.txt
head -c 1000 micro_genome.txt > micro_start.txt
tail -c 1000 micro_genome.txt > micro_end.txt

# Visual comparison
echo "Your start:    $(head -c 50 your_start.txt)"
echo "Your end:      $(tail -c 50 your_end.txt | rev)"
echo "Micro start:   $(head -c 50 micro_start.txt)"  
echo "Micro end:     $(tail -c 50 micro_end.txt | rev)"

# Clean up
rm -f your_genome.txt your_start.txt your_end.txt micro_genome.txt micro_start.txt micro_end.txt 