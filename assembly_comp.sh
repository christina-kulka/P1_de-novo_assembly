#!/bin/bash

# Assembly Evaluation Script
# Compare your assembly vs Microsynth assembly

YOUR_ASSEMBLY="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out/polished_assembly.fasta"
MICROSYNTH_ASSEMBLY="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/00_raw_data_microsynth/B006_results/B006_results/Assembly/B006.fasta"

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