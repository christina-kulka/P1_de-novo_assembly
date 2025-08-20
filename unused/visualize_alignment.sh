#!/bin/bash

# Alignment Visualization Script
echo "=== VISUALIZING ASSEMBLY ALIGNMENT ==="

ALIGNMENT_FILE="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out/alignment.paf"

echo "Your assembly: 173,613 bp"
echo "Microsynth:    138,382 bp" 
echo ""
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

echo "Your assembly covered:  $TOTAL_YOUR_COVERED bp ($(echo "scale=1; $TOTAL_YOUR_COVERED * 100 / 173613" | bc -l)%)"
echo "Microsynth covered:     $TOTAL_MICRO_COVERED bp ($(echo "scale=1; $TOTAL_MICRO_COVERED * 100 / 138382" | bc -l)%)"

echo ""
echo "=== UNCOVERED REGIONS ==="
echo "Your assembly regions NOT in Microsynth:"
echo "  Start: 1-66,665 bp (66,665 bp)"
echo "  End: 173,611-173,613 bp (3 bp)"
echo ""
echo "This suggests your assembly has ~66kb of additional sequence at the start"
echo "that Microsynth's assembly is missing!"



# Extract the extra sequence (first 66,665 bp)
head -c 66665 /home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out/genome_seq.txt > extra_sequence.txt

# Convert to FASTA format
echo ">extra_66kb_sequence" > extra_sequence.fasta
fold -w 80 extra_sequence.txt >> extra_sequence.fasta

# BLAST it
blastn -query extra_sequence.fasta -db nt -remote -outfmt "6 qseqid sseqid pident length evalue stitle" -max_target_seqs 5