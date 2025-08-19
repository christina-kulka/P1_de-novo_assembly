#!/bin/bash

# BLAST Analysis of Extra Sequence
# Check if the additional 66kb is legitimate ORFV sequence

echo "=== ANALYZING EXTRA 66KB SEQUENCE ==="

ASSEMBLY_DIR="/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly/B006/miniasm_out"
ASSEMBLY_FILE="$ASSEMBLY_DIR/polished_assembly.fasta"

cd "$ASSEMBLY_DIR"

echo "Extracting extra sequence (first 66,665 bp)..."

# Extract genome sequence without header
grep -v '^>' "$ASSEMBLY_FILE" | tr -d '\n' > temp_genome.txt

# Extract the extra 66kb sequence
head -c 66665 temp_genome.txt > extra_sequence.txt

# Convert to FASTA format
echo ">extra_66kb_from_your_assembly" > extra_sequence.fasta
fold -w 80 extra_sequence.txt >> extra_sequence.fasta

echo "Extra sequence extracted: $(wc -c < extra_sequence.txt) bp"
echo "First 100 bp: $(head -c 100 extra_sequence.txt)"
echo ""

echo "Running BLAST against NCBI database..."
echo "This may take a few minutes..."

# BLAST the extra sequence
blastn -query extra_sequence.fasta \
       -db nt \
       -remote \
       -outfmt "6 qseqid sseqid pident length evalue stitle" \
       -max_target_seqs 10 \
       -evalue 1e-10 > blast_extra_results.txt

echo ""
echo "=== BLAST RESULTS ==="
echo "Format: Query | Subject | %Identity | Length | E-value | Description"
echo ""

if [ -s blast_extra_results.txt ]; then
    cat blast_extra_results.txt | head -10
    echo ""
    echo "Analysis:"
    
    # Check if results contain ORFV
    if grep -qi "orf\|parapox" blast_extra_results.txt; then
        echo "✅ LEGITIMATE: Extra sequence matches ORFV/Parapoxvirus sequences"
        echo "   Your assembly appears to be more complete than Microsynth's"
    elif grep -qi "virus" blast_extra_results.txt; then
        echo "⚠️  VIRAL: Extra sequence matches other viral sequences"
        echo "   Could be contamination or related virus"
    else
        echo "❌ SUSPICIOUS: Extra sequence doesn't match known viral sequences"
        echo "   Could be assembly artifact or contamination"
    fi
else
    echo "❌ No significant BLAST hits found"
    echo "   Extra sequence may be assembly artifact"
fi

# Clean up temp files
rm -f temp_genome.txt extra_sequence.txt

echo ""
echo "Detailed results saved in: blast_extra_results.txt"
echo "Extra sequence saved in: extra_sequence.fasta"