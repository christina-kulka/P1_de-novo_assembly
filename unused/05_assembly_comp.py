"""
Assembly comparison and ITR extraction script
Usage: python 05_compare_and_extract_itr.py SAMPLE_NAME
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def get_assembly_length(fasta_file):
    """Get total length of assembly"""
    if not os.path.exists(fasta_file):
        return None
    
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
    return total_length

def extract_terminals(fasta_file, extract_size):
    """Extract terminal regions from assembly"""
    if not os.path.exists(fasta_file):
        return None, None
    
    # Get the first (longest) contig
    record = next(SeqIO.parse(fasta_file, "fasta"))
    genome = record.seq
    
    left_terminal = genome[:extract_size]
    right_terminal = genome[-extract_size:]
    
    return left_terminal, right_terminal

def compare_assemblies_and_determine_size(sample):
    """Compare assembly sizes and determine ITR extraction size"""
    
    # Define paths
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    # Your viral assembly (the cleaned one)
    #your_assembly = f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta"
    your_assembly = f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta"


    # Microsynth assembly
    microsynth_assembly = f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta"
    
    print(f"=== ASSEMBLY COMPARISON FOR {sample} ===")
    
    # Get assembly lengths
    your_length = get_assembly_length(your_assembly)
    micro_length = get_assembly_length(microsynth_assembly)
    
    if your_length is None:
        print(f"Error: Your assembly not found: {your_assembly}")
        return None
    
    if micro_length is None:
        print(f"Error: Microsynth assembly not found: {microsynth_assembly}")
        return None
    
    print(f"Your assembly:    {your_length:,} bp")
    print(f"Microsynth:       {micro_length:,} bp")
    
    # Calculate difference
    length_diff = abs(your_length - micro_length)
    print(f"Difference:       {length_diff:,} bp")
    
    # Determine extraction size
    if your_length > micro_length and length_diff > 20000:  # Your assembly much longer
        #extract_size = min(length_diff // 2 + 5000, 25000)  # Half the difference + 5kb, max 25kb
        extract_size = 5000  # Use standard 5kb for Microsynth analysis

        reason = f"Your assembly is {length_diff:,} bp longer - using {extract_size:,} bp"
    elif micro_length > your_length and length_diff > 20000:  # Microsynth much longer
        #extract_size = min(length_diff // 2 + 5000, 25000)
        extract_size = 5000  # Use standard 5kb for Microsynth analysis

        reason = f"Microsynth is {length_diff:,} bp longer - using {extract_size:,} bp"
    else:  # Similar sizes
        extract_size = 5000
        reason = "Similar assembly sizes - using standard 5,000 bp"
    
    print(f"ITR extraction:   {extract_size:,} bp from each end")
    print(f"Reason:           {reason}")
    print()
    
    return your_assembly, microsynth_assembly, extract_size

def simple_similarity(seq1, seq2):
    """Calculate simple similarity percentage between two sequences"""
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100 if len(seq1) > 0 else 0

def find_best_itr_match(left_term, right_term, window_size=100, step_size=50):
    """Find best ITR match using sliding windows"""
    right_rc = right_term.reverse_complement()
    
    best_similarity = 0
    best_left_start = 0
    best_right_start = 0
    best_length = 0
    
    # Try different window sizes
    for win_size in [50, 100, 150, 200, 300]:
        if win_size > len(left_term) or win_size > len(right_rc):
            continue
            
        # Slide window across left terminal
        for left_start in range(0, len(left_term) - win_size + 1, step_size):
            left_seq = left_term[left_start:left_start + win_size]
            
            # Slide window across right terminal (reverse complement)
            for right_start in range(0, len(right_rc) - win_size + 1, step_size):
                right_seq = right_rc[right_start:right_start + win_size]
                
                similarity = simple_similarity(left_seq, right_seq)
                
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_left_start = left_start
                    best_right_start = right_start
                    best_length = win_size
    
    return best_similarity, best_left_start, best_right_start, best_length

def analyze_itr_for_sample(sample):
    """Complete ITR analysis for one sample"""
    
    # Compare assemblies and determine extraction size
    result = compare_assemblies_and_determine_size(sample)
    if result is None:
        return
    
    your_assembly, microsynth_assembly, extract_size = result
    
    print(f"=== ITR ANALYSIS FOR {sample} ===")
    
    # Extract terminals from your assembly
    your_left, your_right = extract_terminals(your_assembly, extract_size)
    if your_left is None:
        print("Error: Could not extract terminals from your assembly")
        return
    
    # Find ITR in your assembly
    print("Analyzing your assembly...")

    # Try searching deeper into the terminal regions
    print("Searching for ITRs in full terminal regions...")
    similarity, left_start, right_start, itr_length = find_best_itr_match(your_left, your_right)
    print(f"Full region search: {similarity:.1f}%")

    # Also try just the outermost 5kb (like a normal genome)
    your_left_5k = your_left[:5000]
    your_right_5k = your_right[-5000:]
    sim_5k, left_5k, right_5k, len_5k = find_best_itr_match(your_left_5k, your_right_5k)
    
    print(f"Outermost 5kb only: {sim_5k:.1f}%")    
    print(f"Best ITR match:")
    print(f"  Similarity:     {similarity:.1f}%")
    print(f"  ITR length:     ~{itr_length:,} bp")
    print(f"  Left position:  {left_start + 1}-{left_start + itr_length}")
    print(f"  Right position: {len(your_right) - right_start - itr_length + 1}-{len(your_right) - right_start}")
    
    # Check for internal ITR patterns (might indicate true genome boundaries)
    print("\nSearching for internal ITR patterns...")
    genome_seq = next(SeqIO.parse(your_assembly, "fasta")).seq
    genome_length = len(genome_seq)

    # Test potential genome boundaries every 5kb
    best_internal_sim = 0
    best_internal_pos = 0

    for pos in range(5000, genome_length - 5000, 5000):
        left_test = genome_seq[pos:pos + 2000]
        right_test = genome_seq[-(pos + 2000):-pos].reverse_complement()
        
        if len(left_test) == len(right_test) == 2000:
            sim = simple_similarity(left_test, right_test)
            if sim > best_internal_sim:
                best_internal_sim = sim
                best_internal_pos = pos

    print(f"Best internal ITR-like pattern: {best_internal_sim:.1f}% at position {best_internal_pos}")

    # Quick comparison with Microsynth if available
    micro_left, micro_right = extract_terminals(microsynth_assembly, 5000)
    if micro_left is not None:
        print("\nComparing with Microsynth assembly (5kb)...")
        micro_sim, _, _, micro_len = find_best_itr_match(micro_left, micro_right)
        print(f"  Microsynth ITR similarity: {micro_sim:.1f}%")
    
    print()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 05_compare_and_extract_itr.py SAMPLE_NAME")
        print("Example: python 05_compare_and_extract_itr.py B006")
        sys.exit(1)
    
    sample = sys.argv[1]
    analyze_itr_for_sample(sample)

if __name__ == "__main__":
    main()