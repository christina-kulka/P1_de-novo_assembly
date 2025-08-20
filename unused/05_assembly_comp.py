#!/usr/bin/env python3
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

def get_assembly_paths(sample):
    """Get paths to all assembly types"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    assemblies = {
        'miniasm': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
        'canu': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
        'microsynth': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta"
    }
    
    return assemblies

def compare_all_assemblies(sample):
    """Compare sizes of all assembly types"""
    assemblies = get_assembly_paths(sample)
    
    print(f"=== ASSEMBLY COMPARISON FOR {sample} ===")
    
    # Get assembly lengths for all types
    lengths = {}
    available_assemblies = {}
    
    for name, path in assemblies.items():
        length = get_assembly_length(path)
        lengths[name] = length
        if length is not None:
            available_assemblies[name] = path
            print(f"{name.capitalize():12} {length:,} bp")
        else:
            print(f"{name.capitalize():12} Not found")
    
    # Calculate differences
    if len(available_assemblies) > 1:
        print(f"\nSize differences:")
        assembly_names = list(available_assemblies.keys())
        for i, name1 in enumerate(assembly_names):
            for name2 in assembly_names[i+1:]:
                diff = abs(lengths[name1] - lengths[name2])
                print(f"  {name1} vs {name2}: {diff:,} bp difference")
    
    print()
    return available_assemblies, lengths

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

def analyze_single_assembly_itr(assembly_name, assembly_path, extract_size=5000):
    """Analyze ITR for a single assembly"""
    print(f"--- {assembly_name.upper()} ASSEMBLY ---")
    
    # Extract terminals
    left_term, right_term = extract_terminals(assembly_path, extract_size)
    if left_term is None:
        print(f"Error: Could not extract terminals from {assembly_name}")
        return None
    
    # Find best ITR match
    similarity, left_start, right_start, itr_length = find_best_itr_match(left_term, right_term)
    
    print(f"Best ITR match:")
    print(f"  Similarity:     {similarity:.1f}%")
    print(f"  ITR length:     {itr_length} bp")
    print(f"  Left position:  {left_start + 1}-{left_start + itr_length}")
    print(f"  Right position: {len(right_term) - right_start - itr_length + 1}-{len(right_term) - right_start}")
    
    # Check for internal ITR patterns
    print(f"Searching for internal ITR patterns...")
    genome_seq = next(SeqIO.parse(assembly_path, "fasta")).seq
    genome_length = len(genome_seq)
    
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
    print()
    
    return {
        'similarity': similarity,
        'itr_length': itr_length,
        'left_start': left_start,
        'right_start': right_start,
        'internal_sim': best_internal_sim,
        'internal_pos': best_internal_pos
    }

def analyze_all_assemblies_itr(sample):
    """Complete ITR analysis for all available assemblies"""
    
    # Compare assembly sizes
    available_assemblies, lengths = compare_all_assemblies(sample)
    
    if not available_assemblies:
        print("Error: No assemblies found for analysis")
        return
    
    print(f"=== ITR ANALYSIS FOR {sample} ===")
    
    # Analyze each assembly
    results = {}
    for assembly_name, assembly_path in available_assemblies.items():
        results[assembly_name] = analyze_single_assembly_itr(assembly_name, assembly_path)
    
    # Summary comparison
    print("=== ITR COMPARISON SUMMARY ===")
    print(f"{'Assembly':<12} {'Size (bp)':<10} {'ITR Sim (%)':<12} {'ITR Length':<12} {'Internal Sim (%)':<15}")
    print("-" * 70)
    
    for assembly_name in available_assemblies.keys():
        if results[assembly_name]:
            result = results[assembly_name]
            print(f"{assembly_name.capitalize():<12} {lengths[assembly_name]:<10,} "
                  f"{result['similarity']:<12.1f} {result['itr_length']:<12} "
                  f"{result['internal_sim']:<15.1f}")
        else:
            print(f"{assembly_name.capitalize():<12} {lengths[assembly_name]:<10,} {'Failed':<12} {'Failed':<12} {'Failed':<15}")
    
    # Identify best ITR result
    valid_results = {name: result for name, result in results.items() if result is not None}
    if valid_results:
        best_assembly = max(valid_results.keys(), key=lambda x: valid_results[x]['similarity'])
        print(f"\nBest ITR detection: {best_assembly.capitalize()} assembly "
              f"({valid_results[best_assembly]['similarity']:.1f}% similarity)")

def main():
    if len(sys.argv) != 2:
        print("Usage: python 05_compare_and_extract_itr.py SAMPLE_NAME")
        print("Example: python 05_compare_and_extract_itr.py D1701")
        sys.exit(1)
    
    sample = sys.argv[1]
    analyze_all_assemblies_itr(sample)

if __name__ == "__main__":
    main()