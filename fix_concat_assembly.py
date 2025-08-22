#!/usr/bin/env python3
"""
Simple concatenation detection and correction script
Detects direct repeats and removes one ITR to create proper inverted structure
Usage: python fix_concatenation.py SAMPLE_NAME ASSEMBLY_TYPE
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO

def get_assembly_path(sample, assembly_type):
    """Get path for specific assembly type"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    paths = {
        'canu': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
        'canu_ultra': f"{base_dir}/03_assembly/{sample}/canu_ultra_output/{sample}.contigs.fasta",
        'miniasm': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta"
    }
    
    return paths.get(assembly_type)

def get_longest_contig(fasta_file):
    """Get the longest contig from assembly"""
    longest_record = None
    max_length = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > max_length:
            max_length = len(record.seq)
            longest_record = record
    
    return longest_record

def detect_terminal_similarity(sequence, test_length=50000):
    """Check if terminals are similar (direct repeat) vs inverted"""
    
    seq_len = len(sequence)
    if seq_len < test_length * 3:
        print(f"Sequence too short ({seq_len:,} bp) for analysis")
        return None
    
    left_term = sequence[:test_length]
    right_term = sequence[-test_length:]
    right_term_rc = right_term.reverse_complement()
    
    # Calculate similarities
    def simple_similarity(seq1, seq2):
        min_len = min(len(seq1), len(seq2))
        seq1, seq2 = seq1[:min_len], seq2[:min_len]
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return (matches / min_len) * 100 if min_len > 0 else 0
    
    direct_similarity = simple_similarity(left_term, right_term)
    inverted_similarity = simple_similarity(left_term, right_term_rc)
    
    print(f"Terminal analysis ({test_length:,} bp):")
    print(f"  Direct similarity: {direct_similarity:.1f}%")
    print(f"  Inverted similarity: {inverted_similarity:.1f}%")
    
    return {
        'direct_similarity': direct_similarity,
        'inverted_similarity': inverted_similarity,
        'is_concatenated': direct_similarity > 70 and direct_similarity > inverted_similarity
    }

def find_best_cut_point(sequence):
    """Find the best position to cut a concatenated assembly"""
    
    seq_len = len(sequence)
    expected_size = 140000
    
    # If not significantly longer than expected, don't cut
    if seq_len < expected_size * 1.2:
        return None
    
    # Test cut points around expected genome size
    best_cut = None
    best_score = 0
    
    search_start = max(100000, seq_len // 3)
    search_end = min(180000, seq_len * 2 // 3)
    
    for cut_point in range(search_start, search_end, 2000):
        truncated_seq = sequence[:cut_point]
        
        # Check if this creates better ITRs
        analysis = detect_terminal_similarity(truncated_seq, 10000)
        
        if analysis:
            # Score based on inverted similarity and size
            size_score = 100 - abs(cut_point - expected_size) / 1000
            itr_score = analysis['inverted_similarity']
            
            # Penalize if still showing direct similarity
            if analysis['direct_similarity'] > 70:
                itr_score *= 0.3
            
            total_score = (size_score + itr_score) / 2
            
            if total_score > best_score:
                best_score = total_score
                best_cut = cut_point
    
    return best_cut if best_score > 40 else None

def fix_assembly(sample, assembly_type):
    """Main function to fix concatenated assembly"""
    
    print(f"=== FIXING {sample} {assembly_type.upper()} ASSEMBLY ===")
    
    # Get assembly path
    assembly_path = get_assembly_path(sample, assembly_type)
    if not assembly_path or not os.path.exists(assembly_path):
        print(f"Assembly not found: {assembly_path}")
        return
    
    # Load sequence
    record = get_longest_contig(assembly_path)
    if not record:
        print("Could not load assembly sequence")
        return
    
    original_seq = record.seq
    seq_len = len(original_seq)
    
    print(f"Original assembly: {seq_len:,} bp")
    
    # Analyze terminal structure
    analysis = detect_terminal_similarity(original_seq)
    
    if not analysis:
        return
    
    if not analysis['is_concatenated']:
        print("Assembly appears to have proper structure - no correction needed")
        return
    
    print("Direct repeat detected - assembly appears concatenated")
    
    # Find optimal cut point
    cut_point = find_best_cut_point(original_seq)
    
    if not cut_point:
        print("Could not find suitable cut point")
        return
    
    # Create corrected sequence
    corrected_seq = original_seq[:cut_point]
    
    print(f"Cutting at position {cut_point:,}")
    print(f"Corrected length: {len(corrected_seq):,} bp")
    
    # Validate correction
    corrected_analysis = detect_terminal_similarity(corrected_seq)
    
    if corrected_analysis:
        print(f"After correction:")
        print(f"  Direct similarity: {corrected_analysis['direct_similarity']:.1f}%")
        print(f"  Inverted similarity: {corrected_analysis['inverted_similarity']:.1f}%")
        
        if corrected_analysis['inverted_similarity'] > analysis['inverted_similarity']:
            print("Correction improved ITR structure")
        else:
            print("Warning: Correction may not have improved structure")
    
    # Save corrected assembly
    output_dir = Path(assembly_path).parent
    output_file = output_dir / f"{sample}.contigs.corrected.fasta"
    
    with open(output_file, 'w') as f:
        f.write(f">{record.id}_corrected\n")
        f.write(str(corrected_seq))
    
    print(f"Corrected assembly saved to: {output_file}")
    
    # Create report
    report_file = output_dir / f"{sample}_correction_report.txt"
    with open(report_file, 'w') as f:
        f.write(f"Concatenation Correction Report - {sample} {assembly_type}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Original file: {assembly_path}\n")
        f.write(f"Original length: {seq_len:,} bp\n")
        f.write(f"Corrected length: {len(corrected_seq):,} bp\n")
        f.write(f"Removed: {seq_len - len(corrected_seq):,} bp\n")
        f.write(f"Cut position: {cut_point:,}\n\n")
        
        f.write("Terminal Analysis:\n")
        f.write(f"Original direct similarity: {analysis['direct_similarity']:.1f}%\n")
        f.write(f"Original inverted similarity: {analysis['inverted_similarity']:.1f}%\n")
        
        if corrected_analysis:
            f.write(f"Corrected direct similarity: {corrected_analysis['direct_similarity']:.1f}%\n")
            f.write(f"Corrected inverted similarity: {corrected_analysis['inverted_similarity']:.1f}%\n")
    
    print(f"Report saved to: {report_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python fix_concatenation.py SAMPLE_NAME ASSEMBLY_TYPE")
        print("Example: python fix_concatenation.py D1701 canu_ultra")
        print("Assembly types: canu, canu_ultra, miniasm")
        sys.exit(1)
    
    sample = sys.argv[1]
    assembly_type = sys.argv[2]
    
    if assembly_type not in ['canu', 'canu_ultra', 'miniasm']:
        print(f"Invalid assembly type: {assembly_type}")
        print("Valid types: canu, canu_ultra, miniasm")
        sys.exit(1)
    
    fix_assembly(sample, assembly_type)

if __name__ == "__main__":
    main()