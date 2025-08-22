#!/usr/bin/env python3
"""
Assembly structural investigation script
Investigates inversions, concatenations, and other structural issues
Usage: python investigate_assembly_issues.py SAMPLE_NAME
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import tempfile
import subprocess

def get_assembly_paths(sample):
    """Get paths to all assembly types"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    assemblies = {
        'miniasm': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
        'canu': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
        'microsynth': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta",
        'canu_ultra': f"{base_dir}/03_assembly/{sample}/canu_ultra_output/{sample}.contigs.fasta"
    }
    
    # Check which assemblies exist
    available = {}
    for name, path in assemblies.items():
        if os.path.exists(path):
            available[name] = path
    
    return available

def get_longest_contig(fasta_file):
    """Get the longest contig from a FASTA file"""
    if not os.path.exists(fasta_file):
        return None
    
    longest_record = None
    max_length = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > max_length:
            max_length = len(record.seq)
            longest_record = record
    
    return longest_record

def load_sequence(fasta_file):
    """Load sequence from FASTA file"""
    record = get_longest_contig(fasta_file)
    return record.seq if record else None

def test_inversion(seq1, seq2, sample_size=10000):
    """Test if seq2 might be inverted relative to seq1"""
    
    # Take samples from both ends
    seq1_start = seq1[:sample_size]
    seq1_end = seq1[-sample_size:]
    seq2_start = seq2[:sample_size]
    seq2_end = seq2[-sample_size:]
    
    # Create reverse complement of seq2
    seq2_rc = seq2.reverse_complement()
    seq2_rc_start = seq2_rc[:sample_size]
    seq2_rc_end = seq2_rc[-sample_size:]
    
    # Test alignments using simple similarity
    def simple_similarity(s1, s2):
        min_len = min(len(s1), len(s2))
        s1, s2 = s1[:min_len], s2[:min_len]
        matches = sum(1 for a, b in zip(s1, s2) if a == b)
        return (matches / min_len) * 100 if min_len > 0 else 0
    
    # Test normal orientation
    normal_start = simple_similarity(seq1_start, seq2_start)
    normal_end = simple_similarity(seq1_end, seq2_end)
    normal_avg = (normal_start + normal_end) / 2
    
    # Test inverted orientation
    inverted_start = simple_similarity(seq1_start, seq2_rc_start)
    inverted_end = simple_similarity(seq1_end, seq2_rc_end)
    inverted_avg = (inverted_start + inverted_end) / 2
    
    # Test cross-alignment (start vs end)
    cross_normal = simple_similarity(seq1_start, seq2_end)
    cross_inverted = simple_similarity(seq1_start, seq2_rc_end)
    
    return {
        'normal_similarity': normal_avg,
        'inverted_similarity': inverted_avg,
        'cross_normal': cross_normal,
        'cross_inverted': cross_inverted,
        'likely_inverted': inverted_avg > normal_avg and inverted_avg > 70
    }

def detect_concatenation(sequence, expected_size=138000, threshold=0.8):
    """Detect potential concatenation by looking for repeated patterns"""
    
    seq_len = len(sequence)
    if seq_len < expected_size * 1.3:  # Less than 30% larger
        return {'is_concatenated': False, 'evidence': 'Size within normal range'}
    
    # Test if sequence could be ~2x concatenated
    if seq_len > expected_size * 1.5:
        mid_point = seq_len // 2
        first_half = sequence[:mid_point]
        second_half = sequence[mid_point:]
        
        # Compare first and second half
        similarity = compare_sequence_similarity(first_half, second_half)
        
        if similarity > threshold * 100:
            return {
                'is_concatenated': True,
                'evidence': f'High similarity ({similarity:.1f}%) between first and second half',
                'suggested_cut': mid_point,
                'type': 'tandem_duplication'
            }
    
    # Test for partial concatenation (overlapping ends)
    overlap_sizes = [1000, 2000, 5000, 10000]
    best_overlap = 0
    best_similarity = 0
    
    for overlap_size in overlap_sizes:
        if overlap_size > seq_len // 4:
            continue
            
        start_seq = sequence[:overlap_size]
        end_seq = sequence[-overlap_size:]
        
        similarity = compare_sequence_similarity(start_seq, end_seq)
        if similarity > best_similarity:
            best_similarity = similarity
            best_overlap = overlap_size
    
    if best_similarity > 80:
        return {
            'is_concatenated': True,
            'evidence': f'High terminal overlap ({best_similarity:.1f}%) over {best_overlap}bp',
            'suggested_cut': seq_len - best_overlap,
            'type': 'circular_overlap'
        }
    
    return {
        'is_concatenated': False,
        'evidence': f'No clear concatenation pattern found (best terminal similarity: {best_similarity:.1f}%)'
    }

def compare_sequence_similarity(seq1, seq2):
    """Compare similarity between two sequences"""
    min_len = min(len(seq1), len(seq2))
    seq1, seq2 = seq1[:min_len], seq2[:min_len]
    
    if min_len == 0:
        return 0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / min_len) * 100

def analyze_terminal_repeats(sequence, max_size=5000):
    """Analyze terminal repeats (potential ITRs)"""
    
    results = {}
    
    for size in [500, 1000, 2000, 3000, max_size]:
        if size > len(sequence) // 3:
            continue
            
        left_term = sequence[:size]
        right_term = sequence[-size:]
        right_term_rc = right_term.reverse_complement()
        
        # Test ITR similarity
        similarity = compare_sequence_similarity(left_term, right_term_rc)
        results[size] = similarity
    
    best_size = max(results.keys(), key=lambda x: results[x])
    best_similarity = results[best_size]
    
    return {
        'best_itr_size': best_size,
        'best_itr_similarity': best_similarity,
        'all_similarities': results,
        'has_good_itrs': best_similarity > 70
    }

def investigate_assembly(sample, assembly_name, assembly_path, reference_seq=None):
    """Complete investigation of a single assembly"""
    
    print(f"\n=== INVESTIGATING {assembly_name.upper()} ASSEMBLY ===")
    
    # Load sequence
    sequence = load_sequence(assembly_path)
    seq_length = len(sequence)
    
    print(f"Assembly length: {seq_length:,} bp")
    
    # Test for concatenation
    concat_result = detect_concatenation(sequence)
    print(f"Concatenation test: {concat_result['evidence']}")
    
    # Test terminal repeats (ITRs)
    itr_result = analyze_terminal_repeats(sequence)
    print(f"Best ITR: {itr_result['best_itr_similarity']:.1f}% similarity over {itr_result['best_itr_size']}bp")
    
    # Test inversion if reference provided
    inversion_result = None
    if reference_seq is not None:
        inversion_result = test_inversion(reference_seq, sequence)
        if inversion_result['likely_inverted']:
            print(f"INVERSION DETECTED: {inversion_result['inverted_similarity']:.1f}% similarity when inverted")
        else:
            print(f"Normal orientation: {inversion_result['normal_similarity']:.1f}% similarity")
    
    return {
        'length': seq_length,
        'concatenation': concat_result,
        'itr': itr_result,
        'inversion': inversion_result,
        'sequence': sequence
    }

def create_investigation_plot(sample, results):
    """Create comprehensive investigation plot"""
    
    n_assemblies = len(results)
    fig, axes = plt.subplots(n_assemblies, 3, figsize=(18, 6 * n_assemblies))
    if n_assemblies == 1:
        axes = axes.reshape(1, -1)
    
    for i, (assembly_name, data) in enumerate(results.items()):
        
        # Plot 1: Length comparison
        ax1 = axes[i, 0]
        expected_length = 138000
        actual_length = data['length']
        
        bars = ax1.bar(['Expected', assembly_name.title()], 
                      [expected_length, actual_length],
                      color=['gray', 'skyblue'])
        
        # Color bar based on issues
        if data['concatenation']['is_concatenated']:
            bars[1].set_color('red')
        elif abs(actual_length - expected_length) > 20000:
            bars[1].set_color('orange')
        else:
            bars[1].set_color('green')
            
        ax1.set_ylabel('Assembly Length (bp)')
        ax1.set_title(f'{assembly_name.title()}: Length Analysis')
        
        # Add text annotation
        diff = actual_length - expected_length
        ax1.text(1, actual_length + 5000, f'{diff:+,} bp', ha='center', fontweight='bold')
        
        # Plot 2: ITR analysis
        ax2 = axes[i, 1]
        itr_data = data['itr']['all_similarities']
        sizes = list(itr_data.keys())
        similarities = list(itr_data.values())
        
        ax2.plot(sizes, similarities, 'o-', linewidth=2, markersize=8)
        ax2.axhline(y=70, color='red', linestyle='--', alpha=0.7, label='Good ITR threshold')
        ax2.set_xlabel('ITR Size (bp)')
        ax2.set_ylabel('Similarity (%)')
        ax2.set_title(f'{assembly_name.title()}: ITR Quality')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Issues summary
        ax3 = axes[i, 2]
        ax3.axis('off')
        
        # Create summary text
        summary_text = f"{assembly_name.upper()} SUMMARY\n\n"
        summary_text += f"Length: {actual_length:,} bp\n"
        summary_text += f"Size difference: {diff:+,} bp\n\n"
        
        # Concatenation
        if data['concatenation']['is_concatenated']:
            summary_text += "⚠️  CONCATENATION DETECTED\n"
            summary_text += f"   {data['concatenation']['evidence']}\n\n"
        else:
            summary_text += "✓  No concatenation detected\n\n"
        
        # ITR quality
        best_itr = data['itr']['best_itr_similarity']
        if best_itr > 80:
            summary_text += f"✓  Good ITRs: {best_itr:.1f}%\n\n"
        elif best_itr > 60:
            summary_text += f"⚠️  Moderate ITRs: {best_itr:.1f}%\n\n"
        else:
            summary_text += f"❌ Poor ITRs: {best_itr:.1f}%\n\n"
        
        # Inversion
        if data['inversion']:
            if data['inversion']['likely_inverted']:
                summary_text += "🔄 INVERSION DETECTED\n"
                summary_text += f"   {data['inversion']['inverted_similarity']:.1f}% when inverted\n\n"
            else:
                summary_text += "✓  Correct orientation\n\n"
        
        # Overall assessment
        issues = []
        if data['concatenation']['is_concatenated']:
            issues.append('concatenation')
        if best_itr < 70:
            issues.append('poor ITRs')
        if data['inversion'] and data['inversion']['likely_inverted']:
            issues.append('inversion')
        
        if not issues:
            summary_text += "Overall: ✓ GOOD ASSEMBLY"
            text_color = 'green'
        elif len(issues) == 1:
            summary_text += f"Overall: ⚠️  MINOR ISSUES\n({', '.join(issues)})"
            text_color = 'orange'
        else:
            summary_text += f"Overall: ❌ MAJOR ISSUES\n({', '.join(issues)})"
            text_color = 'red'
        
        ax3.text(0.05, 0.95, summary_text, transform=ax3.transAxes,
                fontsize=11, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor=text_color, alpha=0.1))
    
    plt.suptitle(f'Assembly Investigation - Sample {sample}', fontsize=16, y=0.98)
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_assembly_investigation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nInvestigation plot saved to: {output_file}")
    
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python investigate_assembly_issues.py SAMPLE_NAME")
        print("Example: python investigate_assembly_issues.py D1701")
        sys.exit(1)
    
    sample = sys.argv[1]
    
    # Get available assemblies
    assemblies = get_assembly_paths(sample)
    
    if not assemblies:
        print(f"No assemblies found for sample {sample}")
        sys.exit(1)
    
    print(f"Investigating assemblies for {sample}: {', '.join(assemblies.keys())}")
    
    # Choose reference assembly (shortest one, likely most accurate)
    reference_assembly = min(assemblies.keys(), key=lambda x: len(load_sequence(assemblies[x])))
    reference_seq = load_sequence(assemblies[reference_assembly])
    print(f"Using {reference_assembly} as reference for inversion testing")
    
    # Investigate each assembly
    results = {}
    for name, path in assemblies.items():
        ref_seq = reference_seq if name != reference_assembly else None
        results[name] = investigate_assembly(sample, name, path, ref_seq)
    
    # Create comprehensive plot
    create_investigation_plot(sample, results)
    
    # Print final recommendations
    print(f"\n=== RECOMMENDATIONS FOR {sample} ===")
    for name, data in results.items():
        issues = []
        fixes = []
        
        if data['concatenation']['is_concatenated']:
            issues.append("concatenated")
            if 'suggested_cut' in data['concatenation']:
                fixes.append(f"trim to {data['concatenation']['suggested_cut']:,} bp")
        
        if data['inversion'] and data['inversion']['likely_inverted']:
            issues.append("inverted")
            fixes.append("reverse complement")
        
        if data['itr']['best_itr_similarity'] < 70:
            issues.append("poor ITRs")
        
        if not issues:
            print(f"{name}: ✓ Use as-is")
        else:
            print(f"{name}: Issues - {', '.join(issues)}")
            if fixes:
                print(f"   Suggested fixes: {', '.join(fixes)}")

if __name__ == "__main__":
    main()