"""
Comprehensive ITR Analysis Script
Usage: python 05_itr_analysis.py SAMPLE_NAME [OPTIONS]
"""

import sys
import os
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def load_config():
    """Load configuration variables"""
    config = {
        'base_dir': '/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly',
        'microsynth_dir': '00_raw_data_microsynth',
        'assembly_dir': '03_assembly',
        'results_dir': '06_itr_analysis'
    }
    return config

def get_assembly_paths(sample, config, assembly_type='flye'):
    """Get paths to different assembly types"""
    base_dir = config['base_dir']
    
    paths = {
        'microsynth': f"{base_dir}/{config['microsynth_dir']}/{sample}_results/{sample}_results/Assembly/{sample}.fasta",
        'miniasm_raw': f"{base_dir}/{config['assembly_dir']}/{sample}/miniasm_out/polished_assembly.fasta",
        'miniasm_viral': f"{base_dir}/{config['assembly_dir']}/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
        'flye': f"{base_dir}/{config['assembly_dir']}/{sample}/flye_out/assembly.fasta",
        'canu': f"{base_dir}/{config['assembly_dir']}/{sample}/canu_out/{sample}.contigs.fasta"
    }
    
    return paths[assembly_type] if assembly_type in paths else None

def get_sequence_info(fasta_file):
    """Get basic sequence information"""
    if not os.path.exists(fasta_file):
        return None
    
    records = list(SeqIO.parse(fasta_file, "fasta"))
    total_length = sum(len(record.seq) for record in records)
    
    return {
        'file': fasta_file,
        'contigs': len(records),
        'total_length': total_length,
        'longest_contig': max(len(record.seq) for record in records),
        'sequence': records[0].seq if records else None  # Use first/longest contig
    }

def calculate_similarity(seq1, seq2):
    """Calculate similarity percentage between two sequences"""
    if len(seq1) != len(seq2):
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    
    if len(seq1) == 0:
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100

def sliding_window_analysis(left_seq, right_seq, window_sizes, step_size=50):
    """Comprehensive sliding window analysis for ITR detection"""
    right_rc = right_seq.reverse_complement()
    results = []
    
    for window_size in window_sizes:
        if window_size > min(len(left_seq), len(right_rc)):
            continue
            
        max_similarity = 0
        best_left_pos = 0
        best_right_pos = 0
        
        # Slide across left sequence
        for left_pos in range(0, len(left_seq) - window_size + 1, step_size):
            left_window = left_seq[left_pos:left_pos + window_size]
            
            # Slide across right sequence (reverse complement)
            for right_pos in range(0, len(right_rc) - window_size + 1, step_size):
                right_window = right_rc[right_pos:right_pos + window_size]
                
                similarity = calculate_similarity(left_window, right_window)
                
                if similarity > max_similarity:
                    max_similarity = similarity
                    best_left_pos = left_pos
                    best_right_pos = right_pos
        
        results.append({
            'window_size': window_size,
            'max_similarity': max_similarity,
            'left_start': best_left_pos + 1,
            'left_end': best_left_pos + window_size,
            'right_start': len(right_seq) - best_right_pos - window_size + 1,
            'right_end': len(right_seq) - best_right_pos,
            'left_seq': str(left_seq[best_left_pos:best_left_pos + window_size]),
            'right_seq': str(right_seq[-(best_right_pos + window_size):-best_right_pos] if best_right_pos > 0 else right_seq[-window_size:])
        })
    
    return results

def analyze_multiple_thresholds(itr_results, thresholds=[60, 70, 80, 90]):
    """Analyze ITR results at multiple similarity thresholds"""
    threshold_results = {}
    
    for threshold in thresholds:
        passing_results = [r for r in itr_results if r['max_similarity'] >= threshold]
        if passing_results:
            best_result = max(passing_results, key=lambda x: x['max_similarity'])
            threshold_results[threshold] = best_result
        else:
            threshold_results[threshold] = None
    
    return threshold_results

def create_output_directory(sample, config, assembly_type):
    """Create output directory for results"""
    output_dir = Path(config['base_dir']) / config['results_dir'] / sample / assembly_type
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir

def save_itr_sequences(sample, itr_results, output_dir):
    """Save ITR sequences to FASTA files"""
    if not itr_results:
        return
    
    best_result = max(itr_results, key=lambda x: x['max_similarity'])
    
    # Save left ITR
    left_file = output_dir / f"{sample}_left_itr.fasta"
    with open(left_file, 'w') as f:
        f.write(f">{sample}_left_ITR_{best_result['window_size']}bp\n")
        f.write(f"{best_result['left_seq']}\n")
    
    # Save right ITR
    right_file = output_dir / f"{sample}_right_itr.fasta"
    with open(right_file, 'w') as f:
        f.write(f">{sample}_right_ITR_{best_result['window_size']}bp\n")
        f.write(f"{best_result['right_seq']}\n")
    
    print(f"ITR sequences saved to {left_file} and {right_file}")

def plot_itr_analysis(sample, itr_results, threshold_results, output_dir):
    """Create visualization plots for ITR analysis"""
    if not itr_results:
        return
    
    # Create dataframe for plotting
    df = pd.DataFrame(itr_results)
    
    # Plot 1: Similarity vs Window Size
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(df['window_size'], df['max_similarity'], 'bo-', linewidth=2, markersize=8)
    plt.xlabel('Window Size (bp)')
    plt.ylabel('Max Similarity (%)')
    plt.title(f'{sample}: ITR Similarity vs Window Size')
    plt.grid(True, alpha=0.3)
    
    # Add threshold lines
    for threshold in [60, 70, 80, 90]:
        plt.axhline(y=threshold, color='red', linestyle='--', alpha=0.5, label=f'{threshold}%')
    
    # Plot 2: ITR positions
    plt.subplot(2, 2, 2)
    best_result = max(itr_results, key=lambda x: x['max_similarity'])
    
    positions = ['Left Start', 'Left End', 'Right Start', 'Right End']
    values = [best_result['left_start'], best_result['left_end'], 
              best_result['right_start'], best_result['right_end']]
    
    plt.bar(positions, values, color=['blue', 'lightblue', 'red', 'lightcoral'])
    plt.ylabel('Genome Position (bp)')
    plt.title(f'{sample}: Best ITR Positions ({best_result["window_size"]}bp, {best_result["max_similarity"]:.1f}%)')
    plt.xticks(rotation=45)
    
    # Plot 3: Threshold summary
    plt.subplot(2, 2, 3)
    threshold_data = []
    for thresh, result in threshold_results.items():
        if result:
            threshold_data.append([thresh, result['max_similarity'], result['window_size']])
    
    if threshold_data:
        thresh_df = pd.DataFrame(threshold_data, columns=['Threshold', 'Similarity', 'Window_Size'])
        plt.scatter(thresh_df['Threshold'], thresh_df['Similarity'], 
                   s=thresh_df['Window_Size']/10, alpha=0.7, c=thresh_df['Window_Size'], cmap='viridis')
        plt.xlabel('Similarity Threshold (%)')
        plt.ylabel('Achieved Similarity (%)')
        plt.title(f'{sample}: ITR Detection at Different Thresholds')
        plt.colorbar(label='Window Size (bp)')
    
    # Plot 4: Summary table
    plt.subplot(2, 2, 4)
    plt.axis('off')
    
    summary_text = f"Sample: {sample}\n\n"
    summary_text += f"Best ITR Found:\n"
    summary_text += f"  Size: {best_result['window_size']} bp\n"
    summary_text += f"  Similarity: {best_result['max_similarity']:.1f}%\n"
    summary_text += f"  Left: {best_result['left_start']}-{best_result['left_end']}\n"
    summary_text += f"  Right: {best_result['right_start']}-{best_result['right_end']}\n\n"
    
    summary_text += "Threshold Analysis:\n"
    for thresh in [90, 80, 70, 60]:
        result = threshold_results.get(thresh)
        if result:
            summary_text += f"  {thresh}%: {result['window_size']}bp ({result['max_similarity']:.1f}%)\n"
        else:
            summary_text += f"  {thresh}%: No ITR found\n"
    
    plt.text(0.1, 0.9, summary_text, transform=plt.gca().transAxes, 
             fontfamily='monospace', fontsize=10, verticalalignment='top')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = output_dir / f"{sample}_itr_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Analysis plot saved to {plot_file}")

def analyze_sample_itr(sample, assembly_type='flye', terminal_size=5000):
    """Complete ITR analysis for a sample"""
    
    config = load_config()
    
    print(f"=== ITR ANALYSIS FOR {sample} ({assembly_type.upper()}) ===")
    
    # Get assembly path
    assembly_path = get_assembly_paths(sample, config, assembly_type)
    if not assembly_path or not os.path.exists(assembly_path):
        print(f"Error: Assembly not found at {assembly_path}")
        return None
    
    # Get sequence information
    seq_info = get_sequence_info(assembly_path)
    if not seq_info or not seq_info['sequence']:
        print(f"Error: Could not read sequence from {assembly_path}")
        return None
    
    print(f"Assembly: {seq_info['contigs']} contig(s), {seq_info['total_length']:,} bp total")
    
    genome_seq = seq_info['sequence']
    
    # Extract terminal regions
    left_terminal = genome_seq[:terminal_size]
    right_terminal = genome_seq[-terminal_size:]
    
    print(f"Analyzing {terminal_size:,} bp from each terminus...")
    
    # Define window sizes to test
    window_sizes = [50, 100, 200, 500, 1000, 1500, 2000, 3000]
    window_sizes = [w for w in window_sizes if w <= terminal_size]
    
    # Perform sliding window analysis
    itr_results = sliding_window_analysis(left_terminal, right_terminal, window_sizes)
    
    if not itr_results:
        print("No ITR analysis results generated")
        return None
    
    # Analyze at multiple thresholds
    threshold_results = analyze_multiple_thresholds(itr_results)
    
    # Create output directory
    output_dir = create_output_directory(sample, config, assembly_type)
    
    # Print results
    print(f"\nITR Analysis Results:")
    print(f"{'Window Size':<12} {'Max Similarity':<15} {'Left Pos':<15} {'Right Pos':<15}")
    print("-" * 65)
    
    for result in sorted(itr_results, key=lambda x: x['max_similarity'], reverse=True):
        print(f"{result['window_size']:<12} {result['max_similarity']:<15.1f} "
              f"{result['left_start']}-{result['left_end']:<8} "
              f"{result['right_start']}-{result['right_end']}")
    
    print(f"\nThreshold Analysis:")
    for threshold in [90, 80, 70, 60]:
        result = threshold_results.get(threshold)
        if result:
            print(f"  {threshold}% threshold: {result['window_size']}bp ITR with {result['max_similarity']:.1f}% similarity")
        else:
            print(f"  {threshold}% threshold: No ITR detected")
    
    # Save results
    save_itr_sequences(sample, itr_results, output_dir)
    plot_itr_analysis(sample, itr_results, threshold_results, output_dir)
    
    # Save detailed results to CSV
    results_df = pd.DataFrame(itr_results)
    csv_file = output_dir / f"{sample}_itr_results.csv"
    results_df.to_csv(csv_file, index=False)
    print(f"Detailed results saved to {csv_file}")
    
    return itr_results, threshold_results

def main():
    parser = argparse.ArgumentParser(description='Comprehensive ITR Analysis')
    parser.add_argument('sample', help='Sample name (e.g., B006, D1701)')
    parser.add_argument('--assembly', choices=['flye', 'miniasm_viral', 'miniasm_raw', 'microsynth', "canu"], 
                       default='flye', help='Assembly type to analyze')
    parser.add_argument('--terminal-size', type=int, default=5000, 
                       help='Size of terminal regions to extract (bp)')
    
    args = parser.parse_args()
    
    analyze_sample_itr(args.sample, args.assembly, args.terminal_size)

if __name__ == "__main__":
    main()