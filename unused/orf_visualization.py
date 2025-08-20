#!/usr/bin/env python3
"""
Linear genome visualization with ORF annotations
Usage: python visualize_genome_orfs.py SAMPLE_NAME
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import numpy as np

def get_assembly_and_annotation_paths(sample):
    """Get paths to assemblies and their corresponding Prokka annotations"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    paths = {
        'miniasm': {
            'assembly': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
            'annotation': f"{base_dir}/05_annotation/{sample}/miniasm_viral/{sample}.gff"
        },
        'canu': {
            'assembly': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
            'annotation': f"{base_dir}/05_annotation/{sample}/canu/{sample}.gff"  # You'll need to run Prokka on Canu
        },
        'microsynth': {
            'assembly': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta",
            'annotation': f"{base_dir}/05_annotation/{sample}/microsynth/{sample}.gff"  # You'll need to run Prokka on Microsynth
        }
    }
    
    # Check which files exist
    available = {}
    for method, files in paths.items():
        if os.path.exists(files['assembly']):
            available[method] = files
            if not os.path.exists(files['annotation']):
                print(f"Warning: Assembly found but annotation missing for {method}: {files['annotation']}")
    
    return available

def parse_gff_file(gff_file):
    """Parse GFF file and extract ORF information"""
    orfs = []
    
    if not os.path.exists(gff_file):
        print(f"Warning: GFF file not found: {gff_file}")
        return orfs
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                orfs.append({
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length': end - start + 1,
                    'locus_tag': attr_dict.get('locus_tag', 'unknown'),
                    'gene': attr_dict.get('gene', ''),
                    'product': attr_dict.get('product', 'hypothetical protein')
                })
    
    return orfs

def get_genome_length(assembly_file):
    """Get genome length from assembly file"""
    if not os.path.exists(assembly_file):
        return 0
    
    record = next(SeqIO.parse(assembly_file, "fasta"))
    return len(record.seq)

def plot_linear_genome(assembly_method, genome_length, orfs, ax, y_position):
    """Plot a linear representation of the genome with ORFs"""
    
    # Draw genome backbone
    genome_line = patches.Rectangle((0, y_position - 0.05), genome_length, 0.1, 
                                   facecolor='lightgray', edgecolor='black', linewidth=1)
    ax.add_patch(genome_line)
    
    # Color scheme for ORFs
    colors = {'forward': '#2E86AB', 'reverse': '#A23B72'}
    
    # Plot ORFs
    for i, orf in enumerate(orfs):
        color = colors['forward'] if orf['strand'] == '+' else colors['reverse']
        
        # ORF height based on strand
        height = 0.15 if orf['strand'] == '+' else -0.15
        y_offset = y_position + (0.1 if orf['strand'] == '+' else -0.25)
        
        # Draw ORF rectangle
        orf_rect = patches.Rectangle((orf['start'], y_offset), orf['length'], height,
                                   facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.7)
        ax.add_patch(orf_rect)
        
        # Add ORF label for larger ORFs
        if orf['length'] > genome_length * 0.01:  # Only label ORFs > 1% of genome length
            label_y = y_offset + height/2
            label_text = orf['locus_tag'].replace(f'{assembly_method}_', '') if orf['locus_tag'] != 'unknown' else f"ORF{i+1}"
            
            ax.text(orf['start'] + orf['length']/2, label_y, label_text,
                   ha='center', va='center', fontsize=6, rotation=90 if orf['length'] < genome_length * 0.05 else 0)
    
    # Add method label
    ax.text(-genome_length * 0.05, y_position, assembly_method.capitalize(), 
           ha='right', va='center', fontweight='bold', fontsize=10)
    
    # Add scale marks
    for pos in range(0, genome_length, 20000):
        ax.axvline(x=pos, ymin=(y_position - 0.3)/3, ymax=(y_position + 0.3)/3, 
                  color='gray', alpha=0.5, linewidth=0.5)
        if pos % 40000 == 0:  # Label every 40kb
            ax.text(pos, y_position - 0.4, f'{pos//1000}kb', ha='center', va='top', fontsize=8)

def create_genome_comparison_plot(sample, available_data):
    """Create comprehensive genome comparison plot"""
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Get genome lengths and prepare data
    plot_data = []
    max_length = 0
    
    for method, paths in available_data.items():
        genome_length = get_genome_length(paths['assembly'])
        orfs = parse_gff_file(paths['annotation'])
        
        if genome_length > 0:
            plot_data.append({
                'method': method,
                'length': genome_length,
                'orfs': orfs,
                'orf_count': len(orfs)
            })
            max_length = max(max_length, genome_length)
    
    if not plot_data:
        print("No valid assembly data found for plotting")
        return
    
    # Sort by genome length
    plot_data.sort(key=lambda x: x['length'], reverse=True)
    
    # Plot each genome
    y_positions = np.arange(len(plot_data), 0, -1)
    
    for i, data in enumerate(plot_data):
        plot_linear_genome(data['method'], data['length'], data['orfs'], ax, y_positions[i])
    
    # Customize plot
    ax.set_xlim(-max_length * 0.1, max_length * 1.05)
    ax.set_ylim(0, len(plot_data) + 1)
    
    ax.set_xlabel('Genome Position (bp)', fontsize=12)
    ax.set_title(f'Linear Genome Comparison - Sample {sample}\nBlue: Forward strand ORFs, Purple: Reverse strand ORFs', 
                fontsize=14, pad=20)
    
    # Remove y-axis ticks
    ax.set_yticks([])
    
    # Add legend
    forward_patch = patches.Patch(color='#2E86AB', alpha=0.7, label='Forward strand (+)')
    reverse_patch = patches.Patch(color='#A23B72', alpha=0.7, label='Reverse strand (-)')
    ax.legend(handles=[forward_patch, reverse_patch], loc='upper right')
    
    # Add statistics text
    stats_text = "Assembly Statistics:\n"
    for data in plot_data:
        stats_text += f"{data['method'].capitalize()}: {data['length']:,} bp, {data['orf_count']} ORFs\n"
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_genome_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Genome comparison plot saved to: {output_file}")
    
    plt.show()

def create_detailed_orf_table(sample, available_data):
    """Create detailed ORF comparison table"""
    print(f"\n=== DETAILED ORF COMPARISON - {sample} ===")
    
    all_orfs = {}
    for method, paths in available_data.items():
        if os.path.exists(paths['annotation']):
            orfs = parse_gff_file(paths['annotation'])
            all_orfs[method] = orfs
    
    if not all_orfs:
        print("No annotation data available")
        return
    
    # Print ORF counts
    print(f"\nORF Counts:")
    for method, orfs in all_orfs.items():
        print(f"  {method.capitalize()}: {len(orfs)} ORFs")
    
    # Show first 10 ORFs for each method
    print(f"\nFirst 10 ORFs per assembly:")
    for method, orfs in all_orfs.items():
        print(f"\n{method.upper()}:")
        print(f"{'#':<3} {'Start':<8} {'End':<8} {'Strand':<7} {'Length':<7} {'Product':<30}")
        print("-" * 70)
        for i, orf in enumerate(orfs[:10]):
            product = orf['product'][:27] + '...' if len(orf['product']) > 30 else orf['product']
            print(f"{i+1:<3} {orf['start']:<8} {orf['end']:<8} {orf['strand']:<7} {orf['length']:<7} {product:<30}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python visualize_genome_orfs.py SAMPLE_NAME")
        print("Example: python visualize_genome_orfs.py D1701")
        print("\nNote: Make sure you have run Prokka on all assemblies first:")
        print("  ./05_annotation_prokka.sh SAMPLE --assembly canu")
        print("  ./05_annotation_prokka.sh SAMPLE --assembly microsynth")
        sys.exit(1)
    
    sample = sys.argv[1]
    
    # Get available data
    available_data = get_assembly_and_annotation_paths(sample)
    
    if not available_data:
        print(f"No assembly files found for sample {sample}")
        sys.exit(1)
    
    print(f"Found assemblies for: {', '.join(available_data.keys())}")
    
    # Create visualizations
    create_genome_comparison_plot(sample, available_data)
    create_detailed_orf_table(sample, available_data)

if __name__ == "__main__":
    main()