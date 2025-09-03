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
        #'miniasm': {
        #    'assembly': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
        #    'annotation': f"{base_dir}/05_annotation/{sample}/miniasm_viral/{sample}.gff"
        #},
        'canu': {
            'assembly': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu/{sample}_orf_annotated.gff"  
        },
        'microsynth': {
            'assembly': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/microsynth/{sample}_orf_annotated.gff" 
        },
        'canu_ultra': {
            'assembly': f"{base_dir}/03_assembly/{sample}/canu_ultra_output/{sample}.contigs.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu_ultra/{sample}_orf_annotated.gff"
        },
        'canu_ultra_trimmed': {
            'assembly': f"{base_dir}/06_trimmed_assembly/{sample}/canu_ultra/{sample}_canu_ultra_trimmed.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu_ultra_trimmed/{sample}_orf_annotated.gff"
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

def parse_gff_file(gff_file, assembly_file=None):
    """Parse GFF file and extract ORF information from longest contig only"""
    orfs = []
    
    if not os.path.exists(gff_file):
        print(f"Warning: GFF file not found: {gff_file}")
        return orfs
    
    # Get the longest contig name if assembly file is provided
    target_contig = None
    if assembly_file and os.path.exists(assembly_file):
        longest_record = get_longest_contig(assembly_file)
        if longest_record:
            target_contig = longest_record.id
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                
                # Skip if we're filtering by contig and this isn't the target
                if target_contig and parts[0] != target_contig:
                    continue
                
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
                
                # Extract original ORF number from inference field
                inference = attr_dict.get('inference', '')
                original_orf = extract_original_orf_from_inference(inference)
                
                orfs.append({
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length': end - start + 1,
                    'locus_tag': attr_dict.get('locus_tag', 'unknown'),
                    'gene': attr_dict.get('gene', ''),
                    'product': attr_dict.get('product', 'hypothetical protein'),
                    'original_orf': original_orf,
                    'contig': parts[0]
                })
    
    return orfs

def extract_original_orf_from_inference(inference_field):
    """Extract original ORF number from inference field"""
    import re
    
    # Look for pattern like "similar to AA sequence:D1701_ORFs_proteins.faa:ORF_134_R_1"
    pattern = r'similar to AA sequence:[^:]+:(ORF_\d+)'
    match = re.search(pattern, inference_field)
    
    if match:
        return match.group(1)
    
    return None

def get_genome_length(assembly_file):
    """Get genome length from assembly file"""
    record = get_longest_contig(assembly_file)
    return len(record.seq) if record else 0

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
        if orf['length'] > genome_length * 0.00005:  # 0.5% instead of 1%
            label_y = y_offset + height/2
            
            # Use original ORF number if available, otherwise use Prokka number
            if orf['original_orf']:
                label_text = orf['original_orf']
            else:
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
    fig, ax = plt.subplots(figsize=(55, 18))
    
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


def create_synteny_plot(sample, available_data):
    """Create synteny plot showing ORF connections between assemblies"""
    
    # Collect all ORF data
    assembly_data = {}
    for method, paths in available_data.items():
        if os.path.exists(paths['annotation']):
            orfs = parse_gff_file(paths['annotation'])
            genome_length = get_genome_length(paths['assembly'])
            assembly_data[method] = {
                'orfs': orfs,
                'length': genome_length
            }
    
    if len(assembly_data) < 2:
        print("Need at least 2 assemblies for synteny plot")
        return
    
    # Set up plot
    fig, ax = plt.subplots(figsize=(55, 18))
    
    # Plot parameters
    assembly_names = list(assembly_data.keys())
    y_positions = {name: i for i, name in enumerate(assembly_names)}
    max_length = max(data['length'] for data in assembly_data.values())
    
    # Draw each assembly
    for method, data in assembly_data.items():
        y_pos = y_positions[method]
        
        # Draw genome backbone
        genome_line = patches.Rectangle((0, y_pos - 0.05), data['length'], 0.1,
                                       facecolor='lightgray', edgecolor='black', linewidth=1)
        ax.add_patch(genome_line)
        
        # Draw ORFs
        for orf in data['orfs']:
            color = '#2E86AB' if orf['strand'] == '+' else '#A23B72'
            height = 0.15 if orf['strand'] == '+' else -0.15
            y_offset = y_pos + (0.1 if orf['strand'] == '+' else -0.25)
            
            orf_rect = patches.Rectangle((orf['start'], y_offset), orf['length'], height,
                                       facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.7)
            ax.add_patch(orf_rect)
        
        # Add assembly label
        ax.text(-max_length * 0.05, y_pos, method.capitalize(),
               ha='right', va='center', fontweight='bold', fontsize=12)
    
    # Find and draw ORF connections
    assembly_list = list(assembly_data.keys())
    connection_colors = ['red', 'blue', 'green', 'orange', 'purple']
    
    for i, assembly1 in enumerate(assembly_list):
        for j, assembly2 in enumerate(assembly_list[i+1:], i+1):
            color = connection_colors[min(i, len(connection_colors)-1)]
            
            # Find matching ORFs between assemblies
            matches = find_orf_matches(assembly_data[assembly1]['orfs'], 
                                     assembly_data[assembly2]['orfs'])
            
            # Draw connection lines
            for match in matches:
                orf1, orf2 = match
                y1 = y_positions[assembly1]
                y2 = y_positions[assembly2]
                
                # Calculate midpoints of ORFs
                x1 = orf1['start'] + orf1['length'] / 2
                x2 = orf2['start'] + orf2['length'] / 2
                
                # Draw connection line
                ax.plot([x1, x2], [y1, y2], color=color, alpha=0.3, linewidth=0.5)
    
    # Customize plot
    ax.set_xlim(-max_length * 0.1, max_length * 1.05)
    ax.set_ylim(-0.5, len(assembly_names) - 0.5)
    ax.set_xlabel('Genome Position (bp)', fontsize=12)
    ax.set_title(f'Assembly Synteny Comparison - Sample {sample}', fontsize=14, pad=20)
    ax.set_yticks([])
    
    # Add scale marks
    for pos in range(0, max_length, 20000):
        for y_pos in range(len(assembly_names)):
            ax.axvline(x=pos, ymin=(y_pos - 0.3)/len(assembly_names), 
                      ymax=(y_pos + 0.3)/len(assembly_names), 
                      color='gray', alpha=0.3, linewidth=0.5)
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_synteny_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Synteny plot saved to: {output_file}")
    
    plt.show()

def find_orf_matches(orfs1, orfs2):
    """Find matching ORFs between two assemblies based on original ORF numbers"""
    matches = []
    
    # Create lookup dict for orfs2
    orfs2_dict = {}
    for orf in orfs2:
        if orf['original_orf']:
            orfs2_dict[orf['original_orf']] = orf
    
    # Find matches
    for orf1 in orfs1:
        if orf1['original_orf'] and orf1['original_orf'] in orfs2_dict:
            matches.append((orf1, orfs2_dict[orf1['original_orf']]))
    
    return matches


def create_orf_presence_matrix(sample, available_data):
    """Create presence/absence matrix showing which ORFs are in each assembly"""
    
    # Collect all ORF data
    assembly_orfs = {}
    for method, paths in available_data.items():
        if os.path.exists(paths['annotation']):
            orfs = parse_gff_file(paths['annotation'])
            assembly_orfs[method] = orfs
    
    if len(assembly_orfs) < 2:
        print("Need at least 2 assemblies for presence/absence matrix")
        return
    
    # Get all unique ORF identifiers
    all_orfs = set()
    for orfs in assembly_orfs.values():
        for orf in orfs:
            if orf['original_orf']:
                all_orfs.add(orf['original_orf'])
    
    all_orfs = sorted(all_orfs)
    assembly_names = list(assembly_orfs.keys())
    
    # Create presence/absence matrix
    matrix = []
    for orf_id in all_orfs:
        row = []
        for assembly in assembly_names:
            # Check if this ORF exists in this assembly
            found = any(orf['original_orf'] == orf_id for orf in assembly_orfs[assembly])
            row.append(1 if found else 0)
        matrix.append(row)
    
    # Convert to numpy array for easier handling
    import numpy as np
    matrix = np.array(matrix)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, max(10, len(all_orfs) * 0.2)))
    
    # Create heatmap
    im = ax.imshow(matrix, cmap='RdYlBu_r', aspect='auto', vmin=0, vmax=1)
    
    # Set ticks and labels
    ax.set_xticks(range(len(assembly_names)))
    ax.set_xticklabels(assembly_names)
    ax.set_yticks(range(len(all_orfs)))
    ax.set_yticklabels(all_orfs, fontsize=8)
    
    # Add text annotations
    for i in range(len(all_orfs)):
        for j in range(len(assembly_names)):
            text = '●' if matrix[i, j] == 1 else '○'
            ax.text(j, i, text, ha="center", va="center", 
                   color='white' if matrix[i, j] == 1 else 'black', fontsize=12)
    
    # Customize plot
    ax.set_xlabel('Assembly Method', fontsize=12)
    ax.set_ylabel('ORF Identifier', fontsize=12)
    ax.set_title(f'ORF Presence/Absence Matrix - Sample {sample}', fontsize=14, pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label('ORF Present', rotation=270, labelpad=15)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['Absent', 'Present'])
    
    # Add statistics
    stats_text = "ORF Statistics:\n"
    total_orfs = len(all_orfs)
    for i, assembly in enumerate(assembly_names):
        present_count = np.sum(matrix[:, i])
        percentage = (present_count / total_orfs) * 100
        stats_text += f"{assembly}: {present_count}/{total_orfs} ({percentage:.1f}%)\n"
    
    # Add shared ORFs statistics
    if len(assembly_names) >= 2:
        shared_all = np.sum(np.all(matrix == 1, axis=1))
        shared_any = np.sum(np.any(matrix == 1, axis=1))
        stats_text += f"\nShared by all: {shared_all}\n"
        stats_text += f"Present in any: {shared_any}"
    
    ax.text(1.15, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_orf_presence_matrix.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"ORF presence matrix saved to: {output_file}")
    
    # Print summary table
    print(f"\n=== ORF PRESENCE SUMMARY - {sample} ===")
    print(f"Total unique ORFs identified: {total_orfs}")
    print(f"ORFs shared by all assemblies: {shared_all}")
    
    # Show which ORFs are missing from each assembly
    for i, assembly in enumerate(assembly_names):
        missing_orfs = [all_orfs[j] for j in range(len(all_orfs)) if matrix[j, i] == 0]
        if missing_orfs:
            print(f"\nMissing from {assembly}: {', '.join(missing_orfs[:10])}")
            if len(missing_orfs) > 10:
                print(f"  ... and {len(missing_orfs) - 10} more")
    
    plt.show()

def create_position_deviation_plot(sample, available_data):
    """Create plot showing positional differences between all assembly pairs"""
    
    # Collect all ORF data with positions
    assembly_orfs = {}
    for method, paths in available_data.items():
        if os.path.exists(paths['annotation']):
            orfs = parse_gff_file(paths['annotation'])
            # Create lookup dict: ORF_ID -> position
            orf_positions = {}
            for orf in orfs:
                if orf['original_orf']:
                    orf_positions[orf['original_orf']] = orf['start']
            assembly_orfs[method] = orf_positions
    
    if len(assembly_orfs) < 2:
        print("Need at least 2 assemblies for position deviation plot")
        return
    
    # Get all ORFs present in at least 2 assemblies
    all_orfs = set()
    for positions in assembly_orfs.values():
        all_orfs.update(positions.keys())
    
    # Filter to ORFs present in multiple assemblies
    shared_orfs = []
    for orf_id in all_orfs:
        assemblies_with_orf = [method for method, positions in assembly_orfs.items() 
                              if orf_id in positions]
        if len(assemblies_with_orf) >= 2:
            shared_orfs.append(orf_id)
    
    # Sort ORFs numerically if possible
    try:
        shared_orfs.sort(key=lambda x: int(x.split('_')[-1]))
    except:
        shared_orfs.sort()
    
    if not shared_orfs:
        print("No shared ORFs found between assemblies")
        return
    
    # Generate all pairwise comparisons
    assembly_names = list(assembly_orfs.keys())
    comparisons = []
    for i, assembly1 in enumerate(assembly_names):
        for assembly2 in assembly_names[i+1:]:
            comparisons.append((assembly1, assembly2))
    
    # Create the plot
    n_comparisons = len(comparisons)
    fig, axes = plt.subplots(n_comparisons, 1, figsize=(14, 5 * n_comparisons))
    if n_comparisons == 1:
        axes = [axes]
    
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
    
    for i, (assembly1, assembly2) in enumerate(comparisons):
        ax = axes[i]
        
        # Get positions for shared ORFs
        positions1 = []
        positions2 = []
        orf_names = []
        
        for orf_id in shared_orfs:
            if (orf_id in assembly_orfs[assembly1] and 
                orf_id in assembly_orfs[assembly2]):
                
                pos1 = assembly_orfs[assembly1][orf_id]
                pos2 = assembly_orfs[assembly2][orf_id]
                
                positions1.append(pos1)
                positions2.append(pos2)
                orf_names.append(orf_id)
        
        if not positions1:
            continue
        
        # Create scatter plot
        color = colors[i % len(colors)]
        ax.scatter(positions1, positions2, alpha=0.7, s=30, color=color)
        
        # Calculate correlation
        correlation = np.corrcoef(positions1, positions2)[0, 1] if len(positions1) > 1 else 0
        
        # Add diagonal line (perfect correlation)
        min_pos = min(min(positions1), min(positions2))
        max_pos = max(max(positions1), max(positions2))
        ax.plot([min_pos, max_pos], [min_pos, max_pos], 'k--', alpha=0.5, label='Perfect agreement')
        
        # Add trend line
        if len(positions1) > 1:
            z = np.polyfit(positions1, positions2, 1)
            p = np.poly1d(z)
            x_trend = np.array([min(positions1), max(positions1)])
            ax.plot(x_trend, p(x_trend), "r-", alpha=0.8, linewidth=2, label='Trend line')
            
            # Determine relationship type
            slope = z[0]
            intercept = z[1]
            
            if correlation < -0.7:
                relationship = "INVERTED"
                rel_color = "red"
            elif 0.8 <= slope <= 1.2 and abs(intercept) < max_pos * 0.1:
                relationship = "GOOD AGREEMENT"
                rel_color = "green"
            elif slope > 1.2:
                relationship = "EXPANDED"
                rel_color = "orange"
            elif slope < 0.8:
                relationship = "COMPRESSED"
                rel_color = "blue"
            else:
                relationship = "SHIFTED"
                rel_color = "purple"
        
        # Customize subplot
        ax.set_xlabel(f'ORF Position in {assembly1} (bp)', fontsize=10)
        ax.set_ylabel(f'ORF Position in {assembly2} (bp)', fontsize=10)
        ax.set_title(f'{assembly1} vs {assembly2}', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add statistics and relationship info
        if len(positions1) > 1:
            stats_text = f"Correlation: {correlation:.3f}\n"
            stats_text += f"Slope: {slope:.3f}\n"
            stats_text += f"Relationship: {relationship}"
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor=rel_color, alpha=0.3),
                   verticalalignment='top', fontsize=10)
            
            # Add interpretation
            if correlation < -0.7:
                interpretation = "Assemblies appear to be inverted relative to each other"
            elif correlation > 0.9 and 0.9 <= slope <= 1.1:
                interpretation = "Assemblies have very similar structure"
            elif slope > 1.5:
                interpretation = f"{assembly2} is expanded ~{slope:.1f}x relative to {assembly1}"
            elif slope < 0.5:
                interpretation = f"{assembly2} is compressed ~{1/slope:.1f}x relative to {assembly1}"
            else:
                interpretation = "Assemblies have moderate structural differences"
            
            ax.text(0.98, 0.02, interpretation, transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   verticalalignment='bottom', horizontalalignment='right', 
                   fontsize=9, style='italic')
    
    plt.suptitle(f'Assembly Structure Comparison - Sample {sample}', fontsize=16, y=0.98)
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_structure_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Structure comparison plot saved to: {output_file}")
    
    # Print summary
    print(f"\n=== STRUCTURAL COMPARISON SUMMARY - {sample} ===")
    for assembly1, assembly2 in comparisons:
        positions1 = [assembly_orfs[assembly1][orf] for orf in shared_orfs 
                     if orf in assembly_orfs[assembly1] and orf in assembly_orfs[assembly2]]
        positions2 = [assembly_orfs[assembly2][orf] for orf in shared_orfs 
                     if orf in assembly_orfs[assembly1] and orf in assembly_orfs[assembly2]]
        
        if len(positions1) > 1:
            correlation = np.corrcoef(positions1, positions2)[0, 1]
            slope = np.polyfit(positions1, positions2, 1)[0]
            
            print(f"\n{assembly1} vs {assembly2}:")
            print(f"  Correlation: {correlation:.3f}")
            print(f"  Slope: {slope:.3f}")
            
            if correlation < -0.7:
                print(f"  → Likely INVERTED assemblies")
            elif correlation > 0.9 and 0.9 <= slope <= 1.1:
                print(f"  → Good structural agreement")
            else:
                print(f"  → Significant structural differences")
    
    plt.show()


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
    #create_synteny_plot(sample, available_data)
    create_orf_presence_matrix(sample, available_data)
    create_position_deviation_plot(sample, available_data)
    create_detailed_orf_table(sample, available_data)

if __name__ == "__main__":
    main()