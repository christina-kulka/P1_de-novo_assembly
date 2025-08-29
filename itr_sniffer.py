#!/usr/bin/env python3
"""
ORF-based ITR analysis script
Identifies ITRs by finding duplicate ORFs at genome termini
Usage: python orf_based_itr_analysis.py SAMPLE_NAME
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import numpy as np
import importlib.util

def get_assembly_paths(sample):
    """Get paths to all assembly types"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    assemblies = {
        'miniasm': {
            'assembly': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/miniasm_viral/{sample}.gff"
        },
        'canu': {
            'assembly': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu/{sample}.gff"
        },
        'canu_ultra': {
            'assembly': f"{base_dir}/03_assembly/{sample}/canu_ultra_output/{sample}.longest_contig.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu_ultra/{sample}.gff"
        },
        'canu_ultra_trimmed': {
            'assembly': f"{base_dir}/06_trimmed_assembly/{sample}/canu_ultra/{sample}_canu_ultra_trimmed.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/canu_ultra_trimmed/{sample}.gff"
        },
        'microsynth': {
            'assembly': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta",
            'annotation': f"{base_dir}/04_annotation/{sample}/microsynth/{sample}.gff"
        }
    }
    
    # Check which files exist
    available = {}
    for method, files in assemblies.items():
        if os.path.exists(files['assembly']) and os.path.exists(files['annotation']):
            available[method] = files
    
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

def extract_original_orf_from_inference(inference_field):
    """Extract original ORF number from inference field"""
    import re
    
    pattern = r'similar to AA sequence:[^:]+:(ORF_\d+)'
    match = re.search(pattern, inference_field)
    
    if match:
        return match.group(1)
    
    return None

def parse_gff_file(gff_file, assembly_file):
    """Parse GFF file and extract ORF information from longest contig only"""
    orfs = []
    
    if not os.path.exists(gff_file):
        print(f"Warning: GFF file not found: {gff_file}")
        return orfs
    
    # Get the longest contig name
    longest_record = get_longest_contig(assembly_file)
    if not longest_record:
        return orfs
    
    target_contig = longest_record.id
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                
                # Skip if this isn't the target contig
                if parts[0] != target_contig:
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

def identify_terminal_orfs(orfs, genome_length, terminal_distance=40000):
    """Identify ORFs within terminal_distance of genome ends"""
    
    left_terminal_orfs = []
    right_terminal_orfs = []
    
    for orf in orfs:
        if orf['original_orf']:  # Only consider ORFs with reference identifiers
            orf_center = (orf['start'] + orf['end']) / 2
            
            # Left terminus
            if orf_center <= terminal_distance:
                left_terminal_orfs.append(orf)
            
            # Right terminus
            elif orf_center >= genome_length - terminal_distance:
                right_terminal_orfs.append(orf)
    
    return left_terminal_orfs, right_terminal_orfs

def find_matching_orfs(left_orfs, right_orfs):
    """Find ORFs that appear in both terminal regions"""
    
    # Create lookup dictionaries
    left_orf_dict = {orf['original_orf']: orf for orf in left_orfs if orf['original_orf']}
    right_orf_dict = {orf['original_orf']: orf for orf in right_orfs if orf['original_orf']}
    
    # Find matching ORF identifiers
    matching_orf_ids = set(left_orf_dict.keys()) & set(right_orf_dict.keys())
    
    matches = []
    for orf_id in matching_orf_ids:
        matches.append({
            'orf_id': orf_id,
            'left_orf': left_orf_dict[orf_id],
            'right_orf': right_orf_dict[orf_id]
        })
    
    return matches

def determine_itr_boundaries(matches, genome_length):
    """Determine ITR boundaries based on matching ORFs"""
    
    if not matches:
        return None
    
    # Find outermost and innermost matching ORFs
    left_positions = [match['left_orf']['start'] for match in matches]
    right_positions = [match['right_orf']['end'] for match in matches]
    
    # ITR boundaries
    left_itr_start = min(left_positions)
    left_itr_end = max(left_positions)
    
    right_itr_start = min(right_positions)
    right_itr_end = max(right_positions)
    
    # Calculate ITR characteristics
    left_itr_length = left_itr_end - left_itr_start
    right_itr_length = right_itr_end - right_itr_start
    
    return {
        'left_start': left_itr_start,
        'left_end': left_itr_end,
        'left_length': left_itr_length,
        'right_start': right_itr_start,
        'right_end': right_itr_end,
        'right_length': right_itr_length,
        'num_shared_orfs': len(matches),
        'shared_orfs': [match['orf_id'] for match in matches]
    }

def analyze_itr_structure(matches):
    """Analyze the structure and organization of ITR ORFs"""
    
    if not matches:
        return {'type': 'no_itrs', 'confidence': 0}
    
    # Sort matches by position
    matches_sorted_left = sorted(matches, key=lambda x: x['left_orf']['start'])
    matches_sorted_right = sorted(matches, key=lambda x: x['right_orf']['start'])
    
    # Check if ORF order is preserved
    left_order = [match['orf_id'] for match in matches_sorted_left]
    right_order = [match['orf_id'] for match in matches_sorted_right]
    right_order_reversed = list(reversed(right_order))
    
    if left_order == right_order:
        itr_type = 'direct_repeat'
    elif left_order == right_order_reversed:
        itr_type = 'inverted_repeat'
    else:
        itr_type = 'mixed_organization'
    
    # Calculate confidence based on number of shared ORFs
    confidence = min(100, len(matches) * 20)  # 20% per shared ORF, max 100%
    
    return {
        'type': itr_type,
        'confidence': confidence,
        'left_order': left_order,
        'right_order': right_order
    }

def analyze_assembly_itrs(sample, assembly_name, assembly_path, annotation_path):
    """Complete ITR analysis for one assembly"""
    
    print(f"\n=== ITR ANALYSIS: {assembly_name.upper()} ===")
    
    # Load genome information
    longest_contig = get_longest_contig(assembly_path)
    if not longest_contig:
        print(f"Could not load assembly: {assembly_path}")
        return None
    
    genome_length = len(longest_contig.seq)
    print(f"Genome length: {genome_length:,} bp")
    
    # Parse ORFs
    orfs = parse_gff_file(annotation_path, assembly_path)
    print(f"Total ORFs: {len(orfs)}")
    
    # Identify terminal ORFs
    left_terminal, right_terminal = identify_terminal_orfs(orfs, genome_length)
    print(f"Left terminal ORFs (40kb): {len(left_terminal)}")
    print(f"Right terminal ORFs (40kb): {len(right_terminal)}")
    
    # Find matching ORFs
    matches = find_matching_orfs(left_terminal, right_terminal)
    print(f"Shared ORFs between termini: {len(matches)}")
    
    if matches:
        shared_orf_names = [match['orf_id'] for match in matches]
        print(f"Shared ORFs: {', '.join(shared_orf_names[:5])}{'...' if len(shared_orf_names) > 5 else ''}")
    
    # Determine ITR boundaries
    itr_boundaries = determine_itr_boundaries(matches, genome_length)
    
    # Analyze ITR structure
    itr_structure = analyze_itr_structure(matches)
    
    if itr_boundaries:
        print(f"ITR type: {itr_structure['type']}")
        print(f"Confidence: {itr_structure['confidence']}%")
        print(f"Left ITR: {itr_boundaries['left_start']:,}-{itr_boundaries['left_end']:,} bp ({itr_boundaries['left_length']:,} bp)")
        print(f"Right ITR: {itr_boundaries['right_start']:,}-{itr_boundaries['right_end']:,} bp ({itr_boundaries['right_length']:,} bp)")
    else:
        print("No ITRs detected")
    
    return {
        'assembly_name': assembly_name,
        'genome_length': genome_length,
        'total_orfs': len(orfs),
        'left_terminal_orfs': left_terminal,
        'right_terminal_orfs': right_terminal,
        'matches': matches,
        'itr_boundaries': itr_boundaries,
        'itr_structure': itr_structure,
        'all_orfs': orfs
    }

def create_itr_visualization(sample, results):
    """Create comprehensive ITR visualization"""
    
    n_assemblies = len(results)
    fig, axes = plt.subplots(n_assemblies, 1, figsize=(16, 6 * n_assemblies))
    if n_assemblies == 1:
        axes = [axes]
    
    colors = {'shared': '#FF6B6B', 'left_only': '#4ECDC4', 'right_only': '#45B7D1', 'other': '#95A5A6'}
    
    for i, (assembly_name, data) in enumerate(results.items()):
        ax = axes[i]
        
        if not data:
            ax.text(0.5, 0.5, f"No data for {assembly_name}", ha='center', va='center', transform=ax.transAxes)
            continue
        
        genome_length = data['genome_length']
        
        # Draw genome backbone
        backbone = patches.Rectangle((0, 0.45), genome_length, 0.1, 
                                   facecolor='lightgray', edgecolor='black', linewidth=1)
        ax.add_patch(backbone)
        
        # Get shared ORF IDs for coloring
        shared_orf_ids = set()
        if data['matches']:
            shared_orf_ids = {match['orf_id'] for match in data['matches']}
        
        # Plot all ORFs
        for orf in data['all_orfs']:
            if orf['original_orf']:
                # Determine color based on ORF type
                if orf['original_orf'] in shared_orf_ids:
                    color = colors['shared']
                    alpha = 0.9
                elif orf['start'] <= 40000:
                    color = colors['left_only']
                    alpha = 0.7
                elif orf['end'] >= genome_length - 40000:
                    color = colors['right_only']
                    alpha = 0.7
                else:
                    color = colors['other']
                    alpha = 0.3
                
                # ORF position and size
                height = 0.15 if orf['strand'] == '+' else -0.15
                y_offset = 0.55 if orf['strand'] == '+' else 0.3
                
                # Draw ORF
                orf_rect = patches.Rectangle((orf['start'], y_offset), orf['length'], height,
                                           facecolor=color, edgecolor='black', linewidth=0.5, alpha=alpha)
                ax.add_patch(orf_rect)
                
                # Label shared ORFs
                if orf['original_orf'] in shared_orf_ids and orf['length'] > genome_length * 0.005:
                    ax.text(orf['start'] + orf['length']/2, y_offset + height/2, orf['original_orf'],
                           ha='center', va='center', fontsize=6, rotation=90, fontweight='bold')
        
        # Draw ITR boundaries if detected
        if data['itr_boundaries']:
            boundaries = data['itr_boundaries']
            
            # Left ITR boundary
            left_boundary = patches.Rectangle((boundaries['left_start'], 0.2), 
                                            boundaries['left_length'], 0.6, 
                                            facecolor='red', alpha=0.2, edgecolor='red', linewidth=2)
            ax.add_patch(left_boundary)
            
            # Right ITR boundary
            right_boundary = patches.Rectangle((boundaries['right_start'], 0.2), 
                                             boundaries['right_length'], 0.6, 
                                             facecolor='red', alpha=0.2, edgecolor='red', linewidth=2)
            ax.add_patch(right_boundary)
            
            # Add ITR labels
            ax.text(boundaries['left_start'] + boundaries['left_length']/2, 0.1, 
                   f"Left ITR\n{boundaries['left_length']:,} bp", ha='center', va='center', 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=8)
            
            ax.text(boundaries['right_start'] + boundaries['right_length']/2, 0.1, 
                   f"Right ITR\n{boundaries['right_length']:,} bp", ha='center', va='center', 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=8)
        
        # Add terminal region markers
        ax.axvline(x=40000, color='blue', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(x=genome_length-40000, color='blue', linestyle='--', alpha=0.5, linewidth=1)
        
        # Customize plot
        ax.set_xlim(-genome_length * 0.05, genome_length * 1.05)
        ax.set_ylim(0, 1)
        ax.set_xlabel('Genome Position (bp)', fontsize=10)
        ax.set_title(f'{assembly_name.title()}: ITR Analysis', fontsize=12)
        ax.set_yticks([])
        
        # Add assembly statistics
        if data['itr_structure']:
            stats_text = f"ITR Type: {data['itr_structure']['type']}\n"
            stats_text += f"Confidence: {data['itr_structure']['confidence']}%\n"
            stats_text += f"Shared ORFs: {len(data['matches'])}"
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
                   verticalalignment='top')
        
        # Add scale marks
        for pos in range(0, genome_length, 20000):
            ax.axvline(x=pos, ymin=0.4, ymax=0.6, color='gray', alpha=0.3, linewidth=0.5)
            if pos % 40000 == 0:
                ax.text(pos, 0.05, f'{pos//1000}kb', ha='center', va='bottom', fontsize=8)
    
    # Add legend
    legend_elements = [
        patches.Patch(color=colors['shared'], label='Shared ORFs (ITR)'),
        patches.Patch(color=colors['left_only'], label='Left terminal only'),
        patches.Patch(color=colors['right_only'], label='Right terminal only'),
        patches.Patch(color=colors['other'], alpha=0.3, label='Other ORFs')
    ]
    
    axes[-1].legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(1, 0))
    
    plt.suptitle(f'ORF-based ITR Analysis - Sample {sample}', fontsize=16, y=0.98)
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/06_itr_analysis/{sample}")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / f"{sample}_orf_itr_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nITR analysis plot saved to: {output_file}")
    
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python orf_based_itr_analysis.py SAMPLE_NAME")
        print("Example: python orf_based_itr_analysis.py D1701")
        sys.exit(1)
    
    sample = sys.argv[1]
    
    # Get available assemblies
    available_data = get_assembly_paths(sample)
    
    if not available_data:
        print(f"No assembly/annotation pairs found for sample {sample}")
        sys.exit(1)
    
    print(f"Analyzing {sample} assemblies: {', '.join(available_data.keys())}")
    
    # Analyze each assembly
    results = {}
    for assembly_name, paths in available_data.items():
        result = analyze_assembly_itrs(sample, assembly_name, paths['assembly'], paths['annotation'])
        results[assembly_name] = result
    
    # Create visualization
    create_itr_visualization(sample, results)
    
    # Summary comparison
    print(f"\n=== ITR COMPARISON SUMMARY - {sample} ===")
    print(f"{'Assembly':<15} {'Shared ORFs':<12} {'ITR Type':<18} {'Confidence':<12} {'ITR Quality'}")
    print("-" * 80)
    
    for assembly_name, data in results.items():
        if data and data['itr_structure']:
            shared_orfs = len(data['matches'])
            itr_type = data['itr_structure']['type']
            confidence = data['itr_structure']['confidence']
            
            if confidence >= 80:
                quality = "Excellent"
            elif confidence >= 60:
                quality = "Good"
            elif confidence >= 40:
                quality = "Moderate"
            else:
                quality = "Poor"
            
            print(f"{assembly_name.capitalize():<15} {shared_orfs:<12} {itr_type:<18} {confidence:<12}% {quality}")
        else:
            print(f"{assembly_name.capitalize():<15} {'N/A':<12} {'N/A':<18} {'N/A':<12} {'Failed'}")

if __name__ == "__main__":
    main()