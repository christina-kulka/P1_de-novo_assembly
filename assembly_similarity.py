#!/usr/bin/env python3
"""
Assembly similarity analysis script
Compares assemblies using minimap2 and analyzes terminal overlaps
Usage: python assembly_similarity_analysis.py SAMPLE_NAME
"""

import sys
import os
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
import pandas as pd

def get_assembly_paths(sample):
    """Get paths to all assembly types"""
    base_dir = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
    
    assemblies = {
        'miniasm': f"{base_dir}/03_assembly/{sample}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta",
        'canu': f"{base_dir}/03_assembly/{sample}/canu_out/{sample}.contigs.fasta",
        'microsynth': f"{base_dir}/00_raw_data_microsynth/{sample}_results/{sample}_results/Assembly/{sample}.fasta"
    }
    
    # Check which assemblies exist
    available = {}
    for name, path in assemblies.items():
        if os.path.exists(path):
            available[name] = path
    
    return available

def get_assembly_stats(fasta_file):
    """Get basic assembly statistics"""
    stats = {'length': 0, 'contigs': 0, 'longest_contig': 0, 'name': ''}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        stats['contigs'] += 1
        contig_len = len(record.seq)
        stats['length'] += contig_len
        if contig_len > stats['longest_contig']:
            stats['longest_contig'] = contig_len
            stats['name'] = record.id
    
    return stats

def run_minimap2_alignment(reference_file, query_file):
    """Run minimap2 alignment between two assemblies"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False) as tmp_file:
            cmd = ['minimap2', reference_file, query_file]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            tmp_file.write(result.stdout)
            paf_file = tmp_file.name
        
        return paf_file
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e}")
        return None

def parse_paf_file(paf_file):
    """Parse PAF alignment file"""
    alignments = []
    
    with open(paf_file, 'r') as f:
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    alignment = {
                        'query_name': fields[0],
                        'query_length': int(fields[1]),
                        'query_start': int(fields[2]),
                        'query_end': int(fields[3]),
                        'strand': fields[4],
                        'target_name': fields[5],
                        'target_length': int(fields[6]),
                        'target_start': int(fields[7]),
                        'target_end': int(fields[8]),
                        'matches': int(fields[9]),
                        'alignment_length': int(fields[10]),
                        'mapping_quality': int(fields[11])
                    }
                    alignments.append(alignment)
    
    return alignments

def analyze_alignment_coverage(alignments, query_length, target_length):
    """Analyze alignment coverage statistics"""
    if not alignments:
        return {'query_coverage': 0, 'target_coverage': 0, 'identity': 0}
    
    # Calculate total aligned bases
    query_aligned = sum(aln['query_end'] - aln['query_start'] for aln in alignments)
    target_aligned = sum(aln['target_end'] - aln['target_start'] for aln in alignments)
    
    # Calculate coverage percentages
    query_coverage = (query_aligned / query_length) * 100 if query_length > 0 else 0
    target_coverage = (target_aligned / target_length) * 100 if target_length > 0 else 0
    
    # Calculate average identity
    total_matches = sum(aln['matches'] for aln in alignments)
    total_alignment = sum(aln['alignment_length'] for aln in alignments)
    identity = (total_matches / total_alignment) * 100 if total_alignment > 0 else 0
    
    return {
        'query_coverage': query_coverage,
        'target_coverage': target_coverage,
        'identity': identity,
        'alignments': len(alignments)
    }

def extract_and_compare_terminals(assembly1_path, assembly2_path, assembly1_name, assembly2_name, terminal_size=10000):
    """Extract terminal regions and compare them using BLAST-like analysis"""
    print(f"\n--- Terminal Region Analysis: {assembly1_name} vs {assembly2_name} ---")
    
    try:
        # Extract sequences
        record1 = next(SeqIO.parse(assembly1_path, "fasta"))
        record2 = next(SeqIO.parse(assembly2_path, "fasta"))
        
        seq1 = record1.seq
        seq2 = record2.seq
        
        # Extract terminals
        seq1_start = seq1[:terminal_size]
        seq1_end = seq1[-terminal_size:]
        seq2_start = seq2[:terminal_size]
        seq2_end = seq2[-terminal_size:]
        
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f1:
            f1.write(f">{assembly1_name}_start\n{seq1_start}\n")
            f1.write(f">{assembly1_name}_end\n{seq1_end}\n")
            temp1 = f1.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f2:
            f2.write(f">{assembly2_name}_start\n{seq2_start}\n")
            f2.write(f">{assembly2_name}_end\n{seq2_end}\n")
            temp2 = f2.name
        
        # Run BLAST comparison
        cmd = [
            'blastn', '-query', temp1, '-subject', temp2,
            '-outfmt', '6 qseqid sseqid qstart qend sstart send length pident evalue bitscore',
            '-word_size', '11'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse BLAST results
        blast_hits = []
        if result.stdout.strip():
            for line in result.stdout.strip().split('\n'):
                fields = line.split('\t')
                if len(fields) >= 10:
                    blast_hits.append({
                        'query': fields[0],
                        'subject': fields[1],
                        'length': int(fields[6]),
                        'identity': float(fields[7]),
                        'evalue': float(fields[8]),
                        'bitscore': float(fields[9])
                    })
        
        # Analyze results
        if blast_hits:
            high_identity_hits = [hit for hit in blast_hits if hit['identity'] >= 90]
            long_hits = [hit for hit in blast_hits if hit['length'] >= 50]
            
            print(f"Terminal comparison results:")
            print(f"  Total BLAST hits: {len(blast_hits)}")
            print(f"  High identity hits (≥90%): {len(high_identity_hits)}")
            print(f"  Long hits (≥50bp): {len(long_hits)}")
            
            if high_identity_hits:
                avg_identity = sum(hit['identity'] for hit in high_identity_hits) / len(high_identity_hits)
                max_length = max(hit['length'] for hit in high_identity_hits)
                print(f"  Average identity of high-identity hits: {avg_identity:.1f}%")
                print(f"  Longest high-identity hit: {max_length}bp")
        else:
            print(f"No significant terminal similarities found")
        
        # Cleanup
        os.unlink(temp1)
        os.unlink(temp2)
        
        return blast_hits
        
    except Exception as e:
        print(f"Error in terminal analysis: {e}")
        return []

def compare_assembly_pair(ref_path, query_path, ref_name, query_name):
    """Compare a pair of assemblies comprehensively"""
    print(f"\n=== COMPARING {ref_name.upper()} vs {query_name.upper()} ===")
    
    # Get assembly statistics
    ref_stats = get_assembly_stats(ref_path)
    query_stats = get_assembly_stats(query_path)
    
    print(f"{ref_name.capitalize()} assembly: {ref_stats['length']:,} bp, {ref_stats['contigs']} contig(s)")
    print(f"{query_name.capitalize()} assembly: {query_stats['length']:,} bp, {query_stats['contigs']} contig(s)")
    print(f"Size difference: {abs(ref_stats['length'] - query_stats['length']):,} bp")
    
    # Run minimap2 alignment
    print(f"\nRunning minimap2 alignment...")
    paf_file = run_minimap2_alignment(ref_path, query_path)
    
    if paf_file:
        alignments = parse_paf_file(paf_file)
        coverage_stats = analyze_alignment_coverage(alignments, query_stats['length'], ref_stats['length'])
        
        print(f"Alignment results:")
        print(f"  {query_name.capitalize()} coverage: {coverage_stats['query_coverage']:.1f}%")
        print(f"  {ref_name.capitalize()} coverage: {coverage_stats['target_coverage']:.1f}%")
        print(f"  Average identity: {coverage_stats['identity']:.1f}%")
        print(f"  Number of alignments: {coverage_stats['alignments']}")
        
        # Show alignment details
        if alignments:
            print(f"\nAlignment segments:")
            for i, aln in enumerate(alignments[:5]):  # Show first 5 alignments
                print(f"  Segment {i+1}: {query_name}[{aln['query_start']}-{aln['query_end']}] -> "
                      f"{ref_name}[{aln['target_start']}-{aln['target_end']}] "
                      f"({aln['alignment_length']}bp, {(aln['matches']/aln['alignment_length']*100):.1f}% identity)")
        
        os.unlink(paf_file)
        
        # Analyze terminal regions
        terminal_hits = extract_and_compare_terminals(query_path, ref_path, query_name, ref_name)
        
        return {
            'query_coverage': coverage_stats['query_coverage'],
            'target_coverage': coverage_stats['target_coverage'],
            'identity': coverage_stats['identity'],
            'alignments': coverage_stats['alignments'],
            'terminal_hits': len(terminal_hits)
        }
    else:
        return None

def analyze_all_assembly_similarities(sample):
    """Compare all available assemblies pairwise"""
    assemblies = get_assembly_paths(sample)
    
    if len(assemblies) < 2:
        print(f"Error: Need at least 2 assemblies for comparison. Found: {list(assemblies.keys())}")
        return
    
    print(f"=== ASSEMBLY SIMILARITY ANALYSIS FOR {sample} ===")
    print(f"Available assemblies: {', '.join(assemblies.keys())}")
    
    # Create comparison matrix
    results = {}
    assembly_names = list(assemblies.keys())
    
    for i, ref_name in enumerate(assembly_names):
        for query_name in assembly_names[i+1:]:
            comparison_key = f"{ref_name}_vs_{query_name}"
            ref_path = assemblies[ref_name]
            query_path = assemblies[query_name]
            
            result = compare_assembly_pair(ref_path, query_path, ref_name, query_name)
            results[comparison_key] = result
    
    # Summary table
    print(f"\n=== SIMILARITY SUMMARY ===")
    print(f"{'Comparison':<25} {'Query Cov':<12} {'Target Cov':<12} {'Identity':<10} {'Alignments':<12}")
    print("-" * 70)
    
    for comparison, result in results.items():
        if result:
            print(f"{comparison:<25} {result['query_coverage']:<12.1f} "
                  f"{result['target_coverage']:<12.1f} {result['identity']:<10.1f} "
                  f"{result['alignments']:<12}")
        else:
            print(f"{comparison:<25} {'Failed':<12} {'Failed':<12} {'Failed':<10} {'Failed':<12}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python assembly_similarity_analysis.py SAMPLE_NAME")
        print("Example: python assembly_similarity_analysis.py B006")
        sys.exit(1)
    
    sample = sys.argv[1]
    analyze_all_assembly_similarities(sample)

if __name__ == "__main__":
    main()