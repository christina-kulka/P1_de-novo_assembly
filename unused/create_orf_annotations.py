#!/usr/bin/env python3

"""
Create ORF annotations from BLAST results.
Takes BLAST output and creates GenBank annotations with proper ORF names.
Usage: python3 create_orf_annotations.py INPUT_FILE BLAST_RESULTS ORF_PROTEIN_FILE OUTPUT_FILE
"""

import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path

def read_blast_results(blast_file, min_identity=98, min_coverage=90):
    """Parse BLAST results and return ORF coordinates."""
    orf_coordinates = []  # Changed to list to allow multiple matches per ORF
    
    print(f"Reading BLAST results from {blast_file}")
    print(f"Filters: ≥{min_identity}% identity, ≥{min_coverage}% coverage")
    
    with open(blast_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 14:
                    orf_name = parts[0]       # query ID (ORF name)
                    subject_id = parts[1]     # subject ID (genome)
                    identity = float(parts[2]) # percent identity
                    length = int(parts[3])    # alignment length
                    q_start = int(parts[6])   # query start
                    q_end = int(parts[7])     # query end  
                    s_start = int(parts[8])   # subject start (genome coordinates)
                    s_end = int(parts[9])     # subject end
                    evalue = float(parts[10]) # E-value
                    bitscore = float(parts[11]) # bit score
                    q_len = int(parts[12])    # query length
                    s_len = int(parts[13])    # subject length
                    
                    # Calculate coverage
                    query_coverage = (abs(q_end - q_start) + 1) / q_len * 100
                    
                    # Strict filtering for complete matches
                    if identity < min_identity:
                        print(f"  Skipping {orf_name} (line {line_num}): {identity:.1f}% identity < {min_identity}%")
                        continue
                    
                    if query_coverage < min_coverage:
                        print(f"  Skipping {orf_name} (line {line_num}): {query_coverage:.1f}% coverage < {min_coverage}%")
                        continue
                    
                    # Determine strand and coordinates
                    if s_start < s_end:
                        strand = 1
                        start = s_start
                        end = s_end
                    else:
                        strand = -1
                        start = s_end
                        end = s_start
                    
                    # Accept all high-quality matches (allowing duplicates)
                    match_data = {
                        'orf_name': orf_name,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'identity': identity,
                        'length': length,
                        'coverage': query_coverage,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'match_id': f"{orf_name}_{start}_{end}"  # Unique ID for each match
                    }
                    
                    orf_coordinates.append(match_data)
                    print(f"  ✓ {orf_name}: {start}-{end} (strand {strand:+d}, {identity:.1f}% id, {query_coverage:.1f}% cov)")
    
    print(f"Found {len(orf_coordinates)} high-quality ORF matches")
    
    # Check for duplicates
    orf_counts = {}
    for match in orf_coordinates:
        orf_name = match['orf_name']
        orf_counts[orf_name] = orf_counts.get(orf_name, 0) + 1
    
    duplicates = {name: count for name, count in orf_counts.items() if count > 1}
    if duplicates:
        print(f"Found duplicate ORFs: {duplicates}")
    
    return orf_coordinates

def read_orf_proteins(protein_file):
    """Read ORF protein sequences."""
    proteins = {}
    
    try:
        for record in SeqIO.parse(protein_file, 'fasta'):
            # Remove stop codon if present
            protein_seq = str(record.seq).rstrip('*')
            proteins[record.id] = protein_seq
        
        print(f"Read {len(proteins)} protein sequences")
        return proteins
        
    except Exception as e:
        print(f"Error reading protein file: {e}")
        return {}

def create_orf_features(orf_matches, orf_proteins, sample_name):
    """Create CDS features for ORF matches."""
    features = []
    
    for match in orf_matches:
        orf_name = match['orf_name']
        
        # Create feature location
        #start_pos = match['start'] - 1  # Convert to 0-based
        start_pos = match['start']  # Keep as 1-based for Biopython
        end_pos = match['end']
        strand = match['strand']
        
        feature_location = FeatureLocation(start_pos, end_pos, strand=strand)
        
        # Get protein sequence
        protein_seq = orf_proteins.get(orf_name, "")
        
        # Create qualifiers with unique locus tag for each match
        qualifiers = {
            'gene': [orf_name],
            'product': [orf_name],
            'locus_tag': [f"{sample_name}_{match['match_id']}"],
            'note': [f"Identity: {match['identity']:.1f}%, Coverage: {match['coverage']:.1f}%"]
        }
        
        if protein_seq:
            qualifiers['translation'] = [protein_seq]
        
        # Create CDS feature
        cds_feature = SeqFeature(
            location=feature_location,
            type="CDS",
            qualifiers=qualifiers
        )
        
        features.append(cds_feature)
    
    print(f"Created {len(features)} CDS features")
    return features

def merge_with_existing(input_file, new_features, output_file):
    """Merge ORF features with existing annotations."""
    try:
        # Read original file
        if input_file.lower().endswith('.dna'):
            try:
                record = SeqIO.read(input_file, 'snapgene')
            except:
                record = SeqIO.read(input_file, 'genbank')
        elif input_file.lower().endswith(('.gbk', '.gb', '.genbank')):
            record = SeqIO.read(input_file, 'genbank')
        else:
            record = SeqIO.read(input_file, 'fasta')
        
        print(f"Original file: {len(record.seq)} bp, {len(record.features)} existing features")
        
        # Add new ORF features
        original_count = len(record.features)
        
        # Check for overlaps and handle them intelligently
        added_count = 0
        replaced_count = 0
        skipped_count = 0
        
        for new_feature in new_features:
            # Check for overlaps with existing CDS features
            overlaps = False
            existing_to_remove = None
            
            for existing in record.features:
                if existing.type == "CDS":
                    # Calculate overlap
                    new_start = int(new_feature.location.start)
                    new_end = int(new_feature.location.end)
                    existing_start = int(existing.location.start)
                    existing_end = int(existing.location.end)
                    
                    overlap_start = max(new_start, existing_start)
                    overlap_end = min(new_end, existing_end)
                    overlap_length = max(0, overlap_end - overlap_start)
                    
                    if overlap_length >= 50:  # Significant overlap
                        existing_product = existing.qualifiers.get('product', [''])[0]
                        new_product = new_feature.qualifiers.get('product', [''])[0]
                        
                        # Replace generic annotations with specific ORFs
                        if existing_product == 'hypothetical protein':
                            print(f"  Replacing '{existing_product}' with '{new_product}'")
                            existing_to_remove = existing
                            break
                        else:
                            print(f"  Skipping {new_product} - overlaps with {existing_product}")
                            overlaps = True
                            break
            
            if existing_to_remove:
                record.features.remove(existing_to_remove)
                record.features.append(new_feature)
                replaced_count += 1
            elif not overlaps:
                record.features.append(new_feature)
                added_count += 1
            else:
                skipped_count += 1
        
        print(f"Added {added_count} ORFs, replaced {replaced_count} generic features, skipped {skipped_count} due to overlaps")
        
        # Write output file
        SeqIO.write(record, output_file, 'genbank')
        print(f"Wrote annotated file: {output_file}")
        print(f"Final feature count: {len(record.features)} (was {original_count})")
        
        return True
        
    except Exception as e:
        print(f"Error processing files: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Create ORF annotations from BLAST results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 create_orf_annotations.py genome.dna blast_results.txt proteins.fasta output.gbk --sample-name B006
        """
    )
    
    parser.add_argument('input_file', help='Input genome file (.dna, .gbk, .fasta)')
    parser.add_argument('blast_results', help='BLAST results file (tabular format)')
    parser.add_argument('protein_file', help='ORF protein sequences (.fasta)')
    parser.add_argument('output_file', help='Output annotated file (.gbk)')
    parser.add_argument('--sample-name', required=True, help='Sample name for locus tags')
    parser.add_argument('--min-identity', type=float, default=95.0,
                       help='Minimum percent identity for BLAST matches (default: 95.0)')
    
    args = parser.parse_args()
    
    # Check input files
    required_files = [args.input_file, args.blast_results, args.protein_file]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
    
    print(f"Creating ORF annotations for sample: {args.sample_name}")
    print(f"Minimum identity threshold: {args.min_identity}%")
    
    # Read BLAST results
    orf_matches = read_blast_results(args.blast_results, args.min_identity, 90)  # 90% coverage required
    if not orf_matches:
        print("No valid ORF matches found")
        sys.exit(1)
    
    # Read protein sequences
    orf_proteins = read_orf_proteins(args.protein_file)
    
    # Create CDS features
    orf_features = create_orf_features(orf_matches, orf_proteins, args.sample_name)
    if not orf_features:
        print("No ORF features created")
        sys.exit(1)
    
    # Merge with existing annotations
    success = merge_with_existing(args.input_file, orf_features, args.output_file)
    
    if success:
        print("\nORF annotation completed successfully!")
        
        # Show summary of what was found
        print(f"\nAnnotated ORFs:")
        for match in sorted(orf_matches, key=lambda x: (x['orf_name'], x['start'])):
            print(f"  {match['orf_name']}: {match['start']}-{match['end']} ({match['identity']:.1f}% identity, {match['coverage']:.1f}% coverage)")
    else:
        print("ORF annotation failed")
        sys.exit(1)

if __name__ == "__main__":
    main()