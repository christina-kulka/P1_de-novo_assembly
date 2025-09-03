#!/usr/bin/env python3

"""
Merge existing annotations with new ORF annotations from Prokka output.
Preserves all original features and adds only ORF-related CDS features.
Usage: python3 merge_annotations.py ORIGINAL_FILE PROKKA_GBK_FILE OUTPUT_FILE
"""

import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse

def read_original_file(original_file):
    """Read the original annotation file with colleague's work."""
    try:
        # Try different formats
        if original_file.lower().endswith('.dna'):
            try:
                record = SeqIO.read(original_file, 'snapgene')
            except:
                record = SeqIO.read(original_file, 'genbank')
        elif original_file.lower().endswith(('.gbk', '.gb', '.genbank')):
            record = SeqIO.read(original_file, 'genbank')
        else:
            record = SeqIO.read(original_file, 'fasta')
        
        print(f"Read original file: {len(record.seq)} bp, {len(record.features)} features")
        return record
    except Exception as e:
        print(f"Error reading original file: {e}")
        return None

def read_prokka_orfs(prokka_file, orf_names_to_include=None):
    """Read ORF annotations from Prokka output."""
    try:
        prokka_record = SeqIO.read(prokka_file, 'genbank')
        orf_features = []
        
        for feature in prokka_record.features:
            if feature.type == "CDS":
                # Get product name
                product = feature.qualifiers.get('product', [''])[0]
                gene = feature.qualifiers.get('gene', [''])[0]
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                
                # Check if this is an ORF we want to include
                is_orf = (
                    product.startswith('ORF_') or 
                    gene.startswith('ORF_') or
                    'ORF_' in product or
                    'protein' in product.lower()
                )
                
                if is_orf:
                    # If specific ORF names provided, filter by them
                    if orf_names_to_include:
                        include_this = any(orf_name in str(feature.qualifiers) for orf_name in orf_names_to_include)
                        if not include_this:
                            continue
                    
                    orf_features.append(feature)
                    print(f"  Found ORF: {product} at {feature.location}")
        
        print(f"Extracted {len(orf_features)} ORF features from Prokka output")
        return orf_features
    
    except Exception as e:
        print(f"Error reading Prokka file: {e}")
        return []

def check_feature_overlap(new_feature, existing_features, min_overlap=50):
    """Check if new feature significantly overlaps with existing features."""
    new_start = int(new_feature.location.start)
    new_end = int(new_feature.location.end)
    
    for existing in existing_features:
        if existing.type in ["CDS", "gene"]:
            existing_start = int(existing.location.start)
            existing_end = int(existing.location.end)
            
            # Calculate overlap
            overlap_start = max(new_start, existing_start)
            overlap_end = min(new_end, existing_end)
            overlap_length = max(0, overlap_end - overlap_start)
            
            if overlap_length >= min_overlap:
                return True, existing
    
    return False, None

def merge_annotations(original_record, orf_features, avoid_overlap=True):
    """Merge ORF features into original record."""
    added_count = 0
    skipped_count = 0
    
    for orf_feature in orf_features:
        if avoid_overlap:
            overlaps, existing_feature = check_feature_overlap(orf_feature, original_record.features)
            if overlaps:
                print(f"  Skipping {orf_feature.qualifiers.get('product', ['unknown'])[0]} - overlaps with existing feature")
                skipped_count += 1
                continue
        
        # Add the ORF feature
        original_record.features.append(orf_feature)
        added_count += 1
    
    print(f"Added {added_count} ORF features, skipped {skipped_count} due to overlaps")
    return original_record

def main():
    parser = argparse.ArgumentParser(
        description='Merge existing annotations with ORF annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic merge
  python3 merge_annotations.py original.dna prokka_output.gbk merged_output.gbk
  
  # Allow overlaps (less conservative)
  python3 merge_annotations.py original.dna prokka_output.gbk merged_output.gbk --allow-overlaps
        """
    )
    
    parser.add_argument('original_file', help='Original annotation file (.dna, .gbk)')
    parser.add_argument('prokka_file', help='Prokka output with ORF annotations (.gbk)')
    parser.add_argument('output_file', help='Output merged file (.gbk)')
    parser.add_argument('--allow-overlaps', action='store_true',
                       help='Allow ORF features to overlap existing features')
    parser.add_argument('--orf-list', help='Optional file with list of ORF names to include')
    
    args = parser.parse_args()
    
    # Check input files exist
    for file_path in [args.original_file, args.prokka_file]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
    
    # Read ORF names if provided
    orf_names = None
    if args.orf_list and os.path.exists(args.orf_list):
        with open(args.orf_list, 'r') as f:
            orf_names = [line.strip() for line in f if line.strip()]
        print(f"Filtering for {len(orf_names)} specific ORFs")
    
    # Read original file
    print("Reading original annotation file...")
    original_record = read_original_file(args.original_file)
    if not original_record:
        sys.exit(1)
    
    # Read ORF features from Prokka
    print("Reading ORF features from Prokka output...")
    orf_features = read_prokka_orfs(args.prokka_file, orf_names)
    if not orf_features:
        print("No ORF features found in Prokka output")
        sys.exit(1)
    
    # Merge annotations
    print("Merging annotations...")
    merged_record = merge_annotations(
        original_record, 
        orf_features, 
        avoid_overlap=not args.allow_overlaps
    )
    
    # Write output
    try:
        SeqIO.write(merged_record, args.output_file, 'genbank')
        print(f"Wrote merged annotations: {args.output_file}")
        print(f"Final feature count: {len(merged_record.features)}")
    except Exception as e:
        print(f"Error writing output: {e}")

if __name__ == "__main__":
    main()