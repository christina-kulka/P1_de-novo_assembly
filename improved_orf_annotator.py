#!/usr/bin/env python3

"""
Improved ORF Annotation Script
Annotates a genome with provided ORF sequences using direct sequence matching
and proper coordinate handling to avoid BLAST coordinate issues.

Usage: python3 improved_orf_annotator.py INPUT_GENOME ORF_FASTA OUTPUT_GBK --sample-name SAMPLE
"""

import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
import re
from pathlib import Path

def find_orf_in_genome(genome_seq, orf_seq, orf_name, min_identity=95):
    """
    Find ORF sequence in genome using direct sequence matching.
    Returns list of matches with coordinates and strand.
    """
    matches = []
    genome_str = str(genome_seq).upper()
    orf_str = str(orf_seq).upper()
    orf_rc = str(orf_seq.reverse_complement()).upper()
    
    # Search forward strand
    for match in re.finditer(re.escape(orf_str), genome_str):
        start = match.start() + 1  # Convert to 1-based coordinates
        end = match.end()
        
        matches.append({
            'orf_name': orf_name,
            'start': start,
            'end': end,
            'strand': 1,
            'identity': 100.0,
            'sequence': orf_str
        })
    
    # Search reverse strand
    for match in re.finditer(re.escape(orf_rc), genome_str):
        start = match.start() + 1  # Convert to 1-based coordinates  
        end = match.end()
        
        matches.append({
            'orf_name': orf_name,
            'start': start,
            'end': end,
            'strand': -1,
            'identity': 100.0,
            'sequence': orf_rc
        })
    
    # If no exact matches, try with some mismatches (optional)
    if not matches and min_identity < 100:
        print(f"  No exact match for {orf_name}, trying with mismatches...")
        matches = find_approximate_matches(genome_str, orf_str, orf_rc, orf_name, min_identity)
    
    return matches

def find_approximate_matches(genome_str, orf_str, orf_rc, orf_name, min_identity):
    """Find approximate matches allowing for some mismatches."""
    matches = []
    orf_len = len(orf_str)
    min_matches = int(orf_len * min_identity / 100)
    
    # This is a simplified approximate matching - for production use, consider using
    # more sophisticated algorithms like local alignment
    for i in range(len(genome_str) - orf_len + 1):
        # Forward strand
        genome_segment = genome_str[i:i + orf_len]
        matches_count = sum(1 for a, b in zip(orf_str, genome_segment) if a == b)
        
        if matches_count >= min_matches:
            identity = (matches_count / orf_len) * 100
            matches.append({
                'orf_name': orf_name,
                'start': i + 1,
                'end': i + orf_len,
                'strand': 1,
                'identity': identity,
                'sequence': genome_segment
            })
        
        # Reverse strand
        matches_count = sum(1 for a, b in zip(orf_rc, genome_segment) if a == b)
        
        if matches_count >= min_matches:
            identity = (matches_count / orf_len) * 100
            matches.append({
                'orf_name': orf_name,
                'start': i + 1,
                'end': i + orf_len,
                'strand': -1,
                'identity': identity,
                'sequence': genome_segment
            })
    
    return matches

def translate_dna_sequence(dna_seq, strand=1):
    """
    Translate DNA sequence to protein, handling strand direction.
    """
    if strand == -1:
        dna_seq = dna_seq.reverse_complement()
    
    # Try all three reading frames and find the one with fewest stop codons
    best_translation = None
    best_stop_count = float('inf')
    best_frame = 0
    
    for frame in range(3):
        try:
            subseq = dna_seq[frame:]
            # Make length divisible by 3
            subseq = subseq[:len(subseq) - (len(subseq) % 3)]
            translation = subseq.translate(to_stop=False)
            stop_count = str(translation).count('*')
            
            if stop_count < best_stop_count:
                best_translation = translation
                best_stop_count = stop_count
                best_frame = frame
        except:
            continue
    
    # Remove terminal stop codon if present
    if best_translation and str(best_translation).endswith('*'):
        best_translation = best_translation[:-1]
    
    return str(best_translation), best_frame

def create_cds_feature(match, genome_seq, sample_name):
    """Create a CDS feature from a match."""
    start = match['start'] - 1  # Convert to 0-based for BioPython
    end = match['end']
    strand = match['strand']
    orf_name = match['orf_name']
    
    # Extract the sequence for translation
    if strand == 1:
        orf_sequence = genome_seq[start:end]
    else:
        orf_sequence = genome_seq[start:end].reverse_complement()
    
    # Translate the sequence
    protein_seq, frame_used = translate_dna_sequence(orf_sequence, 1)  # Already handled strand above
    
    # Create feature location
    feature_location = FeatureLocation(start, end, strand=strand)
    
    # Create qualifiers
    qualifiers = {
        'gene': [orf_name],
        'product': [orf_name],
        'locus_tag': [f"{sample_name}_{orf_name}"],
        'note': [f"Identity: {match['identity']:.1f}%"],
        'translation': [protein_seq]
    }
    
    # Create CDS feature
    cds_feature = SeqFeature(
        location=feature_location,
        type="CDS",
        qualifiers=qualifiers
    )
    
    return cds_feature

def read_orf_sequences(orf_fasta):
    """Read ORF sequences from FASTA file."""
    orfs = {}
    try:
        for record in SeqIO.parse(orf_fasta, 'fasta'):
            orfs[record.id] = record.seq
        print(f"Read {len(orfs)} ORF sequences from {orf_fasta}")
        return orfs
    except Exception as e:
        print(f"Error reading ORF file: {e}")
        return {}

def load_genome_file(input_file):
    """Load genome from various file formats."""
    try:
        # Try different formats in order of preference
        if input_file.lower().endswith(('.fasta', '.fa', '.fas')):
            record = SeqIO.read(input_file, 'fasta')
            print(f"Loaded FASTA file: {len(record.seq)} bp")
            # FASTA files don't have features by default
            record.features = []
        elif input_file.lower().endswith('.dna'):
            try:
                record = SeqIO.read(input_file, 'snapgene')
                print(f"Loaded SnapGene file: {len(record.seq)} bp, {len(record.features)} existing features")
            except:
                record = SeqIO.read(input_file, 'genbank')
                print(f"Loaded GenBank file: {len(record.seq)} bp, {len(record.features)} existing features")
        elif input_file.lower().endswith(('.gbk', '.gb', '.genbank')):
            record = SeqIO.read(input_file, 'genbank')
            print(f"Loaded GenBank file: {len(record.seq)} bp, {len(record.features)} existing features")
        else:
            # Try FASTA first, then GenBank for unknown extensions
            try:
                record = SeqIO.read(input_file, 'fasta')
                print(f"Auto-detected FASTA format: {len(record.seq)} bp")
                record.features = []
            except:
                record = SeqIO.read(input_file, 'genbank')
                print(f"Auto-detected GenBank format: {len(record.seq)} bp, {len(record.features)} existing features")
        
        return record
        
    except Exception as e:
        print(f"Error loading genome file: {e}")
        return None

def check_feature_overlap(new_feature, existing_features, min_overlap=50):
    """Check if new feature overlaps significantly with existing features."""
    new_start = int(new_feature.location.start)
    new_end = int(new_feature.location.end)
    
    for existing in existing_features:
        if existing.type == "CDS":
            existing_start = int(existing.location.start)
            existing_end = int(existing.location.end)
            
            overlap_start = max(new_start, existing_start)
            overlap_end = min(new_end, existing_end)
            overlap_length = max(0, overlap_end - overlap_start)
            
            if overlap_length >= min_overlap:
                existing_product = existing.qualifiers.get('product', ['unknown'])[0]
                return existing, existing_product, overlap_length
    
    return None, None, 0

def choose_best_orf_version(orf_sequences, orf_base_name):
    """
    Choose the best ORF version (F, R, or base) based on available options.
    Returns the ORF name and sequence to use.
    """
    # Possible variations
    candidates = [
        orf_base_name,           # ORF_029
        f"{orf_base_name}_F",    # ORF_029_F  
        f"{orf_base_name}_R",    # ORF_029_R
        f"{orf_base_name}F",     # ORF_029F
        f"{orf_base_name}R"      # ORF_029R
    ]
    
    available = []
    for candidate in candidates:
        if candidate in orf_sequences:
            available.append(candidate)
    
    if not available:
        return None, None
    
    # If only one version available, use it
    if len(available) == 1:
        return available[0], orf_sequences[available[0]]
    
    # If base name available, prefer it
    if orf_base_name in available:
        return orf_base_name, orf_sequences[orf_base_name]
    
    # Otherwise return the first available (F versions usually come first)
    return available[0], orf_sequences[available[0]]

def select_strand_appropriate_orf(orf_sequences, orf_base_name, strand):
    """
    Select the appropriate ORF version based on the strand it matches.
    """
    # Get all available versions
    f_versions = []
    r_versions = []
    base_version = None
    
    candidates = [orf_base_name, f"{orf_base_name}_F", f"{orf_base_name}_R", 
                  f"{orf_base_name}F", f"{orf_base_name}R"]
    
    for candidate in candidates:
        if candidate in orf_sequences:
            if candidate.endswith('_F') or candidate.endswith('F'):
                f_versions.append(candidate)
            elif candidate.endswith('_R') or candidate.endswith('R'):
                r_versions.append(candidate)
            else:
                base_version = candidate
    
    # Choose based on strand
    if strand == 1:  # Forward strand
        if f_versions:
            return f_versions[0], orf_sequences[f_versions[0]]
        elif base_version:
            return base_version, orf_sequences[base_version]
    else:  # Reverse strand
        if r_versions:
            return r_versions[0], orf_sequences[r_versions[0]]
        elif base_version:
            return base_version, orf_sequences[base_version]
    
    # Fallback - return any available version
    all_available = f_versions + r_versions + ([base_version] if base_version else [])
    if all_available:
        return all_available[0], orf_sequences[all_available[0]]
    
    return None, None

def annotate_genome(input_file, orf_fasta, output_file, sample_name, min_identity=95, replace_hypothetical=True):
    """Main annotation function."""
    
    # Load genome
    genome_record = load_genome_file(input_file)
    if not genome_record:
        return False
    
    # Load ORF sequences
    orf_sequences = read_orf_sequences(orf_fasta)
    if not orf_sequences:
        return False
    
    # Group ORFs by base name (remove _F, _R, F, R suffixes)
    orf_groups = {}
    for orf_name in orf_sequences.keys():
        # Extract base name
        base_name = orf_name
        for suffix in ['_F', '_R', 'F', 'R']:
            if base_name.endswith(suffix):
                base_name = base_name[:-len(suffix)]
                break
        
        if base_name not in orf_groups:
            orf_groups[base_name] = []
        orf_groups[base_name].append(orf_name)
    
    print(f"\nGrouped {len(orf_sequences)} ORFs into {len(orf_groups)} base groups")
    
    # Find matches for each base group
    all_matches = []
    processed_positions = set()  # Track positions to avoid duplicates
    
    print(f"\nSearching for ORFs in genome...")
    
    for base_name, orf_variants in orf_groups.items():
        print(f"  Searching for {base_name} variants: {orf_variants}")
        
        # Test all variants and collect matches
        variant_matches = []
        for orf_name in orf_variants:
            orf_seq = orf_sequences[orf_name]
            matches = find_orf_in_genome(genome_record.seq, orf_seq, orf_name, min_identity)
            for match in matches:
                variant_matches.append(match)
        
        if not variant_matches:
            print(f"    ✗ No matches found for {base_name}")
            continue
        
        # For each unique position, choose the best strand-appropriate variant
        position_matches = {}
        for match in variant_matches:
            pos_key = (match['start'], match['end'])
            if pos_key not in position_matches:
                position_matches[pos_key] = []
            position_matches[pos_key].append(match)
        
        for pos_key, matches_at_pos in position_matches.items():
            if pos_key in processed_positions:
                continue
            
            # Choose the best match for this position
            best_match = matches_at_pos[0]
            strand = best_match['strand']
            
            # Try to find strand-appropriate version
            chosen_name, chosen_seq = select_strand_appropriate_orf(orf_sequences, base_name, strand)
            if chosen_name:
                # Update match info with chosen version
                best_match['orf_name'] = chosen_name
                best_match['chosen_base'] = base_name
                all_matches.append(best_match)
                processed_positions.add(pos_key)
                print(f"    ✓ {chosen_name} at {best_match['start']}-{best_match['end']} (strand {strand:+d}, {best_match['identity']:.1f}% identity)")
            else:
                print(f"    ✗ Could not select appropriate version for {base_name}")
    
    if not all_matches:
        print("No ORF matches found!")
        return False
    
    print(f"\nSelected {len(all_matches)} best ORF matches (avoiding duplicates)")
    
    # Create CDS features
    new_features = []
    for match in all_matches:
        feature = create_cds_feature(match, genome_record.seq, sample_name)
        new_features.append(feature)
    
    # Handle overlaps with existing features
    added_count = 0
    replaced_count = 0
    skipped_count = 0
    
    for new_feature in new_features:
        existing_feature, existing_product, overlap = check_feature_overlap(
            new_feature, genome_record.features
        )
        
        new_product = new_feature.qualifiers.get('product', [''])[0]
        
        if existing_feature:
            if replace_hypothetical and 'hypothetical' in existing_product.lower():
                # Replace hypothetical protein with specific ORF
                print(f"  Replacing '{existing_product}' with '{new_product}'")
                genome_record.features.remove(existing_feature)
                genome_record.features.append(new_feature)
                replaced_count += 1
            else:
                print(f"  Skipping {new_product} - overlaps with {existing_product} ({overlap} bp)")
                skipped_count += 1
        else:
            genome_record.features.append(new_feature)
            added_count += 1
    
    print(f"\nFeature summary:")
    print(f"  Added: {added_count} new ORFs")
    print(f"  Replaced: {replaced_count} hypothetical proteins")
    print(f"  Skipped: {skipped_count} due to overlaps")
    print(f"  Total features: {len(genome_record.features)}")
    
    # Set required GenBank annotations if missing
    if not hasattr(genome_record, 'annotations'):
        genome_record.annotations = {}
    
    # Set required fields for GenBank format
    if 'molecule_type' not in genome_record.annotations:
        genome_record.annotations['molecule_type'] = 'DNA'
    
    if 'topology' not in genome_record.annotations:
        genome_record.annotations['topology'] = 'linear'
    
    if 'data_file_division' not in genome_record.annotations:
        genome_record.annotations['data_file_division'] = 'VRL'  # Viral division
    
    if 'date' not in genome_record.annotations:
        from datetime import datetime
        genome_record.annotations['date'] = datetime.now().strftime('%d-%b-%Y').upper()
    
    if 'source' not in genome_record.annotations:
        genome_record.annotations['source'] = f'{sample_name} genome'
    
    if 'organism' not in genome_record.annotations:
        genome_record.annotations['organism'] = 'Parapoxvirus'
    
    # Ensure record has an ID
    if not genome_record.id:
        genome_record.id = sample_name
    
    if not genome_record.name:
        genome_record.name = sample_name
    
    if not genome_record.description:
        genome_record.description = f'{sample_name} genome with ORF annotations'

    # Write output file
    try:
        SeqIO.write(genome_record, output_file, 'genbank')
        print(f"\nSuccessfully wrote annotated genome to: {output_file}")
        return True
    except Exception as e:
        print(f"Error writing output file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Annotate genome with ORF sequences using direct sequence matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 improved_orf_annotator.py genome.fasta orfs.fasta output.gbk --sample-name B006
  python3 improved_orf_annotator.py genome.dna orfs.fasta output.gbk --sample-name B006 --min-identity 98
        """
    )
    
    parser.add_argument('input_genome', help='Input genome file (.fasta, .dna, .gbk)')
    parser.add_argument('orf_fasta', help='FASTA file with ORF sequences to annotate')
    parser.add_argument('output_gbk', help='Output GenBank file with annotations')
    parser.add_argument('--sample-name', required=True, help='Sample name for locus tags')
    parser.add_argument('--min-identity', type=float, default=95.0,
                       help='Minimum percent identity for matches (default: 95.0)')
    parser.add_argument('--no-replace-hypothetical', action='store_true',
                       help='Do not replace hypothetical proteins')
    
    args = parser.parse_args()
    
    # Check input files exist
    if not os.path.exists(args.input_genome):
        print(f"Error: Input genome file not found: {args.input_genome}")
        sys.exit(1)
    
    if not os.path.exists(args.orf_fasta):
        print(f"Error: ORF FASTA file not found: {args.orf_fasta}")
        sys.exit(1)
    
    print(f"=== ORF GENOME ANNOTATION ===")
    print(f"Sample: {args.sample_name}")
    print(f"Input genome: {args.input_genome}")
    print(f"ORF sequences: {args.orf_fasta}")
    print(f"Output: {args.output_gbk}")
    print(f"Min identity: {args.min_identity}%")
    
    # Run annotation
    success = annotate_genome(
        args.input_genome,
        args.orf_fasta,
        args.output_gbk,
        args.sample_name,
        args.min_identity,
        not args.no_replace_hypothetical
    )
    
    if success:
        print("\n✓ Annotation completed successfully!")
    else:
        print("\n✗ Annotation failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()