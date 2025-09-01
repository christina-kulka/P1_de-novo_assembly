#!/usr/bin/env python3

"""
Standalone script to restore ORF names in Prokka output files.
Usage: python3 restore_orf_names.py SAMPLE_NAME ORIGINAL_ORF_FILE
"""

import sys
import re
import os
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 restore_orf_names.py SAMPLE_NAME ORIGINAL_ORF_FILE")
        print("Example: python3 restore_orf_names.py S1-Japan /path/to/S1-Japan_protein.fasta")
        sys.exit(1)
    
    sample_name = sys.argv[1]
    orf_file = sys.argv[2]

    annotation_path = f"/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/04_annotation/{sample_name}/canu_ultra_trimmed/"
    
    # Define expected Prokka output files
    prokka_faa = f"{annotation_path}{sample_name}.faa"
    prokka_gff = f"{annotation_path}{sample_name}.gff"
    prokka_gbk = f"{annotation_path}{sample_name}.gbk"

    # Check if required files exist
    required_files = [orf_file, prokka_faa, prokka_gff]
    missing_files = [f for f in required_files if not os.path.exists(f)]
    
    if missing_files:
        print(f"Error: Missing required files: {missing_files}")
        sys.exit(1)
    
    print(f"Processing sample: {sample_name}")
    print(f"Original ORF file: {orf_file}")
    
    # Read the original ORF protein file
    print("Reading original ORF sequences...")
    orf_sequences = {}
    with open(orf_file, 'r') as f:
        current_orf = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_orf = line[1:].strip()
                orf_sequences[current_orf] = ""
            elif current_orf and line:
                # Remove stop codon for comparison
                orf_sequences[current_orf] += line.rstrip('*')
    
    print(f"Found {len(orf_sequences)} ORF sequences")
    
    # Read Prokka protein sequences
    print("Reading Prokka sequences...")
    prokka_sequences = {}
    with open(prokka_faa, 'r') as f:
        current_locus = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Extract locus tag (first part before space)
                current_locus = line.split()[0][1:]  # Remove > and get first part
                prokka_sequences[current_locus] = ""
            elif current_locus and line:
                prokka_sequences[current_locus] += line
    
    print(f"Found {len(prokka_sequences)} Prokka sequences")
    
    # Create mapping by comparing sequences
    print("Creating sequence mapping...")
    locus_to_orf = {}
    for locus, seq in prokka_sequences.items():
        for orf_name, orf_seq in orf_sequences.items():
            if seq == orf_seq:
                locus_to_orf[locus] = orf_name
                print(f"  Matched {locus} -> {orf_name}")
                break
    
    print(f"Successfully mapped {len(locus_to_orf)} sequences")
    
    # Update GFF file with ORF names
    print("Updating GFF file...")
    updated_gff = f"{annotation_path}{sample_name}.orf_named.gff"
    with open(prokka_gff, 'r') as infile, open(updated_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or '\tCDS\t' not in line:
                outfile.write(line)
                continue
            
            # Extract locus tag from attributes
            match = re.search(r'ID=([^;]+)', line)
            if match:
                locus_tag = match.group(1)
                if locus_tag in locus_to_orf:
                    orf_name = locus_to_orf[locus_tag]
                    # Replace the product name in the attributes
                    line = re.sub(r'product=[^;]+', f'product={orf_name}', line)
                    # Add gene name if not present
                    if 'gene=' not in line:
                        line = re.sub(r'(ID=[^;]+)', f'\\1;gene={orf_name}', line)
            
            outfile.write(line)
    
    print(f"Created: {updated_gff}")
    
    # Update FAA file with ORF names
    print("Updating protein file...")
    updated_faa = f"{annotation_path}{sample_name}.orf_named.faa"
    with open(prokka_faa, 'r') as infile, open(updated_faa, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.split()
                locus_tag = parts[0][1:]  # Remove > and get first part
                if locus_tag in locus_to_orf:
                    orf_name = locus_to_orf[locus_tag]
                    outfile.write(f'>{locus_tag} {orf_name}\n')
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
    
    print(f"Created: {updated_faa}")
    
    # Update GenBank file with ORF names
    if os.path.exists(prokka_gbk):
        print("Updating GenBank file...")
        updated_gbk = f"{annotation_path}{sample_name}.orf_named.gbk"
        last_locus = None
        
        with open(prokka_gbk, 'r') as infile, open(updated_gbk, 'w') as outfile:
            for line in infile:
                # Look for locus_tag lines
                if '/locus_tag=' in line:
                    locus_match = re.search(r'/locus_tag="([^"]+)"', line)
                    if locus_match:
                        last_locus = locus_match.group(1)
                        outfile.write(line)
                        # Add gene name line after locus_tag if we have a mapping
                        if last_locus in locus_to_orf:
                            indent = len(line) - len(line.lstrip())
                            outfile.write(' ' * indent + f'/gene="{locus_to_orf[last_locus]}"\n')
                    else:
                        outfile.write(line)
                # Update product lines for mapped ORFs
                elif '/product=' in line and last_locus and last_locus in locus_to_orf:
                    orf_name = locus_to_orf[last_locus]
                    indent = len(line) - len(line.lstrip())
                    outfile.write(' ' * indent + f'/product="{orf_name}"\n')
                else:
                    outfile.write(line)
        
        print(f"Created: {updated_gbk}")
    else:
        print(f"Warning: {prokka_gbk} not found, skipping GenBank update")
    
    # Create mapping report
    print("Creating mapping report...")
    mapping_file = f"{sample_name}_orf_mapping.txt"
    with open(mapping_file, 'w') as f:
        f.write(f"# ORF Mapping Report for {sample_name}\n")
        f.write(f"# Original ORF file: {orf_file}\n")
        f.write(f"# Total ORFs in original file: {len(orf_sequences)}\n")
        f.write(f"# Total sequences from Prokka: {len(prokka_sequences)}\n")
        f.write(f"# Successfully mapped: {len(locus_to_orf)}\n")
        f.write(f"# Mapping success rate: {len(locus_to_orf)/len(prokka_sequences)*100:.1f}%\n")
        f.write(f"\nProkka_Locus_Tag\tORF_Name\tStatus\n")
        
        for locus in sorted(prokka_sequences.keys()):
            if locus in locus_to_orf:
                f.write(f"{locus}\t{locus_to_orf[locus]}\tMapped\n")
            else:
                f.write(f"{locus}\t-\tNo_match\n")
    
    print(f"Created: {mapping_file}")
    
    print("\n=== SUMMARY ===")
    print(f"✓ Created {len(locus_to_orf)} ORF name mappings")
    print(f"✓ Updated files:")
    print(f"  - {updated_gff} (GFF with ORF names)")
    print(f"  - {updated_faa} (proteins with ORF names)")
    if os.path.exists(f"{sample_name}.orf_named.gbk"):
        print(f"  - {sample_name}.orf_named.gbk (GenBank with ORF names)")
    print(f"  - {mapping_file} (mapping report)")
    
    print(f"\nMapping success: {len(locus_to_orf)}/{len(prokka_sequences)} sequences ({len(locus_to_orf)/len(prokka_sequences)*100:.1f}%)")

if __name__ == "__main__":
    main()