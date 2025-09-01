#!/usr/bin/env python3

"""
Format protein FASTA file for Prokka compatibility.
Usage: python3 format_proteins.py INPUT_FILE OUTPUT_FILE
"""

import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 format_proteins.py INPUT_FILE OUTPUT_FILE")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    annotation_path = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/04_annotation/D1701/canu_ultra_trimmed/"

    print(f"Formatting {input_file} -> {output_file}")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # Extract ORF name (remove > symbol)
                orf_name = line[1:].strip()
                
                # Format for Prokka: >GENE_NAME PRODUCT_DESCRIPTION
                if orf_name.startswith('ORF_'):
                    # Keep the ORF name as the gene name and add a description
                    outfile.write(f'>{orf_name} {orf_name}_protein\n')
                else:
                    # For other names like AcGFP, use as both gene and product
                    outfile.write(f'>{orf_name} {orf_name}_protein\n')
            elif line and not line.startswith('#'):
                # Remove trailing stop codon if present and write sequence
                sequence = line.rstrip('*')
                if sequence:  # Only write non-empty sequences
                    outfile.write(sequence + '\n')
    
    print(f"Formatted protein file created: {output_file}")

if __name__ == "__main__":
    main()