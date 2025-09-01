#!/usr/bin/env python3

import os
from Bio import SeqIO

# Sample list
samples = ["B006", "B021", "B032", "B044", "D1701", "S1-Japan"]

# Base path
base_path = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly/03_assembly"

for sample in samples:
    input_file = f"{base_path}/{sample}/canu_out/{sample}.contigs.fasta"
    output_file = f"{base_path}/{sample}/canu_out/{sample}.longest_contig.fasta"
    
    if os.path.exists(input_file):
        try:
            # Read all sequences
            records = list(SeqIO.parse(input_file, "fasta"))
            
            if records:
                # Find longest sequence
                longest = max(records, key=lambda x: len(x.seq))
                
                # Write longest to new file
                SeqIO.write(longest, output_file, "fasta")
                
                print(f"{sample}: {len(longest.seq)} bp contig saved to {output_file}")
            else:
                print(f"{sample}: No sequences found in file")
                
        except Exception as e:
            print(f"{sample}: Error processing file - {e}")
    else:
        print(f"{sample}: File not found - {input_file}")

print("Done!")