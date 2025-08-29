#!/bin/bash

# Alternative Annotation Script - Creates DNA files with ORF annotations
# Usage: ./05_annotation_alternative.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SAMPLE=$1
ASSEMBLY_TYPE=${2:-"miniasm_viral"}

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    exit 1
fi

# Define assembly file paths (same as your original script)
# Define assembly file paths based on type
case $ASSEMBLY_TYPE in
    "miniasm_viral")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/viral_verification/Prediction_results_fasta/polished_assembly_virus.fasta"
        ;;
    "miniasm_raw")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
        ;;
    "canu")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/canu_ultra_output/${SAMPLE}.longest_contig.fasta"
        ;;
    "canu_ultra_trimmed")
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_1000.fasta" #TODO fix this
        ;;
    "canu_super")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "microsynth")
        ASSEMBLY_FILE="${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"
        ;;
    *)
        echo "Error: Unknown assembly type '$ASSEMBLY_TYPE'"
        echo "Valid types: miniasm_viral, miniasm_raw, canu, microsynth"
        exit 1
        ;;
esac

mkdir -p "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
cd "${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"

if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE with assembly type: $ASSEMBLY_TYPE"

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV_PROKKA

# Option A: Use Artemis/DNAPlotter approach
echo "=== OPTION A: Creating Artemis-compatible files ==="

# 1. Find ORFs using Prodigal (more control than Prokka)
echo "Finding ORFs with Prodigal..."
prodigal -i "$ASSEMBLY_FILE" -o "${SAMPLE}_genes.gff" -a "${SAMPLE}_proteins.faa" -d "${SAMPLE}_genes.fna" -f gff -p single

# 2. BLAST against your ORF database to get mappings
if [ -f "$ORF_PROTEIN_FILE" ]; then
    echo "Mapping genes to reference ORFs..."
    makeblastdb -in "$ORF_PROTEIN_FILE" -dbtype prot -out orf_db -logfile /dev/null
    blastp -query "${SAMPLE}_proteins.faa" -db orf_db -outfmt "6 qseqid sseqid pident evalue" -max_target_seqs 1 -evalue 1e-3 > blast_mapping.txt
    
    # Create high-confidence mapping
    awk '$3 >= 70 {print $1"\t"$2}' blast_mapping.txt > gene_orf_map.txt
    
    # 3. Create custom GFF with ORF names
    echo "Creating ORF-named annotations..."
    awk -F'\t' '
    BEGIN {
        while ((getline line < "gene_orf_map.txt") > 0) {
            split(line, fields, "\t")
            # Extract gene number from prodigal ID (e.g., "1_1" -> "1")
            if (match(fields[1], /_([0-9]+)$/, gene_num)) {
                map[gene_num[1]] = fields[2]
            }
        }
        close("gene_orf_map.txt")
    }
    /^#/ { print; next }
    $3 == "CDS" {
        # Extract gene number from prodigal GFF ID
        if (match($9, /ID=([0-9]+)_/, gene_id)) {
            gene_num = gene_id[1]
            if (gene_num in map) {
                orf_name = map[gene_num]
                # Replace the entire attributes field with ORF info
                $9 = "ID=" orf_name ";Name=" orf_name ";gene=" orf_name ";product=" orf_name
            }
        }
        print
    }
    $3 != "CDS" { print }
    ' "${SAMPLE}_genes.gff" > "${SAMPLE}_orf_annotated.gff"
    
    echo "Created: ${SAMPLE}_orf_annotated.gff"
    
    # 4. Convert to Artemis format (.art) which can be opened in DNA viewers
    if command -v artemis &> /dev/null; then
        echo "Converting to Artemis format..."
        # Artemis can read GFF directly, just copy the assembly
        cp "$ASSEMBLY_FILE" "${SAMPLE}_with_orfs.fasta"
        echo "Use: artemis ${SAMPLE}_with_orfs.fasta +${SAMPLE}_orf_annotated.gff"
    fi
    
    rm -f orf_db.* blast_mapping.txt gene_orf_map.txt
fi

# Option B: Create SnapGene-style DNA file using BioPython
echo ""
echo "=== OPTION B: Creating SnapGene-compatible DNA file ==="

cat > create_dna_file.py << 'EOF'
#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def create_annotated_dna_file(fasta_file, gff_file, output_file, sample_name):
    """Create an annotated DNA file from FASTA and GFF"""
    
    # Read the sequence
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if not sequences:
        print("Error: No sequences found in FASTA file")
        return False
    
    # Use the first/longest sequence
    seq_record = sequences[0]
    if len(sequences) > 1:
        # Find longest sequence
        seq_record = max(sequences, key=len)
        print(f"Using longest sequence: {len(seq_record.seq)} bp")
    
    # Clear existing features and add new ones from GFF
    seq_record.features = []
    seq_record.id = sample_name
    seq_record.name = sample_name
    seq_record.description = f"{sample_name} with ORF annotations"
    
    # Read GFF and add features
    if os.path.exists(gff_file):
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'CDS':
                    start = int(fields[3]) - 1  # Convert to 0-based
                    end = int(fields[4])
                    strand = 1 if fields[6] == '+' else -1
                    
                    # Parse attributes
                    attrs = {}
                    for attr in fields[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attrs[key] = value
                    
                    # Create feature
                    feature = SeqFeature(
                        location=FeatureLocation(start, end, strand=strand),
                        type='CDS',
                        qualifiers={
                            'gene': [attrs.get('gene', attrs.get('Name', 'unknown'))],
                            'product': [attrs.get('product', attrs.get('Name', 'hypothetical protein'))],
                            'locus_tag': [attrs.get('gene', attrs.get('Name', 'unknown'))]
                        }
                    )
                    seq_record.features.append(feature)
    
    print(f"Added {len(seq_record.features)} features")
    
    # Save as GenBank format (can be imported into DNA viewers)
    try:
        SeqIO.write(seq_record, output_file, "genbank")
        print(f"Created annotated file: {output_file}")
        return True
    except Exception as e:
        print(f"Error writing file: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 create_dna_file.py <fasta> <gff> <output> <sample_name>")
        sys.exit(1)
    
    success = create_annotated_dna_file(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    sys.exit(0 if success else 1)
EOF

# Run the DNA file creation
if command -v python3 &> /dev/null && python3 -c "import Bio" 2>/dev/null; then
    if [ -f "${SAMPLE}_orf_annotated.gff" ]; then
        echo "Creating SnapGene-compatible file..."
        python3 create_dna_file.py "$ASSEMBLY_FILE" "${SAMPLE}_orf_annotated.gff" "${SAMPLE}_annotated.gbk" "$SAMPLE"
        
        # Also create a simple tab-delimited feature file for other DNA viewers
        echo -e "Feature\tStart\tEnd\tStrand\tName\tType" > "${SAMPLE}_features.txt"
        awk -F'\t' '
        $3 == "CDS" {
            match($9, /gene=([^;]+)/, gene)
            name = gene[1] ? gene[1] : "unknown"
            strand = ($6 == "+") ? "forward" : "reverse"
            print "CDS\t" $4 "\t" $5 "\t" strand "\t" name "\tCDS"
        }' "${SAMPLE}_orf_annotated.gff" >> "${SAMPLE}_features.txt"
        
        echo "Created feature table: ${SAMPLE}_features.txt"
    fi
    
    rm -f create_dna_file.py
else
    echo "BioPython not available. Install with: conda install -c conda-forge biopython"
fi

# Option C: Simple custom format for easy parsing
echo ""
echo "=== OPTION C: Simple custom annotation format ==="

if [ -f "${SAMPLE}_orf_annotated.gff" ]; then
    # Create a simple DNA annotation file that's easy to parse
    echo "Creating simple annotation format..."
    
    # Create header
    cat > "${SAMPLE}_simple_annotations.txt" << HEADER
# Simple DNA Annotation File
# Sample: $SAMPLE
# Assembly: $ASSEMBLY_TYPE
# Format: ORF_Name	Start	End	Strand	Length	Product
HEADER
    
    # Add annotations
    awk -F'\t' '
    $3 == "CDS" {
        match($9, /gene=([^;]+)/, gene)
        match($9, /product=([^;]+)/, prod)
        name = gene[1] ? gene[1] : "unknown"
        product = prod[1] ? prod[1] : name
        length = $5 - $4 + 1
        strand = ($6 == "+") ? "+" : "-"
        print name "\t" $4 "\t" $5 "\t" strand "\t" length "\t" product
    }' "${SAMPLE}_orf_annotated.gff" >> "${SAMPLE}_simple_annotations.txt"
    
    echo "Created: ${SAMPLE}_simple_annotations.txt"
fi

echo ""
echo "=== SUMMARY ==="
echo "Created annotation files:"
if [ -f "${SAMPLE}_orf_annotated.gff" ]; then
    echo "  ✓ ${SAMPLE}_orf_annotated.gff (GFF3 with ORF names)"
fi
if [ -f "${SAMPLE}_annotated.gbk" ]; then
    echo "  ✓ ${SAMPLE}_annotated.gbk (GenBank format - can import to SnapGene/Benchling)"
fi
if [ -f "${SAMPLE}_features.txt" ]; then
    echo "  ✓ ${SAMPLE}_features.txt (Tab-delimited features)"
fi
if [ -f "${SAMPLE}_simple_annotations.txt" ]; then
    echo "  ✓ ${SAMPLE}_simple_annotations.txt (Simple format)"
fi

echo ""
echo "To use these files:"
echo "  - Import .gbk into SnapGene, Benchling, or other DNA viewers"
echo "  - Use .gff with Artemis, IGV, or other genome browsers"
echo "  - Parse simple_annotations.txt for custom analysis"

echo ""
echo "Alternative annotation pipeline completed!"