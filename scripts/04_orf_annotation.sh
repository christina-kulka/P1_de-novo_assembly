#!/bin/bash

# Manual ORF Annotation Script (Step 5) - Replaces Prokka
# Usage: ./05_annotation_manual.sh SAMPLE_NAME [ASSEMBLY_TYPE]

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

# Parse command line arguments
SAMPLE=$1
ASSEMBLY_TYPE=${2:-"miniasm_viral"}  # Default to miniasm_viral

if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME [ASSEMBLY_TYPE]"
    echo "Example: $0 B006"
    echo "Example: $0 B006 canu"
    echo "Assembly types: miniasm_viral, miniasm_raw, canu, canu_ultra, canu_super, microsynth"
    exit 1
fi

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
        ASSEMBLY_FILE="${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/${SAMPLE}_canu_ultra_trimmed_150.fasta"
        ;;
    "canu_super")
        ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/${SAMPLE}.longest_contig.fasta"
        ;;
    "microsynth")
        ASSEMBLY_FILE="${BASE_DIR}/00_raw_data_microsynth/${SAMPLE}_results/${SAMPLE}_results/Assembly/${SAMPLE}.fasta"
        ;;
    *)
        echo "Error: Unknown assembly type '$ASSEMBLY_TYPE'"
        echo "Valid types: miniasm_viral, miniasm_raw, canu, canu_ultra, canu_super, microsynth"
        exit 1
        ;;
esac

# Define ORF files - we need the DNA sequences, not protein
ORF_DNA_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"

# Create output directory
OUTPUT_DIR="${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/${ASSEMBLY_TYPE}"
mkdir -p "$OUTPUT_DIR"

# Validate input files
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: Assembly file not found for sample $SAMPLE, type $ASSEMBLY_TYPE"
    echo "Expected: $ASSEMBLY_FILE"
    exit 1
fi

if [ ! -f "$ORF_DNA_FILE" ]; then
    echo "Error: ORF DNA file not found: $ORF_DNA_FILE"
    echo "Expected: ${SAMPLE}_ORFs.fasta in ${ORF_DATABASE_DIR}/"
    echo ""
    echo "Note: This script requires ORF DNA sequences, not protein sequences."
    echo "Make sure you have ${SAMPLE}_ORFs.fasta with the DNA sequences of your ORFs."
    exit 1
fi

echo "=== MANUAL ORF ANNOTATION ==="
echo "Processing sample: $SAMPLE"
echo "Assembly type: $ASSEMBLY_TYPE"
echo "Assembly file: $ASSEMBLY_FILE"
echo "ORF DNA file: $ORF_DNA_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if improved annotator exists
if [ ! -f "${SCRIPT_DIR}/improved_orf_annotator.py" ]; then
    echo "Error: improved_orf_annotator.py not found in ${SCRIPT_DIR}"
    echo "Please make sure the improved annotation script is in the same directory."
    exit 1
fi

# Output files
ANNOTATED_GBK="${OUTPUT_DIR}/${SAMPLE}_orf_annotated.gbk"
ANNOTATION_LOG="${OUTPUT_DIR}/${SAMPLE}_annotation.log"

# Run our improved annotation script
echo "Running manual ORF annotation..."
$PYTHON3 "${SCRIPT_DIR}/improved_orf_annotator.py" \
    "$ASSEMBLY_FILE" \
    "$ORF_DNA_FILE" \
    "$ANNOTATED_GBK" \
    --sample-name "$SAMPLE" \
    --min-identity 95.0 \
    2>&1 | tee "$ANNOTATION_LOG"

# Check if annotation succeeded
if [ $? -ne 0 ]; then
    echo "Error: Manual annotation failed"
    exit 1
fi

if [ ! -f "$ANNOTATED_GBK" ]; then
    echo "Error: Output file was not created: $ANNOTATED_GBK"
    exit 1
fi

echo ""
echo "=== ANNOTATION COMPLETED ==="

# Generate additional output files for compatibility with downstream tools
cd "$OUTPUT_DIR"

echo "Generating additional output formats..."

# Extract protein sequences from GenBank file
PROTEIN_FASTA="${SAMPLE}_orf_annotated.faa"
echo "Extracting protein sequences to: $PROTEIN_FASTA"

$PYTHON3 -c "
from Bio import SeqIO
import sys

try:
    record = SeqIO.read('${SAMPLE}_orf_annotated.gbk', 'genbank')
    proteins = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
            product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
            translation = feature.qualifiers.get('translation', [''])[0]
            
            if translation:
                protein_id = f'{locus_tag}|{product}'
                proteins.append(f'>{protein_id}\n{translation}')
    
    with open('$PROTEIN_FASTA', 'w') as f:
        f.write('\n'.join(proteins))
    
    print(f'Extracted {len(proteins)} protein sequences')
    
except Exception as e:
    print(f'Error extracting proteins: {e}')
    sys.exit(1)
"

# Extract nucleotide CDS sequences
NUCLEOTIDE_FASTA="${SAMPLE}_orf_annotated.fna"
echo "Extracting nucleotide CDS sequences to: $NUCLEOTIDE_FASTA"

$PYTHON3 -c "
from Bio import SeqIO
import sys

try:
    record = SeqIO.read('${SAMPLE}_orf_annotated.gbk', 'genbank')
    cds_sequences = []
    
    for feature in record.features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers.get('locus_tag', ['unknown'])[0]
            product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
            
            # Extract the sequence
            cds_seq = feature.location.extract(record.seq)
            
            cds_id = f'{locus_tag}|{product}'
            cds_sequences.append(f'>{cds_id}\n{str(cds_seq)}')
    
    with open('$NUCLEOTIDE_FASTA', 'w') as f:
        f.write('\n'.join(cds_sequences))
    
    print(f'Extracted {len(cds_sequences)} CDS sequences')
    
except Exception as e:
    print(f'Error extracting CDS sequences: {e}')
    sys.exit(1)
"

# Create a simple GFF file
GFF_FILE="${SAMPLE}_orf_annotated.gff"
echo "Creating GFF file: $GFF_FILE"

$PYTHON3 -c "
from Bio import SeqIO
import sys

try:
    record = SeqIO.read('${SAMPLE}_orf_annotated.gbk', 'genbank')
    
    gff_lines = ['##gff-version 3']
    gff_lines.append(f'##sequence-region {record.id} 1 {len(record.seq)}')
    
    for i, feature in enumerate(record.features):
        if feature.type == 'CDS':
            start = int(feature.location.start) + 1  # GFF is 1-based
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'
            
            locus_tag = feature.qualifiers.get('locus_tag', [f'cds_{i+1}'])[0]
            product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
            gene = feature.qualifiers.get('gene', [locus_tag])[0]
            
            # GFF format: seqname source feature start end score strand frame attributes
            attributes = f'ID={locus_tag};Name={gene};product={product}'
            gff_line = f'{record.id}\tmanual\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}'
            gff_lines.append(gff_line)
    
    with open('$GFF_FILE', 'w') as f:
        f.write('\n'.join(gff_lines))
    
    print(f'Created GFF with {len([l for l in gff_lines if l.startswith(record.id)])} features')
    
except Exception as e:
    print(f'Error creating GFF: {e}')
    sys.exit(1)
"

# Generate summary statistics
echo ""
echo "=== ANNOTATION SUMMARY ==="

# Count features
TOTAL_FEATURES=$(grep -c "/gene=" "$ANNOTATED_GBK" 2>/dev/null || echo "0")
echo "Total ORF features annotated: $TOTAL_FEATURES"

# Show genome info
GENOME_SIZE=$(grep -v "^>" "$ASSEMBLY_FILE" | tr -d '\n' | wc -c)
echo "Genome size: $GENOME_SIZE bp"

# Calculate coding density
if [ "$TOTAL_FEATURES" -gt 0 ] && [ "$GENOME_SIZE" -gt 0 ]; then
    CODING_DENSITY=$(python3 -c "print(f'{($TOTAL_FEATURES / $GENOME_SIZE * 1000):.2f}')")
    echo "Gene density: $CODING_DENSITY genes per kb"
fi

echo ""
echo "Sample annotated ORFs:"
grep "/gene=" "$ANNOTATED_GBK" | head -10 | sed 's/.*\/gene="/  - /' | sed 's/".*//' | sort

echo ""
echo "=== OUTPUT FILES ==="
echo "Main files:"
echo "  - ${SAMPLE}_orf_annotated.gbk (GenBank format with annotations)"
echo "  - ${SAMPLE}_orf_annotated.gff (GFF format)"
echo "  - ${SAMPLE}_orf_annotated.faa (Protein sequences)"  
echo "  - ${SAMPLE}_orf_annotated.fna (Nucleotide CDS sequences)"
echo "  - ${SAMPLE}_annotation.log (Detailed annotation log)"

echo ""
echo "Files location: $OUTPUT_DIR"
echo ""
echo "Manual annotation pipeline completed successfully!"
echo ""
echo "Key advantages over Prokka:"
echo "  ✓ No hypothetical proteins added"
echo "  ✓ Uses your specific ORF sequences"  
echo "  ✓ Proper strand selection (F/R variants)"
echo "  ✓ No coordinate system issues"
echo "  ✓ Direct sequence matching for accuracy"