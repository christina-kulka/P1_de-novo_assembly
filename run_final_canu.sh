#!/bin/bash

# Final Canu Assembly Pipeline
# Usage: ./run_final_canu.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Check if sample name is provided
SAMPLE=$1
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    echo ""
    echo "This script runs the complete final Canu assembly pipeline:"
    echo "1. Quality control"
    echo "2. Data preprocessing"
    echo "3. Canu assembly"
    echo "4. ORF annotation (canu_ultra)"
    echo "5. ITR and ORF analysis"
    echo "6. Assembly trimming"
    exit 1
fi

echo "=== FINAL CANU ASSEMBLY PIPELINE ==="
echo "Sample: $SAMPLE"
echo "Assembly type: canu_ultra"
echo ""

# Function to run a step and check for errors
run_step() {
    local script=$1
    local sample=$2
    local assembly_type=$3
    local step_name=$4
    
    echo "Running: $step_name"
    echo "Command: ./$script $sample $assembly_type"
    
    if [ -n "$assembly_type" ]; then
        ./"$script" "$sample" "$assembly_type"
    else
        ./"$script" "$sample"
    fi
    
    if [ $? -ne 0 ]; then
        echo "ERROR: $step_name failed for sample $sample"
        exit 1
    fi
    
    echo "$step_name completed successfully"
    echo ""
}

# Step 1: Quality Control
run_step "01_quality_control.sh" "$SAMPLE" "" "Quality Control"

# Step 2: Data Preprocessing
run_step "02_preprocessing.sh" "$SAMPLE" "" "Data Preprocessing"

# Step 3: Canu Assembly
run_step "03_assembly_canu.sh" "$SAMPLE" "" "Canu Assembly"

# Step 4: ORF Annotation
run_step "04_orf_annotation.sh" "$SAMPLE" "canu_ultra" "ORF Annotation"

# Step 5: ITR and ORF Analysis
run_step "05_itr_orf_eval.sh" "$SAMPLE" "canu_ultra" "ITR and ORF Analysis"

# Step 6: Assembly Trimming
run_step "06_assembly_trimmer.sh" "$SAMPLE" "canu_ultra" "Assembly Trimming"

echo "=== FINAL CANU PIPELINE COMPLETED ==="
echo "Sample: $SAMPLE"
echo ""
echo "Output locations:"
echo "- Quality control: ${QC_OUTPUT_DIR}/${SAMPLE}/"
echo "- Preprocessed reads: ${PREPROC_OUTPUT_DIR}/${SAMPLE}/"
echo "- Assembly: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/super_canu_out/"
echo "- Annotations: ${ANNOTATION_OUTPUT_DIR}/${SAMPLE}/canu_ultra/"
echo "- ITR analysis: ${ITR_OUTPUT_DIR}/${SAMPLE}/canu_ultra/"
echo "- Final trimmed genome: ${FINAL_ASSEMBLY_DIR}/${SAMPLE}/canu_ultra/"
echo ""
echo "Pipeline completed successfully for $SAMPLE"
