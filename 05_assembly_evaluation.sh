#!/bin/bash

# Assembly Evaluation Pipeline (Step 5)
# Execute all evaluation scripts in sequence
# Usage: ./05_assembly_evaluation.sh SAMPLE_NAME

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Check if sample name is provided
SAMPLE=$1
if [ -z "$SAMPLE" ]; then
    echo "Usage: $0 SAMPLE_NAME"
    echo "Example: $0 B006"
    exit 1
fi

# Create output directory for evaluation results
EVAL_OUTPUT_DIR="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/evaluation_results"
mkdir -p "$EVAL_OUTPUT_DIR"

# Set up log file
LOG_FILE="${EVAL_OUTPUT_DIR}/assembly_evaluation_${SAMPLE}_$(date +%Y%m%d_%H%M%S).txt"

# Function to log and display output
log_and_display() {
    echo "$1" | tee -a "$LOG_FILE"
}

# Start logging
log_and_display "=========================================="
log_and_display "ASSEMBLY EVALUATION PIPELINE - STEP 5"
log_and_display "Sample: $SAMPLE"
log_and_display "Started: $(date)"
log_and_display "Log file: $LOG_FILE"
log_and_display "=========================================="
log_and_display ""

# Check if assembly exists
ASSEMBLY_FILE="${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/polished_assembly.fasta"
if [ ! -f "$ASSEMBLY_FILE" ]; then
    log_and_display "❌ Error: Assembly file not found: $ASSEMBLY_FILE"
    log_and_display "   Run the assembly pipeline first (steps 1-4)!"
    exit 1
fi

log_and_display "✅ Assembly file found: $ASSEMBLY_FILE"
log_and_display ""

# Function to run evaluation script and capture output
run_evaluation() {
    local script_name=$1
    local description=$2
    
    log_and_display "=========================================="
    log_and_display "Running: $script_name"
    log_and_display "Description: $description"
    log_and_display "Started: $(date)"
    log_and_display "=========================================="
    log_and_display ""
    
    if [ -f "${SCRIPT_DIR}/${script_name}" ]; then
        # Run the script and capture all output
        log_and_display "Executing ${script_name}..."
        log_and_display "----------------------------------------"
        
        # Run script and capture output to log file
        bash "${SCRIPT_DIR}/${script_name}" "$SAMPLE" 2>&1 | tee -a "$LOG_FILE"
        
        # Check exit status
        if [ ${PIPESTATUS[0]} -eq 0 ]; then
            log_and_display ""
            log_and_display "✅ $script_name completed successfully"
        else
            log_and_display ""
            log_and_display "❌ $script_name failed"
        fi
    else
        log_and_display "❌ Error: Script not found: ${SCRIPT_DIR}/${script_name}"
    fi
    
    log_and_display ""
    log_and_display "Completed: $(date)"
    log_and_display "----------------------------------------"
    log_and_display ""
    log_and_display "Press Enter to continue to next evaluation..."
    read -r
    log_and_display ""
}

# Run all evaluation scripts in sequence
log_and_display "Starting assembly evaluation pipeline..."
log_and_display ""

# 1. ITR Check
run_evaluation "05a_itr_check.sh" "Check for mirror-image sequences (ITRs) in ORFV genome"

# 2. Circularity Check
run_evaluation "05b_circularity_check.sh" "Check for circular genome structure"

# 3. Assembly Comparison
run_evaluation "05c_assembly_comp.sh" "Compare your assembly vs Microsynth assembly"

# 4. Alignment Visualization
run_evaluation "05d_visualize_alignment.sh" "Visualize assembly alignment and coverage"

# 5. BLAST Extra Sequence
run_evaluation "05e_blast_extra.sh" "BLAST analysis of extra 66kb sequence"

log_and_display "=========================================="
log_and_display "ASSEMBLY EVALUATION PIPELINE COMPLETED!"
log_and_display "Completed: $(date)"
log_and_display "=========================================="
log_and_display ""
log_and_display "All evaluation results have been generated for sample: $SAMPLE"
log_and_display ""
log_and_display "Summary of what was analyzed:"
log_and_display "  ✅ ITR sequences (mirror-image repeats)"
log_and_display "  ✅ Circular genome structure"
log_and_display "  ✅ Assembly comparison with Microsynth"
log_and_display "  ✅ Alignment visualization and coverage"
log_and_display "  ✅ BLAST analysis of extra sequence"
log_and_display ""
log_and_display "Check the output files in: ${ASSEMBLY_OUTPUT_DIR}/${SAMPLE}/miniasm_out/"
log_and_display "Evaluation log saved in: $LOG_FILE"
log_and_display ""
log_and_display "Next steps:"
log_and_display "  - Review all evaluation results"
log_and_display "  - Determine if extra sequence is legitimate ORFV"
log_and_display "  - Consider ITR correction if needed (step 4)"
log_and_display "  - Finalize your ORFV genome assembly"
log_and_display ""
log_and_display "=========================================="
log_and_display "Evaluation pipeline completed successfully!"
log_and_display "All results have been saved to: $LOG_FILE"
log_and_display "==========================================" 