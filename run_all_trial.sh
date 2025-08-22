#!/bin/bash

# Batch Analysis Script
# Usage: ./run_batch_analysis.sh

# Define your sample list here
SAMPLES=("B032" "B021" "B044" "S1-Japan" "D1701" "B006")

# Define which assembly types to process
ASSEMBLY_TYPES=("canu" "miniasm_viral" "microsynth")

# Define which analyses to run (comment out lines you don't want)
RUN_QC=true
RUN_PREPROCESSING=true
RUN_CANU=true
RUN_MINIASM=true
RUN_VIRAL_VERIFY=true
RUN_PROKKA=true
RUN_ITR_COMPARISON=true
RUN_ASSEMBLY_SIMILARITY=true

echo "=== BATCH ANALYSIS STARTING ==="
echo "Samples: ${SAMPLES[@]}"
echo "Assembly types: ${ASSEMBLY_TYPES[@]}"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"

    if [ "$RUN_QC" = true ]; then
        echo "  Running QC..."
        ./01_quality_control.sh "$SAMPLE"
    fi

    if [ "$RUN_PREPROCESSING" = true ]; then
        echo "  Running preprocessing..."
        ./02_preprocessing.sh "$SAMPLE"
    fi

    if [ "$RUN_CANU" = true ]; then
        echo "  Running Canu assembly..."
        ./03_assembly_canu.sh "$SAMPLE"
    fi

    if [ "$RUN_MINIASM" = true ]; then
        echo "  Running Miniasm assembly..."
        ./03_assembly_miniasm.sh "$SAMPLE"
    fi

    if [ "$RUN_VIRAL_VERIFY" = true ]; then
        echo "  Running viral verification..."
        ./03b_viralverify.sh "$SAMPLE"
    fi
    
    # Run Prokka annotation for each assembly type
    if [ "$RUN_PROKKA" = true ]; then
        echo "  Running Prokka annotations..."
        for ASSEMBLY_TYPE in "${ASSEMBLY_TYPES[@]}"; do
            echo "    - $ASSEMBLY_TYPE assembly"
            ./05_annotation_prokka.sh "$SAMPLE" "$ASSEMBLY_TYPE"
        done
    fi
    
    echo "  Completed: $SAMPLE"
    echo ""
done

echo "=== BATCH ANALYSIS COMPLETED ==="