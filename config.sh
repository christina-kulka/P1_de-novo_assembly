#!/bin/bash

# Configuration file for P1 de-novo assembly pipeline
# This file contains all paths and variables used across the repository

# Base directory structure
BASE_DIR="$HOME/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"

# Input data directories
RAW_DATA_DIR="${BASE_DIR}/00_raw_data_microsynth"

# Output directories
QC_OUTPUT_DIR="${BASE_DIR}/01_quality_control"
PREPROC_OUTPUT_DIR="${BASE_DIR}/02_preprocessing"
ASSEMBLY_OUTPUT_DIR="${BASE_DIR}/03_assembly"
ANNOTATION_OUTPUT_DIR="${BASE_DIR}/04_annotation"
ITR_OUTPUT_DIR="${BASE_DIR}/05_itr_analysis"
FINAL_ASSEMBLY_DIR="${BASE_DIR}/06_trimmed_assembly"

# Database paths
PFAM_DB="${BASE_DIR}/databases/Pfam-A.hmm.gz"
ORF_DATABASE_DIR="/home/ubuntu/data-volume/001_Raw_Data/Databases/ORF"
PROTEIN_DATABASE_DIR="/home/ubuntu/data-volume/001_Raw_Data/Databases/Proteins"
VIRAL_DB_DIR="${BASE_DIR}/databases/viral"
VIRAL_REF_FASTA="${VIRAL_DB_DIR}/parapoxvirus_references.fasta"

# Tool paths
PORECHOP="/home/ubuntu/miniconda3/bin/porechop"
NANOFILT="/home/ubuntu/miniconda3/envs/viralFlye/bin/NanoFilt"
TANDEMTOOLS="/home/ubuntu/TandemTools"
KRAKEN2="/home/ubuntu/miniconda3/bin/kraken2"
KRAKEN2_BUILD="/home/ubuntu/miniconda3/bin/kraken2-build"
BLASTN="/home/ubuntu/miniconda3/envs/itr_analysis/bin/blastn"
MAKEBLASTDB="/home/ubuntu/miniconda3/envs/itr_analysis/bin/makeblastdb"

# Additional tools (called directly)
MINIMAP2="minimap2"
MINIPOLISH="minipolish" 
VIRALVERIFY="viralverify"
RNAFOLD="RNAfold"
SEQKIT="seqkit"
PYTHON3="python3"

# Conda setup
CONDA_SETUP_PATH="~/miniconda3/etc/profile.d/conda.sh"

# Conda environment names
CONDA_ENV_QC="viralFlye"
CONDA_ENV_FLYE="flye"
CONDA_ENV_MINIASM="miniasm"
CONDA_ENV_PROKKA="prokka"
CONDA_ENV_CANU="canu"
CONDA_ENV_ITR="itr_analysis"


# parameters
GENOME_SIZE="150k"
THREADS="25" 

# Updated filtering parameters for better ITR handling
MIN_READ_LENGTH="3000"
MAX_READ_LENGTH="40000"  # Reduced to avoid concatemers but keep ITR-spanning reads
MIN_QUALITY="13"         # Slightly more permissive
HEAD_CROP="10"
TAIL_CROP="10"

# Updated Canu parameters for viral repeats
CORRECTED_ERROR_RATE="0.16"
MIN_OVERLAP="500"
UTG_REPEATS="20"
CANU_MEMORY="50G"


# ORF analysis paths
ORF_REFERENCE_FILE="${ORF_DATABASE_DIR}/${SAMPLE}_ORFs.fasta"
ORF_PROTEIN_FILE="${PROTEIN_DATABASE_DIR}/${SAMPLE}_protein.fasta"

# ITR analysis output
ITR_ANALYSIS_DIR="${BASE_DIR}/05_itr_analysis"
ANNOTATION_OUTPUT_DIR="${BASE_DIR}/04_annotation"

# Hairpin search parameters
HAIRPIN_MIN_LENGTH="50"
HAIRPIN_MAX_LENGTH="100"
HAIRPIN_SEARCH_REGION="3000"  # bp from contig ends to search

# Analysis parameters
CRS_LENGTH="150"  # Length of CRS region for trimming
MIN_CRS_DISTANCE="10000"  # Minimum distance between CRS positions
TERMINAL_ANALYSIS_SIZE="5000"  # Size for terminal region analysis
EXPECTED_GENOME_SIZE="138000"  # Expected final genome size
OVERLAP_SEARCH_SIZES="1000 2000 5000 10000"  # Overlap sizes for analysis
TERMINAL_DISTANCE="40000"  # Distance from ends to consider terminal ORFs

# Sample list (space-separated)
SAMPLE_LIST="B006 B021 B032 B044 D1701 S1-Japan"

# Annotation parameters
MIN_FEATURE_OVERLAP="50"  # Minimum overlap for feature collision detection in bp

PROXIMITY_THRESHOLD=500

# Temporary directory for processing
TEMP_DIR="/tmp"

# Retrospective annotation paths
PROPOSAL_INPUT_DIR="${BASE_DIR}/proposal"
PROPOSAL_OUTPUT_DIR="${BASE_DIR}/proposal"