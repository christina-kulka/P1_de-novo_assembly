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
ITR_OUTPUT_DIR="${BASE_DIR}/04_final_results"

# Database paths
PFAM_DB="${BASE_DIR}/databases/Pfam-A.hmm.gz"
ORF_DATABASE_DIR="/home/ubuntu/data-volume/001_Raw_Data/Databases/ORF"
VIRAL_DB_DIR="${BASE_DIR}/databases/viral"
VIRAL_REF_FASTA="${VIRAL_DB_DIR}/parapoxvirus_references.fasta"

# Tool paths
PORECHOP="/home/ubuntu/miniconda3/bin/porechop"
NANOFILT="/home/ubuntu/miniconda3/envs/viralFlye/bin/NanoFilt"
TANDEMTOOLS="/home/ubuntu/TandemTools"
KRAKEN2="/home/ubuntu/miniconda3/bin/kraken2"
KRAKEN2_BUILD="/home/ubuntu/miniconda3/bin/kraken2-build"

# Conda environment names
CONDA_ENV_QC="viralFlye"
CONDA_ENV_FLYE="flye"
CONDA_ENV_MINIASM="miniasm"
CONDA_ENV_PROKKA="prokka"
CONDA_ENV_CANU="canu"

# Assembly parameters
GENOME_SIZE="160k"
MIN_OVERLAP="1000"
MIN_READ_LENGTH="1000"
MIN_QUALITY="8"
HEAD_CROP="10"
TAIL_CROP="10"
THREADS="25" 
CANU_MEMORY="50G"
