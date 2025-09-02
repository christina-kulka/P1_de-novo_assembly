# P1 De-novo Assembly Pipeline

This repository contains scripts for de-novo genome assembly of ORFV (Orf virus) samples using Oxford Nanopore long-read sequencing data. The pipeline supports multiple assembly approaches and includes comprehensive analysis tools.

## Samples

The pipeline processes the following samples: B006, B021, B032, B044, D1701, S1-Japan.

## Configuration

All paths and parameters are centralized in `config.sh`. Edit this file to match your system setup before running any scripts.

## Pipeline Overview

The pipeline is organized into three main workflows:

1. **Assembly Trial Workflow** - Testing and comparing different assembly approaches
2. **Final Canu Assembly Workflow** - Production assembly using optimized Canu parameters
3. **Retrospective Annotation** - Annotating existing .dna files with ORF information

---

## 1. Assembly Trial Workflow

This workflow tests multiple assembly approaches to determine the best strategy for each sample.

### Core Assembly Scripts

Run these scripts in order for each sample:

```bash
# Step 1: Quality Control
./01_quality_control.sh SAMPLE_NAME

# Step 2: Data Preprocessing
./02_preprocessing.sh SAMPLE_NAME

# Step 3: Assembly (choose one or more)
./03_assembly_flye.sh SAMPLE_NAME
./03_assembly_miniasm.sh SAMPLE_NAME
./03_assembly_canu.sh SAMPLE_NAME

# Step 4: Viral Verification (for miniasm assemblies)
./03b_viralverify.sh SAMPLE_NAME

# Step 5: ORF Annotation
./04_orf_annotation_prokka.sh SAMPLE_NAME [ASSEMBLY_TYPE]

# Step 6: ITR and ORF Analysis
./05_itr_orf_eval.sh SAMPLE_NAME [ASSEMBLY_TYPE]

# Step 7: Assembly Trimming (if needed)
./06_assembly_trimmer.sh SAMPLE_NAME [ASSEMBLY_TYPE] [MIN_DISTANCE]
```

### Script Descriptions

#### `01_quality_control.sh`
Generates quality control reports for raw Nanopore reads using NanoPlot.
- **Input**: Raw FASTQ files from sequencing
- **Output**: Quality plots and statistics in `01_quality_control/SAMPLE/`
- **Usage**: `./01_quality_control.sh SAMPLE_NAME`

#### `02_preprocessing.sh`
Cleans and filters raw reads for assembly.
- **Steps**: Adapter trimming (Porechop), quality filtering (NanoFilt), outlier removal (seqkit)
- **Output**: Clean reads in `02_preprocessing/SAMPLE/SAMPLE_super_clean_final.fastq`
- **Usage**: `./02_preprocessing.sh SAMPLE_NAME`

#### `03_assembly_flye.sh`
Performs genome assembly using Flye assembler.
- **Input**: Clean reads from preprocessing
- **Output**: Assembly in `03_assembly/SAMPLE/flye_out/assembly.fasta`
- **Usage**: `./03_assembly_flye.sh SAMPLE_NAME`

#### `03_assembly_miniasm.sh`
Performs fast genome assembly using Miniasm + Minipolish.
- **Input**: Raw reads (uses original data for better overlap detection)
- **Output**: Polished assembly in `03_assembly/SAMPLE/miniasm_out/polished_assembly.fasta`
- **Usage**: `./03_assembly_miniasm.sh SAMPLE_NAME`

#### `03_assembly_canu.sh`
Performs high-quality genome assembly using Canu.
- **Input**: Clean reads from preprocessing
- **Output**: Assembly in `03_assembly/SAMPLE/super_canu_out/SAMPLE.contigs.fasta`
- **Usage**: `./03_assembly_canu.sh SAMPLE_NAME`

#### `03b_viralverify.sh`
Verifies that assemblies contain viral sequences using viralVerify.
- **Input**: Miniasm assembly
- **Output**: Viral classification results in `03_assembly/SAMPLE/viral_verification/`
- **Usage**: `./03b_viralverify.sh SAMPLE_NAME`

#### `04_orf_annotation_prokka.sh`
Annotates genome assemblies with ORF information using Prokka.
- **Assembly types**: miniasm_viral, miniasm_raw, canu, canu_ultra, canu_super, microsynth
- **Output**: Annotations in `04_annotation/SAMPLE/ASSEMBLY_TYPE/`
- **Usage**: `./04_orf_annotation_prokka.sh SAMPLE_NAME [ASSEMBLY_TYPE]`

#### `05_itr_orf_eval.sh`
Analyzes assemblies for ITR (Inverted Terminal Repeats) and ORF content.
- **Functions**: Hairpin structure detection, ORF counting, concatemer detection
- **Output**: Analysis reports in `05_itr_analysis/SAMPLE/ASSEMBLY_TYPE/`
- **Usage**: `./05_itr_orf_eval.sh SAMPLE_NAME [ASSEMBLY_TYPE]`

#### `06_assembly_trimmer.sh`
Trims assemblies based on CRS (Concatemer Resolution Sequence) positions.
- **Input**: ITR analysis results with CRS positions
- **Output**: Trimmed genome in `06_trimmed_assembly/SAMPLE/ASSEMBLY_TYPE/`
- **Usage**: `./06_assembly_trimmer.sh SAMPLE_NAME [ASSEMBLY_TYPE] [MIN_DISTANCE]`

### Investigation Scripts (Python)

Located in `py_investigation/` directory:

#### `investigate_assemblies.py`
Compares multiple assemblies for a sample to identify issues like concatenation and inversions.
- **Usage**: `python3 investigate_assemblies.py SAMPLE_NAME`
- **Output**: Comparative analysis plots and reports

#### `itr_sniffer.py`
ORF-based ITR analysis that identifies terminal ORFs and potential ITR regions.
- **Usage**: `python3 itr_sniffer.py SAMPLE_NAME`
- **Output**: ORF distribution plots and ITR predictions

#### `orf_visualization.py`
Visualizes ORF annotations across different assembly methods.
- **Usage**: `python3 orf_visualization.py SAMPLE_NAME`
- **Output**: ORF distribution and comparison plots

#### `assembly_similarity.py`
Analyzes similarity between different assemblies using alignment and terminal comparison.
- **Usage**: `python3 assembly_similarity.py SAMPLE_NAME`
- **Output**: Assembly similarity reports and alignment statistics

#### `get_longest_contig.py`
Extracts the longest contig from multi-contig assemblies.
- **Usage**: `python3 get_longest_contig.py`
- **Output**: Longest contigs saved as separate FASTA files

#### `config_reader.py`
Utility module for Python scripts to read configuration from `config.sh`.
- **Usage**: Import in other Python scripts to access configuration variables

---

## 2. Final Canu Assembly Workflow

For production assemblies using optimized Canu parameters:

```bash
# Use the canu_ultra assembly type for final assemblies
./01_quality_control.sh SAMPLE_NAME
./02_preprocessing.sh SAMPLE_NAME
./03_assembly_canu.sh SAMPLE_NAME
./04_orf_annotation_prokka.sh SAMPLE_NAME canu_ultra
./05_itr_orf_eval.sh SAMPLE_NAME canu_ultra
./06_assembly_trimmer.sh SAMPLE_NAME canu_ultra
```

The `canu_ultra` assembly type uses optimized parameters for viral genome assembly with proper ITR handling.

### Automated Final Pipeline

#### `run_final_canu.sh`
Automated script that runs the complete final Canu assembly pipeline.
- **Usage**: `./run_final_canu.sh SAMPLE_NAME`
- **Function**: Executes all steps from quality control through assembly trimming using canu_ultra
- **Features**: Error checking, progress reporting, and output location summary

---

## 3. Retrospective Annotation Workflow

For annotating existing genome files (e.g., .dna files from colleagues) with ORF information:

### Scripts

#### `prokka_annotation_proposal.sh`
Annotates existing genome files with ORF information using Prokka and custom ORF databases.
- **Input**: Genome file (.dna, .gbk, or .fasta format)
- **Output**: GenBank file with ORF annotations
- **Usage**: `./prokka_annotation_proposal.sh INPUT_FILE SAMPLE_NAME OUTPUT_FILE`
- **Example**: `./prokka_annotation_proposal.sh genome.dna B006 annotated_output.gbk`

#### `merge_annotations.py`
Merges existing annotations with new ORF annotations from Prokka output.
- **Input**: Original annotation file and Prokka output
- **Output**: Merged GenBank file with both original and ORF annotations
- **Usage**: `python3 merge_annotations.py ORIGINAL_FILE PROKKA_GBK_FILE OUTPUT_FILE`
- **Options**: `--allow-overlaps` to permit overlapping features

### Helper Scripts

#### `4a_format_proteins.py`
Formats protein FASTA files for Prokka compatibility.
- **Usage**: `python3 4a_format_proteins.py INPUT_FILE OUTPUT_FILE`
- **Function**: Converts ORF protein sequences to Prokka-compatible format

#### `4b_restore_orf_names.py`
Restores original ORF names in Prokka output files.
- **Usage**: `python3 4b_restore_orf_names.py SAMPLE_NAME ORIGINAL_ORF_FILE [OUTPUT_DIR]`
- **Function**: Maps Prokka locus tags back to original ORF names

---

## Output Directory Structure

```
├── 01_quality_control/          # NanoPlot QC reports
├── 02_preprocessing/            # Cleaned reads
├── 03_assembly/                 # Assembly outputs
│   └── SAMPLE/
│       ├── flye_out/
│       ├── miniasm_out/
│       ├── super_canu_out/
│       └── viral_verification/
├── 04_annotation/               # Prokka annotations
├── 05_itr_analysis/            # ITR and ORF analysis
└── 06_trimmed_assembly/        # Final trimmed genomes
```

## Requirements

- Conda environments as specified in `config.sh`
- Tools: Flye, Canu, Miniasm, Prokka, NanoPlot, viralVerify, RNAfold
- Python packages: BioPython, matplotlib, numpy
- Reference databases: Pfam, custom ORF proteins

## Notes

- Edit `config.sh` before running any scripts
- Assembly types supported: miniasm_viral, miniasm_raw, canu, canu_ultra, canu_super, microsynth
- All scripts include error checking and progress reporting
- Python investigation scripts require the `config_reader.py` module