#!/usr/bin/env python3
"""
Simple config reader for Python scripts to read bash config variables.
"""
import os
import re

def read_config():
    """Read configuration from config.sh file"""
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up one level to find config.sh
    config_path = os.path.join(script_dir, '..', 'config.sh')
    
    config = {}
    
    if not os.path.exists(config_path):
        # Fallback to hardcoded paths if config not found
        config['BASE_DIR'] = "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly"
        return config
    
    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            # Look for variable assignments
            match = re.match(r'^([A-Z_]+)="([^"]*)"', line)
            if match:
                key, value = match.groups()
                # Simple variable substitution for ${VAR} patterns
                while '${' in value:
                    for var in re.findall(r'\$\{([^}]+)\}', value):
                        if var in config:
                            value = value.replace(f'${{{var}}}', config[var])
                        else:
                            break
                    else:
                        continue
                    break
                config[key] = value
    
    return config

def get_assembly_paths(sample, config=None):
    """Get standard assembly paths for a sample"""
    if config is None:
        config = read_config()
    
    base_dir = config.get('BASE_DIR', "/home/ubuntu/data-volume/001_Raw_Data/Whole_Genome_Seq/ORFV_genome_assembly/P1_de-novo_assembly")
    
    return {
        'base_dir': base_dir,
        'assembly_dir': f"{base_dir}/03_assembly",
        'annotation_dir': f"{base_dir}/04_annotation", 
        'itr_analysis_dir': f"{base_dir}/05_itr_analysis",
        'trimmed_assembly_dir': f"{base_dir}/06_trimmed_assembly",
        'microsynth_dir': f"{base_dir}/00_raw_data_microsynth"
    }
