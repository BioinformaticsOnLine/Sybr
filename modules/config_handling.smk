# Configuration handling module

import os
import glob

# Use Snakemake's native configfile directive
# Note: The --configfile parameter from command line will override this
configfile: "run_sybr_config.yaml"  # Changed to match your actual config file

# Now we can access the config through snakemake.config
try:
    # Validate config structure
    required_keys = ["satsuma_alignments", "synteny_results", "scripts", "reference_name", "reference_species", "run_stages"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    # Handle optional scaffolds_file - use "ALL_CHROMOSOMES" as special value
    if "eba_format" in config:
        if "scaffolds_file" not in config["eba_format"]:
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"  # Special value
        elif config["eba_format"]["scaffolds_file"] == "Chromosomes":
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"  # Standardize the value
        
except Exception as e:
    raise ValueError(f"Error loading config: {str(e)}")

# Function to resolve relative paths
def resolve_path(path, base_dir=None):
    """Resolve relative paths to absolute paths relative to workflow directory"""
    if path is None:
        return None
    if isinstance(path, str):
        # If it's already an absolute path, return as is
        if os.path.isabs(path):
            return path
        # If it's a special value, return as is
        if path in ["ALL_CHROMOSOMES", "Chromosomes"]:
            return path
        # Resolve relative path relative to workflow directory
        if base_dir:
            return os.path.join(base_dir, path)
        return os.path.join(workflow.basedir, path)
    return path

# Resolve all paths in config
def resolve_config_paths(config_dict, base_key=None):
    """Recursively resolve all paths in config dictionary"""
    resolved_config = {}
    for key, value in config_dict.items():
        if isinstance(value, dict):
            resolved_config[key] = resolve_config_paths(value, base_key=key)
        elif isinstance(value, str) and any(path_term in key.lower() for path_term in ['dir', 'file', 'path', 'results', 'alignments', 'input']):
            # Only resolve paths for keys that suggest they are file/directory paths
            if base_key == 'eba_tools' and key == 'eba_script_path':
                # Special handling for eba_script_path - relative to workflow basedir
                resolved_config[key] = resolve_path(value, workflow.basedir)
            else:
                resolved_config[key] = resolve_path(value)
        else:
            resolved_config[key] = value
    return resolved_config

# Apply path resolution to config
config = resolve_config_paths(config)

# Global variables accessible in all modules
SATSUMA_ALIGNMENTS = config["satsuma_alignments"]
SYNTENY_RESULTS = config["synteny_results"]
SCRIPTS = config["scripts"]
REFERENCE = config["reference_species"]
RUN_STAGES = config["run_stages"]

# ADD THIS LINE - Define ebrs_dir globally
ebrs_dir = config["eba_format"].get("ebrs_dir", "EBRs")

# Create output directory if it doesn't exist
os.makedirs(SYNTENY_RESULTS, exist_ok=True)

# Auto-discover samples with validation (only if synteny processing is enabled)
if RUN_STAGES.get("synteny_processing", True):
    try:
        # Look for files with pattern: two_names_separated_by_underscore.txt
        sample_files = glob.glob(f"{SATSUMA_ALIGNMENTS}/*_*.txt")
        
        # Extract sample names by removing the .txt extension from filename
        SAMPLES = [os.path.splitext(os.path.basename(f))[0] for f in sample_files]
        SAMPLES = sorted(SAMPLES)  # Ensure consistent order
        
        print(f"Discovered {len(SAMPLES)} samples: {SAMPLES}")
        
    except Exception as e:
        raise ValueError(f"Error discovering samples: {str(e)}")

    # Validate we found samples
    if not SAMPLES:
        raise ValueError(f"No sample files found in {SATSUMA_ALIGNMENTS} - expected files with pattern: name1_name2.txt")
else:
    SAMPLES = []  # Empty list if synteny processing is disabled

# Helper function to get fasta species with configurable directory
def get_fasta_species():
    """Get list of species from fasta files"""
    if not RUN_STAGES.get("chainNet_generation", True):
        return []
    fasta_files = []
    seq_dir = config.get("chainNet", {}).get("seq_dir", "seq")
    for ext in ['.fa', '.fna', '.fasta']:
        fasta_files.extend(glob.glob(f"{seq_dir}/*{ext}"))
    return [os.path.splitext(os.path.basename(f))[0] for f in fasta_files]

# Helper function to get alignment species with configurable directory
def get_alignment_species():
    """Get list of species from alignment files"""
    if not RUN_STAGES.get("chainNet_generation", True):
        return []
    alignments_dir = config.get("chainNet", {}).get("lastZ_alignments", "alignments2")
    axt_files = glob.glob(f"{alignments_dir}/*.axt")
    return [
        os.path.splitext(os.path.basename(f))[0]
        for f in axt_files
        if os.path.splitext(os.path.basename(f))[0] != REFERENCE
    ]

# Function for conditional inputs - must be defined AFTER config is loaded
def get_conditional_inputs():
    inputs = []
    
    # Synteny processing outputs (only if enabled)
    if RUN_STAGES.get("synteny_processing", True) and SAMPLES:
        inputs.extend([
            expand(f"{SYNTENY_RESULTS}/{{sample}}_out/st_input", sample=SAMPLES),
            expand(f"{SYNTENY_RESULTS}/synteny_out/100000/{{sample}}_out/synteny_assign_done", sample=SAMPLES),
            expand(f"{SYNTENY_RESULTS}/synteny_out/300000/{{sample}}_out/synteny_assign_done", sample=SAMPLES),
            expand(f"{SYNTENY_RESULTS}/synteny_out/500000/{{sample}}_out/synteny_assign_done", sample=SAMPLES)
        ])
    
    # EBA analysis outputs (only if enabled)
    if RUN_STAGES.get("eba_analysis", True):
        inputs.append(os.path.join(config["eba_format"]["pre_EBA_dir"], ".done_marker")) 

        inputs.extend([
            os.path.join(config["eba_format"]["pre_EBA_dir"], "resolution_complete.txt"),
            # Updated path for EBA-input
            os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input/input_preparation_done.txt"),
            os.path.join(config["eba_format"]["mshsbs_dir"], "msHSBs_done.txt"),
            os.path.join(ebrs_dir, "all_EBRs/intermediate_files/brkfile3_done.txt"),
            os.path.join(ebrs_dir, "all_EBRs/5KBextenction_brk.txt"),
            # Add EBA working directory outputs
            os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/eba_complete.marker"),
            os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT/300/EBA_OutFiles/final_classify.eba7")
        ])
    
    # Enrichment analysis outputs (only if enabled)
    if RUN_STAGES.get("enrichment_analysis", True):
        msHSBs_dir = config["enrichment"].get("msHSBs_dir", "msHSBs")
        
        # Always include these files
        enrichment_inputs = [
            os.path.join(config["eba_format"]["mshsbs_dir"], "background.txt"),
            os.path.join(ebrs_dir, "all_EBRs/background.txt"),
            os.path.join(ebrs_dir, "EBRs_split_done.txt"),
            os.path.join(ebrs_dir, "EBRs_subdirs_overlap_done.txt"),
            os.path.join(ebrs_dir, "EBR_configs_cloned.txt"),
            os.path.join(ebrs_dir, "all_EBRs/getENRICH_done.txt"),
            expand(os.path.join(ebrs_dir, "{lineage}_lineage/getENRICH_done.txt"),
                   lineage=[os.path.basename(d).replace("_lineage", "") 
                          for d in glob.glob(f"{ebrs_dir}/*_lineage")]),
            os.path.join(config["eba_format"]["mshsbs_dir"], "getENRICH_done.txt"),
            
            os.path.join(msHSBs_dir, "msHSBs_NCBI_genes_overlap.txt"),
            os.path.join(msHSBs_dir, "foreground.txt"),
            os.path.join(msHSBs_dir, "background.txt"),
        ]
        
        # Only include KEGG-related files if resource is "ko"
        if config["getenrich"]["r"] == "ko":
            enrichment_inputs.extend([
                os.path.join(msHSBs_dir, "3kegg_annotationTOgenes.txt"),
                os.path.join(ebrs_dir, "all_EBRs/3kegg_annotationTOgenes.txt")
            ])
        
        inputs.extend(enrichment_inputs)
    
    # Alignment processing outputs (only if enabled)
    if RUN_STAGES.get("chainNet_generation", True):
        # Get alignment output directory from config
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            expand(f"{ALIGN_OUTPUT}/2bit/{{species}}.2bit", species=get_fasta_species()),
            f"{ALIGN_OUTPUT}/clean_axt.done",
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}.chain", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}-chainSplit", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/fasize/{{species}}.size", species=get_fasta_species()),
            # Replace directory with completion marker
            f"{ALIGN_OUTPUT}/chainPreNet/processing_done.txt",
            f"{ALIGN_OUTPUT}/netFiles.done",
            f"{ALIGN_OUTPUT}/netSyntenic.done",
            f"{ALIGN_OUTPUT}/data/{REFERENCE}/reorganized.done"
        ])
    
    # Deschrambler outputs (only if enabled) - CORRECTED
    if RUN_STAGES.get("Ancestor_seq_recunstruction", True):
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            f"{ALIGN_OUTPUT}/deschrambler.done"  # This matches the output of our DESCHRAMBLER rule
        ])
    
    return inputs

# Conditional rule function
def should_run_rule(condition_key):
    """Check if a rule should run based on configuration"""
    return RUN_STAGES.get(condition_key, True)
