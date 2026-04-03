# Configuration handling module

import os
import glob

# Load static paths first, then user config.
# Snakemake does a SHALLOW merge of top-level keys, so dict-valued keys like
# "eba:" that appear in both files would be overwritten entirely.
# We fix this below with an explicit deep-merge after both files are loaded.
configfile: "pipeline_paths.yaml"
configfile: "run_sybr_config.yaml"

# Deep-merge any top-level dicts that are intentionally split across both files.
# pipeline_paths.yaml is loaded first so run_sybr_config.yaml wins on conflict.
def _deep_merge_from_files(cfg, *yaml_files):
    """Re-read each yaml file and deep-merge its dict keys into cfg."""
    import yaml as _yaml
    layers = []
    for path in yaml_files:
        try:
            with open(path) as _f:
                layers.append(_yaml.safe_load(_f) or {})
        except Exception:
            layers.append({})
    for layer in layers:
        for k, v in layer.items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k] = {**cfg[k], **v}   # shallow merge of the sub-dict
    return cfg

config = _deep_merge_from_files(config, "pipeline_paths.yaml", "run_sybr_config.yaml")

# Rebase all "inputs/..." and "outputs/..." paths in the config using
# base_input_dir and base_output_dir from run_sybr_config.yaml.
# This lets users set a single pair of absolute dirs and have every
# downstream path automatically point to the right place.
def _rebase_paths(cfg, in_base, out_base):
    """
    Walk every string value in cfg (recursively).
    - Values starting with "inputs/"  → replace that prefix with in_base
    - Values starting with "outputs/" → replace that prefix with out_base
    Already-absolute paths are left untouched.
    """
    if isinstance(cfg, dict):
        return {k: _rebase_paths(v, in_base, out_base) for k, v in cfg.items()}
    if isinstance(cfg, list):
        return [_rebase_paths(v, in_base, out_base) for v in cfg]
    if isinstance(cfg, str) and not os.path.isabs(cfg):
        if cfg.startswith("inputs/") or cfg == "inputs":
            return os.path.join(in_base, cfg[len("inputs/"):])
        if cfg.startswith("outputs/") or cfg == "outputs":
            return os.path.join(out_base, cfg[len("outputs/"):])
    return cfg

_in_base  = config.get("base_input_dir",  "inputs")
_out_base = config.get("base_output_dir", "outputs")

# Only rebase if explicit base dirs were provided (i.e. not the bare defaults)
if _in_base != "inputs" or _out_base != "outputs":
    config = _rebase_paths(config, _in_base, _out_base)

try:
    # Validate config structure
    required_keys = ["satsuma_alignments", "synteny_results", "scripts", "reference_name", "reference_species", "run_stages"]
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    # Handle optional scaffolds_file
    if "eba_format" in config:
        if "scaffolds_file" not in config["eba_format"]:
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"
        elif config["eba_format"]["scaffolds_file"] == "Chromosomes":
            config["eba_format"]["scaffolds_file"] = "ALL_CHROMOSOMES"
        
except Exception as e:
    raise ValueError(f"Error loading config: {str(e)}")

# Function to resolve relative paths
def resolve_path(path, base_dir=None):
    """Resolve relative paths to absolute paths relative to workflow directory"""
    if path is None:
        return None
    if isinstance(path, str):
        if os.path.isabs(path):
            return path
        if path in ["ALL_CHROMOSOMES", "Chromosomes"]:
            return path
        if base_dir:
            return os.path.join(base_dir, path)
        return os.path.join(workflow.basedir, path)
    return path

def resolve_config_paths(config_dict, base_key=None):
    """Recursively resolve all paths in config dictionary"""
    resolved_config = {}
    for key, value in config_dict.items():
        if isinstance(value, dict):
            resolved_config[key] = resolve_config_paths(value, base_key=key)
        elif isinstance(value, str) and any(path_term in key.lower() for path_term in ['dir', 'file', 'path', 'results', 'alignments', 'input']):
            if base_key == 'eba_tools' and key == 'eba_script_path':
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

# Define ebrs_dir globally
ebrs_dir = config["eba_format"].get("ebrs_dir", "EBRs")

# Get synteny parameters if provided
if "synteny_params" in config:
    WINDOW_SIZES = config["synteny_params"]["window_sizes"]
    STEP_SIZE = config["synteny_params"]["step_size"]
    CUSTOM_RESOLUTIONS = config["synteny_params"].get("custom_resolutions", False)
else:
    WINDOW_SIZES = [100000, 300000, 500000]
    STEP_SIZE = 30000
    CUSTOM_RESOLUTIONS = False

# Create output directory if it doesn't exist
os.makedirs(SYNTENY_RESULTS, exist_ok=True)

# Auto-discover samples
if RUN_STAGES.get("synteny_processing", True):
    try:
        sample_files = glob.glob(f"{SATSUMA_ALIGNMENTS}/*_*.txt")
        SAMPLES = [os.path.splitext(os.path.basename(f))[0] for f in sample_files]
        SAMPLES = sorted(SAMPLES)
        print(f"Discovered {len(SAMPLES)} samples: {SAMPLES}")
    except Exception as e:
        raise ValueError(f"Error discovering samples: {str(e)}")

    if not SAMPLES:
        raise ValueError(f"No sample files found in {SATSUMA_ALIGNMENTS}")
else:
    SAMPLES = []

# Function for conditional inputs - UPDATED to only include actual outputs
def get_conditional_inputs():
    inputs = []
    
    # Synteny processing outputs (only if enabled)
    if RUN_STAGES.get("synteny_processing", True) and SAMPLES:
        inputs.extend([
            expand(f"{SYNTENY_RESULTS}/{{sample}}_out/st_input", sample=SAMPLES),
        ])
        
        # Add outputs for each window size dynamically
        for window_size in WINDOW_SIZES:
            inputs.extend(
                expand(f"{SYNTENY_RESULTS}/synteny_results/{window_size}/{{sample}}_out/synteny_assign_done", sample=SAMPLES)
            )

        # Add EBA format and overlap resolution output markers
        if config.get("out_final"):
            inputs.extend([
                os.path.join(SYNTENY_RESULTS, config["out_final"], ".done_marker"),  # EBA format completion
                os.path.join(SYNTENY_RESULTS, config["out_final"], ".overlap_resolve_done"),  # Overlap resolution completion
                os.path.join(SYNTENY_RESULTS, ".concatination_done"),
                os.path.join(SYNTENY_RESULTS, ".reformatting_done"),
                os.path.join(SYNTENY_RESULTS, ".plots_generated_done"),
                os.path.join(SYNTENY_RESULTS, ".split_by_species_done"),
                os.path.join(SYNTENY_RESULTS, ".linear_synteny_plots_done")
            ])

  

    # EBA analysis outputs (only if enabled)
    if RUN_STAGES.get("eba_analysis", True):
        # Add copy_to_pre_EBA marker file as input when eba_analysis is enabled
        if config.get("eba_format", {}).get("pre_EBA_dir"):
            pre_EBA_dir = os.path.join(workflow.basedir, config["eba_format"]["pre_EBA_dir"])
            inputs.extend([
                os.path.join(pre_EBA_dir, "EBA-input", ".copy_complete"),
                # EBA analysis outputs - only marker files, not directories
                os.path.join(pre_EBA_dir, "EBA_WORKING/eba_complete.marker"),
                # Main EBA result file
                os.path.join(pre_EBA_dir, "EBA_WORKING/EBA_OUT/Merge/Result_Merge.final")
            ])
        

        # Add msHSBs processing outputs
        if config.get("eba_format", {}).get("mshsbs_dir"):
            mshsbs_dir = os.path.join(workflow.basedir, config["eba_format"]["mshsbs_dir"])
            inputs.extend([
                os.path.join(mshsbs_dir, "msHSBs.txt"),
                os.path.join(mshsbs_dir, "msHSBs_done.txt"),
                # msHSBs sequence extraction outputs
                os.path.join(mshsbs_dir, "msHSBs_sequences.fasta"),
                os.path.join(mshsbs_dir, "msHSBs_fasta_done.txt"),
            ])
        
        # Add EBRs processing outputs
        if config.get("eba_format", {}).get("ebrs_dir"):
            ebrs_dir = os.path.join(workflow.basedir, config["eba_format"]["ebrs_dir"])
            inputs.extend([
                os.path.join(ebrs_dir, "EBRs.txt"),
                os.path.join(ebrs_dir, "ebrs_processed.marker"),
                # EBR stats report (produced by ebr_stats_report rule)
                os.path.join(ebrs_dir, "EBRs_stats.txt"),
                os.path.join(ebrs_dir, "EBRs_stats.html"),
            ])

    
    # Enrichment analysis outputs (only if enabled)
    if RUN_STAGES.get("enrichment_analysis", True):
        # Get msHSBs directory from enrichment config, fallback to eba_format if not specified
        msHSBs_dir = config["enrichment"].get("msHSBs_dir", config["eba_format"].get("mshsbs_dir", "msHSBs"))
        msHSBs_dir = os.path.join(workflow.basedir, msHSBs_dir)
        
        enrichment_inputs = [
            os.path.join(msHSBs_dir, "msHSBs_NCBI_genes_overlap.txt"),
            os.path.join(msHSBs_dir, "foreground.txt"),
            os.path.join(msHSBs_dir, "background.txt"),
        ]
        
        if config["getenrich"]["r"] == "ko":
            enrichment_inputs.append(os.path.join(msHSBs_dir, "3kegg_annotationTOgenes.txt"))
        
        # Add EBRs-related outputs if eba_analysis is also enabled
        if RUN_STAGES.get("eba_analysis", True) and config.get("eba_format", {}).get("ebrs_dir"):
            ebrs_dir = os.path.join(workflow.basedir, config["eba_format"]["ebrs_dir"])
            
            # Create the EBRs_Enrichment_results directory path
            all_ebrs_dir = os.path.join(ebrs_dir, "EBRs_Enrichment_results")
            
            enrichment_inputs.extend([
                os.path.join(all_ebrs_dir, "EBRs_NCBI_genes_overlap.txt"),
                os.path.join(all_ebrs_dir, "foreground.txt"),
                os.path.join(all_ebrs_dir, "background.txt"),
                os.path.join(ebrs_dir, "EBRs_split_done.txt"),
                os.path.join(ebrs_dir, "EBRs_subdirs_overlap_done.txt"),
                os.path.join(ebrs_dir, "EBRs_Enrichment_results", "getENRICH_done.txt"),
                os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "getENRICH_done.txt"),
                expand(
                    os.path.join(ebrs_dir, "{lineage}_lineage", "getENRICH_done.txt"),
                    lineage=[
                        os.path.basename(d).replace("_lineage", "")
                        for d in glob.glob(f"{ebrs_dir}/*_lineage")
                    ]
                ),

                
            ])
            
            if config["getenrich"]["r"] == "ko":
                enrichment_inputs.append(os.path.join(all_ebrs_dir, "3kegg_annotationTOgenes.txt"))
        
        inputs.extend(enrichment_inputs)
        

    
    # Alignment processing outputs (only if enabled AND if you have chainNet module)
    if RUN_STAGES.get("chainNet_generation", True):  # Set to True when you have chainNet module
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            expand(f"{ALIGN_OUTPUT}/2bit/{{species}}.2bit", species=get_fasta_species()),
            f"{ALIGN_OUTPUT}/clean_axt.done",
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}.chain", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/chain/{{species}}-chainSplit", species=get_alignment_species()),
            expand(f"{ALIGN_OUTPUT}/fasize/{{species}}.size", species=get_fasta_species()),
            f"{ALIGN_OUTPUT}/chainPreNet/processing_done.txt",
            f"{ALIGN_OUTPUT}/netFiles.done",
            f"{ALIGN_OUTPUT}/netSyntenic.done",
            f"{ALIGN_OUTPUT}/data/{REFERENCE}/reorganized.done"
        ])
    
    # Deschrambler outputs (only if enabled AND if you have deschrambler module)
    if RUN_STAGES.get("Ancestor_seq_recunstruction", True):  # Set to True when you have deschrambler module
        ALIGN_OUTPUT = config["chainNet"].get("output_dir", ".")
        inputs.extend([
            f"{ALIGN_OUTPUT}/deschrambler.done"
        ])
    
    return inputs

# Helper functions for other modules (keep these but they won't be used if modules are disabled)
def get_fasta_species():
    """Get list of species from fasta files"""
    if not RUN_STAGES.get("chainNet_generation", True):
        return []
    fasta_files = []
    seq_dir = config.get("chainNet", {}).get("seq_dir", "seq")
    for ext in ['.fa', '.fna', '.fasta']:
        fasta_files.extend(glob.glob(f"{seq_dir}/*{ext}"))
    return [os.path.splitext(os.path.basename(f))[0] for f in fasta_files]

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

# Conditional rule function
def should_run_rule(condition_key):
    """Check if a rule should run based on configuration"""
    return RUN_STAGES.get(condition_key, True)
