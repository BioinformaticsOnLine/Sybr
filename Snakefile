# Main Snakefile - includes all modules
import os
import glob
import yaml

# Include all modules
include: "modules/config_handling.smk"
include: "modules/synteny_processing.smk"
include: "modules/eba_analysis.smk"
include: "modules/enrichment_analysis.smk"
include: "modules/alignment_processing.smk"
include: "modules/deschrambler.smk"

# The 'all' rule remains in main file
rule all:
    input:
        get_conditional_inputs()
