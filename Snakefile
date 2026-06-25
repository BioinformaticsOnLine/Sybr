# Main Snakefile - includes all modules
import os
import glob
import yaml

# Include all modules
include: "modules/config_handling.smk"
include: "modules/satsuma_alignment.smk"   # de-novo Satsuma alignment (no-op when disabled)
include: "modules/lastz_alignment.smk"     # de-novo LastZ alignment   (no-op when disabled)
include: "modules/synteny_processing.smk"
include: "modules/eba_analysis.smk"
include: "modules/enrichment_analysis.smk"
include: "modules/alignment_processing.smk"
include: "modules/deschrambler.smk"
include: "modules/hgt_overlap.smk"

# The 'all' rule remains in main file
rule all:
    input:
        get_conditional_inputs()
