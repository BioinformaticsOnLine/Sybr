#!/bin/bash

# Enable debug output
set -x

# Check arguments - now expecting 4 parameters with output directory
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <chains_folder> <vaga_size_filename> <size_files_directory> <output_directory>"
    exit 1
fi

chains_folder="$1"
vaga_size_file="$2"
size_files_dir="$3"
output_base="$4"  # Configurable output directory

# Verify inputs
[ ! -d "$chains_folder" ] && echo "Error: Chains folder not found" && exit 1
[ ! -d "$size_files_dir" ] && echo "Error: Size files directory not found" && exit 1

vaga_size_path="$size_files_dir/$vaga_size_file"
[ ! -f "$vaga_size_path" ] && echo "Error: Vaga size file not found" && exit 1

# Create output directory if it doesn't exist
mkdir -p "$output_base"

# Process chain files
process_chain() {
    local chain_file="$1"
    local vaga_size="$2"
    local habrotrocha_size="$3"
    local output_base="$4"
    
    # Create output path maintaining folder structure
    local relative_path="${chain_file#$chains_folder/}"
    local chain_name=$(basename "$chain_file" .chain)
    local output_dir="$output_base/$(dirname "$relative_path")"
    
    # Create output directory if needed
    mkdir -p "$output_dir"
    
    # Define output file paths
    local net_output="$output_dir/${chain_name}T.net"
    local qnet_output="$output_dir/${chain_name}Q.net"  # Changed to append Q
    
    echo "Processing $chain_file"
    echo "Output files: $net_output and $qnet_output"
    echo "Using size files: $vaga_size and $habrotrocha_size"
    
    # Verify input chain exists and isn't empty
    if [ ! -s "$chain_file" ]; then
        echo "Error: Input chain file is empty or missing"
        return 1
    fi
    
    # Run chainNet
    chainNet "$chain_file" "$vaga_size" "$habrotrocha_size" "$net_output" "$qnet_output"
    
    # Verify outputs
    for output_file in "$net_output" "$qnet_output"; do
        if [ ! -f "$output_file" ]; then
            echo "Error: Output file $output_file was not created"
            return 1
        elif [ ! -s "$output_file" ]; then
            echo "Warning: Output file $output_file is empty"
            return 1
        fi
    done
    
    echo "Successfully created $net_output and $qnet_output"
    return 0
}

# Main processing loop
find "$chains_folder" -type f -name "chr*.chain" | while read -r chain_file; do
    subfolder=$(dirname "$chain_file")
    subfolder_name=$(basename "$subfolder")
    habrotrocha_name="${subfolder_name//-chainSplit/}"
    
    # Handle the hobrotrocha/habrotrocha naming inconsistency
    [ "$habrotrocha_name" = "hobrotrocha" ] && habrotrocha_name="habrotrocha"
    
    habrotrocha_size="$size_files_dir/$habrotrocha_name.size"
    
    if [ ! -f "$habrotrocha_size" ]; then
        echo "Skipping $chain_file - size file $habrotrocha_size not found"
        continue
    fi
    
    process_chain "$chain_file" "$vaga_size_path" "$habrotrocha_size" "$output_base" || continue
done

echo "Processing complete"
set +x
