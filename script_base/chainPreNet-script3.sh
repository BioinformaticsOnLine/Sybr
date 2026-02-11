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
output_base="$4"  # Now configurable output directory

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
    local output_path="$output_base/$relative_path"
    local output_dir=$(dirname "$output_path")
    
    # Create output directory if needed
    mkdir -p "$output_dir"
    
    echo "Processing $chain_file -> $output_path"
    echo "Using size files: $vaga_size and $habrotrocha_size"
    
    # Verify input chain exists and isn't empty
    if [ ! -s "$chain_file" ]; then
        echo "Error: Input chain file is empty or missing"
        return 1
    fi
    
    # Run chainPreNet with full path
    chainPreNet "$chain_file" "$vaga_size" "$habrotrocha_size" "$output_path"
    
    # Verify output
    if [ ! -f "$output_path" ]; then
        echo "Error: Output file was not created"
        return 1
    elif [ ! -s "$output_path" ]; then
        echo "Warning: Output file is empty"
        return 1
    fi
    
    echo "Successfully created $output_path"
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
