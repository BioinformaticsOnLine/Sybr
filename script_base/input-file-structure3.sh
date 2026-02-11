#!/bin/bash

# Check if primary species name and output directory were provided
if [ $# -lt 2 ]; then
    echo "Error: Please provide the primary species name and output directory as arguments"
    echo "Usage: $0 <primary_species_name> <output_directory>"
    exit 1
fi

primary_species=$1
output_dir=$2

# Create the new directory structure with primary species as root
mkdir -p "$output_dir/data/$primary_species"

# Get all species names by examining chainPreNet directory and removing -chainSplit
species=()
for dir in "$output_dir/chainPreNet"/*-chainSplit; do
    # Extract just the species name (remove -chainSplit)
    sp=$(basename "$dir" "-chainSplit")
    species+=("$sp")
done

# Process each species
for sp in "${species[@]}"; do
    echo "Processing species: $sp"
    
    # Create the new directory structure for each species
    mkdir -p "$output_dir/data/$primary_species/${sp}/chain"
    mkdir -p "$output_dir/data/$primary_species/${sp}/net"
    
    # Copy chain files (preserve originals)
    cp "$output_dir/chainPreNet/${sp}-chainSplit/chr"*.chain "$output_dir/data/$primary_species/${sp}/chain/" 2>/dev/null
    
    # Copy net files (preserve originals)
    cp "$output_dir/net/${sp}-chainSplit/chr"*.net "$output_dir/data/$primary_species/${sp}/net/" 2>/dev/null
done

echo "Directory reorganization complete. Original files preserved."
echo "New structure created under: $output_dir/data/$primary_species"
