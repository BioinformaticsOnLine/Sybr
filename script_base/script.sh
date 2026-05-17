#!/bin/bash

# Default values
add_column_value="Chromosomes"
scaffolds_file=""
chromosomes_value=""

show_help() {
    echo "Usage: $0 -i <input_parent_dir> -e <ebd_dir> [-s <scaffolds_file> | -c <chromosomes_value>] [-h]"
    echo ""
    echo "Options:"
    echo "  -i    Parent directory containing resolution subdirectories (e.g., 100000, 300000, 500000) (required)"
    echo "  -e    EBA output directory (subdirectories will be created based on input resolutions) (required)"
    echo "  -s    Path to Scaffolds.txt file (use this OR -c, not both)"
    echo "  -c    Use this value for all species (e.g., 'Chromosomes') (use this OR -s, not both)"
    echo "  -h    Show this help message and exit"
}

# Parse command-line arguments
while getopts "i:e:s:c:h" opt; do
    case $opt in
        i) input_parent_dir="$OPTARG" ;;
        e) eba_dir="$OPTARG" ;;
        s) scaffolds_file="$OPTARG" ;;
        c) chromosomes_value="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$input_parent_dir" ] || [ -z "$eba_dir" ]; then
    echo "Error: Arguments -i and -e must be provided."
    show_help
    exit 1
fi

# Check if input directory exists
if [ ! -d "$input_parent_dir" ]; then
    echo "Error: Input directory '$input_parent_dir' not found."
    exit 1
fi

# Check mutual exclusivity
if [ -n "$scaffolds_file" ] && [ -n "$chromosomes_value" ]; then
    echo "Error: Cannot use both -s and -c options. Use one or the other."
    show_help
    exit 1
fi

# Initialize scaffold species array
scaffold_species=()

# Process based on which option was provided
if [ -n "$scaffolds_file" ]; then
    # Read scaffold species names from file
    if [ -f "$scaffolds_file" ]; then
        while IFS= read -r line || [[ -n "$line" ]]; do
            species=$(echo "$line" | tr -d '\r')
            scaffold_species+=("$species")
        done < "$scaffolds_file"
        echo "Using scaffolds file: $scaffolds_file"
        echo "Scaffold species: ${scaffold_species[@]}"
    else
        echo "Error: Scaffolds file not found: $scaffolds_file"
        exit 1
    fi
elif [ -n "$chromosomes_value" ]; then
    # Use the specified value for all species
    add_column_value="$chromosomes_value"
    echo "Using value '$chromosomes_value' for all species"
else
    # Default behavior (should not reach here due to config handling)
    echo "Using default value 'Chromosomes' for all species"
fi

# Function to convert resolution number to K notation
convert_resolution_to_k() {
    local resolution=$1
    local k_resolution
    
    if (( resolution % 1000 == 0 )); then
        k_resolution=$((resolution / 1000))
        echo "${k_resolution}k"
    else
        # If not divisible by 1000, use the number as-is
        echo "${resolution}"
    fi
}

# Function to process a single resolution directory
process_resolution_directory() {
    local input_dir="$1"
    local resolution="$2"  # This is the directory name (e.g., "100000")
    local eba_output_dir="$3"
    
    # Convert resolution to K notation (e.g., 100000 -> 100k)
    local k_resolution=$(convert_resolution_to_k "$resolution")
    
    echo "Processing ${resolution} (${k_resolution}) resolution files..."
    
    for d in "$input_dir"/*_out; do
        if [ -d "$d" ]; then
            local species_name=$(basename "$d" | sed 's/_out$//')
            local second_word=$(echo "$species_name" | cut -d'_' -f2 | tr 'A-Z' 'a-z')
            local output_file="$eba_output_dir/${second_word}_${k_resolution}_final.txt"

            # Determine the value to use
            local current_add_column_value="$add_column_value"
            
            # Only check scaffolds if we're using a scaffolds file
            if [ ${#scaffold_species[@]} -gt 0 ]; then
                local is_scaffold=0
                for scaffold in "${scaffold_species[@]}"; do
                    if [[ "$species_name" == "$scaffold" ]]; then
                        is_scaffold=1
                        break
                    fi
                done

                if [ "$is_scaffold" -eq 1 ]; then
                    current_add_column_value="Scaffolds"
                fi
            fi

            local input_file="$d/blocks_info"
            if [ ! -f "$input_file" ]; then
                echo "  Input file $input_file not found. Skipping..."
                continue
            fi

            # Copy and process the file
            cp "$input_file" "$output_file"
            sed -i '1,2d' "$output_file"

            # Modify only the first column (Adineta_vaga -> Adineta_vaga:resolution)
            awk -v res="$k_resolution" 'BEGIN {FS=OFS="\t"} {$1 = $1 ":" res; print}' "$output_file" > "$output_file.tmp" && mv "$output_file.tmp" "$output_file"

            # Process chromosome names (only in the second column)
            # Your existing chromosome processing commands here...
            # sed -i "s/CP075492.1_Adineta_vaga_breed_AD008_omosome/3/g" "$output_file"
            # etc.

            # Remove last column and add new column with appropriate value
            tmpfile=$(mktemp)
            awk 'BEGIN {FS=OFS="\t"} {NF--; print}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
            awk -v value="$current_add_column_value" 'BEGIN {FS=OFS="\t"} {print $0, value}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
            awk 'BEGIN {FS=OFS="\t"} {if (NF >= 8 && $8 !~ /^[+-]$/) $8 = "+"; print}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
            
            echo "  Created: $output_file"
        fi
    done
}

# Main processing loop
echo "Scanning for resolution directories in: $input_parent_dir"

# Find all directories in the input parent directory that look like resolution directories
# (containing only numbers, likely representing resolution values)
resolution_dirs=()
for dir in "$input_parent_dir"/*; do
    if [ -d "$dir" ]; then
        dir_name=$(basename "$dir")
        # Check if directory name consists only of digits (likely a resolution directory)
        if [[ "$dir_name" =~ ^[0-9]+$ ]]; then
            resolution_dirs+=("$dir_name")
        fi
    fi
done

# Sort resolution directories numerically
IFS=$'\n' resolution_dirs=($(sort -n <<< "${resolution_dirs[*]}"))
unset IFS

if [ ${#resolution_dirs[@]} -eq 0 ]; then
    echo "Error: No resolution directories found in '$input_parent_dir'."
    echo "Expected directories with numeric names (e.g., 100000, 300000, 500000)."
    exit 1
fi

echo "Found resolution directories: ${resolution_dirs[@]}"

# Create output directory if it doesn't exist
mkdir -p "$eba_dir"

# Process each resolution directory
for resolution_dir in "${resolution_dirs[@]}"; do
    # Convert resolution to K notation for output directory name
    k_resolution=$(convert_resolution_to_k "$resolution_dir")
    eba_output_subdir="$eba_dir/$k_resolution"
    
    # Create output subdirectory
    mkdir -p "$eba_output_subdir"
    
    # Process this resolution directory
    process_resolution_directory "$input_parent_dir/$resolution_dir" "$resolution_dir" "$eba_output_subdir"
done

echo ""
echo "Processing complete. Output files are in:"
for resolution_dir in "${resolution_dirs[@]}"; do
    k_resolution=$(convert_resolution_to_k "$resolution_dir")
    echo "  - $eba_dir/$k_resolution"
done
