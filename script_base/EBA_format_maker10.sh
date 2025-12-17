#!/bin/bash

# Default values
add_column_value="Chromosomes"
scaffolds_file=""
chromosomes_value=""

show_help() {
    echo "Usage: $0 -i <input_parent_dir> -e <ebd_dir> [-s <scaffolds_file> | -c <chromosomes_value>] [-h]"
    echo ""
    echo "Options:"
    echo "  -i    Parent directory containing resolution subdirectories (100000, 300000, 500000) (required)"
    echo "  -e    EBA output directory (100/300/500 subdirs will be created) (required)"
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

# Create resolution-specific subdirectories
mkdir -p "$eba_dir"/100 "$eba_dir"/300 "$eba_dir"/500

process_directory() {
    local input_dir="$1"
    local resolution="$2"
    local eba_output_dir="$3"

    for d in "$input_dir"/*_out; do
        if [ -d "$d" ]; then
            local species_name=$(basename "$d" | sed 's/_out$//')
            local second_word=$(echo "$species_name" | cut -d'_' -f2 | tr 'A-Z' 'a-z')
            local resolution_number=$(echo "$resolution" | grep -o '^[0-9]\+')
            local output_file="$eba_output_dir/${second_word}_${resolution_number}_final.txt"

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
                echo "Input file $input_file not found. Skipping..."
                continue
            fi

            cp "$input_file" "$output_file"
            sed -i '1,2d' "$output_file"

            # Modify only the first column (Adineta_vaga -> Adineta_vaga:resolution)
            awk -v res="$resolution" 'BEGIN {FS=OFS="\t"} {$1 = $1 ":" res; print}' "$output_file" > "$output_file.tmp" && mv "$output_file.tmp" "$output_file"

            # Process chromosome names (only in the second column)
            # Your existing chromosome processing commands here...
            # sed -i "s/CP075492.1_Adineta_vaga_breed_AD008_omosome/3/g" "$output_file"
            # etc.

            tmpfile=$(mktemp)
            awk 'BEGIN {FS=OFS="\t"} {NF--; print}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
            awk -v value="$current_add_column_value" 'BEGIN {FS=OFS="\t"} {print $0, value}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
            awk 'BEGIN {FS=OFS="\t"} {if (NF >= 8 && $8 !~ /^[+-]$/) $8 = "+"; print}' "$output_file" > "$tmpfile" && mv "$tmpfile" "$output_file"
        fi
    done
}

# Process each resolution directory
echo "Processing 100k resolution files..."
process_directory "$input_parent_dir/100000" "100k" "$eba_dir/100"

echo "Processing 300k resolution files..."
process_directory "$input_parent_dir/300000" "300k" "$eba_dir/300"

echo "Processing 500k resolution files..."
process_directory "$input_parent_dir/500000" "500k" "$eba_dir/500"

echo "Processing complete. Output files are in:"
echo "  - $eba_dir/100"
echo "  - $eba_dir/300"
echo "  - $eba_dir/500"