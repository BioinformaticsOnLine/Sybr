#!/bin/bash

# Check arguments - now expecting input and output directories
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

input_base="$1"
output_base="$2"

# Create output base directory if it doesn't exist
mkdir -p "$output_base"

# Loop through all species subdirectories
for species_dir in "$input_base"/*/; do
    # Remove trailing slash from path
    species_dir="${species_dir%/}"
    
    # Get just the species directory name
    species_name=$(basename "$species_dir")
    
    # Create corresponding output directory
    output_dir="$output_base/$species_name"
    mkdir -p "$output_dir"
    
    echo "Processing directory: $species_dir"
    echo "Output directory: $output_dir"
    
    # Find all T.net files in the species directory
    for tnet_file in "$species_dir"/*T.net; do
        # Check if file exists and is not empty
        if [ -s "$tnet_file" ]; then
            # Get the base filename without path
            filename=$(basename "$tnet_file")
            
            # Create output filename (remove the T before .net)
            output_filename="${filename/T.net/.net}"
            output_file="$output_dir/$output_filename"
            
            echo "Processing input: $tnet_file"
            echo "Output will be: $output_file"
            
            # Run the netSyntenic command with proper arguments
            netSyntenic "$tnet_file" "$output_file"
            
            # Verify output was created
            if [ -s "$output_file" ]; then
                echo "Successfully created $output_file"
            else
                echo "Warning: Output file is empty - processing may have failed"
                # Remove empty output file
                rm -f "$output_file"
            fi
            
            echo "Finished processing $filename"
        elif [ -f "$tnet_file" ]; then
            echo "Warning: Empty input file $tnet_file - skipping"
        fi
    done
done

echo "All processing complete. Output files are in $output_base"
