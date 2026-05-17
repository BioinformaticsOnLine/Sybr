import os
import sys

# Usage: python split_by_species.py input_data.txt /path/to/output_parent_dir

def split_file_by_species(input_filename, parent_output_dir="."):
    # Dictionary to hold lines for each species
    species_data = {}
    
    with open(input_filename, 'r') as infile:
        for line in infile:
            # Skip empty lines
            if not line.strip():
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
                
            # Split the species field at the first colon
            species_field = parts[3]
            species_name = species_field.split(':', 1)[0]
            
            # Add line to the appropriate species in the dictionary
            if species_name not in species_data:
                species_data[species_name] = []
            species_data[species_name].append(line)
    
    # Create parent directory if it doesn't exist
    os.makedirs(parent_output_dir, exist_ok=True)
    
    # Create species folders and write files
    for species_name, lines in species_data.items():
        # Define folder and file paths
        folder_name = f"{species_name}_lineage"
        folder_path = os.path.join(parent_output_dir, folder_name)
        output_filepath = os.path.join(folder_path, f"{species_name}.txt")
        
        # Create species folder
        os.makedirs(folder_path, exist_ok=True)
        
        # Write the file
        with open(output_filepath, 'w') as outfile:
            outfile.writelines(lines)
        print(f"Created file: {output_filepath} with {len(lines)} entries")

if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file> [parent_output_dir]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Use provided output directory or default to current directory
    parent_output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    
    if os.path.exists(input_file):
        split_file_by_species(input_file, parent_output_dir)
    else:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
