#!/usr/bin/env python3
"""
Script to concatenate all *_nonoverlap.txt files from numeric subdirectories
"""

import os
import sys
import argparse
from pathlib import Path

def combine_nonoverlap_files(parent_dir):
    """
    Combine all *_nonoverlap.txt files from numeric subdirectories
    and save them inside each subdirectory
    """
    parent_path = Path(parent_dir)
    
    if not parent_path.exists():
        print(f"Error: Parent directory '{parent_dir}' does not exist.")
        sys.exit(1)
    
    # Find all subdirectories that start with numeric values
    numeric_subdirs = []
    for item in parent_path.iterdir():
        if item.is_dir() and item.name[0].isdigit():
            numeric_subdirs.append(item)
    
    if not numeric_subdirs:
        print(f"No numeric subdirectories found in '{parent_dir}'")
        sys.exit(1)
    
    # Process each numeric subdirectory
    for subdir in numeric_subdirs:
        print(f"\nProcessing directory: {subdir.name}")
        
        # Find all *_nonoverlap.txt files in this subdirectory
        nonoverlap_files = list(subdir.glob("*_nonoverlap.txt"))
        
        if not nonoverlap_files:
            print(f"  No *_nonoverlap.txt files found in {subdir.name}")
            continue
        
        # Create output filename - save INSIDE the subdirectory
        output_file = subdir / f"all_nonoverlap_combined_{subdir.name}.txt"
        
        # Combine all files WITHOUT headers
        total_lines = 0
        files_count = 0
        
        try:
            with open(output_file, 'w') as out_f:
                for file_path in sorted(nonoverlap_files):
                    print(f"  Reading: {file_path.name}")
                    try:
                        with open(file_path, 'r') as in_f:
                            lines = in_f.readlines()
                            if lines:
                                # Write all lines directly (no headers)
                                out_f.writelines(lines)
                                total_lines += len(lines)
                                files_count += 1
                                
                                # Add newline between files if last line doesn't end with newline
                                if lines[-1].strip() and not lines[-1].endswith('\n'):
                                    out_f.write("\n")
                                    
                    except Exception as e:
                        print(f"  Error reading {file_path}: {e}")
            
            print(f"  Created: {output_file.name}")
            print(f"  Combined {files_count} files")
            print(f"  Total lines: {total_lines}")
            
        except Exception as e:
            print(f"  Error writing to {output_file}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Combine *_nonoverlap.txt files from numeric subdirectories"
    )
    parser.add_argument(
        "parent_dir",
        help="Path to parent directory (e.g., pre-EBA1)"
    )
    
    args = parser.parse_args()
    
    print(f"Starting combination of non-overlap files...")
    print(f"Parent directory: {args.parent_dir}")
    
    combine_nonoverlap_files(args.parent_dir)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
