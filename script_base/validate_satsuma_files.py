#!/usr/bin/env python3
"""
Validate Satsuma alignment file format before running Snakemake pipeline.
"""

import os
import sys
import yaml
import argparse
from pathlib import Path

def validate_satsuma_file(file_path):
    """
    Validate a single Satsuma alignment file.
    
    Expected format: 8 columns separated by tabs
    Columns 1-6: integers (should be non-negative)
    Column 7: float (score/probability, typically between 0-1)
    Column 8: + or - only (strand)
    """
    errors = []
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
        if not lines:
            errors.append(f"File {file_path} is empty")
            return errors
            
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            parts = line.split('\t')
            
            # Check number of columns
            if len(parts) != 8:
                errors.append(f"Line {line_num}: Expected 8 columns, found {len(parts)}")
                continue
            
            # Validate columns 1-6 as integers
            for i in range(6):  # First 6 columns should be integers
                try:
                    val = int(parts[i])
                    if val < 0:
                        errors.append(f"Line {line_num}, column {i+1}: Negative value '{parts[i]}'")
                except ValueError:
                    errors.append(f"Line {line_num}, column {i+1}: Not an integer '{parts[i]}'")
            
            # Validate column 7 as float (score/probability)
            try:
                val = float(parts[6])
                # Optional: check if value is between 0 and 1 (typical for scores)
                # if val < 0 or val > 1:
                #     errors.append(f"Line {line_num}, column 7: Score out of range (0-1): '{parts[6]}'")
            except ValueError:
                errors.append(f"Line {line_num}, column 7: Not a valid number '{parts[6]}'")
            
            # Validate column 8 as + or -
            if parts[7] not in ['+', '-']:
                errors.append(f"Line {line_num}, column 8: Expected '+' or '-', got '{parts[7]}'")
    
    except Exception as e:
        errors.append(f"Error reading file: {str(e)}")
    
    return errors

def validate_filename(filename):
    """
    Validate filename format: firstname_lastname.txt
    For species names, we'll be more flexible
    """
    if not filename.endswith('.txt'):
        return False, "Filename must end with .txt"
    
    basename = filename[:-4]  # Remove .txt
    parts = basename.split('_')
    
    if len(parts) < 2:
        return False, "Filename should be in format: genus_species.txt (e.g., Adineta_vaga.txt)"
    
    # Check if first part and second part are non-empty
    if not parts[0] or not parts[1]:
        return False, "Genus and species names cannot be empty"
    
    # Check for valid species names (allow alphanumeric and optional numbers at end)
    import re
    # Genus: starts with capital letter, followed by lowercase letters
    if not re.match(r'^[A-Z][a-z]+$', parts[0]):
        return False, f"Genus name '{parts[0]}' should start with capital letter followed by lowercase"
    
    # Species: lowercase letters, optionally followed by numbers
    if not re.match(r'^[a-z]+[0-9]*$', parts[1]):
        return False, f"Species name '{parts[1]}' should be lowercase letters optionally followed by numbers"
    
    return True, ""

def validate_all_satsuma_files(satsuma_dir):
    """
    Validate all files in the Satsuma alignments directory.
    """
    all_errors = {}
    valid_files = []
    
    if not os.path.isdir(satsuma_dir):
        return False, {"directory_error": f"Directory not found: {satsuma_dir}"}, []
    
    files = [f for f in os.listdir(satsuma_dir) if f.endswith('.txt')]
    
    if not files:
        return False, {"no_files": "No .txt files found in directory"}, []
    
    for filename in files:
        # Validate filename format
        is_valid_name, name_error = validate_filename(filename)
        if not is_valid_name:
            all_errors[filename] = [f"Invalid filename format: {name_error}"]
            continue
        
        file_path = os.path.join(satsuma_dir, filename)
        errors = validate_satsuma_file(file_path)
        
        if errors:
            all_errors[filename] = errors
        else:
            valid_files.append(filename)
    
    return len(all_errors) == 0, all_errors, valid_files

def load_config(config_file):
    """Load YAML config file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"Error loading config file: {e}")
        sys.exit(1)

def print_summary(satsuma_dir, valid_files):
    """Print a summary of the validation results."""
    print("\n" + "="*60)
    print("SATSUMA ALIGNMENT FILES VALIDATION SUMMARY")
    print("="*60)
    
    total_size = 0
    total_lines = 0
    
    for filename in valid_files:
        file_path = os.path.join(satsuma_dir, filename)
        size_kb = os.path.getsize(file_path) / 1024
        total_size += size_kb
        
        # Count lines
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
            total_lines += line_count
        
        print(f"  ✓ {filename:<30} {line_count:>8} lines  {size_kb:>8.1f} KB")
    
    print("-"*60)
    print(f"  Total: {len(valid_files):>3} files   {total_lines:>8} lines  {total_size:>8.1f} KB")
    print("="*60)

def main():
    parser = argparse.ArgumentParser(description='Validate Satsuma alignment files')
    parser.add_argument('-c', '--config', default='run_sybr_config.yaml',
                       help='Configuration file (default: run_sybr_config.yaml)')
    parser.add_argument('-d', '--satsuma-dir',
                       help='Satsuma alignments directory (overrides config)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output (show all errors)')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Quiet mode (only show summary)')
    
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    
    # Check if synteny_processing is enabled
    if not config.get('run_stages', {}).get('synteny_processing', False):
        if not args.quiet:
            print("Synteny processing is disabled in config, skipping validation.")
        sys.exit(0)
    
    # Get satsuma directory
    if args.satsuma_dir:
        satsuma_dir = args.satsuma_dir
    else:
        satsuma_dir = config.get('satsuma_alignments')
    
    if not satsuma_dir:
        print("ERROR: Satsuma alignments directory not specified in config or command line.")
        sys.exit(1)
    
    # Validate files
    is_valid, errors, valid_files = validate_all_satsuma_files(satsuma_dir)
    
    if is_valid:
        if not args.quiet:
            print(f"✓ All Satsuma alignment files are valid.")
            print_summary(satsuma_dir, valid_files)
        sys.exit(0)
    else:
        if not args.quiet:
            print("✗ Validation failed for the following files:")
            
            # Show first few errors per file
            max_errors_per_file = 10 if not args.verbose else None
            
            for filename, file_errors in errors.items():
                print(f"\n  File: {filename}")
                
                if max_errors_per_file and len(file_errors) > max_errors_per_file:
                    for error in file_errors[:max_errors_per_file]:
                        print(f"    - {error}")
                    print(f"    - ... and {len(file_errors) - max_errors_per_file} more errors")
                else:
                    for error in file_errors:
                        print(f"    - {error}")
                
                # Show sample of the problematic lines
                if not args.verbose and file_errors:
                    try:
                        file_path = os.path.join(satsuma_dir, filename)
                        with open(file_path, 'r') as f:
                            first_line = f.readline().strip()
                            if first_line:
                                print(f"    Sample line: {first_line}")
                    except:
                        pass
        
        # Even if there are errors, show summary of valid files if any
        if valid_files:
            print(f"\n  Note: {len(valid_files)} file(s) passed validation: {', '.join(valid_files)}")
        
        sys.exit(1)

if __name__ == "__main__":
    main()
