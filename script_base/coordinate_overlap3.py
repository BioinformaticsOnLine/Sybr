import pandas as pd
import sys
import os
from glob import glob

def find_overlapping_rows(file1, file2, output_file, chr_col1, start_col1, end_col1, chr_col2, start_col2, end_col2):
    # Read the input files without headers
    df1 = pd.read_csv(file1, sep='\t', header=None)
    df2 = pd.read_csv(file2, sep='\t', header=None)

    # Open the output file
    with open(output_file, 'w') as f_out:
        # Write headers (assuming they are the same for both files)
        headers = [f'file1_{i}' for i in range(len(df1.columns))] + [f'file2_{i}' for i in range(len(df2.columns))]
        f_out.write('\t'.join(headers) + '\n')

        # Iterate over rows in df1
        for _, row1 in df1.iterrows():
            # Filter df2 for matching chromosome
            df2_filtered = df2[df2.iloc[:, chr_col2] == row1.iloc[chr_col1]]

            # Check for overlapping start and end positions
            for _, row2 in df2_filtered.iterrows():
                if (row1.iloc[start_col1] <= row2.iloc[end_col2]) and (row1.iloc[end_col1] >= row2.iloc[start_col2]):
                    # Write the overlapping rows to the output file
                    f_out.write('\t'.join(map(str, row1)) + '\t' + '\t'.join(map(str, row2)) + '\n')

def process_lineage_folders(parent_folder, file2, chr_col1, start_col1, end_col1, chr_col2, start_col2, end_col2):
    # Find all subfolders with '_lineage' suffix
    lineage_folders = glob(os.path.join(parent_folder, '*_lineage'))
    
    if not lineage_folders:
        print(f"No folders with '_lineage' suffix found in {parent_folder}")
        return
    
    for folder in lineage_folders:
        # Find all .txt files in the lineage folder
        txt_files = glob(os.path.join(folder, '*.txt'))
        
        if not txt_files:
            print(f"No .txt files found in {folder}")
            continue
            
        for file1 in txt_files:
            # Create output file path in the same folder
            output_file = os.path.join(folder, 'EBRs_NCBI_genes_overlap.txt')
            
            print(f"Processing {file1} with {file2} -> {output_file}")
            find_overlapping_rows(file1, file2, output_file, 
                                chr_col1, start_col1, end_col1, 
                                chr_col2, start_col2, end_col2)

if __name__ == "__main__":
    # Example usage: python script.py parent_folder file2.txt 0 1 2 0 1 2
    if len(sys.argv) != 9:
        print("Usage: python script.py <parent_folder> <file2> <chr_col1> <start_col1> <end_col1> <chr_col2> <start_col2> <end_col2>")
        print("Note: file2 is a single file that will be compared against all *.txt files in each *_lineage subfolder")
        sys.exit(1)

    parent_folder = sys.argv[1]
    file2 = sys.argv[2]
    chr_col1 = int(sys.argv[3])
    start_col1 = int(sys.argv[4])
    end_col1 = int(sys.argv[5])
    chr_col2 = int(sys.argv[6])
    start_col2 = int(sys.argv[7])
    end_col2 = int(sys.argv[8])

    process_lineage_folders(parent_folder, file2, 
                          chr_col1, start_col1, end_col1, 
                          chr_col2, start_col2, end_col2)
