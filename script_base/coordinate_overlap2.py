import pandas as pd
import sys

def find_overlapping_rows(file1, file2, file3, chr_col1, start_col1, end_col1, chr_col2, start_col2, end_col2):
    # Read the input files without headers
    df1 = pd.read_csv(file1, sep='\t', header=None)
    df2 = pd.read_csv(file2, sep='\t', header=None)

    # Open the output file
    with open(file3, 'w') as f_out:
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

if __name__ == "__main__":
    # Example usage: python script.py file1.txt file2.txt file3.txt 0 1 2 0 1 2
    if len(sys.argv) != 10:
        print("Usage: python script.py <file1> <file2> <file3> <chr_col1> <start_col1> <end_col1> <chr_col2> <start_col2> <end_col2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    file3 = sys.argv[3]
    chr_col1 = int(sys.argv[4])
    start_col1 = int(sys.argv[5])
    end_col1 = int(sys.argv[6])
    chr_col2 = int(sys.argv[7])
    start_col2 = int(sys.argv[8])
    end_col2 = int(sys.argv[9])

    find_overlapping_rows(file1, file2, file3, chr_col1, start_col1, end_col1, chr_col2, start_col2, end_col2)
