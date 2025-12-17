import pandas as pd
import os
import argparse

def resolve_overlaps(input_file, output_file):
    # Resolve relative paths to absolute paths
    input_file = os.path.abspath(input_file)
    output_file = os.path.abspath(output_file)
    
    # Check if the input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The input file does not exist: {input_file}")

    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'start', 'end', 'name'])

    # Sort the DataFrame by chromosome and start
    df.sort_values(by=['chromosome', 'start'], inplace=True)

    # Initialize a list to store the resolved rows
    resolved_rows = []

    # Iterate over the rows
    i = 0
    while i < len(df) - 1:
        current_row = df.iloc[i].copy()
        next_row = df.iloc[i + 1].copy()

        if current_row['chromosome'] == next_row['chromosome']:
            if current_row['end'] < next_row['start']:
                resolved_rows.append(current_row)
                i += 1
            else:
                current_diff = current_row['end'] - current_row['start']
                next_diff = next_row['end'] - next_row['start']

                if current_diff > next_diff:
                    next_row['start'] = current_row['end'] + 1
                    resolved_rows.append(current_row)
                    df.iloc[i + 1] = next_row
                    i += 1
                else:
                    current_row['end'] = next_row['start'] - 1
                    resolved_rows.append(current_row)
                    i += 1
        else:
            resolved_rows.append(current_row)
            i += 1

    # Add the last row
    resolved_rows.append(df.iloc[-1])

    # Write the resolved rows to the output file
    resolved_df = pd.DataFrame(resolved_rows)
    resolved_df.to_csv(output_file, sep='\t', header=False, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Resolve overlapping intervals in a tab-delimited file.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to the output file")

    args = parser.parse_args()

    resolve_overlaps(args.input_file, args.output_file)
