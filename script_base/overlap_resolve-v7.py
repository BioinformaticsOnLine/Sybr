import argparse
import os
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Resolve overlaps in alignment data and report modifications.')
    parser.add_argument('parent_dir', type=str, help='Path to the parent directory.')
    return parser.parse_args()

def read_input_file(input_file):
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ['ref_species', 'ref_chromosome', 'ref_start', 'ref_end',
                  'query_chromosome', 'query_start', 'query_end', 'orientation',
                  'query_species', 'assembly_level']
    return df

def sort_input_file(df):
    df.sort_values(by=['ref_chromosome', 'ref_start'], ascending=[True, True], inplace=True)
    return df

def resolve_overlaps(df):
    df = sort_input_file(df)
    resolved = []
    current_chromosome = None

    for _, row in df.iterrows():
        if not resolved or row['ref_chromosome'] != current_chromosome:
            resolved.append(row)
            current_chromosome = row['ref_chromosome']
        else:
            last = resolved[-1]
            if row['ref_start'] <= last['ref_end']:
                if row['ref_end'] <= last['ref_end']:
                    print(f"Entire overlap detected: Skipping row for block in {row['ref_chromosome']} with start {row['ref_start']} and end {row['ref_end']}")
                    continue
                else:
                    print(f"Overlap detected: Adjusting start from {row['ref_start']} to {last['ref_end'] + 1} for block in {row['ref_chromosome']}")
                    row['ref_start'] = last['ref_end'] + 1
            resolved.append(row)
    return pd.DataFrame(resolved)

def write_output_file(df, output_file):
    df.to_csv(output_file, sep='\t', index=False, header=False)

def main():
    args = parse_arguments()
    parent_dir = args.parent_dir

    # Loop through all subdirectories in the parent directory
    for d in os.listdir(parent_dir):
        sub_dir = os.path.join(parent_dir, d)
        if os.path.isdir(sub_dir):
            for fname in os.listdir(sub_dir):
                if fname.endswith('_final.txt'):
                    input_file = os.path.join(sub_dir, fname)
                    output_file = os.path.join(sub_dir, fname.replace('_final.txt', '_nonoverlap.txt'))

                    if not os.path.isfile(input_file):
                        print(f"Input file {input_file} not found. Skipping...")
                        continue

                    print(f"Processing file: {input_file}")
                    df = read_input_file(input_file)
                    resolved_df = resolve_overlaps(df)
                    write_output_file(resolved_df, output_file)
                    print(f"Written resolved data to: {output_file}")

if __name__ == '__main__':
    main()
