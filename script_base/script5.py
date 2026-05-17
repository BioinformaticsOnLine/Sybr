#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd

THRESHOLD = 0.90


def extract_filtered_dataframe(input_file):
    rows = []

    with open(input_file) as infile:
        header = infile.readline().rstrip("\n").split("\t")

        try:
            chr_idx = header.index("Chromosome")
            brk_idx = header.index("Narrowest_brk")
            prob_idx = header.index("First_Probability")
        except ValueError as e:
            sys.exit(f"Missing required column: {e}")

        for line in infile:
            fields = line.rstrip("\n").split("\t")

            last_entry = fields[prob_idx].split()[-1]

            try:
                species, value_str = last_entry.rsplit(":", 1)
                probability = float(value_str)
            except ValueError:
                continue

            if probability < THRESHOLD:
                continue

            try:
                start, end = fields[brk_idx].split("<->")
                start = int(start)
                end = int(end)
            except ValueError:
                continue

            rows.append([
                fields[chr_idx],
                start,
                end,
                species,
                probability
            ])

    return pd.DataFrame(
        rows,
        columns=["chromosome", "start", "end", "species", "probability"]
    )


def resolve_overlaps(df):
    df = df.sort_values(by=["chromosome", "start"]).reset_index(drop=True)
    resolved_rows = []

    i = 0
    while i < len(df) - 1:
        current = df.iloc[i].copy()
        next_row = df.iloc[i + 1].copy()

        if current["chromosome"] == next_row["chromosome"]:
            if current["end"] < next_row["start"]:
                resolved_rows.append(current)
                i += 1
            else:
                current_len = current["end"] - current["start"]
                next_len = next_row["end"] - next_row["start"]

                if current_len >= next_len:
                    next_row["start"] = current["end"] + 1
                    resolved_rows.append(current)
                    df.iloc[i + 1] = next_row
                    i += 1
                else:
                    current["end"] = next_row["start"] - 1
                    resolved_rows.append(current)
                    i += 1
        else:
            resolved_rows.append(current)
            i += 1

    resolved_rows.append(df.iloc[-1])
    return pd.DataFrame(resolved_rows)


def main(input_file, output_file):
    input_file = os.path.abspath(input_file)
    output_file = os.path.abspath(output_file)

    if not os.path.isfile(input_file):
        sys.exit(f"Input file does not exist: {input_file}")

    df = extract_filtered_dataframe(input_file)
    resolved_df = resolve_overlaps(df)

    resolved_df.to_csv(output_file, sep="\t", header=False, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract EBA results, filter probability ≥0.90, and remove overlaps."
    )
    parser.add_argument("input_file", help="Input TSV file")
    parser.add_argument("output_file", help="Output TSV file")

    args = parser.parse_args()
    main(args.input_file, args.output_file)
