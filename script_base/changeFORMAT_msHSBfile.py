import argparse

def convert_coordinates(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_chromosome = None
        for line in infile:
            line = line.strip()
            if line.startswith("Chromosome"):
                # Extract the chromosome number
                current_chromosome = line.split()[1].strip(':')
            elif line and current_chromosome:
                # Check if the line contains "Common overlapping"
                if "Common overlapping" in line:
                    continue
                # Extract start and end coordinates
                start, end = line.split()[:2]
                # Write the formatted output
                outfile.write(f"{current_chromosome}\t{start}\t{end}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert chromosome coordinates to a tab-separated format")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("output_file", help="Output file path")
    args = parser.parse_args()

    convert_coordinates(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
