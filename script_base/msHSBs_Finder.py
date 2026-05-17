import csv
import argparse
import os

def process_input_file(input_file_path, output_file_path):
    """
    Process the input file to extract specific columns and write them to an output file.
    This function mimics the behavior of an awk command in Python.
    """
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            # Extract columns 2, 3, 4, and 9 (0-based index: 1, 2, 3, 8)
            writer.writerow([row[1], row[2], row[3], row[8]])

def assign_ids(input_file_path, output_file_path):
    """
    Assign unique IDs to names in the fourth column of the input file.
    The IDs are appended as a new third column, and the modified data is written to an output file.
    """
    name_to_id = {}  # Dictionary to map names to unique IDs
    current_id = 1   # Counter for assigning new IDs

    with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            name = row[3]  # Get the name from the fourth column
            if name not in name_to_id:
                # Assign a new ID if the name is not already in the dictionary
                name_to_id[name] = current_id
                current_id += 1
            # Insert the ID as the third column and shift the name to the fourth column
            row.insert(3, name_to_id[name])
            row.pop(4)  # Remove the original name column
            writer.writerow(row)

def find_common_overlapping_intervals(file_path, output_file_path):
    """
    Find common overlapping intervals for each chromosome and write the results to an output file.
    """
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        # Read coordinates as tuples of (chromosome, start, end, ID)
        coordinates = [(int(row[0]), int(row[1]), int(row[2]), int(row[3])) for row in reader]

    # Dictionary to group coordinates by chromosome and ID
    chromosome_to_id_to_coordinates = {}
    for coordinate in coordinates:
        chromosome = coordinate[0]
        id = coordinate[3]
        if chromosome not in chromosome_to_id_to_coordinates:
            chromosome_to_id_to_coordinates[chromosome] = {}
        if id not in chromosome_to_id_to_coordinates[chromosome]:
            chromosome_to_id_to_coordinates[chromosome][id] = []
        chromosome_to_id_to_coordinates[chromosome][id].append((coordinate[1], coordinate[2]))

    with open(output_file_path, 'w') as outfile:
        for chromosome, id_to_coordinates in chromosome_to_id_to_coordinates.items():
            outfile.write(f"Chromosome {chromosome}:\n")
            id_intervals = {}
            for id, coord_list in id_to_coordinates.items():
                coord_list.sort(key=lambda x: x[0])  # Sort intervals by start position
                merged_intervals = [coord_list[0]]
                for current_interval in coord_list[1:]:
                    last_merged_interval = merged_intervals[-1]
                    # Merge overlapping intervals
                    if current_interval[0] <= last_merged_interval[1]:
                        merged_intervals[-1] = (last_merged_interval[0], max(last_merged_interval[1], current_interval[1]))
                    else:
                        merged_intervals.append(current_interval)
                id_intervals[id] = merged_intervals

            # Find common overlapping intervals across all IDs
            common_overlapping_intervals = id_intervals[list(id_intervals.keys())[0]]
            for id, intervals in id_intervals.items():
                if id != list(id_intervals.keys())[0]:
                    new_common_intervals = []
                    for common_interval in common_overlapping_intervals:
                        for interval in intervals:
                            # Check for overlap and update common intervals
                            if interval[0] <= common_interval[1] and interval[1] >= common_interval[0]:
                                new_common_intervals.append((max(common_interval[0], interval[0]), min(common_interval[1], interval[1])))
                    common_overlapping_intervals = new_common_intervals

            outfile.write("Common overlapping intervals:\n")
            for interval in common_overlapping_intervals:
                outfile.write(f"{interval[0]} {interval[1]}\n")
            outfile.write("\n")

def main():
    """
    Main function to parse command-line arguments and execute the processing steps.
    """
    parser = argparse.ArgumentParser(description="Process input file, assign IDs, and find common overlapping intervals")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("output_file", help="Output file path")
    parser.add_argument("-k", "--keep_intermediate", action="store_true", help="Keep intermediate files")
    args = parser.parse_args()

    # Define intermediate file paths
    intermediate_file1 = "intermediate_file1.tsv"
    intermediate_file2 = "intermediate_file2.tsv"

    # Process the input file and assign IDs
    process_input_file(args.input_file, intermediate_file1)
    assign_ids(intermediate_file1, intermediate_file2)
    # Find common overlapping intervals and write to the output file
    find_common_overlapping_intervals(intermediate_file2, args.output_file)

    # Remove intermediate files unless the keep_intermediate flag is set
    if not args.keep_intermediate:
        os.remove(intermediate_file1)
        os.remove(intermediate_file2)

if __name__ == "__main__":
    main()
