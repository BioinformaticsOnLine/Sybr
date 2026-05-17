import os
import json
import argparse
import sys

def create_modified_config(source_config_path, target_folder, project_root):
    """Create a modified config file in the target folder with updated absolute DB paths"""
    try:
        # Load the original config
        with open(source_config_path) as f:
            config = json.load(f)

        # Resolve absolute location of target folder
        target_folder = os.path.abspath(target_folder)

        #
        # Update input files (unchanged)
        #
        config["input_files"].update({
            "background_genes_sb_1": os.path.join(target_folder, "background.txt"),
            "kegg_annotationTOgenes_sb_3": os.path.join(target_folder, "3kegg_annotationTOgenes.txt"),
            "genes_of_interest_sb_2": os.path.join(target_folder, "foreground.txt")
        })

        #
        # ⬇️ Insert absolute database paths using project_root
        #
        if "db" not in config:
            config["db"] = {}

        config["db"].update({
            "term2gene": os.path.join(project_root, "tools/getENRICH/db/term2gene.txt"),
            "term2name": os.path.join(project_root, "tools/getENRICH/db/term2name.txt")
        })

        #
        # Output configuration
        #
        config["output_files"]["outdir"] = os.path.join(target_folder, "enrichment_result")
        config["output_files"]["enrichment_KEGG_results_csv"] = "enrichment_result.csv"

        # Create directory if needed
        os.makedirs(target_folder, exist_ok=True)

        # Save final config file
        target_config_path = os.path.join(target_folder, os.path.basename(source_config_path))
        with open(target_config_path, 'w') as f:
            json.dump(config, f, indent=2)

        print(f"Successfully created config in {target_config_path}")
        return True

    except Exception as e:
        print(f"Error processing {target_folder}: {str(e)}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Clone config file to multiple folders with updated DB paths",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("source_config", help="Path to source config file")
    parser.add_argument("target_folders", nargs="+", help="List of target folders")

    # New argument
    parser.add_argument(
        "--project-root",
        required=True,
        help="Absolute path to Snakemake project root"
    )

    parser.add_argument("--verbose", "-v", action="store_true", help="Show detailed output")

    args = parser.parse_args()

    if not os.path.isfile(args.source_config):
        print(f"Error: Source config file not found: {args.source_config}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.project_root):
        print(f"Error: Project root directory not found: {args.project_root}", file=sys.stderr)
        sys.exit(1)

    success_count = 0

    for folder in args.target_folders:
        if args.verbose:
            print(f"Processing folder: {folder}")

        if create_modified_config(
            args.source_config,
            folder,
            args.project_root
        ):
            success_count += 1

    print(f"\nSummary: Successfully created {success_count} config files out of {len(args.target_folders)}")

    if success_count < len(args.target_folders):
        sys.exit(1)


if __name__ == "__main__":
    main()
