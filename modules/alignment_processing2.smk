# Alignment processing module

# Alignment processing rules (only if enabled)
if should_run_rule("chainNet_generation"):
    # Get paths from config with defaults
    SEQ_DIR = config["chainNet"].get("seq_dir", "seq")
    ALIGNMENTS_DIR = config["chainNet"].get("lastZ_alignments", "alignments2")
    OUTPUT_DIR = config["chainNet"].get("output_dir", ".")
    REFERENCE = config["reference_species"]

    rule fasta_to_2bit:
        input:
            fa=lambda wildcards: next(
                (f"{SEQ_DIR}/{wildcards.species}{ext}"
                 for ext in ['.fa', '.fna', '.fasta']
                 if os.path.exists(f"{SEQ_DIR}/{wildcards.species}{ext}")),
                None
            )
        output:
            f"{OUTPUT_DIR}/2bit/{{species}}.2bit"
        params:
            output_dir = OUTPUT_DIR
        shell:
            "mkdir -p {params.output_dir}/2bit && faToTwoBit {input.fa} {output}"

    rule clean_axt:
        input:
            axt_files=glob.glob(f"{ALIGNMENTS_DIR}/*.axt")
        output:
            f"{OUTPUT_DIR}/clean_axt.done"
        params:
            alignments_dir = ALIGNMENTS_DIR,
            output_dir = OUTPUT_DIR
        shell:
            """
            for file in {input.axt_files}; do
                sed -i '/^#/d' "$file"
            done
            mkdir -p {params.output_dir}
            touch {output}
            """

    rule axt_to_chain:
        input:
            axt = f"{ALIGNMENTS_DIR}/{{species}}.axt",
            target = lambda wildcards: f"{OUTPUT_DIR}/2bit/{REFERENCE}.2bit",
            query = f"{OUTPUT_DIR}/2bit/{{species}}.2bit"
        output:
            f"{OUTPUT_DIR}/chain/{{species}}.chain"
        params:
            output_dir = OUTPUT_DIR
        shell:
            "mkdir -p {params.output_dir}/chain && axtChain -linearGap=medium {input.axt} {input.target} {input.query} {output}"

    rule clean_and_split_chain:
        input:
            chain = f"{OUTPUT_DIR}/chain/{{species}}.chain"
        output:
            directory(f"{OUTPUT_DIR}/chain/{{species}}-chainSplit")
        params:
            output_dir = OUTPUT_DIR
        shell:
            """
            mkdir -p {params.output_dir}/chain
            sed -i '/^#/d' {input.chain}
            chainSplit {output} {input.chain}
            """

    rule generate_chrom_sizes:
        input:
            fa=lambda wildcards: next(
                (f"{SEQ_DIR}/{wildcards.species}{ext}"
                for ext in ['.fa', '.fna', '.fasta']
                if os.path.exists(f"{SEQ_DIR}/{wildcards.species}{ext}")),
                None
            )
        output:
            f"{OUTPUT_DIR}/fasize/{{species}}.size"
        params:
            output_dir = OUTPUT_DIR
        shell:
            "mkdir -p {params.output_dir}/fasize && faSize {input.fa} -detailed > {output}"

    rule chainPreNet_processing:
        input:
            chain_dirs=expand(f"{OUTPUT_DIR}/chain/{{species}}-chainSplit", species=get_alignment_species()),
            all_size_files=expand(f"{OUTPUT_DIR}/fasize/{{species}}.size", species=get_fasta_species()),
            reference_size=f"{OUTPUT_DIR}/fasize/{REFERENCE}.size"
        output:
            directory(f"{OUTPUT_DIR}/chainPreNet"),
            touch(f"{OUTPUT_DIR}/chainPreNet/processing_done.txt")  # Completion marker
        params:
            output_dir = OUTPUT_DIR,
            chains_folder=f"{OUTPUT_DIR}/chain",
            vaga_size_file=f"{REFERENCE}.size",
            size_files_dir=f"{OUTPUT_DIR}/fasize",
            script_path = os.path.join(SCRIPTS, "chainPreNet-script3.sh")
        shell:
            """
            mkdir -p {params.output_dir} && {params.script_path} {params.chains_folder} {params.vaga_size_file} {params.size_files_dir} {params.output_dir}/chainPreNet
            touch {output[1]}  # Create completion marker
            """

    rule chainNet_processing:
        input:
            chainPreNet_dir=f"{OUTPUT_DIR}/chainPreNet",
            reference_size=f"{OUTPUT_DIR}/fasize/{REFERENCE}.size",
            all_size_files=expand(f"{OUTPUT_DIR}/fasize/{{species}}.size", species=get_fasta_species())
        output:
            touch(f"{OUTPUT_DIR}/netFiles.done")
        params:
            output_dir = OUTPUT_DIR,
            chainPreNet_folder=f"{OUTPUT_DIR}/chainPreNet",
            vaga_size_file=f"{REFERENCE}.size",
            size_files_dir=f"{OUTPUT_DIR}/fasize",
            script_path = os.path.join(SCRIPTS, "chainNet-script2.sh")
        shell:
            "mkdir -p {params.output_dir} && {params.script_path} {params.chainPreNet_folder} {params.vaga_size_file} {params.size_files_dir} {params.output_dir}/chainNet"

    rule netSyntenic_processing:
        input:
            f"{OUTPUT_DIR}/netFiles.done",
            chainPreNet_dir=f"{OUTPUT_DIR}/chainPreNet"
        output:
            directory(f"{OUTPUT_DIR}/net"),
            touch(f"{OUTPUT_DIR}/netSyntenic.done")
        params:
            output_dir = OUTPUT_DIR,
            script_path = os.path.join(SCRIPTS, "netSyntenic-script.sh")
        shell:
            """
            mkdir -p {params.output_dir}/net
            {params.script_path} {params.output_dir}/chainNet {params.output_dir}/net
            touch {output[1]}
            """

    rule reorganize_data:
        input:
            f"{OUTPUT_DIR}/netSyntenic.done",
            chainPreNet_dir=f"{OUTPUT_DIR}/chainPreNet",
            net_dir=f"{OUTPUT_DIR}/net"
        output:
            directory(f"{OUTPUT_DIR}/data/{REFERENCE}"),
            f"{OUTPUT_DIR}/data/{REFERENCE}/reorganized.done"
        params:
            output_dir = OUTPUT_DIR,
            primary_species=REFERENCE,
            script_path = os.path.join(SCRIPTS, "input-file-structure3.sh")
        shell:
            """
            mkdir -p {params.output_dir}/data/{params.primary_species}
            {params.script_path} {params.primary_species} {params.output_dir}
            touch {output[1]}
            """
