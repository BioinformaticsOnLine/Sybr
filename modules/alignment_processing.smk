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
            """
            mkdir -p {params.output_dir}/2bit
            echo -e "\e[33;1mConverting\e[0m \e[36m{wildcards.species}\e[0m FASTA to 2bit format" >&2
            faToTwoBit {input.fa} {output}
            echo -e "\e[32;1m✔\e[0m \e[33;1mConverted\e[0m \e[36m{wildcards.species}\e[0m FASTA to 2bit format" >&2
            """

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
            echo -e "\e[33;1mCleaning\e[0m AXT files" >&2
            for file in {input.axt_files}; do
                sed -i '/^#/d' "$file"
            done
            mkdir -p {params.output_dir}
            touch {output}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCleaned\e[0m AXT files" >&2
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
            """
            mkdir -p {params.output_dir}/chain
            echo -e "\e[33;1mConverting\e[0m \e[36m{wildcards.species}\e[0m AXT to chain format" >&2
            axtChain -linearGap=medium {input.axt} {input.target} {input.query} {output}
            echo -e "\e[32;1m✔\e[0m \e[33;1mConverted\e[0m \e[36m{wildcards.species}\e[0m AXT to chain format" >&2
            """

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
            echo -e "\e[33;1mCleaning and splitting\e[0m \e[36m{wildcards.species}\e[0m chain file" >&2
            sed -i '/^#/d' {input.chain}
            chainSplit {output} {input.chain}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCleaned and split\e[0m \e[36m{wildcards.species}\e[0m chain file" >&2
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
            """
            mkdir -p {params.output_dir}/fasize
            echo -e "\e[33;1mGenerating chromosome sizes for\e[0m \e[36m{wildcards.species}\e[0m" >&2
            faSize {input.fa} -detailed > {output}
            echo -e "\e[32;1m✔\e[0m \e[33;1mGenerated chromosome sizes for\e[0m \e[36m{wildcards.species}\e[0m" >&2
            """

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
            mkdir -p {params.output_dir}
            echo -e "\e[33;1mRunning chainPreNet processing\e[0m" >&2
            {params.script_path} {params.chains_folder} {params.vaga_size_file} {params.size_files_dir} {params.output_dir}/chainPreNet
            touch {output[1]}  # Create completion marker
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m chainPreNet processing" >&2
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
            """
            mkdir -p {params.output_dir}
            echo -e "\e[33;1mRunning chainNet processing\e[0m" >&2
            {params.script_path} {params.chainPreNet_folder} {params.vaga_size_file} {params.size_files_dir} {params.output_dir}/chainNet
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m chainNet processing" >&2
            """

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
            echo -e "\e[33;1mRunning netSyntenic processing\e[0m" >&2
            {params.script_path} {params.output_dir}/chainNet {params.output_dir}/net
            touch {output[1]}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m netSyntenic processing" >&2
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
            echo -e "\e[33;1mReorganizing data structure for\e[0m \e[36m{params.primary_species}\e[0m" >&2
            {params.script_path} {params.primary_species} {params.output_dir}
            touch {output[1]}
            echo -e "\e[32;1m✔\e[0m \e[33;1mReorganized data structure for\e[0m \e[36m{params.primary_species}\e[0m" >&2
            """
