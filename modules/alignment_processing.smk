# Alignment processing module

# Alignment processing rules (only if enabled)
if should_run_rule("chainNet_generation"):
    # Get paths from config with defaults
    SEQ_DIR = config["chainNet"].get("seq_dir", "seq")
    ALIGNMENTS_DIR = config["chainNet"].get("lastZ_alignments", "alignments2")
    OUTPUT_DIR = config["chainNet"].get("output_dir", ".")
    REFERENCE = config["reference_species"]

    # FASTA source for chainNet_generation and Ancestor_seq_recunstruction is
    # CHR_FASTA_DIR — files will be placed there by prepare_chr_fasta with the
    # Genus_ prefix stripped and >chr* headers added.
    # e.g.  Genus_sps1.fa (headers >1,>2…) → sps1.fa (headers >chr1,>chr2…)
    _FASTA_DIR = CHR_FASTA_DIR

    # Build species → EXPECTED chr-FASTA path map at parse time.
    # We derive paths from the ORIGINAL fasta_dir (always present) rather than
    # globbing CHR_FASTA_DIR (which doesn't exist yet at parse time).
    # The expected output path in CHR_FASTA_DIR uses the short stem (genus stripped).
    _ORIG_FASTA_DIR_AP = config["satsuma_align"]["fasta_dir"]
    _fasta_map = {}   # key → expected path in CHR_FASTA_DIR
    for _ext in (".fa", ".fna", ".fasta"):
        for _f in glob.glob(os.path.join(_ORIG_FASTA_DIR_AP, f"*{_ext}")):
            _full_stem = os.path.splitext(os.path.basename(_f))[0]   # e.g. "Genus_sps1"
            _short     = _full_stem.split("_", 1)[-1]                # e.g. "sps1"
            _chr_path  = os.path.join(CHR_FASTA_DIR, _short + ".fa") # expected output
            _parts     = _short.lower().split("_")
            _fasta_map[_short.lower()] = _chr_path
            for _part in _parts:
                _fasta_map.setdefault(_part, _chr_path)

    if not _fasta_map:
        raise ValueError(
            f"[alignment_processing] No FASTA files found in original fasta_dir: "
            f"{_ORIG_FASTA_DIR_AP}.\n"
            f"Check satsuma_align.fasta_dir in pipeline_paths.yaml."
        )

    def _fa_for_species(species):
        """Return the expected chr-FASTA path for a species wildcard.

        Built from original fasta_dir stems at parse time, so it works even
        before prepare_chr_fasta has run. The actual file will exist by the
        time fasta_to_2bit / generate_chrom_sizes execute (chr_ready guards it).
        """
        path = _fasta_map.get(species.lower())
        if path is None:
            raise ValueError(
                f"[alignment_processing] Cannot find FASTA for species '{species}'\n"
                f"Expected chr-FASTA dir: {_FASTA_DIR}\n"
                f"Tried key '{species.lower()}' against available short stems:\n"
                f"  {sorted(_fasta_map.keys())}\n"
                f"Ensure the original FASTA filename contains the species name\n"
                f"(e.g.  Genus_sps1.fa  →  sps1.fa  →  key 'sps1')."
            )
        return path

    rule fasta_to_2bit:
        input:
            fa        = lambda wildcards: _fa_for_species(wildcards.species),
            chr_ready = CHR_FASTA_DIR + "/.chr_fasta_done"
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
            axt_files = glob.glob(f"{ALIGNMENTS_DIR}/*.axt"),
            # When de-novo LastZ is active, wait for all .axt files to be
            # produced before cleaning. When using pre-computed files this
            # evaluates to [] and has no effect.
            lastz_done = lambda wildcards: (
                [os.path.join(LASTZ_ALIGNMENTS, ".lastz_alignment_done")]
                if config.get("run_stages", {}).get("run_lastz_alignment", False)
                else []
            )
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
            axt    = f"{ALIGNMENTS_DIR}/{{species}}.axt",
            target = lambda wildcards: f"{OUTPUT_DIR}/2bit/{REFERENCE}.2bit",
            query  = f"{OUTPUT_DIR}/2bit/{{species}}.2bit",
            # When run_lastz_alignment is ON the .axt file is produced by
            # run_lastz; the marker guarantees it exists before we start chaining.
            lastz_done = lambda wildcards: (
                [os.path.join(LASTZ_ALIGNMENTS, ".lastz_alignment_done")]
                if config.get("run_stages", {}).get("run_lastz_alignment", False)
                else []
            )
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
            fa = lambda wildcards: _fa_for_species(wildcards.species)
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
            reference_size=f"{OUTPUT_DIR}/fasize/{REFERENCE}.size",
            # When run_lastz_alignment is ON, get_alignment_species() may return []
            # at parse time (no .axt files yet).  The lastz_alignment_all marker
            # ensures all .axt → chain → chainSplit steps have completed before
            # chainPreNet starts, even if the expand() above resolved to nothing.
            lastz_done = lambda wildcards: (
                [os.path.join(LASTZ_ALIGNMENTS, ".lastz_alignment_done")]
                if config.get("run_stages", {}).get("run_lastz_alignment", False)
                else []
            )
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
