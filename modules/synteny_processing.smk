# Synteny processing module

if should_run_rule("synteny_processing"):
    rule sort:
        input:
            lambda wildcards: f"{SATSUMA_ALIGNMENTS}/{wildcards.sample}.txt"
        output:
            "{synteny_results}/{sample}_out/sorted_satsuma.txt"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort.pl")
        log:
            "{synteny_results}/{sample}_out/logs/sort.log"
#        message:
#            "Sorting {wildcards.sample} synteny results"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {input} ]; then
                echo "Error: Input file {input} not found" >&2
                exit 1
            fi
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi

            # Completion message
#            echo "Completed sorting {wildcards.sample} synteny results" >&2

            """


if should_run_rule("synteny_processing"):
    rule sort1:
        input:
            rules.sort.output
        output:
            "{synteny_results}/{sample}_out/sort1_result"  # Removed temp() to keep file
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1.pl")
        log:
            "{synteny_results}/{sample}_out/logs/sort1.log"
#        message:
#            "Completed secondary sort for {wildcards.sample}"
	
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort1.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule filtera:
        input:
            rules.sort1.output
        output:
            "{synteny_results}/{sample}_out/filtera_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filtera.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filtera.log"
#        message:
#            "Filtering results for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filtera.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule filter2a:
        input:
            rules.filtera.output
        output:
            "{synteny_results}/{sample}_out/filter2a_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter2a.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filter2a.log"
#        message:
#            "Running secondary filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter2a.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule filter:
        input:
            rules.filter2a.output
        output:
            "{synteny_results}/{sample}_out/filter_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filter.log"
#        message:
#            "Running final filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule filter3:
        input:
            rules.filter.output
        output:
            "{synteny_results}/{sample}_out/filter3_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter3.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filter3.log"
#        message:
#            "Running tertiary filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter3.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule filter4new:
        input:
            rules.filter3.output
        output:
            "{synteny_results}/{sample}_out/filter4new_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter4new.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filter4new.log"
#        message:
#            "Running final filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter4new.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """


if should_run_rule("synteny_processing"):
    rule change_reference_1:
        input:
            rules.filter4new.output
        output:
            "{synteny_results}/{sample}_out/chr_merge"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "change_reference_1.pl")
        log:
            "{synteny_results}/{sample}_out/logs/change_reference_1.log"
#        message:
#            "Changing reference for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: change_reference_1.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_sort:
        input:
            rules.change_reference_1.output
        output:
            "{synteny_results}/{sample}_out/chr_sort_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_sort.log"
#        message:
#            "Sorting chromosomes for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_sort1:
        input:
            rules.chr_sort.output
        output:
            "{synteny_results}/{sample}_out/chr_sort1_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_sort1.log"
#        message:
#            "Running secondary chromosome sort for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort1.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_filtera:
        input:
            rules.chr_sort1.output
        output:
            "{synteny_results}/{sample}_out/chr_filtera_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filtera.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_filtera.log"
#        message:
#            "Filtering chromosome data for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filtera.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_filter2a:
        input:
            rules.chr_filtera.output
        output:
            "{synteny_results}/{sample}_out/chr_filter2a_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter2a.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_filter2a.log"
#        message:
#            "Running secondary chromosome filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter2a.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_filter:
        input:
            rules.chr_filter2a.output
        output:
            "{synteny_results}/{sample}_out/chr_filter_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_filter.log"
#        message:
#            "Running final chromosome filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """


if should_run_rule("synteny_processing"):
    rule chr_filter3:
        input:
            rules.chr_filter.output
        output:
            "{synteny_results}/{sample}_out/chr_filter3_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter3.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_filter3.log"
#        message:
#            "Running tertiary chromosome filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter3.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_filter4new:
        input:
            rules.chr_filter3.output
        output:
            "{synteny_results}/{sample}_out/chr_filter4new_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter4new.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_filter4new.log"
#        message:
#            "Running final chromosome filtering for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: filter4new.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule chr_change_reference:
        input:
            rules.chr_filter4new.output
        output:
            "{synteny_results}/{sample}_out/ch1"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "change_reference_1.pl")
        log:
            "{synteny_results}/{sample}_out/logs/chr_change_reference.log"
#        message:
#            "Changing chromosome reference for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: change_reference_1.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule ch1_sort:
        input:
            rules.chr_change_reference.output
        output:
            "{synteny_results}/{sample}_out/ch1_sort_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort.pl")
        log:
            "{synteny_results}/{sample}_out/logs/ch1_sort.log"
#        message:
#            "Sorting chromosome 1 data for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule ch1_sort1:
        input:
            rules.ch1_sort.output
        output:
            "{synteny_results}/{sample}_out/ch1_sort1_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1.pl")
        log:
            "{synteny_results}/{sample}_out/logs/ch1_sort1.log"
#        message:
#            "Running secondary sort on chromosome 1 for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: sort1.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule st_format:
        input:
            rules.ch1_sort1.output
        output:
            "{synteny_results}/{sample}_out/st_input"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "st_format.pl")
        log:
            "{synteny_results}/{sample}_out/logs/st_format.log"
#        message:
#            "Formatting synteny track input for {wildcards.sample}"
        shell:
            """
            mkdir -p {params.outdir}/logs
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi
            
            if ! perl {params.script_path} {input} {output} 2>&1 | tee {log}; then
                echo "Error: st_format.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi
            
            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m sorting & filtering of \e[36;1m{wildcards.sample}\e[0m" >&2

            """

if should_run_rule("synteny_processing"):
    rule synteny_assign_100k:
        input:
            st_input = "{synteny_results}/{sample}_out/st_input"
        output:
            done_file = "{synteny_results}/synteny_out/100000/{sample}_out/synteny_assign_done"
        params:
            script_path = os.path.abspath("script_base/synteny_assign_v3_80_genes.pl"),
            ref_name = config["reference_name"],
            sample_name = lambda wildcards: wildcards.sample,
            target_dir = lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])),
            log_file = lambda wildcards, output: os.path.abspath(os.path.join(os.path.dirname(output[0]), "logs", "synteny_assign.log"))
        log:
            "{synteny_results}/synteny_out/100000/{sample}_out/logs/synteny_assign.log"
#        message:
#            "Running synteny assignment for {wildcards.sample} at 100k resolution"
        shell:
            """
            # STAGE 1: SETUP - Use absolute paths everywhere
            TARGET_DIR="{params.target_dir}"
            LOG_DIR="$(dirname "{params.log_file}")"
            SCRIPT_PATH="{params.script_path}"
            INPUT_FILE="{input.st_input}"
            
            # Create directories with verbose error checking
            if ! mkdir -p "$LOG_DIR"; then
                echo "FATAL: Failed to create log directory: $LOG_DIR" >&2
                exit 1
            fi
            if ! mkdir -p "$TARGET_DIR"; then
                echo "FATAL: Failed to create target directory: $TARGET_DIR" >&2
                exit 1
            fi
            
            # STAGE 2: FILE OPERATIONS
            if ! cp "$INPUT_FILE" "$TARGET_DIR/"; then
                echo "FATAL: Failed to copy input file to target directory" >&2
                exit 1
            fi
            
            # STAGE 3: EXECUTION
            cd "$TARGET_DIR" && \
            if ! perl "$SCRIPT_PATH" {params.ref_name} {params.sample_name} 100000 30000 st_input > "{params.log_file}" 2>&1; then
                echo "----------------------------------------" >&2
                echo "SYNTENY ASSIGNMENT FAILURE DETAILS" >&2
                echo "Working directory: $(pwd)" >&2
                echo "Command: perl $SCRIPT_PATH {params.ref_name} {params.sample_name} 100000 30000 st_input" >&2
                echo "Last 10 lines of log:" >&2
                tail -n 10 "{params.log_file}" >&2
                echo "----------------------------------------" >&2
                exit 1
            fi
            
            # STAGE 4: VERIFICATION
            if [ ! -f "$TARGET_DIR/synteny_assign_done" ]; then
                touch "$TARGET_DIR/synteny_assign_done"
            fi

#            echo "completed synteny tracking of {wildcards.sample}" in 100k resolution >&2
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m synteny tracking of \e[36;1m{wildcards.sample}\e[0m in 100k resolution" >&2
            """


if should_run_rule("synteny_processing"):
    rule synteny_assign_300k:
        input:
            st_input = "{synteny_results}/{sample}_out/st_input"
        output:
            done_file = "{synteny_results}/synteny_out/300000/{sample}_out/synteny_assign_done"
        params:
            script_path = os.path.abspath("script_base/synteny_assign_v3_80_genes.pl"),
            ref_name = config["reference_name"],
            sample_name = lambda wildcards: wildcards.sample,
            target_dir = lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])),
            log_file = lambda wildcards, output: os.path.abspath(os.path.join(os.path.dirname(output[0]), "logs", "synteny_assign.log"))
        log:
            "{synteny_results}/synteny_out/300000/{sample}_out/logs/synteny_assign.log"
#        message:
#            "Running synteny assignment for {wildcards.sample} at 300k resolution"
        shell:
            """
            # STAGE 1: SETUP
            TARGET_DIR="{params.target_dir}"
            LOG_DIR="$(dirname "{params.log_file}")"
            SCRIPT_PATH="{params.script_path}"
            INPUT_FILE="{input.st_input}"
            
            # Create directories
            mkdir -p "$LOG_DIR" && mkdir -p "$TARGET_DIR"
            
            # STAGE 2: FILE OPERATIONS
            cp "$INPUT_FILE" "$TARGET_DIR/"
            
            # STAGE 3: EXECUTION (with 300000 parameter)
            cd "$TARGET_DIR" && \
            perl "$SCRIPT_PATH" {params.ref_name} {params.sample_name} 300000 30000 st_input > "{params.log_file}" 2>&1
            
            # STAGE 4: VERIFICATION
            touch "$TARGET_DIR/synteny_assign_done"

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m synteny tracking of \e[36;1m{wildcards.sample}\e[0m in 300k resolution" >&2

            """


if should_run_rule("synteny_processing"):
    rule synteny_assign_500k:
        input:
            st_input = "{synteny_results}/{sample}_out/st_input"
        output:
            done_file = "{synteny_results}/synteny_out/500000/{sample}_out/synteny_assign_done"
        params:
            script_path = os.path.abspath("script_base/synteny_assign_v3_80_genes.pl"),
            ref_name = config["reference_name"],
            sample_name = lambda wildcards: wildcards.sample,
            target_dir = lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])),
            log_file = lambda wildcards, output: os.path.abspath(os.path.join(os.path.dirname(output[0]), "logs", "synteny_assign.log"))
        log:
            "{synteny_results}/synteny_out/500000/{sample}_out/logs/synteny_assign.log"
#        message:
#            "Running synteny assignment for {wildcards.sample} at 500k resolution"
        shell:
            """
            # STAGE 1: SETUP
            TARGET_DIR="{params.target_dir}"
            LOG_DIR="$(dirname "{params.log_file}")"
            SCRIPT_PATH="{params.script_path}"
            INPUT_FILE="{input.st_input}"
            
            # Create directories
            mkdir -p "$LOG_DIR" && mkdir -p "$TARGET_DIR"
            
            # STAGE 2: FILE OPERATIONS
            cp "$INPUT_FILE" "$TARGET_DIR/"
            
            # STAGE 3: EXECUTION (with 500000 parameter)
            cd "$TARGET_DIR" && \
            perl "$SCRIPT_PATH" {params.ref_name} {params.sample_name} 500000 30000 st_input > "{params.log_file}" 2>&1
            
            # STAGE 4: VERIFICATION
            touch "$TARGET_DIR/synteny_assign_done"

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m synteny tracking of \e[36;1m{wildcards.sample}\e[0m in 500k resolution" >&2

            """

