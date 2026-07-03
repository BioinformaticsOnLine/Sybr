# Synteny processing module

# Import needed at module level
import os
import re as _re

# ── Pin {synteny_results} and {window_size} wildcards ───────────────────────
wildcard_constraints:
    synteny_results = _re.escape(SYNTENY_RESULTS),
    window_size     = r"\d+"

# Get synteny parameters from config at module level
if "synteny_params" in config:
    WINDOW_SIZES = config["synteny_params"].get("window_sizes", [100000, 300000, 500000])
    STEP_SIZE = config["synteny_params"].get("step_size", 30000)
    CUSTOM_RESOLUTIONS = config["synteny_params"].get("custom_resolutions", False)
else:
    WINDOW_SIZES = [100000, 300000, 500000]
    STEP_SIZE = 30000
    CUSTOM_RESOLUTIONS = False

if should_run_rule("synteny_processing"):
    rule sort:
        input:
            alignment = lambda wildcards: f"{SATSUMA_ALIGNMENTS}/{wildcards.sample}.txt",
            # When de-novo alignment is active, wait for all alignments to finish
            # before starting the filter chain. When using pre-computed files
            # this evaluates to [] and has no effect.
            satsuma_done = lambda wildcards: (
                [os.path.join(config["satsuma_align"]["output_dir"], ".satsuma_alignment_done")]
                if config.get("run_stages", {}).get("run_satsuma_alignment", False)
                else []
            )
        output:
            "{synteny_results}/{sample}_out/sorted_satsuma.txt"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/logs/sort.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {input.alignment} ]; then
                echo "Error: Input file {input.alignment} not found" >&2
                exit 1
            fi

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input.alignment} {output} 2>&1 | tee {log}; then
                echo "Error: sort.pl failed for sample {wildcards.sample}" >&2
                exit 1
            fi

            if [ ! -f {output} ]; then
                echo "Error: Output file {output} was not created" >&2
                exit 1
            fi
            """

if should_run_rule("synteny_processing"):
    rule sort1:
        input:
            rules.sort.output
        output:
            "{synteny_results}/{sample}_out/sort1_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/logs/sort1.log"
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
            script_path = os.path.join(SCRIPTS, "filtera_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filtera.log"
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
            script_path = os.path.join(SCRIPTS, "filter2a_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/logs/filter2a.log"
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

# ── From here every rule is per-resolution: outputs live under {window_size}/ ──

if should_run_rule("synteny_processing"):
    rule filter:
        input:
            rules.filter2a.output
        output:
            "{synteny_results}/{sample}_out/{window_size}/filter_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter_NJ.pl"),
            window_size = lambda wildcards: int(wildcards.window_size)
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/filter.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.window_size} 2>&1 | tee {log}; then
                echo "Error: filter_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/filter3_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter3_NJ.pl"),
            window_size = lambda wildcards: int(wildcards.window_size)
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/filter3.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.window_size} 2>&1 | tee {log}; then
                echo "Error: filter3_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/filter4new_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter4new_NJ.pl"),
            step_size = STEP_SIZE
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/filter4new.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.step_size} 2>&1 | tee {log}; then
                echo "Error: filter4new_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/chr_merge"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "change_reference_1_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/change_reference_1.log"
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
            "{synteny_results}/{sample}_out/{window_size}/chr_sort_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_sort.log"
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
            "{synteny_results}/{sample}_out/{window_size}/chr_sort1_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_sort1.log"
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
            "{synteny_results}/{sample}_out/{window_size}/chr_filtera_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filtera_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_filtera.log"
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
            "{synteny_results}/{sample}_out/{window_size}/chr_filter2a_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter2a_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_filter2a.log"
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
            "{synteny_results}/{sample}_out/{window_size}/chr_filter_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter_NJ.pl"),
            window_size = lambda wildcards: int(wildcards.window_size)
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_filter.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.window_size} 2>&1 | tee {log}; then
                echo "Error: filter_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/chr_filter3_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter3_NJ.pl"),
            window_size = lambda wildcards: int(wildcards.window_size)
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_filter3.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.window_size} 2>&1 | tee {log}; then
                echo "Error: filter3_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/chr_filter4new_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "filter4new_NJ.pl"),
            step_size = STEP_SIZE
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_filter4new.log"
        shell:
            """
            mkdir -p {params.outdir}/logs

            if [ ! -f {params.script_path} ]; then
                echo "Error: Script {params.script_path} not found" >&2
                exit 1
            fi

            if ! perl {params.script_path} {input} {output} {params.step_size} 2>&1 | tee {log}; then
                echo "Error: filter4new_NJ.pl failed for sample {wildcards.sample} resolution {wildcards.window_size}" >&2
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
            "{synteny_results}/{sample}_out/{window_size}/ch1"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "change_reference_1_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/chr_change_reference.log"
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
            "{synteny_results}/{sample}_out/{window_size}/ch1_sort_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/ch1_sort.log"
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
            "{synteny_results}/{sample}_out/{window_size}/ch1_sort1_result"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "sort1_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/ch1_sort1.log"
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
            "{synteny_results}/{sample}_out/{window_size}/st_input"
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            script_path = os.path.join(SCRIPTS, "st_format_NJ.pl")
        log:
            "{synteny_results}/{sample}_out/{window_size}/logs/st_format.log"
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

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m sorting & filtering of \e[36;1m{wildcards.sample}\e[0m at resolution \e[36;1m{wildcards.window_size}\e[0m" >&2
            """

# ============================================================================
# DYNAMIC RULES FOR DIFFERENT WINDOW SIZES
# ============================================================================

# Create a function to generate rules for each window size
def create_synteny_assign_rules():
    """Create synteny assignment rules for each window size"""

    for window_size in WINDOW_SIZES:
        window_kb = window_size // 1000

        # Define the rule with a unique name
        rule_name = f"synteny_assign_{window_kb}k"

        # Create a rule dictionary
        rule_dict = {
            'input': lambda wildcards: f"{{synteny_results}}/{{wildcards.sample}}_out/{window_size}/st_input",
            'output': lambda wildcards: f"{{synteny_results}}/synteny_results/{window_size}/{{wildcards.sample}}_out/synteny_assign_done",
            'params': lambda wildcards: {
                'script_path': os.path.abspath("script_base/synteny_assign_v3_80_genes.pl"),
                'ref_name': config["reference_name"],
                'sample_name': wildcards.sample,
                'window_size': window_size,
                'step_size': STEP_SIZE,
                'target_dir': os.path.dirname(f"{{synteny_results}}/synteny_results/{window_size}/{{wildcards.sample}}_out/synteny_assign_done"),
                'log_file': os.path.join(
                    os.path.dirname(f"{{synteny_results}}/synteny_results/{window_size}/{{wildcards.sample}}_out/synteny_assign_done"),
                    "logs", "synteny_assign.log"
                )
            },
            'log': lambda wildcards: f"{{synteny_results}}/synteny_results/{window_size}/{{wildcards.sample}}_out/logs/synteny_assign.log",
            'shell': f"""
                # STAGE 1: SETUP - Use absolute paths everywhere
                TARGET_DIR="{{params.target_dir}}"
                LOG_DIR="$(dirname "{{params.log_file}}")"
                SCRIPT_PATH="{{params.script_path}}"
                INPUT_FILE="{{input}}"

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
                cd "$TARGET_DIR" && \\
                if ! perl "$SCRIPT_PATH" {{params.ref_name}} {{params.sample_name}} {window_size} {STEP_SIZE} st_input > "{{params.log_file}}" 2>&1; then
                    echo "----------------------------------------" >&2
                    echo "SYNTENY ASSIGNMENT FAILURE DETAILS" >&2
                    echo "Working directory: $(pwd)" >&2
                    echo "Command: perl $SCRIPT_PATH {{params.ref_name}} {{params.sample_name}} {window_size} {STEP_SIZE} st_input" >&2
                    echo "Last 10 lines of log:" >&2
                    tail -n 10 "{{params.log_file}}" >&2
                    echo "----------------------------------------" >&2
                    exit 1
                fi

                # STAGE 4: VERIFICATION
                if [ ! -f "$TARGET_DIR/synteny_assign_done" ]; then
                    touch "$TARGET_DIR/synteny_assign_done"
                fi

                echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m synteny tracking of \\e[36;1m{{wildcards.sample}}\\e[0m in {window_size} resolution" >&2
            """
        }

        # Create the rule using a different approach - define it directly
        # We'll use a pattern-based approach with checkpoints instead
        pass



if should_run_rule("synteny_processing"):
    rule synteny_assign:
        input:
            # st_input is now resolution-specific: same window_size feeds into
            # synteny_assign for that resolution
            st_input = "{synteny_results}/{sample}_out/{window_size}/st_input"
        output:
            done_file = "{synteny_results}/synteny_results/{window_size}/{sample}_out/synteny_assign_done"
        params:
            ref_name = config["reference_name"],
            sample_name = lambda wildcards: wildcards.sample,
            window_size = lambda wildcards: int(wildcards.window_size),
            step_size = STEP_SIZE,
            output_dir = lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])),
            log_file = lambda wildcards, output: os.path.abspath(
                os.path.join(os.path.dirname(output[0]), "logs", "synteny_assign.log")
            )

        log:
            "{synteny_results}/synteny_results/{window_size}/{sample}_out/logs/synteny_assign.log"
        shell:
            r"""
            set -euo pipefail

            # STAGE 1: SETUP
            TARGET_DIR="{params.output_dir}"
            LOG_DIR="$(dirname "{params.log_file}")"
            INPUT_FILE="{input.st_input}"

            mkdir -p "$LOG_DIR"
            mkdir -p "$TARGET_DIR"

            # STAGE 2: COPY INPUT
            cp "$INPUT_FILE" "$TARGET_DIR/"

            # STAGE 3: RUN TOOL
            cd "$TARGET_DIR"

            if ! syntenytracker2 {params.ref_name} {params.sample_name} {params.window_size} {params.step_size} st_input > "{params.log_file}" 2>&1; then
                echo "----------------------------------------" >&2
                echo "SYNTENY ASSIGNMENT FAILURE DETAILS" >&2
                echo "Working directory: $(pwd)" >&2
                echo "Command: syntenytracker2 {params.ref_name} {params.sample_name} {params.window_size} {params.step_size} st_input" >&2
                echo "Last 10 lines of log:" >&2
                tail -n 10 "{params.log_file}" >&2
                echo "----------------------------------------" >&2

                # 🔴 SAFE COPY EVEN ON FAILURE
                BLOCK_FILE=$(find "$TARGET_DIR/jobs" -type f -name blocks_info 2>/dev/null | head -n 1 || true)
                if [ -n "$BLOCK_FILE" ]; then
                    cp "$BLOCK_FILE" "$TARGET_DIR/blocks_info"
                    echo "Recovered blocks_info after failure" >&2
                else
                    echo "No blocks_info found after failure" >&2
                fi

                exit 1
            fi

            # STAGE 4: EXTRACT blocks_info (SUCCESS PATH)
            BLOCK_FILE=$(find "$TARGET_DIR/jobs" -type f -name blocks_info 2>/dev/null | head -n 1 || true)

            if [ -n "$BLOCK_FILE" ]; then
                cp "$BLOCK_FILE" "$TARGET_DIR/blocks_info"
                echo "blocks_info copied to $TARGET_DIR" >> "{params.log_file}"
            else
                echo "WARNING: blocks_info not found after successful run" >&2
                find "$TARGET_DIR/jobs" >&2
            fi

            # STAGE 5: CREATE DONE FILE
            touch "{output.done_file}"

            echo -e "\033[32;1m✔\033[0m \033[33;1mCompleted\033[0m synteny tracking of \033[36;1m{wildcards.sample}\033[0m in {params.window_size} resolution" >&2
            """




if should_run_rule("synteny_processing"):
    rule make_eba_format1:
        input:
            # Dynamically create input files for ALL window sizes
            done_files = expand(
                f"{SYNTENY_RESULTS}/synteny_results/{{window_size}}/{{sample}}_out/synteny_assign_done",
                sample=SAMPLES,
                window_size=WINDOW_SIZES  # This should be your list of window sizes
            ) if should_run_rule("synteny_processing") and SAMPLES else [],

            scaffolds = config["eba_format"]["scaffolds_file"] if config["eba_format"]["scaffolds_file"] and config["eba_format"]["scaffolds_file"] != "ALL_CHROMOSOMES" else []

        output:
            eba_dir = directory(os.path.join(SYNTENY_RESULTS, config["out_final"])),
            marker = touch(os.path.join(SYNTENY_RESULTS, config["out_final"], ".done_marker"))

        params:
            synteny_dir = f"{SYNTENY_RESULTS}/synteny_results",
            script_path = os.path.abspath(os.path.join(SCRIPTS, "script.sh")),
            output_dir = os.path.join(SYNTENY_RESULTS, config["out_final"]),
            scaffolds_file = config["eba_format"]["scaffolds_file"]

        log:
            "logs/make_eba_format.log"

        shell:
            """
            # Create parent directory for marker
            mkdir -p {params.output_dir}

            # Verify script exists
            if [ ! -f "{params.script_path}" ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            # Verify synteny_dir exists and has data
            if [ ! -d "{params.synteny_dir}" ] || [ -z "$(ls -A {params.synteny_dir})" ]; then
                echo "Error: Synteny directory is empty or doesn't exist: {params.synteny_dir}" >&2 | tee {log}
                echo "Available windows: $(ls {params.synteny_dir} 2>/dev/null || echo 'none')" >&2 | tee -a {log}
                exit 1
            fi

            # Build the command based on scaffolds_file value
            if [ "{params.scaffolds_file}" = "ALL_CHROMOSOMES" ]; then
                echo "Using Chromosomes for all species (special value: ALL_CHROMOSOMES)" >> {log}
                echo "Processing all window sizes from: {params.synteny_dir}" >> {log}
                echo "Window sizes present: $(ls {params.synteny_dir})" >> {log}

                {params.script_path} \
                    -i "{params.synteny_dir}" \
                    -e "{output.eba_dir}" \
                    -c "Chromosomes" > {log} 2>&1
            else
                # Verify scaffolds file exists if it's a real file path
                if [ ! -f "{input.scaffolds}" ]; then
                    echo "Error: Scaffolds file not found at {input.scaffolds}" >&2 | tee {log}
                    exit 1
                fi
                echo "Using scaffolds file: {input.scaffolds}" >> {log}
                echo "Processing all window sizes from: {params.synteny_dir}" >> {log}

                {params.script_path} \
                    -i "{params.synteny_dir}" \
                    -e "{output.eba_dir}" \
                    -s "{input.scaffolds}" > {log} 2>&1
            fi

            # Verify output
            if [ ! -d "{output.eba_dir}" ] || [ -z "$(ls -A {output.eba_dir})" ]; then
                echo "Error: EBA output directory is empty or doesn't exist: {output.eba_dir}" >&2 | tee {log}
                echo "Last 10 lines of script output:" >&2
                tail -n 10 {log} >&2
                exit 1
            fi

            # Create completion marker
            touch "{output.marker}"

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m creating EBA input format to run EBA" >&2
            echo "Processed all window sizes from: {params.synteny_dir}" >&2
            """


    # Overlap resolution rule - runs AFTER EBA format creation
    rule overlap_resolve1:
        input:
            marker = rules.make_eba_format1.output.marker,  # Wait for make_eba_format to complete
            eba_dir = rules.make_eba_format1.output.eba_dir
        output:
            marker = touch(os.path.join(SYNTENY_RESULTS, config["out_final"], ".overlap_resolve_done"))
        params:
            eba_dir = os.path.join(SYNTENY_RESULTS, config["out_final"]),
            script_path = os.path.abspath(os.path.join(SCRIPTS, "overlap_resolve-v7.py"))
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "overlap_resolve.log")
        shell:
            """
            # Create log directory
            mkdir -p $(dirname {log})

            # Verify script exists
            if [ ! -f "{params.script_path}" ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            # Verify input directory exists
            if [ ! -d "{params.eba_dir}" ]; then
                echo "Error: EBA directory not found at {params.eba_dir}" >&2 | tee {log}
                exit 1
            fi

            # Run the Python script with the directory as argument
            echo "Running overlap resolution on: {params.eba_dir}" >> {log}
            if ! python3 "{params.script_path}" "{params.eba_dir}" 2>&1 | tee -a {log}; then
                echo "Error: overlap_resolve-v7.py failed" >&2
                exit 1
            fi

            echo "Overlap resolution completed successfully" >> {log}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m overlap resolution in EBRs" >&2
            """

        # Rule to concatenate all non-overlap files from numeric subdirectories
    rule concatinate_nonoverlap_files1:
        input:
            marker = rules.overlap_resolve1.output.marker  # Wait for overlap_resolve to complete
        output:
            marker = touch(os.path.join(SYNTENY_RESULTS, ".concatination_done"))
        params:
            eba_dir = os.path.join(SYNTENY_RESULTS, config["out_final"]),
            script_path = os.path.abspath(os.path.join(SCRIPTS, "concatinate_files.py"))
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "concatinate_files.log")
        shell:
            """
            # Create log directory
            mkdir -p $(dirname {log})

            # Verify script exists
            if [ ! -f "{params.script_path}" ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            # Verify input directory exists (output from make_eba_format1)
            if [ ! -d "{params.eba_dir}" ]; then
                echo "Error: EBA directory not found at {params.eba_dir}" >&2 | tee {log}
                exit 1
            fi

            # Verify overlap resolution completed by checking for its marker
            if [ ! -f "{input.marker}" ]; then
                echo "Error: Overlap resolution marker not found: {input.marker}" >&2 | tee {log}
                exit 1
            fi

            # Verify overlap resolution completed (check for numeric subdirectories with non-overlap files)
            echo "Checking for numeric subdirectories with *_nonoverlap.txt files..." >> {log}

            # Count numeric subdirectories
            num_subdirs=$(find "{params.eba_dir}" -maxdepth 1 -type d -name "[0-9]*" | wc -l)
            echo "Found $num_subdirs numeric subdirectories" >> {log}

            if [ "$num_subdirs" -eq 0 ]; then
                echo "Warning: No numeric subdirectories found in {params.eba_dir}" >> {log}
                echo "This might indicate overlap resolution didn't run properly." >> {log}
            fi

            # Run the concatinate_files.py script
            echo "Running concatinate_files.py on: {params.eba_dir}" >> {log}

            if ! python3 "{params.script_path}" "{params.eba_dir}" 2>&1 | tee -a {log}; then
                echo "Error: concatinate_files.py failed" >&2 | tee -a {log}
                exit 1
            fi

            # Simple verification that combined files exist
            echo "Verifying combined files were created..." >> {log}
            combined_files=$(find "{params.eba_dir}" -name "all_nonoverlap_combined_*.txt" | wc -l)
            echo "Found $combined_files combined files" >> {log}

            if [ "$combined_files" -eq 0 ]; then
                echo "Warning: No combined files were created." >> {log}
            fi

            # Create completion marker (in SYNTENY_RESULTS, not inside pre-EBA1)
            touch "{output.marker}"

            echo "File concatenation completed successfully" >> {log}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m concatenation of non-overlap files" >&2
            """




if should_run_rule("synteny_processing"):
    rule reformat_combined_files:
        input:
            # Wait for concatenation to complete
            concatenation_marker = rules.concatinate_nonoverlap_files1.output.marker
        output:
            # Create marker file
            marker = touch(os.path.join(SYNTENY_RESULTS, ".reformatting_done"))
        params:
            eba_dir = os.path.join(SYNTENY_RESULTS, config["out_final"]),
            # Get window sizes from config
            window_sizes = lambda wildcards: WINDOW_SIZES
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "reformat_combined_files.log")
        shell:
            """
            # Create log directory if it doesn't exist
            mkdir -p $(dirname {log})

            echo "Starting reformatting of combined files..." > {log}
            echo "Window sizes: {params.window_sizes}" >> {log}

            # Check if EBA directory exists
            if [ ! -d "{params.eba_dir}" ]; then
                echo "Error: EBA directory not found at {params.eba_dir}" >&2 | tee {log}
                exit 1
            fi

            # Get all numeric subdirectories (window sizes)
            cd "{params.eba_dir}"
            subdirs=$(ls -d [0-9]*/ 2>/dev/null || true)

            if [ -z "$subdirs" ]; then
                echo "Warning: No numeric subdirectories found in {params.eba_dir}" >> {log}
                echo "Cannot reformat files - exiting." >> {log}
                exit 0  # Exit gracefully if no directories
            fi

            echo "Found subdirectories: $subdirs" >> {log}

            # Process each numeric subdirectory
            for subdir in $subdirs; do
                # Remove trailing slash
                res=$(basename "$subdir")

                # Define file paths
                input_file="{params.eba_dir}/$res/all_nonoverlap_combined_${{res}}.txt"
                output_file="{params.eba_dir}/$res/all_nonoverlap_combined_${{res}}_reformatted.txt"

                echo "Processing resolution: $res" >> {log}
                echo "Input file: $input_file" >> {log}
                echo "Output file: $output_file" >> {log}

                # Check if input file exists
                if [ ! -f "$input_file" ]; then
                    echo "Warning: Input file not found: $input_file" >> {log}
                    echo "Skipping reformatting for resolution $res" >> {log}
                    continue
                fi

                # Reformat using awk with tab delimiter
                echo "Reformatting $input_file..." >> {log}
                if ! awk -F"\t" '{{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$1"\t"$9}}' "$input_file" > "$output_file" 2>> {log}; then
                    echo "Error: Failed to reformat file $input_file" >&2 | tee -a {log}
                    exit 1
                fi

                # Verify output file was created and has content
                if [ ! -f "$output_file" ]; then
                    echo "Error: Output file $output_file was not created" >&2 | tee -a {log}
                    exit 1
                fi

                # Count lines
                input_lines=$(wc -l < "$input_file" 2>/dev/null || echo "0")
                output_lines=$(wc -l < "$output_file" 2>/dev/null || echo "0")

                if [ "$input_lines" -eq "$output_lines" ]; then
                    echo "Successfully reformatted $input_lines lines for resolution $res" >> {log}
                else
                    echo "Warning: Line count mismatch for resolution $res" >> {log}
                    echo "Input lines: $input_lines, Output lines: $output_lines" >> {log}
                fi
            done

            # Summary
            echo -e "\n=== Reformatting Summary ===" >> {log}
            reformatted_count=0
            for subdir in $subdirs; do
                res=$(basename "$subdir")
                output_file="{params.eba_dir}/$res/all_nonoverlap_combined_${{res}}_reformatted.txt"
                if [ -f "$output_file" ]; then
                    line_count=$(wc -l < "$output_file" 2>/dev/null || echo "0")
                    echo "✓ $res: $output_file ($line_count lines)" >> {log}
                    reformatted_count=$((reformatted_count + 1))
                else
                    echo "✗ $res: Not reformatted" >> {log}
                fi
            done

            echo "Reformatted $reformatted_count out of $(echo "$subdirs" | wc -w) resolutions" >> {log}

            # Create completion marker
            touch "{output.marker}"

            echo "Reformatting completed successfully" >> {log}
#            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m reformatting of combined files in \\e[36;1m{params.eba_dir}\\e[0m" >&2
            """


def _find_synteny_reference_fasta():
    """
    Locate the reference FASTA in satsuma_align.fasta_dir (inputs/fasta) using
    reference_name from run_sybr_config.yaml. Tries extensions: .fa, .fna, .fasta
    """
    fasta_dir = config.get("satsuma_align", {}).get("fasta_dir", "")
    ref_name  = config.get("reference_name", "")
    for ext in (".fa", ".fna", ".fasta"):
        candidate = os.path.join(fasta_dir, ref_name + ext)
        if os.path.isfile(candidate):
            return candidate
    return ""

# NOTE: this file must live OUTSIDE os.path.join(SYNTENY_RESULTS, config["out_final"]) —
# that path is declared as a directory() output by rule make_eba_format1, and Snakemake
# wipes/recreates directory() outputs before running their owning job, which would delete
# this file (and its log) if it were nested inside.
SYNTENY_CHROM_SIZES = os.path.join(SYNTENY_RESULTS, "chrom_sizes.txt")

# --annot and --sizes are optional for eh_plot-interective8.py, so only wire them in
# when the underlying data actually exists (some input sets — e.g. a custom
# base_input_dir without an enrichment_analysis/ or fasta/<reference_name> — won't have
# them). Otherwise Snakemake would hard-fail the whole rule on a MissingInputException.
_SYNTENY_ANNOT_FILE    = config.get("enrichment", {}).get("annotation_file", "")
_HAS_SYNTENY_ANNOT     = bool(_SYNTENY_ANNOT_FILE) and os.path.isfile(_SYNTENY_ANNOT_FILE)
_HAS_SYNTENY_REF_FASTA = bool(_find_synteny_reference_fasta())

if should_run_rule("synteny_processing") and _HAS_SYNTENY_REF_FASTA:
    rule generate_synteny_chrom_sizes:
        """
        Generate a reference chromosome-sizes file (chr<TAB>size) with seqkit,
        for use as the --sizes argument of eh_plot-interective8.py.
        """
        input:
            ref_fasta = _find_synteny_reference_fasta()
        output:
            sizes = SYNTENY_CHROM_SIZES
        params:
            ref_name = config.get("reference_name", "")
        log:
            os.path.join(SYNTENY_RESULTS, "logs", "generate_synteny_chrom_sizes.log")
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {output.sizes})"
            mkdir -p "$(dirname {log})"

            echo "=== Generating chromosome sizes for synteny plots ===" > {log}
            echo "Reference FASTA : {input.ref_fasta}" >> {log}
            echo "Reference name  : {params.ref_name}"  >> {log}

            if ! command -v seqkit &> /dev/null; then
                echo "Error: seqkit not found. Please install seqkit." >&2 | tee -a {log}
                exit 1
            fi

            seqkit fx2tab --length --name "{input.ref_fasta}" > "{output.sizes}" 2>> {log}

            if [ ! -s "{output.sizes}" ]; then
                echo "Error: chromosome sizes file is empty" >&2 | tee -a {log}
                exit 1
            fi

            echo "Wrote $(wc -l < "{output.sizes}") entries to {output.sizes}" >> {log}
            """

if should_run_rule("synteny_processing"):
    rule generate_synteny_plots:
        input:
            # Wait for reformatting to complete
            reformatting_marker = rules.reformat_combined_files.output.marker,
            annot = [_SYNTENY_ANNOT_FILE] if _HAS_SYNTENY_ANNOT else [],
            sizes = rules.generate_synteny_chrom_sizes.output.sizes if _HAS_SYNTENY_REF_FASTA else [],
            script = os.path.join(workflow.basedir, config["scripts"], "eh_plot-interective8.py")
        output:
            # Create a single marker file for completion
            marker = touch(os.path.join(SYNTENY_RESULTS, ".plots_generated_done"))
        params:
            eba_dir = os.path.join(SYNTENY_RESULTS, config["out_final"])
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "generate_synteny_plots.log")
        shell:
            """
            # Create log directory if it doesn't exist
            mkdir -p $(dirname {log})

            echo "Starting synteny plot generation..." > {log}

            # Check if python3 is available
            if ! command -v python3 &> /dev/null; then
                echo "Error: python3 not found." >&2 | tee -a {log}
                exit 1
            fi

            # Check if the plotting script exists
            if [ ! -f "{input.script}" ]; then
                echo "Error: Plotting script not found at {input.script}" >&2 | tee -a {log}
                exit 1
            fi

            annot_file="{input.annot}"
            sizes_file="{input.sizes}"
            echo "Using annotation file: ${{annot_file:-(none provided)}}" >> {log}
            echo "Using chromosome sizes file: ${{sizes_file:-(none provided)}}" >> {log}

            # --annot / --sizes are optional in eh_plot-interective8.py; only pass
            # them through when the corresponding file is actually available
            plot_extra_args=()
            if [ -n "$annot_file" ]; then
                plot_extra_args+=(--annot "$annot_file")
            fi
            if [ -n "$sizes_file" ]; then
                plot_extra_args+=(--sizes "$sizes_file")
            fi

            # Check if EBA directory exists
            eba_dir="{params.eba_dir}"
            if [ ! -d "$eba_dir" ]; then
                echo "Error: EBA directory not found at $eba_dir" >&2 | tee {log}
                exit 1
            fi

            # Get all numeric subdirectories (window sizes)
            cd "$eba_dir"
            subdirs=$(ls -d [0-9]*/ 2>/dev/null || true)

            if [ -z "$subdirs" ]; then
                echo "Warning: No numeric subdirectories found in $eba_dir" >> {log}
                echo "Cannot generate plots - exiting." >> {log}
                touch "{output.marker}"
                exit 0  # Exit gracefully if no directories
            fi

            echo "Found subdirectories: $subdirs" >> {log}

            # Process each numeric subdirectory (resolution)
            total_plots=0
            for subdir in $subdirs; do
                # Remove trailing slash
                res=$(basename "$subdir")

                res_dir="$eba_dir/$res"
                echo "Processing resolution $res in directory: $res_dir" >> {log}

                if [ ! -d "$res_dir" ]; then
                    echo "Warning: Directory $res_dir not found. Skipping..." >> {log}
                    continue
                fi

                # Check if reformatted file exists
                reformatted_file="$res_dir/all_nonoverlap_combined_${{res}}_reformatted.txt"
                if [ ! -f "$reformatted_file" ]; then
                    echo "Warning: Reformatted file not found: $reformatted_file" >> {log}
                    echo "Skipping resolution $res..." >> {log}
                    continue
                fi

                # Get unique chromosome numbers from first column
                echo "Extracting unique chromosome numbers from $reformatted_file..." >> {log}
                chromosomes=$(awk '{{print $1}}' "$reformatted_file" | sort -n | uniq)

                if [ -z "$chromosomes" ]; then
                    echo "Warning: No chromosomes found in $reformatted_file" >> {log}
                    echo "Skipping resolution $res..." >> {log}
                    continue
                fi

                echo "Found chromosomes for resolution $res: $chromosomes" >> {log}

                # Count total lines in reformatted file
                total_lines=$(wc -l < "$reformatted_file" 2>/dev/null || echo "0")
                echo "Total data points in $reformatted_file: $total_lines" >> {log}

                # Build a comma-separated chromosome list and generate a single
                # multi-tab plot for the whole resolution
                chr_list=$(echo "$chromosomes" | paste -sd, -)
                echo "Generating combined plot for chromosomes $chr_list in resolution $res..." >> {log}

                res_plots=0
                plot_prefix="synteny_plot_${{res}}"

                if ! python3 "{input.script}" "$reformatted_file" \\
                    --chr "$chr_list" \\
                    "${{plot_extra_args[@]}}" \\
                    --output "$res_dir/$plot_prefix" >> {log} 2>&1; then
                    echo "Error: Failed to generate combined plot for resolution $res" >&2 | tee -a {log}
                else
                    # Verify plot was created
                    plot_file="$res_dir/${{plot_prefix}}.html"
                    if [ -f "$plot_file" ]; then
                        echo "Successfully created: $plot_file" >> {log}
                        res_plots=1
                        total_plots=$((total_plots + 1))
                    else
                        echo "Warning: Plot file not created: $plot_file" >> {log}
                    fi
                fi

                echo "Generated $res_plots plot(s) for resolution $res" >> {log}

                # Create completion marker for this resolution
                marker_file="$res_dir/.plots_${{res}}_done"
                touch "$marker_file"
                echo "Created marker file: $marker_file" >> {log}
            done

            # Summary
            echo -e "\n=== Plot Generation Summary ===" >> {log}
            echo "Total plots generated: $total_plots" >> {log}
            echo "Directories processed: $(echo "$subdirs" | wc -w)" >> {log}

            # List all generated plot files
            echo -e "\nGenerated plot files:" >> {log}
            find "$eba_dir" -name "synteny_plot_*.html" -type f | while read plot; do
                echo "  ✓ $plot" >> {log}
            done || echo "  No plot files found" >> {log}

            # Create completion marker
            touch "{output.marker}"

            echo "Plot generation completed successfully" >> {log}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m synteny plot generation" >&2
            """




if should_run_rule("synteny_processing"):
    rule split_by_species:
        input:
            # Wait for reformatting to complete
            reformatting_marker = rules.reformat_combined_files.output.marker
        output:
            # Create a single marker file for completion
            marker = touch(os.path.join(SYNTENY_RESULTS, ".split_by_species_done"))
        params:
            eba_dir = os.path.join(SYNTENY_RESULTS, config["out_final"])
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "split_by_species.log")
        shell:
            """
            # Create log directory if it doesn't exist
            mkdir -p $(dirname {log})

            echo "Starting to split combined files by species..." > {log}

            # Function to split a file by last column (species name)
            split_by_species() {{
                local input_file="$1"
                local output_dir="$2"
                local resolution="$3"

                echo "Processing file: $input_file" >> {log}

                if [ ! -f "$input_file" ]; then
                    echo "Warning: Input file not found: $input_file" >> {log}
                    return 1
                fi

                # Get all unique species names from last column (column 9)
                echo "Extracting unique species names from $input_file..." >> {log}
                species_list=$(awk -F'\t' '{{print $NF}}' "$input_file" | sort -u)

                if [ -z "$species_list" ]; then
                    echo "Warning: No species names found in $input_file" >> {log}
                    return 1
                fi

                echo "Found species: $species_list" >> {log}

                # Count total lines in input file
                total_lines=$(wc -l < "$input_file" 2>/dev/null || echo "0")
                echo "Total lines in input file: $total_lines" >> {log}

                # Create output directory if it doesn't exist
                mkdir -p "$output_dir"

                # Split the file for each species
                species_count=0
                for species in $species_list; do
                    # Clean the species name for filename (replace special characters with underscores)
                    clean_species=$(echo "$species" | sed 's/[[:space:]:\/]/_/g')
                    output_file="$output_dir/${{clean_species}}_${{resolution}}_reformatted_split.txt"

                    echo "Extracting data for species: $species -> $output_file" >> {log}

                    # Extract rows where last column matches the species
                    awk -F'\t' -v species="$species" '$NF == species {{print $0}}' "$input_file" > "$output_file"

                    # Count lines in output file
                    line_count=$(wc -l < "$output_file" 2>/dev/null || echo "0")

                    if [ "$line_count" -gt 0 ]; then
                        echo "  Created $output_file with $line_count lines" >> {log}
                        species_count=$((species_count + 1))
                    else
                        echo "  Warning: $output_file is empty or not created" >> {log}
                        # Remove empty file
                        rm -f "$output_file"
                    fi
                done

                echo "Created $species_count species-specific files for resolution $resolution" >> {log}

                return 0
            }}

            # Check if EBA directory exists
            eba_dir="{params.eba_dir}"
            if [ ! -d "$eba_dir" ]; then
                echo "Error: EBA directory not found at $eba_dir" >&2 | tee {log}
                exit 1
            fi

            # Get all numeric subdirectories (window sizes)
            cd "$eba_dir"
            subdirs=$(ls -d [0-9]*/ 2>/dev/null || true)

            if [ -z "$subdirs" ]; then
                echo "Warning: No numeric subdirectories found in $eba_dir" >> {log}
                echo "Cannot split files - exiting." >> {log}
                touch "{output.marker}"
                exit 0  # Exit gracefully if no directories
            fi

            echo "Found subdirectories: $subdirs" >> {log}

            # Process each numeric subdirectory (resolution)
            total_species_files=0
            for subdir in $subdirs; do
                # Remove trailing slash
                res=$(basename "$subdir")

                res_dir="$eba_dir/$res"
                echo "Processing resolution $res in directory: $res_dir" >> {log}

                if [ ! -d "$res_dir" ]; then
                    echo "Warning: Directory $res_dir not found. Skipping..." >> {log}
                    continue
                fi

                # Check if reformatted file exists
                reformatted_file="$res_dir/all_nonoverlap_combined_${{res}}_reformatted.txt"
                if [ ! -f "$reformatted_file" ]; then
                    echo "Warning: Reformatted file not found: $reformatted_file" >> {log}
                    echo "Skipping resolution $res..." >> {log}
                    continue
                fi

                # Create a subdirectory for split files
                split_dir="$res_dir/split_by_species"
                echo "Creating split directory: $split_dir" >> {log}
                mkdir -p "$split_dir"

                # Split the file
                echo "Splitting $reformatted_file by species..." >> {log}
                if split_by_species "$reformatted_file" "$split_dir" "$res"; then
                    echo "Successfully split file for resolution $res" >> {log}

                    # Count files created for this resolution
                    res_files=$(ls "$split_dir"/*.txt 2>/dev/null | wc -l)
                    total_species_files=$((total_species_files + res_files))
                    echo "Created $res_files species files for resolution $res" >> {log}

                    # Create completion marker for this resolution
                    marker_file="$res_dir/.split_by_species_${{res}}_done"
                    touch "$marker_file"
                    echo "Created marker file: $marker_file" >> {log}

                    # List the created files
                    echo "Files created in $split_dir:" >> {log}
                    ls -la "$split_dir"/*.txt 2>/dev/null | head -20 >> {log} || echo "No files found" >> {log}
                    if [ "$res_files" -gt 20 ]; then
                        echo "... and $((res_files - 20)) more files" >> {log}
                    fi
                else
                    echo "Error: Failed to split file for resolution $res" >&2 | tee -a {log}
                    # Continue with other resolutions instead of exiting
                    echo "Continuing with other resolutions..." >> {log}
                fi
            done

            # Summary
            echo -e "\n=== Species Splitting Summary ===" >> {log}
            echo "Total species-specific files created: $total_species_files" >> {log}
            echo "Resolutions processed: $(echo "$subdirs" | wc -w)" >> {log}

            # List all split directories
            echo -e "\nSplit directories created:" >> {log}
            find "$eba_dir" -name "split_by_species" -type d | while read split_dir; do
                file_count=$(ls "$split_dir"/*.txt 2>/dev/null | wc -l)
                echo "  ✓ $split_dir ($file_count files)" >> {log}
            done || echo "  No split directories found" >> {log}

            # Create completion marker
            touch "{output.marker}"

            echo "Species splitting completed successfully" >> {log}
            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m splitting EBRs by lineages" >&2
            """


if should_run_rule("synteny_processing"):
    rule linear_synteny_plots:
        """
        Generate linear synteny plots using syntenyPlotteR (draw.linear) for every
        resolution directory produced by the split_by_species rule.
        One combined PDF (all query species vs reference) is written per resolution
        into the split_by_species/ subfolder.
        """
        input:
            split_done = rules.split_by_species.output.marker
        output:
            marker = touch(os.path.join(SYNTENY_RESULTS, ".linear_synteny_plots_done"))
        params:
            eba_dir      = os.path.join(SYNTENY_RESULTS, config["out_final"]),
            lengths_file = config["sequence_lengths_file"],
            script_path  = os.path.join(SCRIPTS, "linear_synteny.R"),
        log:
            os.path.join(SYNTENY_RESULTS, config["out_final"], "logs", "linear_synteny_plots.log")
        shell:
            """
            mkdir -p "$(dirname {log})"
            echo "=== Linear synteny plot generation ===" > {log}
            echo "EBA dir : {params.eba_dir}"    >> {log}
            echo "Lengths : {params.lengths_file}" >> {log}
            echo "R script: {params.script_path}"  >> {log}

            if [ ! -f "{params.script_path}" ]; then
                echo "ERROR: R script not found at {params.script_path}" | tee -a {log} >&2
                exit 1
            fi

            if [ ! -f "{params.lengths_file}" ]; then
                echo "ERROR: Lengths file not found at {params.lengths_file}" | tee -a {log} >&2
                exit 1
            fi

            total_plots=0
            error_count=0

            # Iterate over each resolution subdirectory (100k, 300k, 500k …)
            for res_dir in "{params.eba_dir}"/*/; do
                res=$(basename "$res_dir")
                split_dir="${{res_dir%/}}/split_by_species"

                if [ ! -d "$split_dir" ]; then
                    echo "No split_by_species dir in $res_dir — skipping" >> {log}
                    continue
                fi

                txt_files=$(ls "$split_dir"/*.txt 2>/dev/null | wc -l)
                if [ "$txt_files" -eq 0 ]; then
                    echo "No .txt files in $split_dir — skipping" >> {log}
                    continue
                fi

                echo "Processing resolution $res ($txt_files alignment files) …" >> {log}

                # Pass the full lengths file directly — the R script handles col swapping
                # and resolution-suffix stripping internally.
                if Rscript "{params.script_path}" \
                        "{params.lengths_file}" \
                        "$split_dir" \
                        >> {log} 2>&1; then
                    new_pdfs=$(find "$split_dir" -name "*.pdf" | wc -l)
                    if [ "$new_pdfs" -gt 0 ]; then
                        echo "  ✓ $new_pdfs PDF(s) written to $split_dir" >> {log}
                        total_plots=$((total_plots + new_pdfs))
                    else
                        echo "  ✗ Rscript succeeded but no PDFs found in $split_dir" >> {log}
                        error_count=$((error_count + 1))
                    fi
                else
                    echo "  ✗ Rscript failed for resolution $res (see log above)" >> {log}
                    error_count=$((error_count + 1))
                fi
            done

            echo "" >> {log}
            echo "=== Summary ===" >> {log}
            echo "Total PDFs generated : $total_plots" >> {log}
            echo "Resolutions with errors: $error_count" >> {log}

            if [ "$error_count" -gt 0 ]; then
                echo "WARNING: $error_count resolution(s) failed — check log for details" | tee -a {log} >&2
            fi

            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m linear synteny plot generation" >&2
            """
