# Enrichment analysis module

# Enrichment analysis rules (only if enabled)
if should_run_rule("enrichment_analysis"):
    rule msHSBs_background:
        input:
            msHSBs = rules.generate_msHSBs.output.final if should_run_rule("eba_analysis") else config["enrichment"].get("msHSBs_file", ""),
            annotation = config["enrichment"]["annotation_file"],
            kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
            script = os.path.join(SCRIPTS, "coordinate_overlap2.py")  # ADD THIS LINE - missing script input
        output:
            overlap = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "msHSBs_NCBI_genes_overlap.txt"),
            foreground = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "foreground.txt"),
            background = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "background.txt"),
            copied_kegg = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "3kegg_annotationTOgenes.txt") if config["getenrich"]["r"] == "ko" else []
        params:
            outdir = config["enrichment"].get("msHSBs_dir", "msHSBs"),
            resource = config["getenrich"]["r"]  # Get resource type
        log:
            "logs/msHSBs_background.log"
#        message:
#            "Generating foreground and background gene lists for msHSBs enrichment analysis"
        shell:
            """
            mkdir -p {params.outdir}
            
            # Step 1: Run overlap script
            python3 {input.script} {input.msHSBs} {input.annotation} {output.overlap} 0 1 2 3 1 2 > {log} 2>&1

            # Step 2: Remove header from overlap file
            tail -n +2 {output.overlap} > {params.outdir}/temp && mv {params.outdir}/temp {output.overlap}

            # Step 3: Create foreground list
            awk -F '\\t' '{{print $8}}' {output.overlap} > {output.foreground}
            sed -i '1i gene' {output.foreground}

            # Step 4: Only copy KEGG annotation and generate background list if resource is "ko"
            if [ "{params.resource}" = "ko" ]; then
                cp {input.kegg} {output.copied_kegg}
                awk -F '\\t' '{{print $2}}' {output.copied_kegg} | sort | uniq > {output.background}
                sed -i '/gene/d' {output.background}
                sed -i '1i gene' {output.background}
            else
                # For other resources, create background from annotation file
                awk -F '\\t' '{{print $2}}' {input.annotation} | sort | uniq > {output.background}
                sed -i '/gene/d' {output.background}
                sed -i '1i gene' {output.background}
            fi
            """

    if should_run_rule("enrichment_analysis"):
        rule EBRs_background:
            input:
                EBRs = rules.process_breakpoints.output.ebrs_file if should_run_rule("eba_analysis") else [],
                annotation = config["enrichment"]["annotation_file"],
                kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
                script = os.path.join(SCRIPTS, "coordinate_overlap2.py")  # ADD THIS LINE
            output:
                overlap = os.path.join(ebrs_dir, "all_EBRs/EBRs_NCBI_genes_overlap.txt"),
                foreground = os.path.join(ebrs_dir, "all_EBRs/foreground.txt"),
                background = os.path.join(ebrs_dir, "all_EBRs/background.txt"),
                copied_kegg = os.path.join(ebrs_dir, "all_EBRs/3kegg_annotationTOgenes.txt") if config["getenrich"]["r"] == "ko" else []
            params:
                outdir = os.path.join(ebrs_dir, "all_EBRs"),
                resource = config["getenrich"]["r"]  # Get resource type
            log:
                "logs/EBRs_background.log"
#            message:
#                "Generating foreground and background gene lists for EBRs enrichment analysis"
            shell:
                """
                mkdir -p {params.outdir}
                
                # Step 1: Run overlap script
                python3 {input.script} {input.EBRs} {input.annotation} {output.overlap} 0 1 2 3 1 2 > {log} 2>&1

                # Step 2: Remove header from overlap file
                tail -n +2 {output.overlap} > {params.outdir}/temp && mv {params.outdir}/temp {output.overlap}

                # Step 3: Create foreground list (using $9 for EBRs instead of $8 for msHSBs)
                awk -F '\\t' '{{print $9}}' {output.overlap} > {output.foreground}
                sed -i '1i gene' {output.foreground}

                # Step 4: Only copy KEGG annotation and generate background list if resource is "ko"
                if [ "{params.resource}" = "ko" ]; then
                    cp {input.kegg} {output.copied_kegg}
                    awk -F '\\t' '{{print $2}}' {output.copied_kegg} | sort | uniq > {output.background}
                    sed -i '/gene/d' {output.background}
                    sed -i '1i gene' {output.background}
                else
                    # For other resources, create background from annotation file
                    awk -F '\\t' '{{print $2}}' {input.annotation} | sort | uniq > {output.background}
                    sed -i '/gene/d' {output.background}
                    sed -i '1i gene' {output.background}
                fi
                """
    rule split_EBRs:
        input:
            ebrs_file = rules.process_breakpoints.output.ebrs_file
        output:
            done = touch(os.path.join(ebrs_dir, "EBRs_split_done.txt"))
        params:
            outdir = ebrs_dir,
            script_path = os.path.join(SCRIPTS, "EBRs_file_splitter3.py")
        log:
            "logs/split_EBRs.log"
#        message:
#            "Splitting EBRs file into individual chromosome files"
        shell:
            """
            mkdir -p {params.outdir}
            mkdir -p $(dirname {log})
            
            if [ ! -f {params.script_path} ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            python3 {params.script_path} {input.ebrs_file} {params.outdir} > {log} 2>&1
            touch {output.done}
            """

    if should_run_rule("enrichment_analysis"):
        rule EBRs_subdirs_overlap:
            input:
                brkfile3_done = os.path.join(ebrs_dir, "all_EBRs/intermediate_files/brkfile3_done.txt"),
                split_done = os.path.join(ebrs_dir, "EBRs_split_done.txt"),
                annotation = config["enrichment"]["annotation_file"],
                kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
                script = os.path.join(SCRIPTS, "coordinate_overlap3.py")
            output:
                done = touch(os.path.join(ebrs_dir, "EBRs_subdirs_overlap_done.txt"))
            params:
                ebrs_dir = ebrs_dir,
                annotation = config["enrichment"]["annotation_file"],
                resource = config["getenrich"]["r"],  # Get resource type
                columns = "0 1 2 3 1 2",
                log_dir = "logs"
            log:
                "logs/EBRs_subdirs_overlap.log"
#            message:
#                "Running coordinate overlap analysis for EBRs in all lineage subdirectories"
            shell:
                """
                # Create log directory if needed
                mkdir -p {params.log_dir}

                # Verify script exists
                if [ ! -f "{input.script}" ]; then
                    echo "Error: Script not found at {input.script}" >&2 | tee {log}
                    exit 1
                fi

                # First run the original coordinate overlap script
                echo "Running coordinate overlap script on EBRs directory" >> {log}
                if ! python3 {input.script} {params.ebrs_dir} {params.annotation} {params.columns} >> {log} 2>&1; then
                    echo "Error: coordinate_overlap3.py failed" >&2 | tee -a {log}
                    exit 1
                fi

                # Then process each lineage subdirectory
                echo "Processing lineage subdirectories" >> {log}
                for lineage_dir in {params.ebrs_dir}/*_lineage; do
                    if [ -d "$lineage_dir" ]; then
                        echo "Processing $lineage_dir" >> {log}
                        
                        if [ -f "$lineage_dir/EBRs_NCBI_genes_overlap.txt" ]; then
                            # Step 2: Remove header from overlap file
                            if ! tail -n +2 "$lineage_dir/EBRs_NCBI_genes_overlap.txt" > "$lineage_dir/temp" 2>> {log} || \
                            ! mv "$lineage_dir/temp" "$lineage_dir/EBRs_NCBI_genes_overlap.txt" 2>> {log}; then
                                echo "Error processing $lineage_dir/EBRs_NCBI_genes_overlap.txt" >&2 | tee -a {log}
                                exit 1
                            fi
                            
                            # Step 3: Create foreground list
                            if ! awk -F '\\t' '{{print $9}}' "$lineage_dir/EBRs_NCBI_genes_overlap.txt" > "$lineage_dir/foreground.txt" 2>> {log} || \
                            ! sed -i '1i gene' "$lineage_dir/foreground.txt" 2>> {log}; then
                                echo "Error creating foreground for $lineage_dir" >&2 | tee -a {log}
                                exit 1
                            fi
                            
                            # Step 4: Handle background generation based on resource type
                            if [ "{params.resource}" = "ko" ]; then
                                # Copy KEGG annotation and generate background list for KEGG
                                if ! cp {input.kegg} "$lineage_dir/3kegg_annotationTOgenes.txt" 2>> {log} || \
                                ! awk -F '\\t' '{{print $2}}' "$lineage_dir/3kegg_annotationTOgenes.txt" | sort | uniq > "$lineage_dir/background.txt" 2>> {log} || \
                                ! sed -i '/gene/d' "$lineage_dir/background.txt" 2>> {log} || \
                                ! sed -i '1i gene' "$lineage_dir/background.txt" 2>> {log}; then
                                    echo "Error creating background for $lineage_dir" >&2 | tee -a {log}
                                    exit 1
                                fi
                            else
                                # For other resources, create background from annotation file
                                if ! awk -F '\\t' '{{print $2}}' {params.annotation} | sort | uniq > "$lineage_dir/background.txt" 2>> {log} || \
                                ! sed -i '/gene/d' "$lineage_dir/background.txt" 2>> {log} || \
                                ! sed -i '1i gene' "$lineage_dir/background.txt" 2>> {log}; then
                                    echo "Error creating background for $lineage_dir" >&2 | tee -a {log}
                                    exit 1
                                fi
                            fi
                        else
                            echo "Warning: No overlap file found in $lineage_dir" >> {log}
                        fi
                    fi
                done

                touch {output.done}
                """
        rule clone_EBR_configs:
            input:
                overlap_done = os.path.join(ebrs_dir, "EBRs_subdirs_overlap_done.txt"),
                config_template = "tools/getENRICH/config2.json",
                all_ebrs_file = os.path.join(ebrs_dir, "all_EBRs/5KBextenction_brk.txt"),
                mshsbs_file = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "msHSBs_done.txt"),
                script = os.path.join(SCRIPTS, "clone_config.py")
            output:
                all_ebrs_config = os.path.join(ebrs_dir, "all_EBRs", "config2.json"),
                ms_config = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "config2.json"),
                lineage_configs = expand(
                    os.path.join(ebrs_dir, "{lineage}_lineage", "config2.json"),
                lineage=[os.path.basename(d).replace("_lineage", "") 
                         for d in glob.glob(f"{ebrs_dir}/*_lineage")]
                ),
                done = touch(os.path.join(ebrs_dir, "EBR_configs_cloned.txt"))
            params:
                ebrs_dir = ebrs_dir,
                mshsbs_dir = config["enrichment"].get("msHSBs_dir", "msHSBs"),
                log_dir = "logs"
            log:
                "logs/clone_EBR_configs.log"
#            message:
#                "Cloning config files for EBRs and msHSBs enrichment analysis"
            shell:
               """
        	mkdir -p {params.log_dir}

           	if [ ! -f "{input.script}" ]; then
        	    echo "Error: Script not found at {input.script}" >&2 | tee {log}
		    exit 1
           	fi

           	if [ ! -d "{params.ebrs_dir}" ]; then
        	    echo "Error: EBRs directory not found at {params.ebrs_dir}" >&2 | tee {log}
                    exit 1
            	fi

           	if [ ! -d "{params.mshsbs_dir}" ]; then
  		    echo "Error: msHSBs directory not found at {params.mshsbs_dir}" >&2 | tee {log}
                    exit 1
            	fi

            	if [ ! -f "{input.config_template}" ]; then
                    echo "Error: Config template not found at {input.config_template}" >&2 | tee {log}
                    exit 1
            	fi

            	if [ ! -f "{input.all_ebrs_file}" ]; then
                    echo "Error: EBRs file not found at {input.all_ebrs_file}" >&2 | tee {log}
                    exit 1
            	fi

            	if [ ! -f "{input.mshsbs_file}" ]; then
                    echo "Error: msHSBs file not found at {input.mshsbs_file}" >&2 | tee {log}
                    exit 1
           	fi

            	echo "Starting config cloning process" >> {log}

            	# Project root = Snakefile directory
            	PROJECT_ROOT=$(pwd)

            	if ! python3 {input.script} {input.config_template} \
            	    {params.ebrs_dir}/all_EBRs \
                    {params.ebrs_dir}/*_lineage \
                    {params.mshsbs_dir} \
                    --project-root $PROJECT_ROOT >> {log} 2>&1; then
                    echo "Error: clone_config.py failed" >&2 | tee -a {log}
                    exit 1
            	fi

            	touch {output.done}
            	echo "Config cloning completed successfully" >> {log}
            	"""



    if should_run_rule("enrichment_analysis"):
        rule run_getENRICH:
            input:
                configs_ready = os.path.join(ebrs_dir, "EBR_configs_cloned.txt"),
                # Explicitly list all config files to ensure they exist
                all_ebrs_config = os.path.join(ebrs_dir, "all_EBRs", "config2.json"),
                ms_config = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "config2.json"),
                lineage_configs = expand(os.path.join(ebrs_dir, "{lineage}_lineage", "config2.json"),
                                        lineage=[os.path.basename(d).replace("_lineage", "") 
                                                for d in glob.glob(f"{ebrs_dir}/*_lineage")]),
                getENRICH = "tools/getENRICH/getENRICH/getENRICH",
                # Conditionally require KEGG file only if resource is "ko"
                kegg_file = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else []
            output:
                ebrs_all_done = os.path.join(ebrs_dir, "all_EBRs", "getENRICH_done.txt"),
                ebrs_lineage_done = expand(os.path.join(ebrs_dir, "{lineage}_lineage", "getENRICH_done.txt"), 
                                        lineage=[os.path.basename(d).replace("_lineage", "") 
                                                for d in glob.glob(f"{ebrs_dir}/*_lineage")]),
                mshsbs_done = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "getENRICH_done.txt")
            params:
                resource = config["getenrich"]["r"],  # Get resource from config
                getENRICH_dir = os.path.abspath("tools/getENRICH/getENRICH"),
                # Store absolute paths for configs
                all_ebrs_config_abs = lambda wc: os.path.abspath(os.path.join(ebrs_dir, "all_EBRs", "config2.json")),
                ms_config_abs = lambda wc: os.path.abspath(os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "config2.json")),
                # Fix lineage paths to include _lineage suffix
                lineage_configs_abs = lambda wc: [os.path.abspath(os.path.join(ebrs_dir, f"{lineage}_lineage", "config2.json")) 
                                                for lineage in [os.path.basename(d).replace("_lineage", "") 
                                                                for d in glob.glob(f"{ebrs_dir}/*_lineage")]],
                log_dir = "logs"
            log:
                all_log = "logs/getENRICH_all.log",
                lineage_log = "logs/getENRICH_lineages.log",
                ms_log = "logs/getENRICH_msHSBs.log"
            threads: 40  # Allow parallel execution
    #        message:
    #            "Running getENRICH analysis for all EBRs, lineages and msHSBs in parallel"
            shell:
                """
                # Create log directory if needed
                mkdir -p {params.log_dir}
                
                # Create log files explicitly to avoid "No such file" errors
                touch {log.all_log}
                touch {log.lineage_log}
                touch {log.ms_log}

                # Verify getENRICH executable exists
                if [ ! -f "{input.getENRICH}" ]; then
                    echo "Error: getENRICH executable not found at {input.getENRICH}" >&2 | tee -a {log.all_log}
                    exit 1
                fi

                # Verify all config files exist
                if [ ! -f "{input.all_ebrs_config}" ]; then
                    echo "Error: Config file not found at {input.all_ebrs_config}" >&2 | tee -a {log.all_log}
                    exit 1
                fi
                
                if [ ! -f "{input.ms_config}" ]; then
                    echo "Error: Config file not found at {input.ms_config}" >&2 | tee -a {log.all_log}
                    exit 1
                fi

                # Function to run getENRICH with error handling
                run_analysis() {{
                    config=$1
                    outfile=$2
                    logfile=$3
                    analysis_type=$4
                    
                    echo "Starting $analysis_type analysis" >> "$logfile"
                    echo "Working directory: $(pwd)" >> "$logfile"
                    echo "Config path: $config" >> "$logfile"
                    echo "Resource: {params.resource}" >> "$logfile"
                    
                    # Use absolute path for config and run from the correct directory
                    # Continue even if getENRICH returns non-zero exit code (some tools do this)
                    (cd {params.getENRICH_dir} && \
                        ./getENRICH -c "$config" -r {params.resource} -k -l -j -i -v -a 2>&1 | tee -a "$logfile") || true
                    
                    # Check if enrichment results were actually generated
                    result_dir=$(dirname "$config")/enrichment_result
                    if [ -d "$result_dir" ] && [ -n "$(ls -A "$result_dir" 2>/dev/null)" ]; then
                        echo "$analysis_type enrichment results successfully generated in $result_dir" >> "$logfile"
                        echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m pathway enrichment analysis of \e[36m$analysis_type\e[0m" >&2
                        touch "$outfile"
                        return 0
                    else
                        echo "Warning: No enrichment results found for $analysis_type in $result_dir" >> "$logfile"
                        echo -e "\e[31;1m✗\e[0m \e[33;1mFailed\e[0m pathway enrichment analysis of $analysis_type" >&2
                        return 1
                    fi
                }}

                # Export function so it can be used in subshells
                export -f run_analysis

                # Run analyses in parallel using background processes
                echo "Running all analyses in parallel..." | tee -a {log.all_log} {log.lineage_log} {log.ms_log}
                
                # Track all background process PIDs
                declare -A pids
                declare -A analysis_types
                
                # Run all_EBRs analysis in background
                echo -e "\e[33;1mStarting\e[0m pathway enrichment analysis of all EBRs" >&2
                run_analysis "{params.all_ebrs_config_abs}" "{output.ebrs_all_done}" "{log.all_log}" "all EBRs" &
                pids["all_ebrs"]=$!
                
                # Run lineage analyses in background
                echo -e "\e[33;1mStarting\e[0m pathway enrichment analysis of all lineage specific EBRs" >&2
                
                # Count total lineage analyses
                lineage_count=0
                for config_path in {params.lineage_configs_abs}; do
                    if [ -f "$config_path" ]; then
                        lineage_count=$((lineage_count + 1))
                    fi
                done
                
                echo "Found $lineage_count lineage analyses to run" >> {log.lineage_log}
                
                # Start each lineage analysis
                for config_path in {params.lineage_configs_abs}; do
                    if [ ! -f "$config_path" ]; then
                        echo "Warning: Config file $config_path not found, skipping" | tee -a {log.lineage_log}
                        continue
                    fi
                    lineage_dir=$(dirname "$config_path")
                    lineage_name=$(basename "$lineage_dir" | sed 's/_lineage$//')
                    outfile="$lineage_dir/getENRICH_done.txt"
                    echo "Starting analysis for lineage: $lineage_name" | tee -a {log.lineage_log}
                    run_analysis "$config_path" "$outfile" "{log.lineage_log}" "$lineage_name lineage EBRs" &
                    pids["lineage_$lineage_name"]=$!
                done
                
                # Run msHSBs analysis in background
                echo -e "\e[33;1mStarting\e[0m pathway enrichment analysis of all msHSBs" >&2
                run_analysis "{params.ms_config_abs}" "{output.mshsbs_done}" "{log.ms_log}" "msHSBs" &
                pids["mshsbs"]=$!
                
                # Function to check if a process is still running
                is_running() {{
                    local pid=$1
                    kill -0 "$pid" 2>/dev/null
                }}
                
                # Wait for all background processes to complete and track results
                echo "Waiting for all analyses to complete..." | tee -a {log.all_log} {log.lineage_log} {log.ms_log}
                
                declare -A results
                for key in "${{!pids[@]}}"; do
                    pid=${{pids[$key]}}
                    wait "$pid"
                    results[$key]=$?
                    if [ ${{results[$key]}} -eq 0 ]; then
                        echo "Analysis $key completed successfully" | tee -a {log.all_log}
                    else
                        echo "Analysis $key completed with warnings or errors" | tee -a {log.all_log}
                    fi
                done
                
                # Final verification that all output files were created
                failed_analyses=0
                
                if [ ! -f "{output.ebrs_all_done}" ]; then
                    echo "Error: all_EBRs completion file not found" >&2 | tee -a {log.all_log}
                    failed_analyses=$((failed_analyses + 1))
                fi
                
                for lineage in {params.lineage_configs_abs}; do
                    if [ -f "$lineage" ]; then
                        lineage_dir=$(dirname "$lineage")
                        outfile="$lineage_dir/getENRICH_done.txt"
                        if [ ! -f "$outfile" ]; then
                            lineage_name=$(basename "$lineage_dir" | sed 's/_lineage$//')
                            echo "Warning: Lineage $lineage_name completion file not found" | tee -a {log.lineage_log}
                        fi
                    fi
                done
                
                if [ ! -f "{output.mshsbs_done}" ]; then
                    echo "Error: msHSBs completion file not found" >&2 | tee -a {log.ms_log}
                    failed_analyses=$((failed_analyses + 1))
                fi

                if [ $failed_analyses -eq 0 ]; then
                    echo "All getENRICH analyses completed successfully" | tee -a {log.all_log}
                    echo -e "\e[32;1m✔\e[0m All pathway enrichment analyses completed successfully" >&2
                else
                    echo "Some analyses completed with issues. Check logs for details." | tee -a {log.all_log}
                    echo -e "\e[33;1m⚠\e[0m Some pathway enrichment analyses completed with issues" >&2
                fi
                """
