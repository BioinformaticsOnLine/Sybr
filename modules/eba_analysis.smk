if should_run_rule("eba_analysis"):
    rule copy_to_pre_EBA:
        input:
            # Lambda functions in input section MUST accept wildcards parameter
            marker = lambda wildcards: rules.overlap_resolve1.output.marker if should_run_rule("synteny_processing") else [],
            eba_dir = lambda wildcards: rules.make_eba_format1.output.eba_dir if should_run_rule("synteny_processing") else []
        output:
            marker = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input", ".copy_complete")
        params:
            source_dir = lambda wildcards: os.path.join(SYNTENY_RESULTS, config["out_final"]),
            dest_dir = lambda wildcards: os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input")
        log:
            os.path.join(config["eba_format"]["pre_EBA_dir"], "logs", "copy_to_pre_EBA.log")
        shell:
            """
            # Create directories
            mkdir -p {params.dest_dir}
            mkdir -p $(dirname {log})
            
            echo "Source directory: {params.source_dir}" > {log}
            echo "Destination directory: {params.dest_dir}" >> {log}
            
            # Check if source exists
            if [ ! -d "{params.source_dir}" ]; then
                echo "Error: Source directory not found: {params.source_dir}" >&2 | tee -a {log}
                exit 1
            fi
            
            # Get numeric subdirectories (those starting with digits)
            cd "{params.source_dir}"
            num_dirs=$(ls -d [0-9]*k 2>/dev/null || true)
            
            if [ -z "$num_dirs" ]; then
                echo "No numeric directories found in {params.source_dir}" >> {log}
                echo "Available directories:" >> {log}
                ls -la >> {log}
            fi
            
            # Copy numeric subdirectories and rename them (remove 'k' suffix)
            for dir in $num_dirs; do
                if [ -n "$dir" ] && [ -d "$dir" ]; then
                    # Remove 'k' from directory name (e.g., 100k -> 100)
                    new_dir=$(echo "$dir" | sed 's/k$//')
                    
                    echo "Processing directory: $dir -> $new_dir" >> {log}
                    
                    # Create destination directory
                    mkdir -p "{params.dest_dir}/$new_dir"
                    
                    # Copy only *_nonoverlap.txt files
                    find "$dir" -name "*_nonoverlap.txt" -type f | while read file; do
                        filename=$(basename "$file")
                        echo "Copying: $filename to {params.dest_dir}/$new_dir/" >> {log}
                        cp "$file" "{params.dest_dir}/$new_dir/"
                    done
                    
                    # Verify files were copied
                    file_count=$(find "{params.dest_dir}/$new_dir" -name "*_nonoverlap.txt" -type f | wc -l)
                    echo "Copied $file_count files to {params.dest_dir}/$new_dir/" >> {log}
                fi
            done
            
            # Verify at least some files were copied
            total_files=$(find "{params.dest_dir}" -name "*_nonoverlap.txt" -type f | wc -l)
            if [ "$total_files" -eq 0 ]; then
                echo "Error: No files were copied to destination!" >&2 | tee -a {log}
                echo "Available source directories:" >&2
                ls -la "{params.source_dir}/" >&2
                exit 1
            else
                echo "Successfully copied $total_files files to {params.dest_dir}" >> {log}
                # Create marker file
                touch "{output}"
            fi
            
#            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m copying to pre-EBA directory: \e[36;1m{params.dest_dir}\e[0m" >&2
            """

if should_run_rule("eba_analysis"):
    rule run_eba:
        input:
            # Wait for copy_to_pre_EBA to complete
            eba_input_dir = lambda wildcards: os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input", ".copy_complete"),
            classification_file = lambda wildcards: config["eba"]["c"],
            eba_script = lambda wildcards: os.path.join(workflow.basedir, config["eba_tools"]["eba_script_path"]),
            # Use chr_size_file from config
            chr_size_file = lambda wildcards: os.path.join(workflow.basedir, config["eba_format"]["chr_size_file"])
        output:
            eba_working_dir = directory(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING")),
            eba_out_dir = directory(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT")),
            # More flexible output location - main merged results
            final_classify = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT/Merge/Result_Merge.final"),
            marker = touch(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/eba_complete.marker"))
        params:
            n_threads = lambda wildcards: config["eba"]["n"],
            ref_name = lambda wildcards: config["eba"]["r"],
            p_param = lambda wildcards: config["eba"].get("p", 300),  # Added from eba config
            t_param = lambda wildcards: config["eba"].get("t", 20),   # Added from eba config
            k_param = lambda wildcards: config["eba"].get("k", True), # Added from eba config
            classification_abs = lambda wildcards: os.path.abspath(os.path.join(workflow.basedir, config["eba"]["c"])),
            log_file_abs = lambda wildcards: os.path.abspath("logs/run_eba.log"),
            eba_input_abs = lambda wildcards: os.path.abspath(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input")),
            working_dir = lambda wildcards: os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING"),
            copy_dest_dir = lambda wildcards: os.path.join(workflow.basedir, config["eba_format"]["copy_destination_dir"]),
            chr_size_abs = lambda wildcards: os.path.abspath(os.path.join(workflow.basedir, config["eba_format"]["chr_size_file"]))
        log:
            "logs/run_eba.log"
        shell:
            """
            # Create log directory
            mkdir -p "$(dirname {params.log_file_abs})"
            
            # Create and set up the separate working directory
            mkdir -p {params.working_dir}
            cd {params.working_dir} || exit 1
            
            echo "=== EBA Analysis Started ===" | tee -a {params.log_file_abs}
            echo "Setting up EBA working directory in: $(pwd)" | tee -a {params.log_file_abs}
            echo "EBA input directory: {params.eba_input_abs}" | tee -a {params.log_file_abs}
            
            # Copy required files to working directory
            # 1. Copy the entire EBA-input structure from pre-EBA location
            echo "Copying EBA input files from pre-EBA directory..." | tee -a {params.log_file_abs}
            # Copy the actual directory (not the marker file)
            cp -r {params.eba_input_abs} . || {{
                echo "Failed to copy EBA-input directory" | tee -a {params.log_file_abs}
                exit 1
            }}
            
            # 2. Copy classification file
            echo "Copying classification file from: {params.classification_abs}" | tee -a {params.log_file_abs}
            cp {params.classification_abs} . || {{
                echo "Failed to copy classification file" | tee -a {params.log_file_abs}
                exit 1
            }}
            
            # 3. Copy chr_size.txt from config location
            echo "Copying chr_size.txt from: {params.chr_size_abs}" | tee -a {params.log_file_abs}
            if [ -f "{params.chr_size_abs}" ]; then
                cp "{params.chr_size_abs}" . || {{
                    echo "Failed to copy chr_size.txt" | tee -a {params.log_file_abs}
                    exit 1
                }}
                echo "Successfully copied chr_size.txt" | tee -a {params.log_file_abs}
            else
                echo "ERROR: chr_size.txt not found at {params.chr_size_abs}" >&2 | tee -a {params.log_file_abs}
                echo "Please check the path in config file: eba_format.chr_size_file" >&2 | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify the EBA script exists
            if [ ! -f "{input.eba_script}" ]; then
                echo "ERROR: EBA script not found at {input.eba_script}" >&2 | tee -a {params.log_file_abs}
                echo "Please check the path in config file: eba_tools.eba_script_path" >&2 | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify input directory was copied correctly
            local_input_dir="EBA-input"
            if [ ! -d "$local_input_dir" ]; then
                echo "ERROR: EBA input directory not found at $local_input_dir" >&2 | tee -a {params.log_file_abs}
                echo "Available files:" >&2 | tee -a {params.log_file_abs}
                ls -la | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Check contents of EBA-input directory
            echo "=== Contents of EBA-input directory ===" | tee -a {params.log_file_abs}
            ls -la "$local_input_dir/" | tee -a {params.log_file_abs}
            
            # Check if there are resolution directories (100, 300, 500, etc.)
            echo "=== Resolution directories found ===" | tee -a {params.log_file_abs}
            find "$local_input_dir" -maxdepth 1 -type d -name "[0-9]*" | tee -a {params.log_file_abs}
            
            # Check if nonoverlap files exist in resolution directories
            echo "=== Checking for nonoverlap files ===" | tee -a {params.log_file_abs}
            for res_dir in "$local_input_dir"/*/; do
                if [ -d "$res_dir" ]; then
                    echo "Directory: $res_dir" | tee -a {params.log_file_abs}
                    nonoverlap_count=$(find "$res_dir" -name "*_nonoverlap.txt" -type f | wc -l)
                    echo "  Found $nonoverlap_count *_nonoverlap.txt files" | tee -a {params.log_file_abs}
                fi
            done
            
            # Run EBA in the working directory
            echo "=== Starting EBA analysis ===" | tee -a {params.log_file_abs}
            echo "Working directory: $(pwd)" | tee -a {params.log_file_abs}
            echo "Local input directory: $local_input_dir" | tee -a {params.log_file_abs}
            echo "EBA script: {input.eba_script}" | tee -a {params.log_file_abs}
            echo "Reference species: {params.ref_name}" | tee -a {params.log_file_abs}
            echo "Threads: {params.n_threads}" | tee -a {params.log_file_abs}
            
            # Build the EBA command based on config parameters
            EBA_CMD="perl {input.eba_script} \\
                -n {params.n_threads} \\
                -d ./$local_input_dir/ \\
                -r {params.ref_name} \\
                -p {params.p_param} \\
                -t {params.t_param} \\
                -c ./$(basename {params.classification_abs})"
            
            # Add -k flag if k_param is True
            if [ "{params.k_param}" = "True" ] || [ "{params.k_param}" = "true" ] || [ "{params.k_param}" = "TRUE" ]; then
                EBA_CMD="$EBA_CMD -k"
                echo "Using -k flag (keep temp files)" | tee -a {params.log_file_abs}
            fi
            
            echo "Running command: $EBA_CMD" | tee -a {params.log_file_abs}
            echo "Full command details:" | tee -a {params.log_file_abs}
            echo "  -n {params.n_threads} (threads)" | tee -a {params.log_file_abs}
            echo "  -d ./$local_input_dir/ (input directory)" | tee -a {params.log_file_abs}
            echo "  -r {params.ref_name} (reference)" | tee -a {params.log_file_abs}
            echo "  -p {params.p_param} (p parameter)" | tee -a {params.log_file_abs}
            echo "  -t {params.t_param} (t parameter)" | tee -a {params.log_file_abs}
            echo "  -c ./$(basename {params.classification_abs}) (classification)" | tee -a {params.log_file_abs}
            
            # Execute EBA
            echo "Executing EBA..." | tee -a {params.log_file_abs}
            if ! $EBA_CMD >> {params.log_file_abs} 2>&1; then
                echo "ERROR: EBA script failed with exit code $?" | tee -a {params.log_file_abs}
                echo "=== Last 30 lines of EBA output ===" | tee -a {params.log_file_abs}
                tail -n 30 {params.log_file_abs} | tee -a {params.log_file_abs}
                echo "=== Checking intermediate files ===" | tee -a {params.log_file_abs}
                ls -la | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify the expected output file was created
            echo "=== Verifying EBA output ===" | tee -a {params.log_file_abs}
            
            # Check for result files in multiple possible locations
            # 1. First check the Merge directory (main results)
            if [ -f "EBA_OUT/Merge/Result_Merge.final" ]; then
                echo "SUCCESS: Main result file found: EBA_OUT/Merge/Result_Merge.final" | tee -a {params.log_file_abs}
                FINAL_OUTPUT="EBA_OUT/Merge/Result_Merge.final"
            elif [ -f "EBA_OUT/Merge/final_classify.final" ]; then
                echo "SUCCESS: Alternative result file found: EBA_OUT/Merge/final_classify.final" | tee -a {params.log_file_abs}
                FINAL_OUTPUT="EBA_OUT/Merge/final_classify.final"
            # 2. Check resolution-specific directories using the p_param value
            elif [ -d "EBA_OUT/{params.p_param}" ] && [ -f "EBA_OUT/{params.p_param}/EBA_OutFiles/final_classify.eba7" ]; then
                echo "SUCCESS: Resolution-specific result found: EBA_OUT/{params.p_param}/EBA_OutFiles/final_classify.eba7" | tee -a {params.log_file_abs}
                FINAL_OUTPUT="EBA_OUT/{params.p_param}/EBA_OutFiles/final_classify.eba7"
            else
                echo "ERROR: No expected output files were created" | tee -a {params.log_file_abs}
                echo "=== Directory structure ===" | tee -a {params.log_file_abs}
                
                # Check what was actually created
                echo "Working directory contents:" | tee -a {params.log_file_abs}
                ls -la | tee -a {params.log_file_abs}
                
                if [ -d "EBA_OUT" ]; then
                    echo "EBA_OUT directory contents:" | tee -a {params.log_file_abs}
                    ls -la "EBA_OUT/" | tee -a {params.log_file_abs}
                    
                    # Check Merge directory
                    if [ -d "EBA_OUT/Merge" ]; then
                        echo "EBA_OUT/Merge directory contents:" | tee -a {params.log_file_abs}
                        ls -la "EBA_OUT/Merge/" | tee -a {params.log_file_abs}
                    fi
                    
                    # Check all resolution directories
                    echo "Checking all resolution directories:" | tee -a {params.log_file_abs}
                    for res_dir in "EBA_OUT"/*/; do
                        if [ -d "$res_dir" ]; then
                            echo "Directory: $res_dir" | tee -a {params.log_file_abs}
                            ls -la "$res_dir/" | tee -a {params.log_file_abs}
                            
                            # Check for EBA_OutFiles subdirectory
                            if [ -d "$res_dir/EBA_OutFiles" ]; then
                                echo "  EBA_OutFiles subdirectory:" | tee -a {params.log_file_abs}
                                ls -la "$res_dir/EBA_OutFiles/" | tee -a {params.log_file_abs}
                            fi
                        fi
                    done
                fi
                
                # Also check for any .final files anywhere
                echo "=== Searching for any .final files ===" | tee -a {params.log_file_abs}
                find . -name "*.final" -type f | tee -a {params.log_file_abs}
                
                echo "=== Searching for any .eba7 files ===" | tee -a {params.log_file_abs}
                find . -name "*.eba7" -type f | tee -a {params.log_file_abs}
                
                exit 1
            fi
            
            # Get the actual output file
            ACTUAL_OUTPUT="$FINAL_OUTPUT"
            echo "Actual output file: $ACTUAL_OUTPUT" | tee -a {params.log_file_abs}
            echo "File size: $(du -h "$ACTUAL_OUTPUT" | cut -f1)" | tee -a {params.log_file_abs}
            echo "Line count: $(wc -l < "$ACTUAL_OUTPUT")" | tee -a {params.log_file_abs}
            
            # Verify the file exists at the expected location
            # If it's not at the exact expected path, copy it there
            if [ ! -f "{output.final_classify}" ]; then
                echo "Copying output file to expected location" | tee -a {params.log_file_abs}
                mkdir -p "$(dirname {output.final_classify})" || true
                if [ -f "$ACTUAL_OUTPUT" ]; then
                    cp "$ACTUAL_OUTPUT" "{output.final_classify}"
                    echo "Copied: $ACTUAL_OUTPUT -> {output.final_classify}" | tee -a {params.log_file_abs}
                else
                    echo "ERROR: Actual output file not found: $ACTUAL_OUTPUT" | tee -a {params.log_file_abs}
                    exit 1
                fi
            fi
            
            echo "SUCCESS: Output file verified: {output.final_classify}" | tee -a {params.log_file_abs}
            echo "File size: $(du -h "{output.final_classify}" | cut -f1)" | tee -a {params.log_file_abs}
            echo "Line count: $(wc -l < "{output.final_classify}")" | tee -a {params.log_file_abs}
            
            # Create completion marker
            touch {output.marker}
            echo "=== EBA analysis completed successfully ===" | tee -a {params.log_file_abs}
            
            # Copy results to destination directory
            echo "=== Copying results to destination directory ===" | tee -a {params.log_file_abs}
            echo "Destination: {params.copy_dest_dir}" | tee -a {params.log_file_abs}
            
            mkdir -p {params.copy_dest_dir} || {{
                echo "Failed to create destination directory" | tee -a {params.log_file_abs}
                exit 1
            }}
            
            # Copy the entire EBA_OUT directory
            if [ -d "EBA_OUT" ]; then
                echo "Copying EBA_OUT directory..." | tee -a {params.log_file_abs}
                cp -r "EBA_OUT" "{params.copy_dest_dir}/" || {{
                    echo "Failed to copy EBA_OUT directory" | tee -a {params.log_file_abs}
                    exit 1
                }}
                echo "Copied EBA_OUT directory successfully" | tee -a {params.log_file_abs}
            else
                echo "WARNING: EBA_OUT directory not found" | tee -a {params.log_file_abs}
            fi
            
            # Copy any other important files that were generated
            echo "Copying additional output files..." | tee -a {params.log_file_abs}
            important_files=("betaScore" "gaps_brks.stats" "species.sps" "sps.txt")
            for file in "${{important_files[@]}}"; do
                if [ -f "$file" ]; then
                    cp -v "$file" "{params.copy_dest_dir}/" | tee -a {params.log_file_abs}
                else
                    echo "NOTE: $file not found in working directory" | tee -a {params.log_file_abs}
                fi
            done
            
            # List what was copied
            echo "=== Files copied to {params.copy_dest_dir} ===" | tee -a {params.log_file_abs}
            ls -la "{params.copy_dest_dir}/" | tee -a {params.log_file_abs}
            
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m EBRs finding in dataset" >&2
            echo -e "\e[34;1mResults available at: {params.copy_dest_dir}\e[0m" >&2
            """

if should_run_rule("eba_analysis"):
    rule process_ebrs:
        input:
            eba_result = rules.run_eba.output.final_classify,
            eba_complete = rules.run_eba.output.marker,
            script5 = os.path.join(workflow.basedir, config["scripts"], "script5.py"),
            anti_script = os.path.join(workflow.basedir, config["scripts"], "antioverlaping2.py")
        output:
            ebrs_file = os.path.join(config["eba_format"]["ebrs_dir"], "EBRs.txt"),
            marker = touch(os.path.join(config["eba_format"]["ebrs_dir"], "ebrs_processed.marker"))
        params:
            ebrs_dir = lambda wildcards: config["eba_format"]["ebrs_dir"],
            intermediate_dir = lambda wildcards: config["eba_format"]["ebrs_dir"],
            log_dir = "logs"
        log:
            "logs/process_ebrs.log"
        shell:
            """
            mkdir -p {params.log_dir}
            mkdir -p {params.ebrs_dir}

            echo "=== Processing EBA results to extract EBRs ===" | tee -a {log}

            if [ ! -s "{input.eba_result}" ]; then
                echo "ERROR: Input file {input.eba_result} is empty or missing" | tee -a {log}
                exit 1
            fi

            echo "Running script5.py to extract raw EBRs..." | tee -a {log}
            python3 {input.script5} {input.eba_result} {params.intermediate_dir}/brkfile2 >> {log} 2>&1

            # ------------------------------------------------------------------
            # STEP 1 — Extend EBRs by ±5000 bp (your reference logic)
            # ------------------------------------------------------------------
            echo "Extending EBR coordinates by ±5000 bp..." | tee -a {log}
            awk -F '\\t' '{{print $1"\\t"$2-5000"\\t"$3+5000"\\t"$4}}' \
                {params.intermediate_dir}/brkfile2 \
                > {params.intermediate_dir}/brkfile3

            echo "5kb extension complete → brkfile3" | tee -a {log}

            # ------------------------------------------------------------------
            # STEP 2 — Anti-overlap cleanup (your reference logic)
            # ------------------------------------------------------------------
            echo "Running anti-overlap script..." | tee -a {log}
            python3 {input.anti_script} \
                {params.intermediate_dir}/brkfile3 \
                {output.ebrs_file} >> {log} 2>&1

            echo "Anti-overlap complete → final EBRs.txt" | tee -a {log}

            # Create completion marker
            touch {output.marker}

            echo -e "\\e[32;1m✔\\e[0m \\e[33;1mCompleted\\e[0m EBRs extraction, 5kb extension, and cleanup" >&2
            echo -e "\\e[34;1mFinal file: {output.ebrs_file}\\e[0m" >&2
            """




if should_run_rule("eba_analysis"):
    rule generate_msHSBs:
        input:
            # Use the copy_complete marker file to ensure EBA-input is ready
            eba_input_ready = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input", ".copy_complete"),
            script1 = os.path.join(SCRIPTS, "msHSBs_Finder.py"),
            script2 = os.path.join(SCRIPTS, "changeFORMAT_msHSBfile.py")
        output:
            final = os.path.join(config["eba_format"]["mshsbs_dir"], "msHSBs.txt"),
            done = touch(os.path.join(config["eba_format"]["mshsbs_dir"], "msHSBs_done.txt"))
        params:
            outdir = config["eba_format"]["mshsbs_dir"],
            eba_input_path = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input"),
            resolution = lambda wildcards: config["eba"]["p"]  # Get resolution from config
        log:
            "logs/generate_msHSBs.log"
        shell:
            """
            # Create output directory
            mkdir -p {params.outdir}
            
            # Check if the resolution directory exists
            RESOLUTION_DIR="{params.eba_input_path}/{params.resolution}"
            if [ ! -d "$RESOLUTION_DIR" ]; then
                echo "ERROR: Resolution directory $RESOLUTION_DIR does not exist" >&2
                echo "Available resolution directories:" >&2
                ls -d {params.eba_input_path}/*/ 2>/dev/null | grep -o '[0-9]*/' | sed 's|/$||' >&2 || true
                exit 1
            fi
            
            # Check for files in the resolution directory
            echo "Looking for files in $RESOLUTION_DIR" >&2
            FILE_COUNT=$(find "$RESOLUTION_DIR" -maxdepth 1 -type f -name "*_nonoverlap.txt" | wc -l)
            if [ "$FILE_COUNT" -eq 0 ]; then
                echo "ERROR: No *_nonoverlap.txt files found in $RESOLUTION_DIR" >&2
                echo "Files in directory:" >&2
                ls -la "$RESOLUTION_DIR/" >&2 || true
                exit 1
            fi
            
            echo "Found $FILE_COUNT *_nonoverlap.txt files in {params.resolution}k resolution directory" >&2
            
            # Combine all resolution files
            echo "Combining files from {params.resolution}k resolution..." | tee -a {log}
            cat "$RESOLUTION_DIR"/*_nonoverlap.txt > {params.outdir}/all_align.txt
            
            # Check if combined file was created
            if [ ! -s "{params.outdir}/all_align.txt" ]; then
                echo "ERROR: Combined file is empty" >&2
                echo "Combined file size: $(wc -c < "{params.outdir}/all_align.txt") bytes" >&2
                exit 1
            fi
            
            echo "Combined file created: {params.outdir}/all_align.txt" | tee -a {log}
            echo "File size: $(du -h "{params.outdir}/all_align.txt" | cut -f1)" | tee -a {log}
            echo "Line count: $(wc -l < "{params.outdir}/all_align.txt")" | tee -a {log}
            
            # Run msHSBs finder
            echo "Running msHSBs finder..." | tee -a {log}
            python3 {input.script1} \
                {params.outdir}/all_align.txt \
                {params.outdir}/inter-msHSBs.txt >> {log} 2>&1
                
            # Check if msHSBs finder succeeded
            if [ ! -f "{params.outdir}/inter-msHSBs.txt" ]; then
                echo "ERROR: msHSBs finder did not create output file" >&2
                echo "=== Last 20 lines of log ===" >&2
                tail -n 20 {log} >&2
                exit 1
            fi
            
            echo "Intermediate msHSBs file created: {params.outdir}/inter-msHSBs.txt" | tee -a {log}
            echo "File size: $(du -h "{params.outdir}/inter-msHSBs.txt" | cut -f1)" | tee -a {log}
            
            # Format the msHSBs file
            echo "Formatting msHSBs file..." | tee -a {log}
            python3 {input.script2} \
                {params.outdir}/inter-msHSBs.txt \
                {output.final} >> {log} 2>&1
            
            # Check if final file was created
            if [ ! -f "{output.final}" ]; then
                echo "ERROR: Formatted msHSBs file was not created" >&2
                echo "=== Last 20 lines of log ===" >&2
                tail -n 20 {log} >&2
                exit 1
            fi
            
            echo "Final msHSBs file created: {output.final}" | tee -a {log}
            echo "File size: $(du -h "{output.final}" | cut -f1)" | tee -a {log}
            echo "Line count: $(wc -l < "{output.final}")" | tee -a {log}
            
            # Show first few lines for verification
            echo "=== First 3 lines of msHSBs output ===" | tee -a {log}
            head -n 3 "{output.final}" | tee -a {log}

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m msHSBs detection in the dataset" >&2
            echo -e "\e[34;1mUsed resolution: {params.resolution}k from config\e[0m" >&2
            """
