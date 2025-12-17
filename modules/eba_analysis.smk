# EBA analysis module

if should_run_rule("eba_analysis"):
    rule make_eba_format:
        input:
            hundred_k = expand(f"{SYNTENY_RESULTS}/synteny_out/100000/{{sample}}_out/synteny_assign_done", sample=SAMPLES) if should_run_rule("synteny_processing") and SAMPLES else [],
            three_hundred_k = expand(f"{SYNTENY_RESULTS}/synteny_out/300000/{{sample}}_out/synteny_assign_done", sample=SAMPLES) if should_run_rule("synteny_processing") and SAMPLES else [],
            five_hundred_k = expand(f"{SYNTENY_RESULTS}/synteny_out/500000/{{sample}}_out/synteny_assign_done", sample=SAMPLES) if should_run_rule("synteny_processing") and SAMPLES else [],
            scaffolds = config["eba_format"]["scaffolds_file"] if config["eba_format"]["scaffolds_file"] and config["eba_format"]["scaffolds_file"] != "ALL_CHROMOSOMES" else []  # Only require file if it's a real path
        output:
            eba_dir = directory(config["eba_format"]["pre_EBA_dir"]),
            marker = touch(os.path.join(config["eba_format"]["pre_EBA_dir"], ".done_marker"))
        params:
            synteny_dir = f"{SYNTENY_RESULTS}/synteny_out",
            script_path = os.path.abspath(os.path.join(SCRIPTS, "EBA_format_maker10.sh")),
            output_dir = config["eba_format"]["pre_EBA_dir"],
            scaffolds_file = config["eba_format"]["scaffolds_file"]  # Pass the actual value
        log:
            "logs/make_eba_format.log"
#        message:
#            "Creating EBA format files from all synteny resolutions"
        shell:
            """
            # Create parent directory for marker
            mkdir -p {params.output_dir}
            
            # Verify script exists
            if [ ! -f "{params.script_path}" ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            # Build the command based on scaffolds_file value
            if [ "{params.scaffolds_file}" = "ALL_CHROMOSOMES" ]; then
                echo "Using Chromosomes for all species (special value: ALL_CHROMOSOMES)" >> {log}
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
                {params.script_path} \
                    -i "{params.synteny_dir}" \
                    -e "{output.eba_dir}" \
                    -s "{input.scaffolds}" > {log} 2>&1
            fi

            # Verify output
            if [ ! -d "{output.eba_dir}" ] || [ -z "$(ls -A {output.eba_dir})" ]; then
                echo "Error: EBA output failed" >&2 | tee {log}
                exit 1
            fi

            # Create completion marker
            touch "{output.marker}"
            """


if should_run_rule("eba_analysis"):
    rule resolve_overlaps:
        input:
            marker = os.path.join(config["eba_format"]["pre_EBA_dir"], ".done_marker")
        output:
            resolution_file = os.path.join(config["eba_format"]["pre_EBA_dir"], "resolution_complete.txt")
        params:
            script_path = os.path.abspath(os.path.join(SCRIPTS, "overlap_resolve-v7.py")),
            log_dir = "logs",
            eba_dir = config["eba_format"]["pre_EBA_dir"]
        log:
            "logs/resolve_overlaps.log"
#        message:
#            "Resolving overlaps in EBA format files"
        shell:
            """
            # Create log directory if needed
            mkdir -p {params.log_dir}

            # Verify script exists
            if [ ! -f "{params.script_path}" ]; then
                echo "Error: Script not found at {params.script_path}" >&2 | tee {log}
                exit 1
            fi

            # Execute overlap resolution on the directory
            python3 {params.script_path} {params.eba_dir} > {log} 2>&1

            # Create completion marker
            touch {output.resolution_file}
            """


if should_run_rule("eba_analysis"):
    rule prepare_eba_input:
        input:
            resolution_complete = os.path.join(config["eba_format"]["pre_EBA_dir"], "resolution_complete.txt")
        output:
            eba_input_dir = directory(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input")),  # Now in pre-EBA directory
            completion_marker = touch(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input/input_preparation_done.txt"))
        params:
            log_dir = "logs",
            source_dir = config["eba_format"]["pre_EBA_dir"]
        log:
            "logs/prepare_eba_input.log"
#        message:
#            "Preparing EBA input directory structure in pre-EBA location"
        shell:
            """
            # Create log directory if needed
            mkdir -p {params.log_dir}

            # Create main output directory in pre-EBA location
            mkdir -p {output.eba_input_dir}

            # Copy each resolution folder's contents
            for folder in 100 300 500; do
                echo "Processing $folder resolution" >> {log}
                mkdir -p {output.eba_input_dir}/$folder
                
                # Verify source directory exists
                if [ ! -d "{params.source_dir}/$folder" ]; then
                    echo "Error: Source directory {params.source_dir}/$folder not found" >> {log}
                    exit 1
                fi
                
                # Copy files (excluding _final.txt files)
                rsync -a --exclude='*_final.txt' "{params.source_dir}/$folder/" "{output.eba_input_dir}/$folder/" >> {log} 2>&1
                
                # Verify files were copied
                if [ -z "$(ls -A {output.eba_input_dir}/$folder)" ]; then
                    echo "Error: No files copied to {output.eba_input_dir}/$folder" >> {log}
                    exit 1
                fi
            done
            
            # Create completion marker
            touch {output.completion_marker}
            """

if should_run_rule("eba_analysis"):
    rule run_eba:
        input:
            eba_input_dir = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input"),  # Updated path
            classification_file = config["eba"]["c"],
            eba_script = os.path.join(workflow.basedir, config["eba_tools"]["eba_script_path"]),
            chr_size_file = os.path.join(os.path.dirname(config["eba_format"]["eba_input_dir"]), "chr_size.txt")
        output:
            eba_working_dir = directory(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING")),
            eba_out_dir = directory(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT")),
            final_classify = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT/300/EBA_OutFiles/final_classify.eba7"),
            marker = touch(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/eba_complete.marker"))
        params:
            n_threads = config["eba"]["n"],
            ref_name = config["eba"]["r"],
            classification_abs = lambda wildcards: os.path.abspath(config["eba"]["c"]),
            log_file_abs = lambda wildcards: os.path.abspath("logs/run_eba.log"),
            eba_input_abs = lambda wildcards: os.path.abspath(os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input")),  # Updated
            working_dir = lambda wildcards: os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING"),
            copy_dest_dir = config["eba_format"]["copy_destination_dir"],
            chr_size_abs = lambda wildcards: os.path.abspath(os.path.join(os.path.dirname(config["eba_format"]["eba_input_dir"]), "chr_size.txt"))
        log:
            "logs/run_eba.log"
        shell:
            """
            # Create log directory
            mkdir -p "$(dirname {params.log_file_abs})"
            
            # Create and set up the separate working directory
            mkdir -p {params.working_dir}
            cd {params.working_dir} || exit 1
            
            echo "Setting up EBA working directory in: $(pwd)" | tee -a {params.log_file_abs}
            echo "EBA input directory: {params.eba_input_abs}" | tee -a {params.log_file_abs}
            
            # Copy required files to working directory
            # 1. Copy the entire EBA-input structure from pre-EBA location
            echo "Copying EBA input files from pre-EBA directory..." | tee -a {params.log_file_abs}
            cp -r {params.eba_input_abs} . || exit 1
            
            # 2. Copy classification file
            echo "Copying classification file..." | tee -a {params.log_file_abs}
            cp {params.classification_abs} . || exit 1
            
            # 3. Copy chr_size.txt from the correct location
            echo "Copying chr_size.txt from: {params.chr_size_abs}" | tee -a {params.log_file_abs}
            if [ -f {input.chr_size_file} ]; then
                cp {input.chr_size_file} . || exit 1
                echo "Successfully copied chr_size.txt" | tee -a {params.log_file_abs}
            else
                echo "Error: chr_size.txt not found at {input.chr_size_file}" >&2 | tee -a {params.log_file_abs}
                echo "Available files in EBA directory:" >&2 | tee -a {params.log_file_abs}
                ls -la "$(dirname {input.chr_size_file})" >&2 | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify the EBA script exists
            if [ ! -f {input.eba_script} ]; then
                echo "Error: EBA script not found at {input.eba_script}" >&2 | tee -a {params.log_file_abs}
                echo "Current directory: $(pwd)" >&2 | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify input directory was copied correctly
            local_input_dir="$(basename {params.eba_input_abs})"
            if [ ! -d "$local_input_dir" ]; then
                echo "Error: EBA input directory not found at $local_input_dir" >&2 | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Run EBA in the working directory
            echo "Starting EBA analysis in working directory..." | tee -a {params.log_file_abs}
            echo "Working directory: $(pwd)" | tee -a {params.log_file_abs}
            echo "Local input directory: $local_input_dir" | tee -a {params.log_file_abs}
            echo "EBA script path: {input.eba_script}" | tee -a {params.log_file_abs}
            
            if ! perl {input.eba_script} \
                -n {params.n_threads} \
                -d "./$local_input_dir/" \
                -r {params.ref_name} \
                -p 300 \
                -t 20 \
                -c "./$(basename {params.classification_abs})" \
                -k >> {params.log_file_abs} 2>&1; then
                echo "EBA script failed with exit code $?" | tee -a {params.log_file_abs}
                echo "Last 20 lines of log:" | tee -a {params.log_file_abs}
                tail -n 20 {params.log_file_abs} | tee -a {params.log_file_abs}
                exit 1
            fi
            
            # Verify the expected output file was created
            if [ ! -f {output.final_classify} ]; then
                echo "Error: Expected output file {output.final_classify} was not created" | tee -a {params.log_file_abs}
                echo "Contents of EBA_WORKING directory:" | tee -a {params.log_file_abs}
                ls -la {params.working_dir} | tee -a {params.log_file_abs}
                if [ -d "EBA_OUT" ]; then
                    echo "Contents of EBA_OUT directory:" | tee -a {params.log_file_abs}
                    ls -la "EBA_OUT" | tee -a {params.log_file_abs}
                    if [ -d "EBA_OUT/300" ]; then
                        echo "Contents of EBA_OUT/300 directory:" | tee -a {params.log_file_abs}
                        ls -la "EBA_OUT/300" | tee -a {params.log_file_abs}
                        if [ -d "EBA_OUT/300/EBA_OutFiles" ]; then
                            echo "Contents of EBA_OUT/300/EBA_OutFiles directory:" | tee -a {params.log_file_abs}
                            ls -la "EBA_OUT/300/EBA_OutFiles" | tee -a {params.log_file_abs}
                        fi
                    fi
                fi
                exit 1
            fi
            
            # Create completion marker
            touch {output.marker}
            echo "EBA analysis completed successfully in working directory" | tee -a {params.log_file_abs}
            
            # Copy results to destination directory
            echo "Copying results to destination directory: {params.copy_dest_dir}" | tee -a {params.log_file_abs}
            mkdir -p {params.copy_dest_dir} || {{ echo "Failed to create destination directory" | tee -a {params.log_file_abs}; exit 1; }}
            
            # Copy the entire working directory or specific results
            cp -rv "EBA_OUT" "{params.copy_dest_dir}/" | tee -a {params.log_file_abs}
            
            # Copy any other important files that were generated
            for file in betaScore gaps_brks.stats species.sps sps.txt; do
                if [ -f "$file" ]; then
                    cp -v "$file" "{params.copy_dest_dir}/" | tee -a {params.log_file_abs}
                else
                    echo "Warning: $file not found in working directory" | tee -a {params.log_file_abs}
                fi
            done
            
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m EBRs finding in dataset" >&2
            """
    
if should_run_rule("eba_analysis"):
    rule process_breakpoints:
        input:
            final_classify = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/EBA_OUT/300/EBA_OutFiles/final_classify.eba7"),  # Updated path
            eba_marker = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA_WORKING/eba_complete.marker")  # Updated path
        output:
            done = os.path.join(config["eba_format"]["ebrs_dir"], "all_EBRs/intermediate_files/brkfile3_done.txt"),
            ebrs_file = os.path.join(config["eba_format"]["ebrs_dir"], "all_EBRs/5KBextenction_brk.txt")
        params:
            script_path = os.path.join(SCRIPTS, "antioverlaping2.py"),
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            intermediate_dir = os.path.join(config["eba_format"]["ebrs_dir"], "all_EBRs/intermediate_files"),
            ebrs_base_dir = config["eba_format"]["ebrs_dir"]
        log:
            "logs/process_breakpoints.log"
        shell:
            """
            mkdir -p {params.intermediate_dir}
            mkdir -p {params.ebrs_base_dir}/all_EBRs
            mkdir -p $(dirname {log})

            awk -F'\\t' -v col1="Chromosome" -v col2="Narrowest_brk" 'NR==1{{for(i=1;i<=NF;i++){{if($i==col1)c1=i; if($i==col2)c2=i}}}} c1&&c2&&NR>1{{print $c1 "\\t" $c2 "\\t" $(c2+3)}}' {input.final_classify} > {params.intermediate_dir}/brkfile
            awk -F '\\t' '$3 ~ /:0\\.9/ {{print}}' {params.intermediate_dir}/brkfile > {params.intermediate_dir}/brkfile2
            sed -i 's/<->/\\t/g' {params.intermediate_dir}/brkfile2
            awk -F '\\t' '{{print $1"\\t"$2-5000"\\t"$3+5000"\\t"$4}}' {params.intermediate_dir}/brkfile2 > {params.intermediate_dir}/brkfile3

            python3 {params.script_path} {params.intermediate_dir}/brkfile3 {output.ebrs_file} > {log} 2>&1
            touch {output.done}

            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m 5K extension of EBRs on both side" >&2
            """


if should_run_rule("eba_analysis"):
    rule generate_msHSBs:
        input:
            eba_input_dir = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input"),  # Updated to new location
            script1 = os.path.join(SCRIPTS, "msHSBs_Finder.py"),
            script2 = os.path.join(SCRIPTS, "changeFORMAT_msHSBfile.py")
        output:
            final = os.path.join(config["eba_format"]["mshsbs_dir"], "msHSBs.txt"),
            done = touch(os.path.join(config["eba_format"]["mshsbs_dir"], "msHSBs_done.txt"))
        params:
            outdir = config["eba_format"]["mshsbs_dir"],
            eba_input_path = os.path.join(config["eba_format"]["pre_EBA_dir"], "EBA-input")  # Updated path
        log:
            "logs/generate_msHSBs.log"
#        message:
#            "Generating microsynteny HSBs (msHSBs) from EBA input"
        shell:
            """
            mkdir -p {params.outdir}
            
            # Combine all 300k resolution files using the new config path
            cat {params.eba_input_path}/300/* > {params.outdir}/all_align.txt
            
            # Run msHSBs finder
            python3 {input.script1} \
                {params.outdir}/all_align.txt \
                {params.outdir}/inter-msHSBs.txt > {log} 2>&1
                
            # Format the msHSBs file
            python3 {input.script2} \
                {params.outdir}/inter-msHSBs.txt \
                {output.final} >> {log} 2>&1


            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m msHSBs detection in the dataset" >&2

            """

