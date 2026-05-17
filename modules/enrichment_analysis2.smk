# Enrichment analysis module

# Enrichment analysis rules (only if enabled)
if should_run_rule("enrichment_analysis"):
    rule msHSBs_background:
        input:
            msHSBs = rules.generate_msHSBs.output.final if should_run_rule("eba_analysis") else config["enrichment"].get("msHSBs_file", ""),
            annotation = config["enrichment"]["annotation_file"],
#            kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
            script = os.path.join(workflow.basedir, config["scripts"], "coordinate_overlap2.py")  # Fixed script path
        output:
            overlap = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "msHSBs_NCBI_genes_overlap.txt"),
            foreground = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "foreground.txt"),
            background = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "background.txt"),
#            copied_kegg = os.path.join(config["enrichment"].get("msHSBs_dir", "msHSBs"), "3kegg_annotationTOgenes.txt") if config["getenrich"]["r"] == "ko" else []
        params:
            outdir = config["enrichment"].get("msHSBs_dir", "msHSBs"),
            resource = config["getenrich"]["r"]  # Get resource type
        log:
            "logs/msHSBs_background.log"
        shell:
            """
            # Create log directory first
            mkdir -p "$(dirname {log})"
            
            echo "=== Generating foreground and background gene lists ===" | tee -a {log}
            
            # Create output directory
            mkdir -p {params.outdir}
            
            echo "Input msHSBs: {input.msHSBs}" | tee -a {log}
            echo "Input annotation: {input.annotation}" | tee -a {log}
            echo "Resource type: {params.resource}" | tee -a {log}
            echo "Script: {input.script}" | tee -a {log}
            
            # Verify input files exist
            if [ ! -f "{input.msHSBs}" ]; then
                echo "ERROR: msHSBs file not found: {input.msHSBs}" | tee -a {log}
                exit 1
            fi
            
            if [ ! -f "{input.annotation}" ]; then
                echo "ERROR: Annotation file not found: {input.annotation}" | tee -a {log}
                exit 1
            fi
            
            # Step 1: Run overlap script
            echo "Running coordinate overlap script..." | tee -a {log}
            python3 {input.script} {input.msHSBs} {input.annotation} {output.overlap} 0 1 2 3 1 2 >> {log} 2>&1
            
            if [ ! -f "{output.overlap}" ]; then
                echo "ERROR: Overlap script did not create output file" | tee -a {log}
                echo "=== Last 20 lines of log ===" | tee -a {log}
                tail -n 20 {log} | tee -a {log}
                exit 1
            fi
            
            echo "Overlap file created: {output.overlap}" | tee -a {log}
            echo "File size: $(du -h "{output.overlap}" | cut -f1)" | tee -a {log}
            echo "Line count: $(wc -l < "{output.overlap}")" | tee -a {log}
            
            # Step 2: Remove header from overlap file
            echo "Removing header from overlap file..." | tee -a {log}
            if [ -f "{output.overlap}" ]; then
                tail -n +2 "{output.overlap}" > "{params.outdir}/temp" && mv "{params.outdir}/temp" "{output.overlap}"
                echo "Header removed" | tee -a {log}
            else
                echo "ERROR: Overlap file missing before header removal" | tee -a {log}
                exit 1
            fi
            
            # Step 3: Create foreground list
            echo "Creating foreground gene list..." | tee -a {log}
            awk -F '\\t' '{{print $8}}' "{output.overlap}" > "{output.foreground}"
            sed -i '1i gene' "{output.foreground}"
            echo "Foreground file created: {output.foreground}" | tee -a {log}
            echo "Foreground line count: $(wc -l < "{output.foreground}")" | tee -a {log}
            
            # Step 4: Create background list from annotation file
            echo "Creating background from annotation file..." | tee -a {log}
            awk -F '\\t' '{{print $5}}' "{input.annotation}" | sort | uniq > "{output.background}"
            sed -i '/gene/d' "{output.background}"
            sed -i '1i gene' "{output.background}"
            echo "Background from annotation created: {output.background}" | tee -a {log}
            
            echo "Background line count: $(wc -l < "{output.background}")" | tee -a {log}
            
            # Show first few lines of outputs for verification
            echo "=== First 3 lines of foreground ===" | tee -a {log}
            head -n 3 "{output.foreground}" | tee -a {log}
            echo "=== First 3 lines of background ===" | tee -a {log}
            head -n 3 "{output.background}" | tee -a {log}
            
            echo "=== Foreground/background generation completed ===" | tee -a {log}
            echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m foreground/background generation for enrichment analysis in msHSBs" >&2
            """


if should_run_rule("enrichment_analysis"):
    rule EBRs_background:
        input:
            EBRs = rules.process_ebrs.output.ebrs_file if should_run_rule("eba_analysis") else [],
            annotation = config["enrichment"]["annotation_file"],
#            kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
            script = os.path.join(SCRIPTS, "coordinate_overlap2.py")  # ADD THIS LINE
        output:
            overlap = os.path.join(ebrs_dir, "EBRs_Enrichment_results/EBRs_NCBI_genes_overlap.txt"),
            foreground = os.path.join(ebrs_dir, "EBRs_Enrichment_results/foreground.txt"),
            background = os.path.join(ebrs_dir, "EBRs_Enrichment_results/background.txt"),
#            copied_kegg = os.path.join(ebrs_dir, "EBRs_Enrichment_results/3kegg_annotationTOgenes.txt") if config["getenrich"]["r"] == "ko" else []
        params:
            outdir = os.path.join(ebrs_dir, "EBRs_Enrichment_results"),
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

            # Step 4: Create background list from annotation file
            awk -F '\\t' '{{print $5}}' {input.annotation} | sort | uniq > {output.background}
            sed -i '/gene/d' {output.background}
            sed -i '1i gene' {output.background}
            """

if should_run_rule("enrichment_analysis"):
    rule split_EBRs:
        input:
            ebrs_file = rules.process_ebrs.output.ebrs_file
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
            split_done = os.path.join(ebrs_dir, "EBRs_split_done.txt"),
            annotation = config["enrichment"]["annotation_file"],
#            kegg = config["enrichment"]["kegg_file"] if config["getenrich"]["r"] == "ko" else [],  # Conditionally require KEGG
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
                        
                        # Step 4: Create background list from annotation file
                        if ! awk -F '\\t' '{{print $5}}' {params.annotation} | sort | uniq > "$lineage_dir/background.txt" 2>> {log} || \
                            ! sed -i '/gene/d' "$lineage_dir/background.txt" 2>> {log} || \
                            ! sed -i '1i gene' "$lineage_dir/background.txt" 2>> {log}; then
                                echo "Error creating background for $lineage_dir" >&2 | tee -a {log}
                                exit 1
                            fi
                    else
                        echo "Warning: No overlap file found in $lineage_dir" >> {log}
                    fi
                fi
            done

            touch {output.done}
            """
