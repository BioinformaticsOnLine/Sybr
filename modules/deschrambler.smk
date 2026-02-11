# Deschrambler module

# Deschrambler rules (only if enabled)
if should_run_rule("Ancestor_seq_recunstruction"):
    rule run_deschrambler:
        input:
            reorg_done = f"{OUTPUT_DIR}/data/{REFERENCE}/reorganized.done",
            tree_file = config["deschrambler"]["tree_file"],
            species_file = config["deschrambler"]["species_file"],
            makefile_sfs = os.path.join(workflow.basedir, config["deschrambler"]["pipeline_tool_dir"], "Makefile.SFs"),
            deschrambler_script = os.path.join(workflow.basedir, config["deschrambler"]["pipeline_tool_dir"], "DESCHRAMBLER.pl")
        output:
            f"{OUTPUT_DIR}/deschrambler.done"
        params:
            output_dir = OUTPUT_DIR,
            reference = REFERENCE,
            pipeline_tool_dir = os.path.join(workflow.basedir, config["deschrambler"]["pipeline_tool_dir"]),
            data_dir = f"{OUTPUT_DIR}/data",
            species_file = config["deschrambler"]["species_file"],
            tree_file = config["deschrambler"]["tree_file"]
        resources:
            mem_mb=16000,
            time_min=120
        shell:
            """
            echo -e "\e[33;1mStarting ancestral sequence reconstruction using DESCHRAMBLER\e[0m" >&2
            
            # Change to the main output directory
            cd {params.output_dir}
            
            echo -e "\e[33;1mCreating configuration files for DESCHRAMBLER\e[0m" >&2
            
            # Read species information from file
            SPECIES_INFO=$(cat {params.species_file})
            
            # Create config.SFs file
            cat > config.SFs << EOF
# path of a net files 
>netdir
{params.data_dir}

# path of a chain files 
>chaindir
{params.data_dir}

# species-name tag1 tag2
# tag1 (0: ref-species, 1: descendents, 2: outgroup)
# tag2 (1: chromosome assembly, 0: others)
>species
$SPECIES_INFO

# block resolution (bp): DO NOT CHANGE
>resolution
<resolutionwillbechanged>
EOF

            echo -e "\e[32;1m✓\e[0m Created config.SFs file" >&2

            # Create params.txt file
            cat > params.txt << EOF
# Reference species
REFSPC={params.reference}

# Output directory
OUTPUTDIR=APCFs.300K

# Block resolution (bp)
RESOLUTION=300000

# Newick tree file
TREEFILE=tree.txt

# Minimum adjacency scores
MINADJSCR=0.0001

# Config and make files for syntenic fragment construction
CONFIGSFSFILE=config.SFs
MAKESFSFILE=Makefile.SFs
EOF

            echo -e "\e[32;1m✓\e[0m Created params.txt file" >&2

            # Copy tree file
            cp {params.tree_file} tree.txt
            echo -e "\e[32;1m✓\e[0m Copied tree file: $(basename {params.tree_file})" >&2
            
            # Copy Makefile.SFs from pipeline tools directory
            cp {input.makefile_sfs} .
            echo -e "\e[32;1m✓\e[0m Copied Makefile.SFs" >&2
            
            # Debug: Check if files exist
            echo -e "\e[33;1mVerifying input files in output directory:\e[0m" >&2
            echo "Files in {params.output_dir}:" >&2
            ls -la >&2
            
            # Run DESCHRAMBLER from pipeline tools directory
            echo -e "\e[33;1mRunning DESCHRAMBLER for ancestral sequence reconstruction...\e[0m" >&2
            echo -e "\e[33;1mThis may take several minutes...\e[0m" >&2
            
            perl {input.deschrambler_script} params.txt
            
            DESCHRAMBLER_EXIT_CODE=$?
            
            # After DESCHRAMBLER completes, check if APCFs.300K directory was created
            if [ $DESCHRAMBLER_EXIT_CODE -eq 0 ] && [ -d "APCFs.300K" ]; then
                echo -e "\e[32;1m✔\e[0m \e[33;1mCompleted\e[0m ancestral sequence reconstruction using DESCHRAMBLER" >&2
                echo -e "\e[32;1m✓\e[0m APCFs.300K directory created successfully" >&2
            else
                echo -e "\e[31;1m✗\e[0m \e[33;1mFailed ancestral sequence reconstruction using DESCHRAMBLER\e[0m" >&2
                if [ ! -d "APCFs.300K" ]; then
                    echo -e "\e[31;1m✗\e[0m APCFs.300K directory was not created" >&2
                fi
                echo -e "\e[31;1m✗\e[0m DESCHRAMBLER exited with code: $DESCHRAMBLER_EXIT_CODE" >&2
                exit 1
            fi
            
            # Mark completion
            touch {output}
            echo -e "\e[32;1m✓\e[0m Created completion marker: $(basename {output})" >&2
            """
