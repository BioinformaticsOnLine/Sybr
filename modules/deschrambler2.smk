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
            # Change to the main output directory
            cd {params.output_dir}
            
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

            # Copy tree file
            cp {params.tree_file} tree.txt
            
            # Copy Makefile.SFs from pipeline tools directory
            cp {input.makefile_sfs} .
            
            # Debug: Check if files exist
            echo "Checking files in output directory:"
            ls -la
            
            # Run DESCHRAMBLER from pipeline tools directory
            echo "Running DESCHRAMBLER..."
            perl {input.deschrambler_script} params.txt
            
            # After DESCHRAMBLER completes, check if APCFs.300K directory was created
            if [ -d "APCFs.300K" ]; then
                echo "DESCHRAMBLER completed successfully. APCFs.300K directory created."
            else
                echo "WARNING: APCFs.300K directory was not created. DESCHRAMBLER may have failed."
                exit 1
            fi
            
            # Mark completion
            touch {output}
            """
