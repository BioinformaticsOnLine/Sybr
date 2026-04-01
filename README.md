# Sybr

<p align="center">
  <img src="https://github.com/user-attachments/assets/c6400623-4d1f-41ad-a2b2-95245da9c983" alt="Sybr Logo" width="250">
</p>

**Sybr** is a powerful bioinformatics tool meticulously designed for the discovery of synteny blocks, the precise identification of evolutionary breakpoints, and robust ancestral genome reconstruction. It serves as an essential resource for researchers in comparative genomics and evolutionary biology, enabling deeper insights into genomic architecture and evolutionary relationships.

## ✨ Features

*   **Synteny Block Discovery:** Efficient algorithms to identify conserved genomic regions across multiple species.
*   **Evolutionary Breakpoint Identification:** Pinpoint the precise locations where genomic rearrangements have occurred using the EBA module.
*   **Ancestral Genome Reconstruction:** Algorithms to infer the genomic organization of common ancestors with DESCHRAMBLER.
*   **Enrichment Analysis:** Perform functional enrichment on identified regions using getENRICH.
*   **Modular Workflow:** Flexible `config.yaml` to run specific stages of the analysis.

## 🚀 Getting Started
### Prerequisite
##### 1. Sample input data
link for sample data https://figshare.com/s/efacec3b0f7589cfcc32

##### 2. Conda Installation

This project requires **Conda** to manage dependencies. If you already have Conda (Miniconda or Anaconda) installed, skip to this.

If Conda is not installed on your system, follow this step-by-step guide for Ubuntu:

👉 https://medium.com/@mustafa_kamal/a-step-by-step-guide-to-installing-conda-in-ubuntu-and-creating-an-environment-d4e49a73fc46

---

##### Check if Conda is Already Installed

Open a terminal and run:
```bash
conda --version
```

-  If you see something like `conda 24.x.x` — **Conda is already installed, skip the section below.**
-  If you see `command not found` — **follow the installation steps below.**

---

### Installation
##### 1. Clone the Repository
Clone the project from GitHub using the following command:

```bash
git clone https://github.com/your-username/sybr.git
```
##### 2. Set File Permissions
Grant the necessary permissions to all files and directories:
```bash
chmod -R 777 Sybr
```
##### 3. Navigate to Project Directory
Move into the project folder:
```bash
cd Sybr
```
##### 4. Create Conda Environment
Install all required dependencies by creating the Conda environment from the provided YAML file:
```bash
conda env create -f install_sybr_dependence.yml
```
##### 5. Activate Conda Environmanet
```bash
conda activate sybr
```
---

### Usage

The link for sample data given above. To use this sample data, downloa this data and move inside the Sybr folder. Sybr's workflow is controlled via a `run_sybr_config.yaml` file, allowing you to selectively run different analysis stages. User can check the Documentation for detailed understanding about config settings.

##### 1. Input Folder Structure
- the inputs folder has a fixed structure. for every selected module, there are subfoldes. and these subfolders contains specific input files with specific format, for reerance user can check teh link for example inputs folder. for detailed info of input files given in the documentation. 
```
./inputs/
├── Ancestor_seq_recunstruction
│   ├── LastZ_alignments   
│   │   ├── Species1.axt  --------------------------------------------------------------------------------------------------------------------------          
│   │   ├── Species2.axt         | in LastZ_alignments subfolder alignemnt files contains with parameters:
│   │   ├── Species3.axt         | lastz reference.fa[multiple] query.fa C=0 E=30 H=2000 K=2200 L=2200 O=400 Y=3400 --format=axt --output=query.axt
│   │   ├── Species4.axt         | file name should same as used in species_info.txt,tree.txt and in fasta files in seq folder.
│   │   └── Species5.axt  --------------------------------------------------------------------------------------------------------------------------
│   ├── seq
│   │   ├── Species1.fa   --------------------------------------------------------------------------------------------------------------------------
│   │   ├── Species2.fa          | in seq subfolder all the query and referance fasta files should present with same name as in LastZ_alignments.
│   │   ├── Species3.fa          |
│   │   ├── Species4.fa          |
│   │   ├── Species5.fa          |
│   │   └── refSpecies.fa --------------------------------------------------------------------------------------------------------------------------
│   ├── species_info.txt  --------------------------------------------------------------------------------------------------------------------------
│   └── tree.txt                 | species_info.txt file contains information of query species, referance species and genome assembly level. tree.txt
|                                | files conatine newick tree with same names as in species_info.txt and files in seq, LastZ_alignments subfolder. details explained in documentation
|                         --------------------------------------------------------------------------------------------------------------------------
|
|                         --------------------------------------------------------------------------------------------------------------------------
├── eba_analysis                 | eba_analysis subfolder should conain two files chr_size.txt and classification.eba which conatins information of chromosome size and classification. 
│   ├── chr_size.txt             | detailed information present in documentation.
│   └── classification.eba       |
|                         --------------------------------------------------------------------------------------------------------------------------
|
|                
├── enrichment_analysis   --------------------------------------------------------------------------------------------------------------------------
│   ├── 3kegg_annotationTOgenes.txt | enrichment_analysis subfolder should contains two files 3kegg_annotationTOgenes.txt and protein_annotation.tsv which contains inormation of
│   └── protein_annotation.tsv      | kegg annotation and coordinates of protein coding genes.
|                         --------------------------------------------------------------------------------------------------------------------------
|                        
└── synteny_processing    --------------------------------------------------------------------------------------------------------------------------
    ├── all_sequence_lengths.txt    | in synteny_processing subfilder conatine two files all_sequence_lengths.txt and Scaffolds.txt and one subfolder Satsuma_alignments.
    ├── Satsuma_alignments          | all_sequence_lengths.txt contains length of all the sequences of all the genomes, this file can be created with the help of genome_length_maker.sh
    │   ├── Genus_species1.txt      | script available in the tool. format of this file is expalined in documantation.
    │   ├── Genus_species2.txt      | Scaffolds.txt files conatins teh list of species which are not chrmosome leve genome assembly.
    │   ├── Genus_species3.txt      | Satsuma_alignments subfolder contains the alignment file with single referance generated by Satsuma tool.
    │   ├── Genus_species4.txt      |
    │   └── Genus_species5.txt      |
    └── Scaffolds.txt               |
                        ---------------------------------------------------------------------------------------------------------------------------
```


##### 1. Config Setting for Sybr
- In run_stages section, user can choose the pipeline modules to run. in frount of each module mane, type **true** for activate the module and **false** for deactivate the module
- 
```bash
# ─────────────────────────────────────────────
#  run_sybr_config.yaml  —  User-facing config
#  Edit this file to control pipeline behaviour
# ─────────────────────────────────────────────

# ── Base I/O directories ────────────────────
base_input_dir:  "/home/ajay/ajay_bhatia/my_writings/sybr/sybr_test4/inputs"
base_output_dir: "/home/ajay/ajay_bhatia/my_writings/sybr/sybr_test4/outputs"

# ── Pipeline stages to run ──────────────────
run_stages:
  synteny_processing: true
  eba_analysis: true
  enrichment_analysis: true
  chainNet_generation: true
  Ancestor_seq_recunstruction: true

# ── Species / reference names ───────────────
reference_name: "Adineta_vaga"
reference_species: "vaga"

# ── EBA parameters (user-facing) ───────────────────────────────
eba:
  n: 5           # number of EBA iterations
  r: "Adineta_vaga"  # reference species name
  p: 300         # primary resolution parameter

# ── Enrichment KEGG options ─────────────────
getenrich:
  r: "ko"  # KEGG code, or model-organism code

```

##### 2. Sybr help command to explor all the options
```bash
./sybr.sh -h


┏━┓╻ ╻┏┓ ┏━┓
┗━┓┗┳┛┣┻┓┣┳┛
┗━┛ ╹ ┗━┛╹┗╸
============

========================================
    Snakemake Pipeline
========================================

Usage: ./sybr.sh [OPTIONS]

Options:
  -c, --config FILE      Configuration file (default: run_sybr_config.yaml)
  -j, --cores N          Number of cores (default: all available)
  -t, --target RULE      Target rule (default: all)
  -l, --log FILE         Log output to file
  -u, --unlock           Unlock working directory
  -n, --dry-run          Dry run (simulate pipeline)
  -k, --keep-going       Keep going on independent job failures
  -v, --verbose          Verbose Snakemake output
  -s, --skip-validation  Skip input validation
  -w, --window-sizes     Comma-separated window sizes in bp (e.g., 100000,300000,500000)
  -p, --step-size        Step size in bp for synteny assignment (default: 30000)
  -h, --help             Show this help

Examples:
  ./sybr.sh -c config.yaml -j 8                         # Run with default settings
  ./sybr.sh --window-sizes 200000,400000 --step-size 50000  # Custom window sizes and step size
  ./sybr.sh --window-sizes 100000                        # Single window size
  ./sybr.sh --step-size 25000                           # Custom step size with default windows

Note: Window sizes and step size only affect synteny_assign rules. If not specified,
      defaults are: window sizes = 100000,300000,500000 and step size = 30000.
```
##### 3. Basic Sybr Commanda
```bash
./sybr.sh -c run_sybr_config.yaml -j 8
```
or
```bash
./sybr.sh -j 8
```
##### 4. Custom window-sizes Sybr Commands
```bash
./sybr.sh -w 200000,400000,500000 -j 8
```


## 📖 Documentation

A comprehensive user manual and API documentation can be found in the `docs/` directory or at [link_to_read_the_docs_if_applicable].

## 🤝 Contributing

We welcome contributions to Sybr! If you'd like to contribute, please read our [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to submit pull requests, report issues, and suggest new features.

## 🐛 Reporting Issues

Encountered a bug or have a feature request? Please open an issue on our [GitHub Issues page](https://github.com/BioinformaticsOnLine/sybr/issues).

## 📄 License

Sybr is released under the [LICENSE_NAME] License. See the [LICENSE](LICENSE) file for more details.

## 💬 Contact

For questions or support, please contact [BioinformaticsOnLine@gmail.com] or open a discussion on GitHub.
