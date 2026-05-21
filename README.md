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

<img width="3352" height="1895" alt="Copy of Input folder (4)" src="https://github.com/user-attachments/assets/69d4d8cc-7dc1-4266-ada8-dce818350edc" />


## 🚀 Getting Started
### Prerequisite
##### 1. Sample input data
###### Fast track: 
Example input files with pre-computed alignments:https://doi.org/10.6084/m9.figshare.32315682

###### Slow track:
Example input files without pre-computed alignments:https://doi.org/10.6084/m9.figshare.32315892

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
-  If you see `command not found` — **install Conda.**

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
Install all required dependencies by creating the Conda environment from the provided YML file:
```bash
conda env create -f install_sybr_dependence.yml
```
##### 5. Activate Conda Environmanet
```bash
conda activate sybr
```
---

### Usage

The link for sample data given above. To use this sample data, downloa this data. Sybr's workflow is controlled via a `run_sybr_config.yaml` file, allowing you to selectively run different analysis stages. User can check the Documentation for detailed understanding about config settings.


##### 1. Input Folder Structure
```
inputs/
│
├── Ancestor_seq_recunstruction/          # Stages ④ chainNet + ⑤ Ancestor reconstruction
│   ├── LastZ_alignments/
│   │   ├── habrotrocha.axt               # one .axt per non-reference species
│   │   ├── ricciae.axt
│   │   ├── roseola.axt
│   │   ├── rotatoria.axt
│   │   └── vagaN.axt
│   ├── seq/
│   │   ├── habrotrocha.fa                # query + reference FASTA (.fa / .fasta / .fna)
│   │   ├── ricciae.fa
│   │   ├── roseola.fa
│   │   ├── rotatoria.fa
│   │   ├── vaga.fa
│   │   ├── vagaN.fa
│   │   └── refSpecies.fa                 # reference genome must also be present
│   ├── species_info.txt                  # fixed filename — see format below
│   └── tree.txt                          # Newick tree — same names as all other files
│
├── eba_analysis/                         # Stage ② EBA analysis
│   ├── chr_size.txt                      # fixed filename — reference chromosomes only
│   ├── classification.eba                # fixed filename — must contain lineage= entry
│   └── reference.fasta                   # reference genome FASTA 
│
├── enrichment_analysis/                  # Stage ③ Enrichment analysis
│   ├── 3kegg_annotationTOgenes.txt       # fixed filename — required when getenrich.r: "ko"
│   └── protein_annotation.tsv            # fixed filename — 5-column TSV, no header
│
└── synteny_processing/                   # Stage ① Synteny processing
    ├── all_sequence_lengths.txt          # fixed filename — generate with genome_length_maker.sh
    ├── Satsuma_alignments/
    │   ├── Adineta_ricciae.txt           # one .txt per query species (Satsuma output)
    │   ├── Adineta_vaga2.txt
    │   ├── Habrotrocha_rosa.txt
    │   ├── Philodina_roseola.txt
    │   └── Rotaria_rotatoria.txt
    └── Scaffolds.txt                     # scaffold-level species list (or set ALL_CHROMOSOMES)
```

> **Fixed filename** — must use the exact name shown; the pipeline looks for it by name.  
> **Variable** — any filename is accepted; only the extension matters. Any number of files allowed.

---

##### Input Files format

###### `LastZ_alignments/*.axt`
Run LastZ with these recommended parameters for non-vertebrate species:
```bash
lastz reference.fa[multiple] query.fa \
    C=0 E=30 H=2000 K=2200 L=2200 O=400 Y=3400 \
    --format=axt --output=SpeciesName.axt
```
Use the **HoxD55** scoring matrix for distant/non-vertebrate comparisons.  
The stem of each `.axt` filename (without extension) **must match** the names used in `seq/`, `species_info.txt`, and `tree.txt`.

###### `seq/*.fa`
Standard FASTA format. Accepts `.fa`, `.fasta`, `.fna`.  
Both query and reference genomes must be present. The format of header should be **>chr1, >chr2, >chr3, ...** for chromosome level genome assembly and **>scaf1, >scaf2, scaf3, ...** for scaffold level genome assemply.

###### `species_info.txt`
Three space/tab-separated fields per line, no header:

| Field | Values | Meaning |
|-------|--------|---------|
| Species name | string | must match `seq/`, `LastZ_alignments/`, and `tree.txt` exactly |
| Role | `0` / `1` / `2` | `0` = reference · `1` = descendant · `2` = outgroup |
| Assembly level | `1` / `0` | `1` = chromosome-scale · `0` = scaffold-level |

###### `tree.txt`
Newick format, single line. Must start with `(` and end with `;`.

###### `chr_size.txt`
Two-column TSV: `chromosome_name  size_bp`. Reference species only. Integer sizes.

###### `all_sequence_lengths.txt`
Three-column TSV: `sequence_name  length_bp  species`. All sequences from all genomes.  
Use the included `genome_length_maker.sh` script to generate this file.

###### `Satsuma_alignments/*.txt`
Eight-column TSV output from Satsuma, no header:
`query_chr  q_start  q_end  ref_chr  r_start  r_end  score  strand`  
Strand must be `+` or `-`.


---
##### 2. Config Setting for Sybr
- In run_stages section, user can choose the pipeline modules to run. in frount of each module mane, type **true** for activate the module and **false** for deactivate the module. 
- Use need to provide path for **base_input_dir** ( path of structured input folder) and 
**base_output_dir** (path of output folder).
- In **reference_species** provide the species name of referance only and in **referance_name** and **r:** provide genus and species name of referance.
- for **r:** in **getenrich** provide **KEGG Organism code**. user can check the availability of KEGG organims code in https://www.genome.jp/kegg/tables/br08606.html. if KEGG organism code is not availabe the user can use **ko**. for this option user need to provide kegg annotation file. to get the information about preparation of kegg annotation file refer to https://getenrich.igib.res.in/assets/files/getENRICH-documentation.pdf 
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

##### 3. Sybr help command to explor all the options
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
  -P, --paths FILE       Static paths file (default: pipeline_paths.yaml)
  -j, --cores N          Number of cores (default: all available)
  -t, --target RULE      Target rule (default: all)
  -l, --log FILE         Log output to file
  -u, --unlock           Unlock working directory
  -n, --dry-run          Dry run (simulate pipeline)
  -k, --keep-going       Keep going on independent job failures
  -v, --verbose          Verbose Snakemake output
  -s, --skip-validation  Skip input validation
  -C, --clean            Remove intermediate files after successful completion
  -w, --window-sizes     Comma-separated window sizes in bp (e.g., 100000,300000,500000)
  -p, --step-size        Step size in bp for synteny assignment (default: 30000)
  -h, --help             Show this help

Examples:
  ./sybr.sh -c config.yaml -j 8                         # Run with default settings
  ./sybr.sh --window-sizes 200000,400000 --step-size 50000  # Custom window sizes and step size
  ./sybr.sh --window-sizes 100000                        # Single window size
  ./sybr.sh --step-size 25000                           # Custom step size with default windows
  ./sybr.sh --clean                                      # Delete intermediate files after run

Note: Window sizes and step size only affect synteny_assign rules. If not specified,
      defaults are: window sizes = 100000,300000,500000 and step size = 30000.
```
##### 4. Basic Sybr Commanda
```bash
./sybr.sh -c run_sybr_config.yaml -j 8
```
or
```bash
./sybr.sh -j 8
```
or with clean **-C** flag
```
./sybr.sh -j 8 -C
```
##### 5. Custom window-sizes Sybr Commands
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
