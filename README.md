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

### Prerequisites

Sybr integrates several powerful modules (getENRICH, EBA, DESCHRAMBLER), each with its own set of dependencies. Ensure you have the following installed:

#### For `getENRICH` (R dependencies)

```R
install.packages(c("jsonlite", "dplyr", "tidyverse", "ggplot2", "ontologyIndex", "plotly"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "UpSetR", "pheatmap", "visNetwork", "enrichplot", "ComplexHeatmap", "circlize", "pathview"))
```

#### For `EBA` (Perl dependencies)

EBA primarily relies on Perl core modules and its own internal libraries.
*   **Perl Core Modules:** `strict`, `warnings`, `English`, `FileHandle`, `Getopt::Long`, `Pod::Usage`, `File::Path`, `Cwd`, `File::Find`, `File::Copy`, `File::Basename`.
*   **Perl CPAN Modules:**
    *   `Math::Round` (install via `cpan Math::Round` or your package manager)
    *   `List::Compare` (install via `cpan List::Compare` or your package manager)
*   **EBA Internal Libraries:** Located in `EBALib/` and `bin/EBALib/` within the EBA tool directory.

#### For `DESCHRAMBLER` (Python & Perl)

DESCHRAMBLER requires a dedicated Conda environment and specific Perl modules.

1.  **Conda Environment Setup:**
    ```bash
    conda create -n deschrambler python=3.7
    conda activate deschrambler
    ```
2.  **BioPerl Installation:**
    ```bash
    conda install -c bioconda perl-bioperl
    ```
3.  **Find BioPerl Path (adjust `~/miniforge3` to your conda base path):**
    ```bash
    find ~/miniforge3 -name "TreeIO.pm" 2>/dev/null
    # Example output: /home/user/miniforge3/envs/deschrambler/lib/perl5/site_perl/5.22.0/Bio/TreeIO.pm
    ```
4.  **Create Symbolic Link for BioPerl in DESCHRAMBLER:**
    ```bash
    # Navigate to your DESCHRAMBLER installation directory first
    cd /path/to/DESCHRAMBLER_tool_directory 
    rm -rf lib/perl/Bio
    ln -s /home/user/miniforge3/envs/deschrambler/lib/perl5/site_perl/5.22.0/Bio lib/perl/Bio
    ```
    *(**Important:** Replace `/home/user/miniforge3/envs/deschrambler/lib/perl5/site_perl/5.22.0/Bio` with the actual path found in step 3, up to `/Bio`)*
5.  **Verify BioPerl Installation:**
    ```bash
    perl -Ilib/perl -MBio::TreeIO -e 'print "SUCCESS: BioPerl working\n"'
    ```
6.  **Other DESCHRAMBLER dependencies:**
    *   `Perl 5.22+`
    *   `make`
    *   `build-essential`

#### Snakemake Workflow Dependencies

The overall Sybr workflow relies on Snakemake and several UCSC tools, installed via Bioconda within a Python 3.11 environment.

```yaml
# Recommended conda environment for the main Sybr workflow
python=3.11
bioconda::ucsc-fatotwobit
bioconda::ucsc-axtchain
bioconda::ucsc-chainsplit
bioconda::ucsc-fasize
bioconda::ucsc-chainprenet
bioconda::ucsc-chainnet
bioconda::ucsc-netsyntenic
bioconda::snakemake=9.9.0
```
**(It's highly recommended to create a dedicated `environment.yml` for this main workflow.)**

### Installation

To get Sybr set up, clone the repository. After cloning, you will need to place the `getENRICH`, `EBA3.0`, and `DESCHRAMBLER` tools into your `tools/` directory as specified in the usage sections.

```bash
git clone https://github.com/your-username/sybr.git
cd sybr

# Ensure sub-tools are placed in the 'tools/' directory:
# sybr/
# ├── tools/
# │   ├── getENRICH/
# │   ├── EBA3.0/
# │   └── DESCHRAMBLER/
# └── ... (Sybr's main files)
```
**(You'll likely need to provide specific instructions or scripts on how to obtain/install `getENRICH`, `EBA3.0`, and `DESCHRAMBLER` and place them correctly, as these seem to be external tools Sybr orchestrates.)**

### Usage

Sybr's workflow is controlled via a `config.yaml` file, allowing you to selectively run different analysis stages. Below are examples for the three main functionalities.

#### 1. Find Evolutionary Breakpoints (EBRs) and Multi-species Homologous Synteny Blocks (msHSBs)

This stage leverages `EBA3.0` to identify genomic rearrangements and conserved syntenic regions.

**`config.yaml` configuration:**
```yaml
run_stages:
  synteny_processing: true
  eba_analysis: true
  enrichment_analysis: false
  alignment_processing: false
  deschrambler: false
```

**Inputs:**
*   `alignments/`: Directory containing alignment files generated by the `SatsumaSynteny2` tool.
*   `Scaffolds.txt`: A plain text file listing the names of organisms with scaffold-level assembly.
*   `Classification.eba`: (Describe purpose if not self-explanatory, e.g., a pre-computed classification file)
*   `Reference name`: Specified in the `config.yaml` (e.g., `Adineta vaga`).
*   `EBA3.0/`: The EBA tool directory, expected to be located in `sybr/tools/EBA3.0/`.

**Outputs:**
*   `EBRs/`: Output directory containing identified Evolutionary Breakpoints.
*   `msHSBs/`: Output directory containing Multi-species Homologous Synteny Blocks.

#### 2. Perform Enrichment Analysis of EBRs and msHSBs

This stage uses `getENRICH` to functionally characterize the genomic regions identified in the previous step.

**`config.yaml` configuration:**
```yaml
run_stages:
  synteny_processing: true
  eba_analysis: true
  enrichment_analysis: true
  alignment_processing: false
  deschrambler: false
```

**Inputs:**
*   `alignments/`: (Same as above)
*   `Scaffolds.txt`: (Same as above)
*   `Classification.eba`: (Same as above)
*   `Reference name`: Specified in the `config.yaml` (e.g., `Adineta vaga`).
*   `getENRICH_input/`: Directory containing KEGG and NCBI annotation files required by `getENRICH`.
*   `getENRICH/`: The getENRICH tool directory, expected to be located in `sybr/tools/getENRICH/`.

**Outputs:**
*   `EBRs/`: Output directory with enrichment analysis results for EBRs.
*   `msHSBs/`: Output directory with enrichment analysis results for msHSBs.

#### 3. Run DESCHRAMBLER for Ancestral Reconstruction

This stage utilizes `DESCHRAMBLER` to perform ancestral genome reconstruction.

**`config.yaml` configuration:**
```yaml
run_stages:
  synteny_processing: false
  eba_analysis: false
  enrichment_analysis: false
  alignment_processing: true
  deschrambler: true
```

**Inputs:**
*   `alignments2/`: Directory containing alignment files generated by the `LastZ` tool.
*   `DESCHRAMBLER-inputs/`: Directory containing `config.SFs`, `params.txt`, and `tree.txt` for DESCHRAMBLER.
*   `reference/`: Fasta file of the reference genome.
*   `seq/`: Directory containing fasta files of query genomes.
*   `DESCHRAMBLER/`: The DESCHRAMBLER tool directory, expected to be located in `sybr/tools/DESCHRAMBLER/`.

**Outputs:**
*   `APCFs.300K/`: Output directory containing ancestral reconstruction results.

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
