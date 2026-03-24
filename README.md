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
link for sample data https://figshare.com/s/49a5e76634a3683362f5

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

##### 1. Sybr help command to explor all the options
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
##### 2. Basic Sybr Commanda
```bash
./sybr.sh -c run_sybr_config.yaml -j 8
```
or
```bash
./sybr.sh -j 8
```
##### 3. Custom window-sizes Sybr Commands
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
