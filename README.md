# ImmunoTyper2 Paper Methods

This repository contains the methods and workflows used in the ImmunoTyper2 paper. The experiments are organized into two main divisions.

## Setup Instructions

1. Create and activate the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate immunotyper-methods
   ```

2. Install ImmunoTyper2 in the conda environment:
   ```bash
   pip install immunotyper-2.0-py3-none-any.whl
   ```

## Repository Structure

### 1. 1KGP_Trios
Analysis of 1000 Genomes Project trio samples using ImmunoTyper2.

#### Steps to Reproduce:
1. Fetch metadata:
   ```bash
   nextflow run fetch_trio_metadata.nf
   ```

2. Download samples:
   ```bash
   # Set your Globus target UUID
   export TARGET_UUID=<Your target UUID>
   ./globus_download.sh
   ```

3. Run ImmunoTyper2:
   ```bash
   nextflow run run_immunotyper.nf
   ```

#### Pre-computed Results
The ImmunoTyper2 outputs for 1KGP_Trios are available on Zenodo:
- [Download Link](https://zenodo.org/records/14455863/files/immunotyper-1kgp-trios-output.tar.gz?download=1)
- Extract to: `1KGP_Trios/immunotyper-output`

### 2. ImmunoTyper-assembly-benchmarking
Benchmarking experiments using genome assemblies.

## Requirements
- Nextflow
- Conda
- Globus CLI (for 1KGP data download)
- ImmunoTyper2