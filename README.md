# Oral-Gut Microbiome PhD Analysis Scripts

This repository contains scripts and configuration files used for analysing and processing microbiome data in my PhD research, focusing specifically on the oral-gut microbiome interactions in inflammatory bowel disease.

# 1. Meta-analysis of the oral microbiota in IBD

## Primer Scripts

These shell scripts are written to remove primer sequences from 16S rRNA sequencing reads before running the nf-core ampliseq pipeline. The cutadapt package is used.

### `Ab_primer.sh`

* **Purpose:** Removes primer sequences from 16S rRNA sequencing reads from the Ab 2022 study.
* **Usage:**

```bash
./Ab_primer.sh
```

### `Mo_primer.sh`

* **Purpose:** Removes primer sequences from 16S rRNA sequencing reads from the Mo 2022 study.
* **Usage:**

```bash
./Mo_primer.sh
```

### `Xi_primer.sh`

* **Purpose:** Removes primer sequences from 16S rRNA sequencing reads from the Xi 2022 study.
* **Usage:**

```bash
./Xi_primer.sh
```

## Ampliseq nf-core pipeline

### `nfcore_ampliseq_command.txt`

* The command used for running the [`nf-core/ampliseq`](https://nf-co.re/ampliseq) pipeline via Nextflow on Eddie (high performance compute cluster). 

## R Analysis

### `oral_meta-analysis_phyloseq.R`

* An R script for conducting meta-analysis of oral microbiome datasets using Phyloseq. Includes data import, preprocessing, normalisation, and diversity analyses.


# 2. Taxonomic profiling the oral microbiota in IBD using shotgun metagenomics

These configuration files and scripts are related to running bioinformatics pipelines using Nextflow on Azure Batch Cloud.

## Taxprofiler nf-core pipeline

### `nfcore_taxprofiler_command.txt`

* The command used for running the [`nf-core/taxprofiler`](https://nf-co.re/taxprofiler) pipeline via Nextflow on Eddie (high performance compute cluster). 

## Azure batch configuation

### `azurebatch.config`

* Contains configuration settings for running Nextflow workflows using Azure Batch, including computational resource specifications.

---

### Usage and Requirements

Ensure you have the necessary permissions and dependencies installed (R packages, Nextflow, Azure CLI) before running these scripts.

For further details or questions, please contact:

* Robert Whelan
* robert.whelan@ed.ac.uk

© 2025 Robert Whelan
