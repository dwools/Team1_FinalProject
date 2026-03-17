# KEGNI (Knowledge-Enhanced Gene Network Inference)

## Overview
KEGNI is a multi-step bioinformatics pipeline that infers cell-type-specific gene regulatory networks. It leverages R for robust data preprocessing and biological validation, alongside a Python-based Masked Autoencoder (MAE) deep learning model to learn latent gene representations from single-cell RNA-seq data and prior knowledge graphs.

The pipeline assumes the project is housed within a specific directory structure: `Team1_FinalProject/morphic-mini-challenge/KEGNI`.

---

## Prerequisites

### R Environment
You must have R installed with the following packages:
* **Data Manipulation & Utilities:** `readr`, `dplyr`, `tibble`, `stringr`, `tidyr`, `tools`, `data.table`, `plyr`
* **Parallel Processing & Workflow:** `pacman`, `foreach`, `parallel`, `docstring`, `this.path`
* **Visualization:** `ggpubr`
* **Bioinformatics (Bioconductor):** `decoupleR`, `viper`

### Python Environment
You must have **Python 3.11** (I used **Python 3.11.9**) accessible via `python`, `python3`, or `py` in your system PATH. **Later versions of Python are incompatible with DGL**.
Required libraries include:
* `torch` (torch 2.3.0+cpu)
* `pandas` (pandas 3.0.1)
* `numpy` (numpy 1.26.4)
* `dgl` (dgl 2.2.1)

---

## Pipeline Architecture

The workflow is orchestrated sequentially by `main.R`:

1.  **Step 1: Input Generation (`STEP1_make_kegni_inputs.R`)**
    * Filters global expression matrices and metadata by matching cell barcodes.
    * Normalizes gene symbols to uppercase and collapses duplicates by sum.
    * Subsets KEGG and TRRUST Knowledge Graphs (KG) into isolated, cell-type-specific directories.
2.  **Step 2: Model Training (`train.py`)**
    * Ingests the processed expression CSVs and merged KG TSVs via a `KGEdataloader`.
    * Trains a graph-based Masked Autoencoder (MAE) to generate latent `scg_embedding.csv` files for each cell type.
3.  **Step 3: Edge Extraction (`STEP3_extract_kegni_tf_targets.R`)**
    * Reads the output embeddings and mathematically reconstructs the network.
    * Calculates edge weights utilizing $z = \tanh(z)$ followed by a similarity scoring matrix.
    * Filters and ranks interactions specifically from Transcription Factors (TFs) to target genes.
4.  **Step 4: Postprocessing (`STEP4_postprocess-CHOOSE-kegni.R`)**
    * Aggregates the ranked TF-target edges across all processed cell types into a single compressed CSV.
5.  **Step 5: Biological Validation (`STEP5_apply-viper-evaluation_kegni.R`)**
    * Applies the VIPER algorithm to compute TF activity scores.
    * Compares inferred network activity between Wild-Type (WT) cells and TF-Knockout (KO) cells using a one-sided Wilcoxon test.
    * Generates violin plots showing the condition distributions.

---

## Usage

The entire pipeline can be run using the orchestrator script. It automatically resolves paths, allowing you to run it from any directory.

### Basic Command
```bash
Rscript main.R
