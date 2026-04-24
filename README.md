# SPARC: Single Cell Polyploidy analysis via RNA Characterization
![SPARC Graphical Abstract](results/pc_plots_v4/SPARC_graphical_abstract.png)
A computational pipeline for identifying polyploid cancer cells in single-cell RNA-seq data.
<small>Graphical abstract created via BioRender (License under results/pc_plots_v4)</small>
---

## Overview

SPARC combines fast CNV inference from scRNA-seq count data (via [CopyKat](https://github.com/navinlabcode/copykat)) with a suite of 5 ML/DL classifiers to predict the probability that each cell is a polyploid cancer (PC) cell.

---

## Installation
 
### 1. Clone the repository
 
```bash
git clone https://github.com/dsaha0295/SPARC.git  
cd SPARC
```
 
### 2. Install R dependencies
 
Requires **R ≥ 4.1**. Install the following R packages:
 
```r
# CRAN packages for single cell processing and data wrangling
install.packages(c("Seurat", "optparse", "tidyverse", "Matrix"))
 
# CopyKat for CNV inference 
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("navinlabcode/copykat")
```
 
### 3. Install Python dependencies
 
Requires **Python ≥ 3.8**. We would recommend creating a virtual environment e.g Conda:
 
```bash
pip install pandas numpy scikit-learn joblib
```
 
### 4. WDL pipeline dependencies
 
To run the WDL pipeline you'll also need:
 
- [Cromwell](https://github.com/broadinstitute/cromwell/releases) (v86 recommended) — download the `.jar` file
- Java 11+ (e.g. `openjdk:11`)
- Docker (for containerized execution)
---
 

## Pipeline Steps

Run the following scripts **in order**:

1. `sparc_cnv_inference.R` — CNV inference and feature extraction
2. `sparc_run_models.py` — Run ML/DL classifiers and output per-cell predictions

---

## Usage

### Step 1 — CNV Estimation (`sparc_cnv_inference.R`)

```bash
Rscript sparc_cnv_inference.R \
  --seurat        /path/to/object.rds \
  --epi_count     4000 \
  --tme_count     4000 \
  --seed          42 \
  --celltype_col  coarse_ano \
  --epi_label     "Epi_Neuroendo" \
  --out_dir       /path/to/output/
```

> **Note:** Your Seurat object metadata must include a column named `cell` with unique barcodes for each cell.

#### Required Arguments

| Argument | Description |
|---|---|
| `--seurat` | Path to input Seurat `.rds` file |
| `--epi_count` | Number of epithelial cells to downsample |
| `--tme_count` | Number of TME / reference cells to downsample |
| `--seed` | Random seed for reproducibility |
| `--celltype_col` | Metadata column used to identify cell types |
| `--epi_label` | Value in `--celltype_col` that marks epithelial/tumor cells |
| `--out_dir` | Output directory for all results |

#### Optional Arguments

| Argument | Default | Description |
|---|---|---|
| `--spikein` | *(none)* | Path to an external spike-in Seurat `.rds` used as the diploid reference. If omitted, non-epithelial cells from the main object are used instead. |
| `--spikein_malignancy_col` | `malignancy` | Metadata column in spike-in to filter on |
| `--spikein_malignancy_val` | `Non-Malignant` | Value to keep in the malignancy column |
| `--spikein_id_col` | — | ID column in spike-in object |
| `--spikein_id_val` | — | ID value to filter spike-in on |
| `--mt_col` | `percent.mt` | Metadata column for mitochondrial percentage |
| `--ncores` | `8` | Number of cores passed to CopyKat |

---

### Step 2 — Run Classifiers (`sparc_run_models.py`)

```bash
sparc_run_models.py \
  --out               paccs_predictions.csv \
  --model_dir         models/ \
  --cell_id_col       cell \
  input.csv
```

The `input.csv` is the output from Step 1.

#### Arguments

| Argument | Default | Description |
|---|---|---|
| `input_csv` | *(required)* | Feature CSV output from `sparc_cnv_inference.R` |
| `--out` | `paccs_predictions.csv` | Output CSV file for predictions |
| `--model_dir` | `models/` | Directory containing saved `.pkl` models and `scaler.pkl` |
| `--features` | See below | Feature columns to use for inference |
| `--drop_na` | `True` | Drop rows with missing values before inference |
| `--cell_id_col` | *(none)* | Column to use as a cell identifier in the output |

**Default features:** `Percent.MT`, `S.Score`, `cnv_total`, `nCount_RNA`, `nFeature_RNA`, `G2M.Score`

---

## Running the WDL Pipeline

### From the Command Line

```bash
java -Xms4g -Xmx12g \
  -Dconfig.file=cromwell_compute1.conf \
  -jar cromwell-86.jar \
  run sparc.wdl \
  --inputs sparc.json
```

### From an HPC (e.g LSF batch script for Compute1 at WashU)

```bash
bsub \
  -J sparc \
  -G compute-christophermaher \
  -g /saha.d/max100 \
  -q general \
  -n 4 \
  -R 'select[mem>16000] rusage[mem=16000]' \
  -M 16000000 \
  -oo Logs/sparc.out \
  -eo Logs/sparc.err \
  -a 'docker(openjdk:11.0.11-jdk-slim)' \
  /usr/local/openjdk-11/bin/java -Xms4g -Xmx12g \
    -Dconfig.file=cromwell_compute1.conf \
    -jar cromwell-86.jar \
    run sparc.wdl \
    --inputs sparc.json
```

See `submit_wdl.sh` for a ready-to-use submission script.

---

## Repository Structure

```
├── src/
│   └── pipeline/                   # WDL pipeline, JSONs, Dockerfiles, configs
|       └──sparc_cnv_inference.R    # Step 1: CNV inference
│       └──sparc_run_models.py      # Step 2: ML/DL classification
├── data/                           # Example CSV files for replication
└── models/                         # Saved .pkl model files
```

Additional processed data generated from cell line experiments during this study are available from the corresponding authors upon reasonable request. Published data from the scRNA-seq mCRPC cohort used in this analysis was accessed through the Gene Expression Omnibus repository (Accession ID GSE264573). 

---

## Citation

TBD