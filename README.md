# Dock8 scRNAseq Analysis Code — Chen et al.

This repository contains all the code used in the analysis for the manuscript:

**"DOCK8-STAT3 Complexes Regulate Mucosal Immune Tolerance"**  
_Chen et al._

- **Manuscript:** [Link] (DOI: [..])
- **Contact:** Qian Chen — (email?) / Talal Chatila -  talal.chatila@childrens.harvard.ed
---

## Overview

This repository collects the scripts for processing and analyzing single-cell RNA sequencing (scRNAseq) data as described in the manuscript.  
Preliminary analyses were performed using [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) (10x Genomics).
The analysis workflow was performed using [Seurat](https://satijalab.org/seurat/) v5, with additional supporting R packages.

---

## Analysis Workflow

1. **Quality Control and Processing (Seurat)**
    - See: `script/1_QC_and_processing.R`
    - Import filtered matrices from CellRanger, quality control, filtering, normalization, integration, and clustering.

2. **Imputation and Differential Gene Expression (Seurat)**
    - See: `script/2_Imputation_DE_genes.R`
    - Imputation and Identification of differentially expressed genes.

3. **Visualization and Figure Generation (Seurat)**
    - See: `script/3_Plots.R`
    - Generation of all figures and plots included in the manuscript.


---

## Data Availability

- **Raw sequencing data:** Available at GEO/SRA: [GSE303265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303265)
- **Processed Seurat object:** Available on [Zenodo]() (DOI: )

---

## Requirements

All R package dependencies and the R version used for this analysis are specified in the `renv.lock` file included in this repository.

To recreate the exact software environment, install the [renv](https://rstudio.github.io/renv/) package in R and run:

```r
install.packages("renv")
renv::restore()
```

This will install the correct R packages and versions as used for the original analysis.

> **Note:**  
> The analysis was performed using R version 4.4.1  
> For details about all package versions, please refer to the `renv.lock` file.

---

## Usage

1. Clone the repository:
    ```bash
    git clone https://github.com/tuo-utente/nome-repo.git
    cd nome-repo
    ```

2. Follow the scripts in the order described above.  
   Please refer to comments within each script for details on requirements and parameters.

**Note:**  
For Cellranger usage, refer to the [10x Genomics documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest).

---

## Issues & Contact

If you encounter any problems or have questions, please open an issue in this repository or contact Qian Chen at [....](....) and Dr. Olga Lanzetta at [olga.lanzetta@edu.unige.it](mailto:olga.lanzetta@edu.unige.it)

---

