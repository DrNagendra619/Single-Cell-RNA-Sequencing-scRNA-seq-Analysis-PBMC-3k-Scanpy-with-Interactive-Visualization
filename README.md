# Single-Cell-RNA-Sequencing-scRNA-seq-Analysis-PBMC-3k-Scanpy-with-Interactive-Visualization
Comprehensive Single Cell RNA Sequencing (scRNA seq) Analysis of PBMC 3k using Scanpy with Interactive Visualization
# Comprehensive Single-Cell RNA Sequencing (scRNA-seq) Analysis of PBMC 3k

This Jupyter Notebook provides a complete, end-to-end pipeline for analyzing **Peripheral Blood Mononuclear Cells (PBMC) 3k scRNA-seq data** using the `scanpy` ecosystem. It emphasizes standard quality control (QC), filtering, normalization, dimensionality reduction, clustering, and marker gene identification, with a focus on **interactive visualization** using `plotly` for enhanced data exploration.

---

## 🚀 Overview

The script processes raw gene expression data to identify distinct cell populations within the PBMC 3k dataset. It leverages `scanpy` for core single-cell analysis and `plotly` for generating dynamic, zoomable plots that are ideal for in-depth inspection of QC metrics, highly variable genes, and UMAP clustering results.

---

## ✨ Key Features

* **Standardized QC:** Calculates key quality metrics including mitochondrial percentage (`pct_counts_mt`), total UMI counts (`total_counts`), and number of genes detected (`n_genes_by_counts`).
* **Rigorous Filtering:** Implements filtering thresholds to remove low-quality cells and uninformative genes:
    * Minimum 200 genes per cell.
    * Maximum 2500 genes per cell (to exclude potential doublets).
    * Maximum 5% mitochondrial counts per cell (to exclude stressed/dying cells).
    * Minimum 3 cells per gene.
* **Dimension Reduction:** Performs Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) to visualize data structure, retaining **40 principal components** based on variance ratio insights.
* **Clustering:** Identifies distinct cell populations using the **Leiden clustering algorithm** (found 6 clusters with `resolution=0.5`).
* **Marker Gene Analysis:** Ranks the most differentially expressed genes for each cluster using the **Wilcoxon test**.
* **Interactive Visualization:** Generates dynamic Plotly scatter and box plots for hands-on exploration of QC and clustering results.

---

## 🛠️ Requirements

The script is a Python Jupyter Notebook designed for single-cell analysis.

### Required Python Libraries

Install the necessary libraries using `pip`:

```bash
pip install scanpy leidenalg plotly pandas matplotlib seaborn
Key Dependencies Used:

    scanpy (sc)

    leidenalg

    plotly.express (px)

    pandas

    matplotlib.pyplot (plt)

    seaborn (sns)

    anndata (ad)

🚀 Usage
1. Setup and Environment

The notebook begins by installing all required packages and setting general plot parameters.
Python

# Installation & Imports
!pip install scanpy leidenalg plotly pandas matplotlib seaborn

import scanpy as sc
import anndata as ad
# ... other imports

2. Data Loading

The script automatically downloads and loads the raw PBMC 3k dataset (approximately 5.6MB) for immediate analysis:
Python

# Load Data
adata = sc.datasets.pbmc3k()

3. Execution

Execute the notebook cells sequentially to run the full analysis workflow. The main steps are:

    QC Calculation: sc.pp.calculate_qc_metrics is used to quantify the quality of cells.

    Filtering: Cells are filtered based on stringent QC thresholds (min_genes=200, n_genes_by_counts < 2500, pct_counts_mt < 5, min_cells=3).

    Preprocessing: Data is normalized and log-transformed.

    Highly Variable Genes: sc.pp.highly_variable_genes is run to select the top 2000 genes for downstream analysis.

    Dimensionality Reduction: PCA and UMAP embeddings are computed.

    Clustering & Marker Genes: Leiden clustering identifies cell populations, and sc.tl.rank_genes_groups finds the distinguishing marker genes.

4. Interpretation (Interactive Plots)

The notebook provides several interactive visualizations (using plotly.express) at critical steps. These plots allow users to zoom, pan, and hover over individual cells or genes to inspect their values, which is essential for validating the analysis:

    Interactive QC Plot: Visualizes library size vs. detected genes, colored by mitochondrial percentage.

    Interactive HVG Plot: Shows mean expression vs. dispersion, highlighting the 2000 highly variable genes selected for analysis.

    Interactive 3D PCA Plot: Displays cells in the top 3 PCA dimensions, useful for confirming major population differences.

    Interactive UMAP Plot: Shows the final 6 cell clusters identified by the Leiden algorithm in 2D space.

    Interactive Marker Box Plots: Displays the scaled expression of key marker genes across the identified clusters.
