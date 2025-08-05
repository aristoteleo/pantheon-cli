# Getting Started with Scanpy

## Overview

Scanpy is a scalable toolkit for analyzing single-cell gene expression data. It includes preprocessing, visualization, clustering, trajectory inference and differential expression testing. The Python-based implementation efficiently deals with datasets of more than one million cells.

## Installation

```bash
# Basic installation
pip install scanpy

# With optional dependencies
pip install scanpy[leiden]  # For leiden clustering
pip install scanpy[rapids]  # For GPU acceleration
```

## Quick Start

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set up scanpy
sc.settings.verbosity = 3  # Verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')
```

## Basic Workflow

### 1. Load Data

```python
# Read 10X data
adata = sc.read_10x_h5('path/to/data.h5')

# Read H5AD file
adata = sc.read_h5ad('path/to/data.h5ad')

# Read CSV/TSV
adata = sc.read_csv('path/to/counts.csv')
```

### 2. Basic Information

```python
# Show data structure
print(adata)

# Basic info
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
```

### 3. Quality Control

```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Add mitochondrial gene info
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Basic QC plots
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```

### 4. Filtering

```python
# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter based on QC metrics
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

### 5. Normalization and Scaling

```python
# Save raw counts
adata.raw = adata

# Normalize every cell to 10,000 UMI
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Scale data
sc.pp.scale(adata, max_value=10)
```

### 6. Dimensionality Reduction

```python
# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP
sc.tl.umap(adata)

# t-SNE (optional)
sc.tl.tsne(adata)
```

### 7. Clustering

```python
# Leiden clustering
sc.tl.leiden(adata)

# Louvain clustering (alternative)
sc.tl.louvain(adata)

# Visualize
sc.pl.umap(adata, color=['leiden'])
```

### 8. Finding Marker Genes

```python
# Find markers for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

### 9. Save Results

```python
# Save processed data
adata.write('processed_data.h5ad')

# Save figures
sc.pl.umap(adata, color=['leiden'], save='_leiden.pdf')
```

## Next Steps

- Explore [preprocessing options](basics/quality_control.md)
- Learn about [advanced clustering](basics/clustering.md)
- Understand [trajectory inference](advanced/trajectories.md)
- Master [visualization techniques](visualization/basic_plots.md)

## Useful Resources

- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [Tutorials repository](https://github.com/scverse/scanpy-tutorials)
- [Single cell best practices](https://www.sc-best-practices.org/)