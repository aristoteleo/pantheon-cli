# Quality Control in Single-Cell Analysis

Quality control (QC) is a crucial first step in single-cell RNA-seq analysis. This tutorial covers how to identify and filter low-quality cells and genes.

## Overview

Poor quality cells can arise from:
- Cell damage during dissociation
- Failure in library preparation
- Low RNA capture efficiency
- Ambient RNA contamination

## QC Metrics

### 1. Calculate Basic Metrics

```python
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load data
adata = sc.read_10x_h5('data.h5')

# Basic metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# This adds to adata.obs:
# - n_genes_by_counts: number of genes expressed
# - total_counts: total UMI counts
# - pct_counts_in_top_20_genes: percentage of counts in top 20 genes
```

### 2. Mitochondrial Gene Percentage

High mitochondrial gene expression often indicates dying cells:

```python
# Identify mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Human
# adata.var['mt'] = adata.var_names.str.startswith('mt-')  # Mouse

# Calculate mitochondrial QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# This adds pct_counts_mt to adata.obs
```

### 3. Ribosomal Gene Percentage

```python
# Identify ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

# Calculate ribosomal QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
```

## Visualization

### 1. Violin Plots

```python
# Key QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Additional metrics
sc.pl.violin(adata, ['pct_counts_ribo', 'pct_counts_in_top_20_genes'],
             jitter=0.4, multi_panel=True)
```

### 2. Scatter Plots

```python
# Relationship between metrics
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

### 3. Histograms

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].hist(adata.obs['n_genes_by_counts'], bins=50)
axes[0].set_xlabel('Number of genes')
axes[0].set_ylabel('Number of cells')

axes[1].hist(adata.obs['total_counts'], bins=50)
axes[1].set_xlabel('Total counts')

axes[2].hist(adata.obs['pct_counts_mt'], bins=50)
axes[2].set_xlabel('Mitochondrial gene percentage')

plt.tight_layout()
plt.show()
```

## Filtering Strategies

### 1. Basic Filtering

```python
# Before filtering
print(f"Before filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# Filter cells with too few or too many genes
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]

# Filter based on mitochondrial percentage
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Filter genes present in too few cells
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
```

### 2. MAD-based Filtering

More robust filtering using Median Absolute Deviation:

```python
def mad_based_outlier(data, threshold=3):
    """Identify outliers using MAD."""
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    modified_z_scores = 0.6745 * (data - median) / mad
    return np.abs(modified_z_scores) > threshold

# Apply MAD filtering
adata.obs['outlier_n_genes'] = mad_based_outlier(adata.obs['n_genes_by_counts'])
adata.obs['outlier_total_counts'] = mad_based_outlier(adata.obs['total_counts'])
adata.obs['outlier_mt'] = mad_based_outlier(adata.obs['pct_counts_mt'])

# Remove outliers
adata = adata[~(adata.obs['outlier_n_genes'] | 
                adata.obs['outlier_total_counts'] | 
                adata.obs['outlier_mt']), :]
```

### 3. Sample-specific Thresholds

For datasets with multiple samples:

```python
# If you have batch information
for batch in adata.obs['batch'].unique():
    batch_mask = adata.obs['batch'] == batch
    
    # Calculate thresholds per batch
    mt_threshold = adata.obs.loc[batch_mask, 'pct_counts_mt'].quantile(0.95)
    
    # Apply batch-specific filtering
    adata = adata[~(batch_mask & (adata.obs['pct_counts_mt'] > mt_threshold)), :]
```

## Doublet Detection

```python
# Using Scrublet
import scrublet as scr

# Run scrublet
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Add to adata
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublets'] = predicted_doublets

# Remove predicted doublets
adata = adata[~adata.obs['predicted_doublets'], :]
```

## QC Report

Create a comprehensive QC report:

```python
def generate_qc_report(adata):
    """Generate QC statistics."""
    qc_stats = pd.DataFrame(index=['Before_QC', 'After_QC'])
    
    # Add your statistics
    qc_stats.loc['Before_QC', 'n_cells'] = initial_n_cells
    qc_stats.loc['After_QC', 'n_cells'] = adata.n_obs
    
    # Calculate percentages kept
    qc_stats['percent_kept'] = (qc_stats.loc['After_QC'] / 
                                 qc_stats.loc['Before_QC'] * 100)
    
    return qc_stats

# Generate report
qc_report = generate_qc_report(adata)
print(qc_report)
```

## Best Practices

1. **Visualize before filtering**: Always plot QC metrics first
2. **Document thresholds**: Record filtering criteria used
3. **Iterative approach**: Start conservative, refine based on downstream results
4. **Sample-specific**: Consider batch effects in QC
5. **Save raw data**: Keep unfiltered data for reference

## Common Issues

- **Too stringent filtering**: Losing rare cell types
- **Too lenient filtering**: Including low-quality cells
- **Batch effects**: Different samples may need different thresholds
- **Cell type bias**: Some cell types naturally have different QC metrics

## Next Steps

After QC, proceed to:
- [Normalization](normalization.md)
- [Finding highly variable genes](normalization.md#highly-variable-genes)
- [Dimensionality reduction](../getting_started.md#6-dimensionality-reduction)