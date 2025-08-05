# Getting Started with CellRank

## Overview

CellRank is a modular framework to study cellular dynamics based on Markov state modeling of multi-view single-cell data. It computes initial and terminal states, fate probabilities, and driver genes using various input data modalities including RNA velocity, pseudotime, developmental potential, and experimental time points.

## Installation

```bash
# Basic installation
pip install cellrank

# Development version
pip install git+https://github.com/theislab/cellrank@main
```

## Quick Start

```python
import cellrank as cr
import scanpy as sc
import scvelo as scv

# CellRank builds on top of scanpy and scvelo
import numpy as np
import pandas as pd
```

## Basic Workflow

### 1. Prepare Data

CellRank expects an AnnData object with computed embeddings:

```python
# Load preprocessed data
adata = sc.read_h5ad('preprocessed_data.h5ad')

# If using RNA velocity, compute it first
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
```

### 2. Compute Transition Matrix

Choose a kernel based on your data:

```python
# RNA velocity kernel
from cellrank.kernels import VelocityKernel
vk = VelocityKernel(adata)
vk.compute_transition_matrix()

# Pseudotime kernel
from cellrank.kernels import PseudotimeKernel
pk = PseudotimeKernel(adata, time_key='dpt_pseudotime')
pk.compute_transition_matrix()

# Combine kernels
from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata)
ck.compute_transition_matrix()

# Combined kernel (80% velocity + 20% connectivity)
combined_kernel = 0.8 * vk + 0.2 * ck
```

### 3. Identify Initial and Terminal States

```python
# Initialize estimator
from cellrank.estimators import GPCCA
g = GPCCA(combined_kernel)

# Compute macrostates
g.compute_macrostates(n_states=5)
g.plot_macrostates(which='all')

# Identify terminal states
g.compute_terminal_states()
g.plot_terminal_states()

# Identify initial states
g.compute_initial_states()
g.plot_initial_states()
```

### 4. Compute Fate Probabilities

```python
# Compute absorption probabilities (fate probabilities)
g.compute_fate_probabilities()

# Visualize fate probabilities
g.plot_fate_probabilities()
```

### 5. Identify Driver Genes

```python
# Compute lineage drivers
g.compute_lineage_drivers(lineages=['Lineage_1', 'Lineage_2'])

# Plot top driver genes
g.plot_lineage_drivers()
```

### 6. Analyze Gene Expression Trends

```python
# Define model
from cellrank.models import GAM
model = GAM(adata)

# Fit model for specific genes
g.plot_gene_trends(
    ['Gene1', 'Gene2', 'Gene3'],
    data_key='X',
    ncols=3
)
```

## Key Concepts

### Kernels
- **VelocityKernel**: Uses RNA velocity
- **PseudotimeKernel**: Uses pseudotime ordering
- **CytoTRACEKernel**: Uses developmental potential
- **RealTimeKernel**: Uses experimental time points

### Estimators
- **GPCCA**: Generalized Perron Cluster Cluster Analysis
- **CFLARE**: Clustering and Filtering of Left and Right Eigenvectors

### Models
- **GAM**: Generalized Additive Models
- **SKLearnModel**: Scikit-learn compatible models

## Example: Complete Analysis

```python
# Load data
adata = sc.read_h5ad('my_data.h5ad')

# Compute velocity if needed
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Create kernel
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

# Run analysis
g = cr.estimators.GPCCA(vk)
g.compute_macrostates(n_states=10)
g.compute_terminal_states()
g.compute_fate_probabilities()

# Save results
g.plot_fate_probabilities(save='fate_probabilities.pdf')
adata.write('cellrank_results.h5ad')
```

## Next Steps

- Learn about different [kernels](kernels/velocity_kernel.md)
- Understand [fate mapping](fate_mapping/initial_terminal_states.md)
- Explore [visualization options](visualization/plotting_states.md)
- Read about [theoretical concepts](../concepts/markov_chains.md)