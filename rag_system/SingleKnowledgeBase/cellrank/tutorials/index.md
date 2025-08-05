# CellRank Tutorials

Welcome to CellRank tutorials! Learn how to perform fate mapping analysis on single-cell data.

## Getting Started

```{toctree}
:maxdepth: 1

getting_started
```

## Tutorial Categories

### Kernels

Learn about different ways to compute transition matrices:

```{toctree}
:maxdepth: 2

kernels/velocity_kernel
kernels/pseudotime_kernel
kernels/combining_kernels
```

### Fate Mapping

Master the core CellRank workflow:

```{toctree}
:maxdepth: 2

fate_mapping/initial_terminal_states
fate_mapping/fate_probabilities
fate_mapping/driver_genes
```

### Visualization

Create informative plots:

```{toctree}
:maxdepth: 2

visualization/plotting_states
visualization/gene_trends
```

## Conceptual Guides

Understand the theory:

```{toctree}
:maxdepth: 1

../concepts/markov_chains
../concepts/gpcca_explained
```

## Quick Reference

### Import CellRank
```python
import cellrank as cr
import scanpy as sc
import scvelo as scv
```

### Common Workflow
```python
# 1. Create kernel
kernel = cr.kernels.VelocityKernel(adata)
kernel.compute_transition_matrix()

# 2. Initialize estimator
g = cr.estimators.GPCCA(kernel)

# 3. Compute states
g.compute_macrostates()
g.compute_terminal_states()
g.compute_initial_states()

# 4. Compute fate probabilities
g.compute_fate_probabilities()

# 5. Find driver genes
g.compute_lineage_drivers()
```

## Prerequisites

Before using CellRank, you should:
1. Have preprocessed single-cell data (using Scanpy)
2. Have computed dimensionality reduction (PCA, UMAP)
3. (Optional) Have computed RNA velocity (using scVelo)

## Choosing the Right Kernel

| Data Type | Recommended Kernel |
|-----------|-------------------|
| RNA velocity | VelocityKernel |
| Time series | RealTimeKernel |
| Snapshot data | PseudotimeKernel |
| Developmental potential | CytoTRACEKernel |
| Multiple modalities | Combined kernels |