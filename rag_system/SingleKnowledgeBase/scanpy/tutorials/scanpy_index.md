# Scanpy Tutorials

Welcome to the Scanpy tutorials! These guides will help you master single-cell analysis.

## Getting Started

```{toctree}
:maxdepth: 1

getting_started
```

## Tutorial Categories

### Basic Workflows

Learn fundamental single-cell analysis steps:

```{toctree}
:maxdepth: 2

basics/data_loading
basics/quality_control
basics/normalization
basics/clustering
basics/integration
```

### Visualization

Master plotting techniques:

```{toctree}
:maxdepth: 2

visualization/basic_plots
visualization/advanced_plots
```

### Advanced Topics

Explore sophisticated analyses:

```{toctree}
:maxdepth: 2

advanced/trajectories
advanced/experimental
```

## Learning Path

1. **New to single-cell?** Start with [Getting Started](getting_started.md)
2. **Ready for analysis?** Follow the Basic Workflows in order
3. **Need specific plots?** Check Visualization guides
4. **Advanced user?** Explore trajectories and experimental features

## Quick Reference

### Import Scanpy
```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
```

### Common Operations
- **Read data**: `sc.read_10x_h5()`, `sc.read_h5ad()`
- **QC metrics**: `sc.pp.calculate_qc_metrics()`
- **Filter**: `sc.pp.filter_cells()`, `sc.pp.filter_genes()`
- **Normalize**: `sc.pp.normalize_total()`, `sc.pp.log1p()`
- **Find variable genes**: `sc.pp.highly_variable_genes()`
- **PCA**: `sc.tl.pca()`
- **Neighbors**: `sc.pp.neighbors()`
- **UMAP**: `sc.tl.umap()`
- **Cluster**: `sc.tl.leiden()`
- **Find markers**: `sc.tl.rank_genes_groups()`