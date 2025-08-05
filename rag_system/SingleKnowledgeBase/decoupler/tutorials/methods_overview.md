# Decoupler Methods Overview

## Introduction

Decoupler implements an ensemble of computational methods to infer biological activities from omics data. Each method has its own strengths and is suitable for different types of analyses.

## Available Methods

The methods module (`dc.mt`) contains various enrichment analysis approaches:

### Statistical Methods
- **ULM (Univariate Linear Model)**: Fast linear modeling approach
- **MLM (Multivariate Linear Model)**: Accounts for gene correlations
- **VIPER**: Virtual Inference of Protein activity by Enriched Regulon
- **GSEA (Gene Set Enrichment Analysis)**: Classic enrichment method
- **AUCell**: Area Under the Curve cell-based scoring
- **ORA (Over-Representation Analysis)**: Hypergeometric test-based
- **GSVA (Gene Set Variation Analysis)**: Sample-wise enrichment
- **Consensus**: Combines multiple method results

### Usage Pattern

All methods follow a similar API pattern:

```python
import decoupler as dc

# Run a specific method
results = dc.run_[method_name](
    mat,           # Expression matrix
    net,           # Network/gene sets
    source='source',
    target='target',
    weight='weight',
    verbose=True
)
```

### Choosing the Right Method

- **For speed**: Use ULM or ORA
- **For accuracy**: Use MLM or consensus approaches
- **For single-cell data**: Consider AUCell or ULM
- **For bulk RNA-seq**: GSEA, GSVA, or VIPER work well
- **For robustness**: Use consensus of multiple methods

## Data Resources

Decoupler integrates with Omnipath (`dc.op`) to provide access to:
- Transcription factor targets
- Pathway databases
- Protein-protein interactions
- Metabolic networks

## Visualization

Use the plotting module (`dc.pl`) to visualize:
- Enrichment scores
- Activity heatmaps
- Dotplots
- Network visualizations

## Further Reading

- See the [API documentation](../api/methods/mt.md) for detailed method parameters
- Check [benchmarking results](../api/methods/bm.md) for method comparisons
- Explore [preprocessing options](../api/analysis/pp.md) for data preparation