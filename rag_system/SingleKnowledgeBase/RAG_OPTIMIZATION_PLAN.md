# RAG System Optimization Plan for Bioinformatics Code Generation

## Current Challenges for RAG

1. **Mixed Content**: Tutorials contain explanatory text mixed with code
2. **Inconsistent Code Formatting**: Different styles across packages
3. **Missing Context**: Code snippets lack metadata (dependencies, versions, prerequisites)
4. **Poor Chunk Boundaries**: Long documents make retrieval less precise
5. **No Semantic Tagging**: Hard to find specific functionality

## Proposed RAG-Optimized Structure

### 1. Separate Code from Documentation

```
SingleKnowledgeBase_RAG/
├── code_snippets/              # Pure code examples
│   ├── by_task/               # Task-oriented organization
│   ├── by_package/            # Package-specific code
│   └── workflows/             # Complete workflows
├── context/                   # Explanatory content
├── metadata/                  # Rich metadata for retrieval
└── index/                     # Pre-built indices
```

### 2. Code Snippet Structure

Each code snippet should be a separate file with metadata header:

```markdown
---
task: quality_control
packages: [scanpy]
prerequisites: [data_loaded]
outputs: [filtered_adata]
keywords: [qc, filter, mitochondrial]
complexity: beginner
---

# Calculate QC metrics
import scanpy as sc

# Calculate mitochondrial gene percentage
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells
adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.n_genes_by_counts > 200, :]
```

### 3. Task-Based Organization

```
code_snippets/
├── by_task/
│   ├── data_loading/
│   │   ├── load_10x_data.md
│   │   ├── load_h5ad.md
│   │   └── load_csv.md
│   ├── quality_control/
│   │   ├── calculate_qc_metrics.md
│   │   ├── filter_cells.md
│   │   └── identify_doublets.md
│   ├── normalization/
│   │   ├── normalize_total.md
│   │   ├── log_transform.md
│   │   └── scale_data.md
│   ├── clustering/
│   │   ├── leiden_clustering.md
│   │   ├── louvain_clustering.md
│   │   └── hierarchical_clustering.md
│   └── trajectory_inference/
│       ├── velocity_kernel.md
│       ├── pseudotime_kernel.md
│       └── fate_probabilities.md
```

### 4. Workflow Templates

Complete analysis workflows with decision points:

```
workflows/
├── standard_scrnaseq/
│   ├── 01_load_and_qc.md
│   ├── 02_normalization.md
│   ├── 03_feature_selection.md
│   ├── 04_dimensionality_reduction.md
│   ├── 05_clustering.md
│   └── 06_differential_expression.md
├── trajectory_analysis/
├── spatial_analysis/
└── perturbation_analysis/
```

### 5. Metadata System

#### Global Metadata Index
```json
{
  "code_snippets": {
    "calculate_qc_metrics": {
      "file": "code_snippets/by_task/quality_control/calculate_qc_metrics.md",
      "task": "quality_control",
      "packages": ["scanpy"],
      "prerequisites": ["adata_loaded"],
      "outputs": ["qc_metrics_calculated"],
      "related_snippets": ["filter_cells", "plot_qc_metrics"],
      "common_errors": ["missing_gene_names", "empty_adata"],
      "performance": "O(n_cells * n_genes)"
    }
  }
}
```

### 6. Enhanced Code Structure

#### Before (Current):
```markdown
# Quality Control Tutorial
Quality control is important...
[100 lines of explanation]
```python
sc.pp.calculate_qc_metrics(adata)
```
[More explanation]
```

#### After (RAG-Optimized):
```markdown
---
id: qc_001
title: Calculate QC Metrics
category: preprocessing/quality_control
---

## Code
```python
def calculate_qc_metrics(adata, species='human'):
    """
    Calculate quality control metrics for single-cell data.
    
    Parameters:
    - adata: AnnData object
    - species: 'human' or 'mouse' (affects mitochondrial gene detection)
    
    Returns:
    - adata with QC metrics in .obs
    """
    import scanpy as sc
    
    # Identify mitochondrial genes
    mt_prefix = 'MT-' if species == 'human' else 'mt-'
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    
    # Calculate metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    
    return adata
```

## Context
- Use after: data_loading
- Use before: filtering
- Common threshold: 5% mitochondrial

## Variations
- Mouse data: use 'mt-' prefix
- Include ribosomal: add 'ribo' calculation
```

### 7. Semantic Search Tags

Create semantic tags for better retrieval:

```yaml
semantic_tags:
  quality_control:
    - filter_low_quality_cells
    - remove_dying_cells
    - mitochondrial_percentage
    - cell_viability
    - technical_artifacts
    
  normalization:
    - count_depth_correction
    - library_size_normalization
    - log_transformation
    - scaling_centering
```

### 8. Error Handling Patterns

```
error_patterns/
├── common_errors/
│   ├── empty_adata.md
│   ├── missing_raw_counts.md
│   └── incompatible_versions.md
└── debugging_guides/
    ├── memory_issues.md
    └── performance_optimization.md
```

### 9. Package Compatibility Matrix

```markdown
# Compatibility Matrix

| Operation | Scanpy | CellRank | Dynamo | OmicVerse |
|-----------|--------|----------|---------|-----------|
| Load 10X  | ✓ 2.0+ | ✓ 1.5+   | ✓ 1.0+  | ✓ 1.2+    |
| QC Filter | ✓ 1.8+ | -        | ✓ 1.1+  | ✓ 1.0+    |
```

### 10. Implementation Steps

1. **Extract Code Snippets**
   - Parse existing tutorials
   - Extract code blocks
   - Add metadata headers

2. **Create Task Taxonomy**
   - Define standard task names
   - Map snippets to tasks
   - Build dependency graph

3. **Build Search Index**
   - Create embeddings for code
   - Index by task/package/keyword
   - Enable semantic search

4. **Version Control**
   - Track package versions
   - Note breaking changes
   - Maintain compatibility info

## Benefits for RAG

1. **Precise Retrieval**: Smaller, focused chunks
2. **Better Context**: Metadata provides usage info
3. **Code Reusability**: Standardized snippets
4. **Error Prevention**: Common issues documented
5. **Workflow Guidance**: Complete examples available

## Example Query Handling

**User Query**: "How to filter low quality cells in scanpy?"

**RAG Retrieves**:
1. Exact code snippet for QC filtering
2. Prerequisites (data must be loaded)
3. Common thresholds
4. Related operations (plotting QC metrics)
5. Potential errors and solutions