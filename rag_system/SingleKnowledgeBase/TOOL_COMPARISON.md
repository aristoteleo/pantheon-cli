# Single-Cell Analysis Tools Comparison Guide

## Overview

This guide helps you choose the right tool for your single-cell analysis needs.

## Core Analysis Packages

### Scanpy vs OmicVerse

| Feature | Scanpy | OmicVerse |
|---------|--------|-----------|
| **Purpose** | Core single-cell toolkit | Comprehensive omics platform |
| **Strengths** | - Industry standard<br>- Extensive documentation<br>- Large community | - All-in-one solution<br>- Advanced methods<br>- GPU acceleration |
| **Best for** | Standard scRNA-seq analysis | Complex multi-omics projects |
| **Learning curve** | Moderate | Steeper |
| **GPU support** | Limited | Extensive |

### When to use each:
- **Scanpy**: Standard single-cell RNA-seq analysis, reproducible workflows
- **OmicVerse**: Advanced analyses, GPU acceleration needed, multi-omics integration

## Trajectory Analysis

### CellRank vs Dynamo

| Feature | CellRank | Dynamo |
|---------|----------|--------|
| **Focus** | Fate mapping | RNA dynamics |
| **Methods** | Markov chains, GPCCA | Vector fields, metabolic labeling |
| **Inputs** | RNA velocity, pseudotime | RNA velocity, labeling data |
| **Outputs** | Fate probabilities, terminal states | Continuous dynamics, regulation |
| **Visualization** | Excellent | Good |

### When to use each:
- **CellRank**: Identifying cell fates, terminal states, lineage probabilities
- **Dynamo**: Understanding continuous dynamics, gene regulation, metabolic labeling data

## Functional Analysis

### Decoupler vs PertPy

| Feature | Decoupler | PertPy |
|---------|-----------|--------|
| **Purpose** | Enrichment analysis | Perturbation analysis |
| **Methods** | Multiple statistical methods | CRISPR/drug response tools |
| **Strengths** | Method consensus, benchmarking | Perturbation-specific tools |
| **Use cases** | Pathway analysis, TF activity | CRISPR screens, drug studies |

### When to use each:
- **Decoupler**: General enrichment analysis, pathway activity inference
- **PertPy**: Analyzing perturbation experiments, CRISPR screens

## Analysis Workflows

### Typical Pipeline Combinations

#### Standard scRNA-seq Analysis
```
Scanpy → Decoupler → CellRank
- Basic processing → Enrichment → Trajectory
```

#### Advanced Multi-Omics
```
OmicVerse → Dynamo → Decoupler
- Integrated analysis → Dynamics → Function
```

#### Perturbation Studies
```
Scanpy → PertPy → Decoupler
- Preprocessing → Perturbation analysis → Enrichment
```

## Feature Matrix

| Feature | Scanpy | CellRank | Dynamo | PertPy | Decoupler | OmicVerse |
|---------|--------|----------|--------|--------|-----------|-----------|
| **Preprocessing** | ✅ | ❌ | ✅ | ✅ | ❌ | ✅ |
| **Clustering** | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ |
| **Trajectory** | ⚠️ | ✅ | ✅ | ❌ | ❌ | ✅ |
| **Velocity** | ❌ | ✅ | ✅ | ❌ | ❌ | ✅ |
| **Enrichment** | ⚠️ | ❌ | ❌ | ❌ | ✅ | ✅ |
| **Perturbations** | ❌ | ❌ | ❌ | ✅ | ❌ | ⚠️ |
| **Spatial** | ⚠️ | ❌ | ❌ | ❌ | ❌ | ✅ |
| **Multi-omics** | ❌ | ❌ | ⚠️ | ❌ | ❌ | ✅ |
| **GPU support** | ⚠️ | ❌ | ⚠️ | ❌ | ❌ | ✅ |

Legend: ✅ Full support | ⚠️ Limited support | ❌ Not supported

## Decision Tree

```
What type of analysis?
├── Standard scRNA-seq → Scanpy
├── Need GPU acceleration → OmicVerse
├── Trajectory/Fate mapping
│   ├── Focus on endpoints → CellRank
│   └── Focus on dynamics → Dynamo
├── Functional analysis
│   ├── Enrichment/Pathways → Decoupler
│   └── Perturbations → PertPy
└── Multi-omics/Spatial → OmicVerse
```

## Integration Capabilities

### Data Flow Between Tools

```
[Raw Data]
    ↓
[Scanpy/OmicVerse] ← Preprocessing
    ↓
[AnnData Object] ← Shared data format
    ↓
┌─────────────┬──────────────┬──────────────┐
│  CellRank   │   Dynamo     │  Decoupler   │
│ (Fate map)  │ (Dynamics)   │ (Enrichment) │
└─────────────┴──────────────┴──────────────┘
```

## Performance Considerations

| Tool | Memory Usage | Speed | Scalability |
|------|--------------|-------|-------------|
| Scanpy | Medium | Fast | Good |
| CellRank | High | Medium | Medium |
| Dynamo | High | Slow | Limited |
| PertPy | Medium | Fast | Good |
| Decoupler | Low | Very Fast | Excellent |
| OmicVerse | Variable | Fast (GPU) | Excellent |

## Recommendations by Dataset Size

- **< 10k cells**: Any tool works well
- **10k-100k cells**: Scanpy, Decoupler, OmicVerse
- **100k-1M cells**: OmicVerse (GPU), Scanpy (with sampling)
- **> 1M cells**: OmicVerse with GPU acceleration

## Community and Support

| Tool | Documentation | Community | Updates |
|------|---------------|-----------|---------|
| Scanpy | Excellent | Very Large | Regular |
| CellRank | Good | Growing | Regular |
| Dynamo | Good | Medium | Moderate |
| PertPy | Good | Growing | Regular |
| Decoupler | Good | Medium | Regular |
| OmicVerse | Extensive | Growing | Frequent |