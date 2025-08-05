# Single-Cell Analysis Knowledge Base

A comprehensive collection of tutorials and documentation for single-cell bioinformatics tools.

## 📚 Repository Overview

This knowledge base contains organized tutorials from 6 major single-cell analysis packages:

| Package | Description | Tutorials |
|---------|-------------|-----------|
| **[Scanpy](scanpy/)** | Core single-cell analysis toolkit | 27 files |
| **[CellRank](cellrank/)** | Trajectory inference and fate mapping | 10 files |
| **[Dynamo](dynamo/)** | RNA velocity and dynamics | 24 files |
| **[PertPy](pertpy/)** | Perturbation analysis | 17 files |
| **[Decoupler](decoupler/)** | Enrichment analysis | 15 files |
| **[OmicVerse](omicverse/)** | Comprehensive omics toolkit | 63 files |

**Total: 156 tutorial files**

## 🚀 Quick Start Guides

### For Beginners
1. **[Scanpy - Getting Started](scanpy/tutorials/getting_started.md)** - Basic single-cell analysis workflow
2. **[OmicVerse - Preprocessing](omicverse/omicverse_guide/docs/Tutorials-single/t_preprocess.md)** - Data preprocessing tutorial
3. **[Scanpy - Quality Control](scanpy/tutorials/basics/quality_control.md)** - QC best practices

### By Analysis Type

#### Single-Cell RNA-seq
- **Basic Analysis**: [Scanpy tutorials](scanpy/tutorials/)
- **Advanced Methods**: [OmicVerse single-cell](omicverse/omicverse_guide/docs/Tutorials-single/)
- **Quality Control**: [Scanpy QC guide](scanpy/tutorials/basics/quality_control.md)

#### Trajectory Inference
- **RNA Velocity**: [CellRank velocity kernel](cellrank/tutorials/getting_started.md)
- **Pseudotime**: [Dynamo tutorials](dynamo/tutorials/)
- **Fate Mapping**: [CellRank fate probabilities](cellrank/tutorials/)

#### Perturbation Analysis
- **CRISPR screens**: [PertPy tools](pertpy/tutorials/)
- **Drug response**: [OmicVerse scDrug](omicverse/omicverse_guide/docs/Tutorials-single/t_scdrug.md)

#### Spatial Transcriptomics
- **Spatial Analysis**: [OmicVerse spatial tutorials](omicverse/omicverse_guide/docs/Tutorials-space/)
- **Integration**: [Single to spatial mapping](omicverse/omicverse_guide/docs/Tutorials-bulk2single/t_single2spatial.md)

#### Bulk RNA-seq
- **Differential Expression**: [OmicVerse bulk tutorials](omicverse/omicverse_guide/docs/Tutorials-bulk/)
- **Enrichment Analysis**: [Decoupler methods](decoupler/tutorials/)

## 📖 Learning Paths

### Path 1: Complete Single-Cell Analysis
1. [Data Loading and QC](scanpy/tutorials/getting_started.md)
2. [Preprocessing](omicverse/omicverse_guide/docs/Tutorials-single/t_preprocess.md)
3. [Clustering](omicverse/omicverse_guide/docs/Tutorials-single/t_cluster.md)
4. [Cell Type Annotation](omicverse/omicverse_guide/docs/Tutorials-single/t_cellanno.md)
5. [Differential Expression](omicverse/omicverse_guide/docs/Tutorials-single/t_deg_single.md)

### Path 2: Trajectory Analysis
1. [RNA Velocity Basics](cellrank/tutorials/getting_started.md)
2. [Trajectory Inference](omicverse/omicverse_guide/docs/Tutorials-single/t_traj.md)
3. [Fate Mapping](cellrank/tutorials/fate_mapping/)
4. [Gene Dynamics](dynamo/tutorials/index_conventional.md)

### Path 3: Multi-Omics Integration
1. [MOFA Analysis](omicverse/omicverse_guide/docs/Tutorials-single/t_mofa.md)
2. [Batch Correction](omicverse/omicverse_guide/docs/Tutorials-single/t_single_batch.md)
3. [Multi-Modal Integration](pertpy/tutorials/)

## 🛠️ Tool Selection Guide

### Which tool should I use?

| Task | Recommended Tool | Tutorial |
|------|------------------|----------|
| Basic single-cell analysis | Scanpy | [Getting Started](scanpy/tutorials/getting_started.md) |
| RNA velocity | CellRank + scVelo | [Velocity Kernel](cellrank/tutorials/getting_started.md) |
| Trajectory inference | Dynamo or CellRank | [Dynamo](dynamo/tutorials/) / [CellRank](cellrank/tutorials/) |
| Cell annotation | OmicVerse | [Cell Annotation](omicverse/omicverse_guide/docs/Tutorials-single/t_cellanno.md) |
| Enrichment analysis | Decoupler | [Methods Overview](decoupler/tutorials/methods_overview.md) |
| Perturbation analysis | PertPy | [Tools](pertpy/tutorials/tools.md) |
| Spatial analysis | OmicVerse | [Spatial Tutorials](omicverse/omicverse_guide/docs/Tutorials-space/) |

## 📊 Package Ecosystem

```
Single-Cell Analysis Ecosystem
├── Data Processing
│   ├── Scanpy (core toolkit)
│   └── OmicVerse (comprehensive)
├── Trajectory Analysis
│   ├── CellRank (fate mapping)
│   └── Dynamo (dynamics)
├── Functional Analysis
│   ├── Decoupler (enrichment)
│   └── PertPy (perturbations)
└── Specialized Methods
    └── OmicVerse (multi-omics)
```

## 🔗 External Resources

- [scverse](https://scverse.org/) - Single-cell analysis ecosystem
- [Single-cell best practices](https://www.sc-best-practices.org/) - Comprehensive guide
- [CZI cellxgene](https://cellxgene.cziscience.com/) - Data portal

## 📝 Contributing

This knowledge base is a compilation of documentation from open-source projects. Please refer to each project's repository for contributing guidelines:

- [Scanpy](https://github.com/scverse/scanpy)
- [CellRank](https://github.com/theislab/cellrank)
- [Dynamo](https://github.com/aristoteleo/dynamo-release)
- [PertPy](https://github.com/theislab/pertpy)
- [Decoupler](https://github.com/saezlab/decoupler-py)
- [OmicVerse](https://github.com/Starlitnightly/omicverse)

## 📄 License

Each package maintains its own license. Please refer to individual repositories for license information.