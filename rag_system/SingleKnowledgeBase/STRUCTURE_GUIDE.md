# SingleKnowledgeBase Structure Guide

## 📁 Overview

The SingleKnowledgeBase is organized with a clear hierarchy that makes it easy to find what you need. Here's a detailed breakdown:

```
SingleKnowledgeBase/
├── README.md                    # Main entry point - overview of all tools
├── TOOL_COMPARISON.md          # Comprehensive comparison of all tools
├── STRUCTURE_GUIDE.md          # This file - explains the organization
│
├── scanpy/                     # Core single-cell analysis toolkit
├── cellrank/                   # Trajectory inference and fate mapping
├── dynamo/                     # RNA velocity and dynamics
├── pertpy/                     # Perturbation analysis
├── decoupler/                  # Enrichment analysis
└── omicverse/                  # Comprehensive omics toolkit
```

## 📚 Detailed Structure by Tool

### 1. Scanpy Directory (`scanpy/`)
The foundation for single-cell analysis:

```
scanpy/
├── README.md                   # Package overview
├── api/                        # API reference documentation
│   ├── index.md               # API overview
│   ├── preprocessing/         # Data preprocessing functions
│   │   └── pp.md             # sc.pp.* functions
│   ├── tools/                 # Analysis tools
│   │   └── tl.md             # sc.tl.* functions
│   ├── plotting/              # Visualization
│   │   └── pl.md             # sc.pl.* functions
│   ├── io/                    # Input/output
│   │   └── io.md             # Reading/writing data
│   └── datasets/              # Built-in datasets
│       └── datasets.md        # Available datasets
├── tutorials/                  # Step-by-step guides
│   ├── index.md               # Tutorial overview
│   ├── getting_started.md     # Quick start guide
│   ├── basics/                # Fundamental concepts
│   │   ├── quality_control.md # QC best practices
│   │   ├── normalization.md   # Data normalization
│   │   ├── clustering.md      # Cell clustering
│   │   └── integration.md     # Batch integration
│   ├── visualization/         # Plotting tutorials
│   └── advanced/              # Advanced topics
└── how-to/                    # Quick task guides
```

**Key files for beginners:**
- `tutorials/getting_started.md` - Start here!
- `tutorials/basics/quality_control.md` - Essential QC steps

### 2. CellRank Directory (`cellrank/`)
For trajectory inference and fate mapping:

```
cellrank/
├── README.md                   # Package overview
├── installation.md            # Installation guide
├── api/                       # API documentation
│   ├── kernels/              # Transition matrix computation
│   │   └── overview.md       # Available kernels
│   ├── estimators/           # State identification
│   │   └── overview.md       # GPCCA, CFLARE
│   ├── models/               # Gene expression models
│   └── plotting/             # Visualization functions
├── tutorials/                 # Learning materials
│   ├── getting_started.md    # Quick introduction
│   ├── kernels/              # Kernel tutorials
│   ├── fate_mapping/         # Core workflow
│   └── visualization/        # Plotting guides
└── concepts/                 # Theoretical background
    ├── markov_chains.md      # Mathematical foundation
    └── gpcca_explained.md    # Algorithm details
```

**Key concepts:**
- Kernels compute transition matrices from different data types
- Estimators identify initial/terminal states
- Fate probabilities show cell differentiation paths

### 3. Dynamo Directory (`dynamo/`)
For studying RNA dynamics:

```
dynamo/
├── api/                       # API reference
├── introduction/              # Conceptual introductions
│   ├── index_cellfate.md     # Cell fate concepts
│   ├── index_geo.md          # Geometric analysis
│   ├── index_silico.md       # In silico experiments
│   └── index_time.md         # Time-resolved analysis
├── tutorials/                 # Detailed tutorials
│   ├── index_preprocessing.md # Data preparation
│   ├── index_conventional.md  # Standard velocity
│   ├── index_labeling.md     # Metabolic labeling
│   ├── index_differential_geometry.md # Advanced math
│   └── index_multivelo.md   # Multi-modal velocity
└── user_guide/               # Comprehensive guides
```

**Unique features:**
- Metabolic labeling data analysis
- Vector field reconstruction
- In silico perturbations

### 4. PertPy Directory (`pertpy/`)
For perturbation experiments:

```
pertpy/
├── about/                     # Package information
│   ├── background.md         # Scientific background
│   └── cite.md              # Citation info
├── api/                      # API documentation
│   ├── tools_index.md       # Analysis tools
│   ├── preprocessing_index.md # Data prep
│   └── metadata_index.md    # Metadata handling
├── tutorials/                # Usage guides
│   ├── tools.md             # Using analysis tools
│   ├── preprocessing.md     # Preparing perturbation data
│   └── metadata.md          # Managing metadata
└── usecases.md              # Real-world examples
```

**Best for:**
- CRISPR screen analysis
- Drug response studies
- Perturbation experiments

### 5. Decoupler Directory (`decoupler/`)
For enrichment analysis:

```
decoupler/
├── api/                      # API organization
│   ├── methods/             # Statistical methods
│   │   ├── mt.md           # All enrichment methods
│   │   └── bm.md           # Benchmarking tools
│   ├── data/               # Data resources
│   │   ├── ds.md           # Data structures
│   │   └── op.md           # Omnipath integration
│   └── analysis/           # Analysis tools
├── tutorials/              # Learning resources
│   ├── getting_started.md  # Quick start
│   └── methods_overview.md # Method comparison
└── references.md           # Citations
```

**Key features:**
- Multiple enrichment methods
- Consensus approaches
- Omnipath database integration

### 6. OmicVerse Directory (`omicverse/`)
The most comprehensive toolkit:

```
omicverse/
├── omicverse_guide/
│   └── docs/
│       ├── Tutorials-single/      # Single-cell (33 tutorials!)
│       │   ├── t_preprocess.md   # Preprocessing
│       │   ├── t_cluster.md      # Clustering
│       │   ├── t_cellanno.md     # Cell annotation
│       │   ├── t_traj.md         # Trajectories
│       │   └── ...               # Many more!
│       ├── Tutorials-bulk/        # Bulk RNA-seq (6 tutorials)
│       ├── Tutorials-space/       # Spatial (11 tutorials)
│       ├── Tutorials-plotting/    # Visualization (3 tutorials)
│       └── Tutorials-bulk2single/ # Integration (3 tutorials)
└── sample/                        # Example code
```

**Extensive coverage:**
- 63 total tutorials
- GPU acceleration
- Multi-omics integration
- Spatial transcriptomics

## 🗺️ Navigation Tips

### Finding What You Need

1. **New to single-cell?** 
   - Start: `scanpy/tutorials/getting_started.md`
   - Then: `scanpy/tutorials/basics/quality_control.md`

2. **Need specific analysis?**
   - Trajectories: Check `cellrank/` or `dynamo/`
   - Enrichment: Go to `decoupler/`
   - Perturbations: Use `pertpy/`
   - Advanced/GPU: Try `omicverse/`

3. **Looking for examples?**
   - Code examples: `omicverse/sample/`
   - Workflows: `*/tutorials/getting_started.md`
   - API usage: `*/api/` directories

### File Naming Conventions

- `index.md` - Overview/table of contents
- `getting_started.md` - Quick start guides
- `t_*.md` - Tutorial files (in OmicVerse)
- `*_index.md` - Category overviews
- API files named by module (e.g., `pp.md` for preprocessing)

## 🔍 Quick Reference Paths

Most commonly needed files:

```bash
# Getting started
SingleKnowledgeBase/scanpy/tutorials/getting_started.md
SingleKnowledgeBase/README.md

# Quality control
SingleKnowledgeBase/scanpy/tutorials/basics/quality_control.md

# Trajectories
SingleKnowledgeBase/cellrank/tutorials/getting_started.md
SingleKnowledgeBase/omicverse/omicverse_guide/docs/Tutorials-single/t_traj.md

# Spatial analysis
SingleKnowledgeBase/omicverse/omicverse_guide/docs/Tutorials-space/

# Tool comparison
SingleKnowledgeBase/TOOL_COMPARISON.md
```

## 📊 Statistics

- **Total files**: 158 markdown files
- **Largest section**: OmicVerse (63 files)
- **Most tutorials**: OmicVerse single-cell (33 tutorials)
- **Best documented**: Scanpy and OmicVerse
- **Most specialized**: PertPy (perturbations only)