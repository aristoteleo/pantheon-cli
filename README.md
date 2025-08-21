# Pantheon CLI

<!-- Banner -->
<p align="center">
  <img src="assets/pantheon_banner_tri_dark.svg#gh-dark-mode-only" alt="Pantheon ASCII (dark)" />
  <img src="assets/pantheon_banner_tri_light.svg#gh-light-mode-only" alt="Pantheon ASCII (light)" />
</p>

<div align="center">

***We're not just building another CLI tool.  
We're defining how scientists interact with data in the AI era.***

**The first fully open-source, infinitely extensible scientific "vibe analysis" framework**

</div>

<div align="center">

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Status: Beta](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![AI-Native](https://img.shields.io/badge/AI-Native-purple.svg)]()

</div>


## `1` [What is Pantheon-CLI?](#1-what-is-pantheon-cli)

Pantheon-CLI is the **first fully open-source "vibe analysis" framework** built specifically for scientific research. We are defining a new way for scientists to interact with data in the AI era.

### **PhD-Level Scientific Assistant**
- Pantheon-CLI is the first command-line intelligent agent assistant for complex real-world scientific analysis, capable of handling PhD-level single-cell and spatial genomics tasks. This is not just a tool—**it's an AI scientist on your research team**.

### **Mixed Programming** 
- With Pantheon-CLI, you can work within the same environment to:
  - Write Python code on the first line
  - Use natural language descriptions on the next line  
  - Even mix in R/Julia languages

All scientists need only focus on the analysis itself, without switching between different tools and environments.

## `2` [Quick Start](#2-quick-start)

### Experience Mixed Programming
```bash
# Start Pantheon CLI
pantheon-cli

# Perfect integration of natural language and code:
> Generate 10 random numbers and calculate their mean and standard deviation
# System automatically generates and executes Python code

> Load single-cell RNA-seq data, the file is data.h5ad
# Automatically uses scanpy to load data

> Now use Seurat for cell clustering analysis
# Seamlessly switches to R language environment

> Use Julia to solve optimization problem, then visualize results in Python
```

### Installation

#### Simple Installation (Recommended)
```bash
pip install pantheon-cli
```

#### Development Installation
```bash
# Install from source (recommended for development)
git clone https://github.com/aristoteleo/pantheon-cli.git
cd Pantheon-cli
pip install -e .

# Make sure dependencies are installed
pip install pantheon-agents pantheon-toolsets
```

#### Verify Installation
```bash
pantheon-cli --version
```

**Note**: Pantheon-CLI requires both `pantheon-agents` and `pantheon-toolsets` to be installed. These provide the core agent functionality and distributed toolsets respectively.

### Basic Usage

#### First Launch
```bash
# Start Pantheon-CLI
pantheon-cli
```

The system will prompt you to configure an API key or select a local model. For quick experience, you can configure an OpenAI or Anthropic API key.

#### API Key Configuration
```bash
# Once the CLI is running, setup your API keys:
/api-key list  # List current API keys

# Set API keys (examples):
/api-key openai sk-your-key-here
/api-key anthropic sk-your-key-here
/api-key deepseek sk-your-key-here
```

#### Launch Options
```bash
# Start with default settings
pantheon-cli

# Start with different model
pantheon-cli --model claude-sonnet-4-20250514

# Start without RAG database
pantheon-cli --disable_rag

# Start with custom workspace
pantheon-cli --workspace /path/to/project

# Start with external toolsets
pantheon-cli --disable_ext False --ext_dir ./ext_toolsets
```

### With RAG Database

If you have a RAG database prepared:

```bash
pantheon-cli --rag_db path/to/rag/database
```

Default RAG database location: `tmp/sc_cli_tools_rag/single-cell-cli-tools`.

**Note that, if a default RAG database is not found, the CLI will automatically run with RAG functionality disabled.**

## `3` [RAG System Setup](#3-rag-system-setup)

To use the RAG knowledge base, build it from the provided configuration:

```bash
python -m pantheon.toolsets.utils.rag build \
    pantheon/cli/rag_system_config.yaml \
    tmp/pantheon_cli_tools_rag
```

This creates a vector database at `tmp/pantheon_cli_tools_rag/pantheon-cli-tools` with genomics tools documentation.


### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--rag_db` | Path to RAG database | `tmp/pantheon_cli_tools_rag/pantheon-cli-tools` |
| `--model` | AI model to use | Loaded from config or `gpt-4.1` |
| `--agent_name` | Name of the agent | `general_bot` |
| `--workspace` | Working directory | Current directory |
| `--instructions` | Custom instructions | Built-in instructions |
| `--disable_rag` | Disable RAG toolset | `False` |
| `--disable_web` | Disable web toolset | `False` |
| `--disable_notebook` | Disable notebook toolset | `False` |
| `--disable_r` | Disable R interpreter toolset | `False` |
| `--disable_julia` | Disable Julia interpreter toolset | `False` |
| `--disable_code_validator` | Disable code validation toolset | `False` |
| `--disable_bio` | Disable bio analysis toolsets | `False` |
| `--disable_ext` | Disable external toolsets loader | `True` |
| `--ext_toolsets` | Comma-separated list of external toolsets to load | All available |
| `--ext_dir` | Directory containing external toolsets | `./ext_toolsets` |

## `4` [Core Features](#4-core-features)

### **AI-Driven Scientific Intelligent Agent**
Built-in intelligent agent designed specifically for scientific computing, capable of handling various complex data analysis tasks. The intelligent agent not only executes commands but also:
- **Understands scientific context**: Knows what type of analysis you're doing
- **Recommends best methods**: Automatically selects appropriate algorithms and parameters  
- **Explains analysis results**: Provides professional biological interpretations
- **PhD-level domain knowledge**: Context-aware professional advice

### **Hybrid Programming Paradigm**
Seamlessly switch between multiple programming approaches within the same environment:
- **Variable persistence**: Python/R/Julia variables directly shared in memory
- **Natural language-driven**: Fluidly transforms thoughts into code execution
- **Multi-language support**: Python, R, Julia in the same session
- **Tool integration**: Access to comprehensive scientific computing ecosystems

### **Open Source & Privacy-First**
- **Fully Open Source**: Transparent, auditable source code
- **Data Privacy Protection**: All computation performed locally
- **Local model support**: Can be used completely offline
- **Zero data uploads**: Research data never leaves your control
- **Infinitely Extensible**: Based on Python ecosystem

## `5` [Available Tools](#5-available-tools)

### Core Tools (Always Enabled)
- **Shell**: System commands and genomics tools with auto-installer
- **Python**: Data analysis and visualization (pandas, matplotlib, scanpy)
- **R**: Statistical analysis and Seurat single-cell workflows with sample data
- **Julia**: High-performance scientific computing (DataFrames.jl, Plots.jl, DifferentialEquations.jl)
- **File Editor**: Read, edit, and create files with diffs
- **Code Search**: Find files (glob), search content (grep), list directories (ls)
- **Code Validation**: Verify Python code, commands, function calls, and detect common errors
- **Todo**: Claude Code-style task management with smart task breakdown and auto-progression
- **Generator**: AI-powered external toolset creation for any domain
- **Bio Tools**: Comprehensive bioinformatics analysis pipelines (ATAC-seq, RNA-seq, etc.)

### Optional Tools
- **RAG**: Vector-based knowledge search with built-in scientific knowledge
- **Web**: Intelligent web operations with automatic URL intent analysis
- **Notebook**: Jupyter notebook editing (no execution)

## `6` [Configuration Files](#6-configuration-files)

Pantheon CLI supports project-specific configuration files similar to Claude Code's `CLAUDE.md`:

- **`PANTHEON.md`**: Project-wide configuration, commands, and guidelines (safe to commit)
- **`PANTHEON.local.md`**: Personal preferences and local settings (add to `.gitignore`)

These files are automatically discovered in your current directory or any parent directory and integrated into the AI assistant's context.

**Example `PANTHEON.md`:**
```markdown
# My Project

## Commands
- Run analysis: `python scripts/analyze.py`
- Quick data load: `%adata = sc.read_h5ad('data.h5ad')`

## Guidelines  
- Use scanpy for Python analysis
- Use Seurat for R analysis
```

See [`CONFIG_FILES.md`](CONFIG_FILES.md) for detailed documentation and examples.





## `7` [Architecture](#7-architecture)

Pantheon-CLI is built as a standalone package that depends on:

- **pantheon-agents**: Core agent functionality and reasoning
- **pantheon-toolsets**: Distributed toolsets for various tasks
- Clean separation of concerns with modular design
- Enterprise-grade distributed architecture

### Package Structure

```
Pantheon-cli/
├── pantheon_cli/              # Main package (renamed to avoid conflicts)
│   ├── __init__.py           # Entry point with cli_main()
│   ├── cli/                  # CLI implementation
│   │   ├── core.py          # Main CLI logic with toolset integration
│   │   └── manager/         # API key and model management
│   └── repl/                # REPL implementation  
│       ├── core.py          # REPL core with updated imports
│       ├── ui.py            # User interface and tool call display
│       └── bio_handler.py   # Bio command handling
├── pyproject.toml           # Package configuration
└── README.md               # This file
```


## `8` [Requirements](#8-requirements)

- Python 3.10+
- Required packages: `fire`, `rich`, `pantheon-agents`, `pantheon-toolsets`, `hypha_rpc`, `pandas`
- Optional: R for statistical analysis, Julia for high-performance computing

## `9` [Real Application Scenarios](#9-real-application-scenarios)

### Biomedical Research
```bash
> I have a 10x single-cell dataset and want to analyze T cell differentiation trajectories

> Load spatial transcriptomics data and identify gene expression patterns in tissue structures

> Integrate multi-omics data to find disease-related biomarkers
```

### Data Science Analysis
```bash
> Perform time series analysis and forecasting on this sales data

> Build a machine learning model to predict customer churn

> Use deep learning to analyze these medical images
```

### Teaching and Learning
```bash
> Explain the principles of principal component analysis and demonstrate with code

> Compare the performance of different clustering algorithms on this dataset

> Show how to perform statistical analysis for A/B testing
```

## `10` [Why Choose Pantheon-CLI?](#10-why-choose-pantheon-cli)

### Problems It Solves

**Tool Fragmentation**
- *Traditional*: Switch between multiple tools, data passed through file save/load
- *Pantheon-CLI*: Hybrid programming with variable persistence, direct memory sharing

**High Technical Barriers** 
- *Traditional*: Master multiple programming languages and complex tool chains
- *Pantheon-CLI*: "Vibe analysis" that understands research intent and auto-matches tools

**Lack of Intelligent Guidance**
- *Traditional*: Manual method selection, parameter adjustment, result interpretation  
- *Pantheon-CLI*: PhD-level domain knowledge with context-aware professional advice

**Data Privacy Concerns**
- *Traditional*: Many AI tools require uploading sensitive data to cloud
- *Pantheon-CLI*: Completely offline scientific computing AI with local models




