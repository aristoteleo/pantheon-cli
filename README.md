# Pantheon CLI

<!-- Banner -->
<p align="center">
  <img src="assets/pantheon_ascii_gradient.svg#gh-dark-mode-only" width="600" alt="Pantheon ASCII (dark)" />
  <img src="assets/pantheon_ascii_gradient_light.svg#gh-light-mode-only" width="600" alt="Pantheon ASCII (light)" />
</p>

<p align="center">
  <img src="assets/pantheon_cli_output_dark.svg#gh-dark-mode-only" width="600" alt="Pantheon CLI output (dark)" />
  <img src="assets/pantheon_cli_output_light.svg#gh-light-mode-only" width="600" alt="Pantheon CLI output (light)" />
</p>




## Features

- **Paradigm Shift**: The first CLI that truly integrates natural language with programming languages (Python, R, Julia), enabling scientists to think and analyze without switching contexts.
- **Open Source**: The first Python-Based CLI.
- **Science-First**: Purpose-built for scientific workflows with comprehensive genomics pipelines, not a generic tool trying to do everything
- **Distributed by Design**: Built on enterprise-grade distributed architecture, not a monolithic application that breaks at scale  
- **AI-Native**: Every component is designed from the ground up for AI interaction, not traditional software with AI bolted on
- **Future-Proof**: MCP protocol support and extensible architecture means it grows with the ecosystem, not against it


## Quick Start

### Experience Mixed Programming
```bash
# Start Pantheon CLI
pantheon-cli

# Now you can seamlessly mix natural language, Python, R, and Julia:
> Analyze single-cell data using Seurat, then create a Python visualization
> 加载PBMC数据，用R做预处理，然后用Python画UMAP图
> Load my_data.csv, run differential equations in Julia, plot with R
> Solve optimization problem with Julia, visualize results in Python
```

### Installation

```bash
# Install from source (recommended for development)
cd Pantheon-cli
pip install -e .

# Make sure dependencies are installed
pip install pantheon-agents pantheon-toolsets
```

**Note**: Pantheon-CLI requires both `pantheon-agents` and `pantheon-toolsets` to be installed. These provide the core agent functionality and distributed toolsets respectively.

### Basic Usage

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

Default RAG database location: `tmp/sc_cli_tools_rag/single-cell-cli-tools`

## RAG System Setup

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

## Available Tools

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
- **RAG**: Vector-based knowledge search (requires database)
- **Web**: Intelligent web operations with automatic URL intent analysis
- **Notebook**: Jupyter notebook editing (no execution)





## Architecture

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

**Key Design Decisions:**
- Package renamed from `pantheon` to `pantheon_cli` to avoid import conflicts
- All relative imports converted to absolute imports (`from pantheon.agent import Agent`)
- Uses local REPL instead of `agent.chat()` to avoid import issues
- Graceful fallback for missing toolsets with error handling

## Requirements

- Python 3.10+
- Required packages: `fire`, `rich`, `pantheon-agents`, `pantheon-toolsets`, `hypha_rpc`, `pandas`
- Optional: R for statistical analysis, Julia for high-performance computing

## Tips

1. **RAG Database**: If not found, the CLI will automatically disable RAG functionality
2. **Memory Optimization**: Use `--disable-*` flags to reduce memory usage for unused tools
3. **Workspace Management**: Defaults to current directory, change with `--workspace`
4. **Custom Instructions**: Completely customize agent behavior and specialization
5. **Web Operations**: The AI automatically analyzes URL intent - just paste URLs naturally
6. **Bio Analysis**: Use `/bio list` to see all available analysis tools
7. **External Toolsets**: Generate custom tools for any domain with `generate_toolset`
8. **Todo Management**: Let the AI break down complex tasks automatically with `add_todo`
9. **Code Validation**: Always validate generated code with built-in validation tools
10. **Julia Computing**: Use Julia for high-performance numerical computing and scientific packages
11. **Multilingual**: Works seamlessly in both English and Chinese contexts


