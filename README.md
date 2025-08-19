# Pantheon CLI

<!-- Banner -->
<p align="center">
  <img src="assets/pantheon_banner_tri_dark.svg#gh-dark-mode-only" alt="Pantheon ASCII (dark)" />
  <img src="assets/pantheon_banner_tri_light.svg#gh-light-mode-only" alt="Pantheon ASCII (light)" />
</p>

<div align="center">

***We're not just building another CLI tool.  
We're re-defining how scientists will interact with data in the AI era.***

</div>

<div align="center">

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Status: Beta](https://img.shields.io/badge/Status-Beta-orange.svg)]()
[![AI-Native](https://img.shields.io/badge/AI-Native-purple.svg)]()

</div>


## Quick Start

### Experience Mixed Programming
```bash
# Start Pantheon CLI
pantheon-cli

# Now you can seamlessly mix natural language, Python, R, and Julia:
> Solve optimization problem with Julia, visualize results in Python
> Analyze single-cell data using Seurat, then create a Python visualization
> Load my_data.csv, fit the data with differential equations in Julia, plot with R
```

### Installation

```bash
# Install from source (recommended for development)
git clone https://github.com/aristoteleo/pantheon-cli.git
cd Pantheon-cli
pip install -e .

# Make sure dependencies are installed
pip install pantheon-agents pantheon-toolsets
```

**Note**: Pantheon-CLI requires both `pantheon-agents` and `pantheon-toolsets` to be installed. These provide the core agent functionality and distributed toolsets respectively. The complete list of features implemented in `pantheon-agents` and `pantheon-toolsets` will be introduced elsewhere.

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

# Once the cli is running, you will need to setup your API keys and model preferences.
/api-key list  # List current API keys
# The following guidance will then show up to help you set your API keys:
💡 Usage:
  /api-key list - Show this status
  /api-key <provider> <key> - Set API key
  Examples:
    /api-key openai sk-... - Set OpenAI key
    /api-key anthropic sk-... - Set Anthropic key
    /api-key google ai... - Set Google key
    /api-key deepseek sk-... - Set DeepSeek key
    /api-key qwen sk-... - Set Qwen key
    /api-key kimi sk-... - Set Kimi key
    /api-key grok sk-... - Set Grok key

# Ensure you have proper internet and API connection. You can now use the CLI for various tasks, enjoy!
```

### With RAG Database

If you have a RAG database prepared:

```bash
pantheon-cli --rag_db path/to/rag/database
```

Default RAG database location: `tmp/sc_cli_tools_rag/single-cell-cli-tools`.

**Note that, if a default RAG database is not found, the CLI will automatically run with RAG functionality disabled.**

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

## Configuration Files

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
- Package renamed from `pantheon-cli` to `pantheon_cli` to avoid import conflicts
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


