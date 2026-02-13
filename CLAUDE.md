# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```bash
# Install in development mode
pip install -e .

# Install with notebook support
pip install -e '.[notebook]'

# Install with dev dependencies (pytest)
pip install -e '.[dev]'

# Run tests
pytest

# Run a single test
pytest tests/test_file.py::test_function

# Start the CLI
pantheon-cli

# Start with specific model
pantheon-cli --model claude-sonnet-4-20250514

# Start without RAG
pantheon-cli --disable_rag

# Build RAG database (requires SCRAPER_API_KEY env var)
pantheon-cli --build-rag [output_dir] [--rag-config config.yaml]
```

## Architecture Overview

Pantheon-CLI is a scientific AI assistant built on a **persistent interpreter paradigm** - Python/R/Julia interpreters maintain session state across all tool calls, avoiding redundant file I/O for large datasets.

### Entry Point Flow

```
pantheon_cli/__init__.py::cli_main()
  → cli/core.py::cli() [Fire CLI]
    → cli/core.py::main()
      → Agent initialization + Toolset registration
      → repl/core.py::Repl.run()
```

### Key Components

| Component | Path | Purpose |
|-----------|------|---------|
| **CLI Core** | `pantheon_cli/cli/core.py` | Main entry, toolset init, agent creation, DEFAULT_INSTRUCTIONS |
| **REPL** | `pantheon_cli/repl/core.py` | Interactive loop, command parsing, tool call display |
| **REPL UI** | `pantheon_cli/repl/ui.py` | Rich console formatting |
| **API Key Manager** | `pantheon_cli/cli/manager/api_key_manager.py` | Multi-provider API key handling |
| **Model Manager** | `pantheon_cli/cli/manager/model_manager.py` | Model switching |
| **Config Loader** | `pantheon_cli/utils/config_loader.py` | Loads PANTHEON.md / PANTHEON.local.md |
| **Bio Handler** | `pantheon_cli/repl/bio_handler.py` | `/bio` command processing |
| **Bio Prompts** | `pantheon_cli/cli/prompt/` | 13+ specialized analysis workflows (ATAC, RNA-seq, spatial, etc.) |
| **Dev Loop Mode** | `pantheon_cli/cli/modes/dev_loop.py` | Autonomous plan → code → review workflow |

### Dependencies

- `pantheon-agents`: Agent framework for reasoning/planning
- `pantheon-toolsets`: Toolset implementations (Python/R/Julia interpreters, file ops, web, bio tools)

### REPL Command Prefixes

- `!command` → Shell execution
- `%python_code` → Direct Python execution (stateful)
- `>r_code` → Direct R execution (stateful)
- `]julia_code` → Direct Julia execution (stateful)
- `/api-key` → Configure API keys
- `/bio <tool> <command>` → Bio analysis tools
- `/clear`, `/exit` → Session management

## Configuration

### API Key Priority

1. Environment variables (`OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, etc.)
2. Local config: `.pantheon_config.json` (current directory)
3. Global config: `~/.pantheon/config.json`

### Project Config Files

- `PANTHEON.md` - Project-wide config (safe to commit)
- `PANTHEON.local.md` - Personal config (add to .gitignore)

These are auto-discovered in current + parent directories and appended to agent instructions.

## Dev Loop Mode

Autonomous plan → code → review cycle:

```bash
pantheon-cli --mode devloop --dev_goal "your goal" --max_iters 10
```

Uses different models for each phase (plan: gpt-5, code: gpt-4.1, review: gpt-5).

## Toolsets

Core toolsets are always enabled: Shell, Python, R, Julia, File Editor, Code Search, Code Validation, Todo, Generator, Bio Tools.

Optional toolsets can be disabled via flags: `--disable_rag`, `--disable_web`, `--disable_dr`, `--disable_notebook`, `--disable_r`, `--disable_julia`, `--disable_bio`, etc.

## Testing

```bash
# Pytest with async support
pytest

# Test config in pyproject.toml uses:
#   asyncio_mode = "auto"
#   asyncio_default_fixture_loop_scope = "function"
```
