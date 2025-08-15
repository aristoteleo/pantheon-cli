# Pantheon CLI Configuration Files

Pantheon CLI supports local configuration files similar to Claude Code's `CLAUDE.md` functionality. These files allow you to provide project-specific context, custom instructions, and preferences to the AI assistant.

## Configuration Files

### `PANTHEON.md`
- **Purpose**: Project-wide configuration and documentation
- **Location**: Place in your project root directory
- **Content**: Project information, commands, workflows, and team guidelines
- **Sharing**: Safe to commit to version control and share with team

### `PANTHEON.local.md`  
- **Purpose**: Personal preferences and local environment settings
- **Location**: Place in your project root directory
- **Content**: Personal shortcuts, local paths, development preferences
- **Sharing**: Should be added to `.gitignore` - keep it local and personal

## How It Works

1. **Automatic Discovery**: Pantheon CLI automatically searches for these files in:
   - Current working directory
   - All parent directories (walking up the directory tree)

2. **Loading Order**: 
   - `PANTHEON.md` is loaded first (project configuration)
   - `PANTHEON.local.md` is loaded second (local overrides)

3. **Integration**: The content is automatically appended to the AI assistant's system prompt

4. **Graceful Fallback**: If no configuration files exist, Pantheon CLI works normally without any errors

## File Examples

### Example `PANTHEON.md`
```markdown
# PANTHEON.md

## Project: Single-Cell RNA Analysis

### Commands
- Run analysis: `python scripts/analyze.py`
- Generate plots: `python scripts/visualize.py`

### Data Paths
- Raw data: `data/raw/`
- Processed: `data/processed/`

### Analysis Guidelines
- Use scanpy for Python analysis
- Use Seurat for R analysis
- Save results in `results/` directory
```

### Example `PANTHEON.local.md`
```markdown
# PANTHEON.local.md

## Personal Settings
- Preferred tool: Python with scanpy
- Local data: `/Users/myname/data/`
- Backup location: `/Users/myname/Dropbox/results/`

## Personal Shortcuts
- Quick load: `%adata = sc.read_h5ad('my_data.h5ad')`
```

## Best Practices

### What to Include in `PANTHEON.md`
- Project overview and goals
- Standard commands and workflows  
- File organization and naming conventions
- Team coding standards
- Environment setup instructions
- Common troubleshooting tips

### What to Include in `PANTHEON.local.md`
- Personal development preferences
- Local file paths and directories
- Personal shortcuts and helper functions
- Individual analysis notes
- Local environment configuration

## Advanced Usage

### Multi-Level Configuration
Configuration files are searched in directory hierarchy:
```
/project/
├── PANTHEON.md              # Project root config
└── analysis/
    ├── PANTHEON.local.md     # Analysis-specific local config  
    └── experiment1/
        └── [current directory]
```

The CLI will find and load both files automatically.

### Integration with Direct Execution
Configuration can include examples of direct execution commands:

```markdown
## Quick Commands
- Load data: `%adata = sc.read_h5ad('data.h5ad')`
- Basic QC: `%sc.pp.calculate_qc_metrics(adata)`
- Plot genes: `%sc.pl.violin(adata, ['n_genes'])`
- R analysis: `>obj <- readRDS('seurat_obj.rds')`
- Shell tasks: `!ls data/`
```

## Version Control

### Recommended `.gitignore` entries:
```gitignore
# Personal configuration - keep local
PANTHEON.local.md

# Project configuration - safe to share
# PANTHEON.md  # <- Don't ignore this one
```

## Migration from Other Tools

### From Claude Code
If you have existing `CLAUDE.md` files, you can:
1. Copy `CLAUDE.md` to `PANTHEON.md`
2. Adapt any Claude-specific instructions to Pantheon CLI context
3. Add Pantheon-specific commands and workflows

### From Other Documentation
Convert existing project documentation:
1. Extract key workflows and put them in `PANTHEON.md`
2. Add personal notes to `PANTHEON.local.md`
3. Include direct execution examples using `!`, `%`, `>`, `]` prefixes

## Troubleshooting

### Configuration Not Loading
- Check file encoding (must be UTF-8)
- Verify file permissions (must be readable)
- Ensure files are in search path (current or parent directories)

### Conflicting Instructions
- `PANTHEON.local.md` takes precedence over `PANTHEON.md`
- Later sections override earlier sections
- Personal preferences in local file override project defaults

### Debugging
To see which configuration files were loaded:
- Look for startup message: `📋 Loaded configuration from: ...`
- No message means no configuration files were found (this is normal)

## Example Templates

See the included example files:
- `PANTHEON.md.example` - Project configuration template
- `PANTHEON.local.md.example` - Personal configuration template

Copy these files and remove the `.example` extension to start using them.