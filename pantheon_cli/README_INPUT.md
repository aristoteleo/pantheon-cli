# Pantheon CLI: Natural Language to Tool Mapping Guide

This guide helps you understand which tools to call based on your natural language requests in the Pantheon CLI. Simply describe what you want to do, and find the corresponding tool function.

## How to Use This Guide

When you have a task in mind, find it in the **"What You Want to Say"** sections below, then use the corresponding **Tool Call** in the CLI.

---

## Task Management & Planning

### What You Want to Say → Tool to Call

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "I need to track my project progress" | `add_todo("Complete data analysis project")` | Creates a new task with automatic breakdown into subtasks |
| "Show me what I'm working on" | `show_todos()` | Displays all current tasks with status and progress |
| "What should I do next?" | `execute_current_task()` | Provides intelligent guidance for the current active task |
| "I finished this step" | `mark_task_done("Completed data preprocessing")` | Marks current task as complete and moves to next task |
| "Clear all my tasks" | `clear_all_todos()` | Removes all tasks to start fresh |
| "Remove completed tasks" | `clear_completed_todos()` | Cleans up finished tasks while keeping active ones |
| "Update task status" | `update_todo_status("task_id", "completed")` | Changes status of specific task |
| "Mark specific task complete" | `complete_todo("task_id")` | Marks specific todo item as completed |
| "Start working on task" | `start_todo("task_id")` | Sets specific task as in progress |
| "Delete a task" | `remove_todo("task_id")` | Removes specific todo from list |
| "Get next task" | `get_next_todo()` | Returns the next pending task to work on |
| "Work on next item" | `work_on_next_todo()` | Starts working on the next pending task |

---

## Bioinformatics & ATAC-seq Analysis

**For comprehensive bioinformatics and ATAC-seq analysis tools, see the detailed guide:**

📖 **[Bioinformatics Tools Guide](README_BIOINFORMATICS.md)**

**Quick Overview:**
- **ATAC-seq Analysis**: `/bio atac init` (enter mode), `/bio atac upstream <folder>` (run analysis)
- **Single-Cell ATAC**: `/bio scatac init`, `/bio scatac upstream <folder>` (Cell Ranger ATAC toolset)
- **Single-Cell RNA**: `/bio scrna init`, `/bio scrna analysis <file>` (omicverse-powered analysis)
- **Core Functions**: Species detection, genome setup, quality control, alignment, peak calling, cell annotation
- **Key Tools**: `bio_atac_scan_folder()`, `bio_scatac_install()`, `bio_scrna_analysis()`, `bio_scrna_annotate()`

The bioinformatics guide includes detailed mappings for:
- Project setup and species detection
- Genome resource management
- Quality control and preprocessing  
- Alignment and BAM processing
- Peak calling and analysis
- Visualization and reporting
- Single-cell ATAC-seq workflows (Cell Ranger ATAC)
- Single-cell RNA-seq analysis (omicverse integration)
- Cell type annotation and subtype analysis
- Tool installation and dependency management

---

## Bioinformatics Analysis Tools

### 🧬 GeneAgent - Original 7-Step Iterative Verification System

GeneAgent is a specialized biological analysis tool that uses a unique 7-step iterative verification methodology. It combines LLM reasoning with real biological database APIs to produce scientifically validated gene set analyses.

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Analyze genes TP53,BRCA1,EGFR using GeneAgent" | `bio GeneAgent TP53,BRCA1,EGFR` | Run full 7-step iterative verification for gene set |
| "Perform enrichment with genes ERBB2,ERBB4,FGFR2,FGFR4,HRAS,KRAS" | `bio GeneAgent ERBB2,ERBB4,FGFR2,FGFR4,HRAS,KRAS --analysis_type enrichment` | Enrichment analysis with biological verification |
| "Run comprehensive GeneAgent analysis on MYC,JUN,FOS" | `bio GeneAgent MYC,JUN,FOS --analysis_type comprehensive` | Complete analysis with all biological data sources |
| "Study protein interactions for immune genes CD4,CD8A,IL2" | `bio GeneAgent CD4,CD8A,IL2 --analysis_type interactions` | Protein interaction network analysis with validation |
| "GeneAgent analysis of cytokines with summary output" | `bio GeneAgent IFNG,TNF,IL6 --output_format summary` | Concise evidence-based summary analysis |
| "Verify biological claims for tumor suppressor genes" | `bio GeneAgent TP53,RB1,BRCA1 --save_results true` | Scientific validation with saved results |

#### GeneAgent 7-Step Process:
1. **Generate baseline analysis** (LLM) - Initial biological process analysis
2. **Extract testable biological claims** (LLM) - Generate verifiable statements  
3. **Verify claims using 8 real biological APIs** - NCBI, g:Profiler, Enrichr, STRING, etc.
4. **Modify analysis based on verification evidence** (LLM) - Evidence-driven refinement
5. **Generate new analysis claims** (LLM) - Second round claim generation
6. **Second verification round with biological APIs** - Additional evidence gathering
7. **Final evidence-based synthesis** (LLM) - Scientifically validated conclusions

#### Real Biological Data Sources:
- **NCBI E-utilities** - Gene function summaries and annotations
- **g:Profiler** - GO enrichment analysis and functional profiling  
- **Enrichr** - Pathway analysis (KEGG, Reactome, MSigDB databases)
- **STRING/BioGRID** - Protein-protein interaction networks
- **DisGeNET** - Disease-gene associations and phenotypes
- **InterPro** - Protein domain and family classifications
- **CORUM** - Protein complex membership data
- **PubMed** - Literature evidence and citations

#### Analysis Types:
- `comprehensive`: Full 7-step verification with all biological data sources
- `functional`: Focus on biological functions and cellular processes
- `enrichment`: GO/KEGG enrichment analysis with database verification
- `interactions`: Protein interaction networks with experimental validation  
- `clinical`: Disease associations and therapeutic target analysis
- `custom`: Custom questions answered with biological evidence

#### What Makes GeneAgent Special:
- ✅ **Evidence-Based**: Every biological claim verified against real databases
- ✅ **Iterative Refinement**: Analysis improved through multiple verification rounds  
- ✅ **Scientifically Validated**: Final results are publication-quality
- ✅ **No Hallucinations**: Zero fabricated biological information
- ✅ **Real APIs**: Uses 8 authoritative biological databases
- ✅ **Original Methodology**: Implements the complete research paper algorithm

### ATAC-seq Analysis

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Initialize ATAC-seq project" | `/bio atac init` | Enter ATAC-seq analysis mode and clear existing todos |
| "Analyze ATAC data in folder" | `/bio atac upstream ./data` | Run complete upstream ATAC-seq analysis pipeline |
| "Scan folder for ATAC data" | `bio_atac_scan_folder("./data")` | Scan directory for ATAC-seq FASTQ files |
| "Check ATAC dependencies" | `bio_atac_check_dependencies()` | Verify ATAC analysis tool installation |
| "Setup genome resources" | `bio_atac_setup_genome_resources("human", "hg38")` | Download genome, GTF, and blacklist files |
| "Auto-align FASTQ files" | `bio_atac_auto_align_fastq("./data")` | Automated alignment pipeline with BWA-MEM2 |
| "Call peaks with MACS2" | `bio_atac_call_peaks_macs2("./bam_files")` | Peak calling using MACS2 algorithm |
| "Generate QC report" | `bio_atac_generate_atac_qc_report("./results")` | Create comprehensive quality control report |

### Single-Cell ATAC-seq Analysis

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Initialize scATAC project" | `/bio scatac init` | Enter single-cell ATAC-seq analysis mode |
| "Install Cell Ranger ATAC" | `/bio scatac install` | Download and install cellranger-atac v2.2.0 |
| "Analyze 10X Chromium ATAC data" | `/bio scatac upstream ./10x_data` | Run cellranger-atac upstream analysis pipeline |
| "Process single sample" | `/bio scatac count Sample1` | Run cellranger-atac count for single sample |

### Single-Cell RNA-seq Analysis

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Initialize scRNA project" | `/bio scrna init` | Enter single-cell RNA-seq analysis mode |
| "Load and analyze scRNA data" | `/bio scrna analysis ./data.h5ad` | Comprehensive scRNA-seq analysis with omicverse |
| "Run quality control" | `/bio scrna qc ./data.h5ad` | Quality control analysis with metrics and filtering |
| "Preprocess scRNA data" | `/bio scrna preprocess ./data.h5ad` | Normalization and preprocessing pipeline |
| "Annotate cell types" | `/bio scrna annotate ./data.h5ad` | Cell type annotation using CellOntologyMapper |
| "Perform subtype analysis" | `/bio scrna subtype ./data.h5ad` | Advanced cell subtype identification and analysis |

### General Bio Tools Management

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "List all bio tools" | `bio list` | Display all available bioinformatics analysis tools |
| "Get tool information" | `bio info atac` | Get detailed information about specific bio tool |
| "Get bio help" | `bio help scrna` | Get help and usage examples for bio tools |

---

## Direct Code Execution (REPL-Style Commands)

### Quick Execution Prefixes

| **Prefix** | **Language** | **Example** | **Description** |
|------------|--------------|-------------|-----------------|
| `!` | **Shell/Bash** | `!ls -la` | Execute shell commands directly in REPL |
| `%` | **Python** | `%print("Hello World")` | Execute Python code directly with instant output |
| `>` | **R** | `>summary(mtcars)` | Execute R statistical code directly |
| `]` | **Julia** | `]println("Hello Julia")` | Execute Julia high-performance computing code |

### Direct Execution Examples

```bash
# Shell commands - instant system operations
!pwd                           # Show current directory
!ls -la *.py                  # List Python files
!git status                   # Check git status
!pip install pandas          # Install Python packages

# Python code - data analysis and scripting
%import pandas as pd         # Import libraries
%df = pd.read_csv("data.csv") # Load data
%df.head()                   # Show first rows
%x = [1,2,3]; sum(x)         # Quick calculations

# R statistical computing - data science
>library(ggplot2)            # Load R packages
>data(mtcars)                # Load dataset
>summary(mtcars)             # Statistical summary
>plot(mtcars$mpg, mtcars$hp) # Create plots

# Julia high-performance computing - numerical analysis
]using LinearAlgebra         # Load Julia packages
]A = [1 2; 3 4]             # Matrix operations
]det(A)                     # Calculate determinant
]using Plots; plot([1,2,3]) # Create visualizations
```

### Key Features

- **Shared Interpreter Sessions**: All prefixed commands use the same interpreter session as the AI agent
- **REPL-Style Output**: Only shows final results and print statements, suppressing intermediate assignments
- **Instant Execution**: No need to call functions - just type the prefix and your code
- **Language Switching**: Switch between languages seamlessly in the same session
- **Persistent State**: Variables and imports persist across commands within each language

---

## Python Code Execution

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Run this Python code" | `run_python("print('Hello World')")` | Executes Python code with full package support |
| "Quick Python execution" | `%print("Hello World")` | **NEW**: Direct Python execution with `%` prefix |
| "Run code in specific interpreter" | `run_code_in_interpreter("interpreter_id", "print('test')")` | Executes code in specific Python interpreter session |
| "Create new Python interpreter" | `new_interpreter("my_python_session")` | Creates new isolated Python interpreter session |
| "Delete Python interpreter" | `delete_interpreter("interpreter_id")` | Removes Python interpreter session |

---

## R Statistical Computing

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Run R statistical analysis" | `run_r("summary(mtcars)")` | Executes R code for statistical computing |
| "Quick R execution" | `>summary(mtcars)` | **NEW**: Direct R execution with `>` prefix |
| "Run R code in specific session" | `run_code_in_interpreter("r_session", "plot(iris)")` | Executes R code in specific interpreter session |
| "Create new R interpreter" | `new_interpreter("my_r_session")` | Creates new isolated R interpreter session |
| "Delete R interpreter" | `delete_interpreter("r_session_id")` | Removes R interpreter session |
| "Get R interpreter output" | `get_interpreter_output("r_session_id")` | Retrieves output from R interpreter session |

---

## Julia High-Performance Computing

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Run Julia scientific computing" | `julia("using LinearAlgebra; det([1 2; 3 4])")` | Executes Julia code for high-performance numerical analysis |
| "Quick Julia execution" | `]println("Hello Julia")` | **NEW**: Direct Julia execution with `]` prefix |
| "Create Julia environment" | `new_interpreter("julia_session")` | Creates new isolated Julia interpreter session |
| "Run Julia code in session" | `run_code_in_interpreter("julia_session", "using Plots; plot([1,2,3])")` | Executes Julia code in specific interpreter session |
| "Install Julia package" | `julia_install_package("DataFrames")` | Installs Julia packages for scientific computing |
| "Load sample data in Julia" | `load_sample_data("iris")` | Loads sample datasets for analysis |
| "Quick Julia analysis" | `quick_analysis(df)` | Performs automated data analysis workflow |
| "Auto-save Julia plots" | `auto_save_plot()` | Automatically saves generated plots |
| "Delete Julia interpreter" | `delete_interpreter("julia_session_id")` | Removes Julia interpreter session |

---

## File Operations & Management

### Basic File Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "List files in directory" | `list_files("/path/to/directory")` | Lists all files and directories in specified path |
| "Read this file" | `read_file("/path/to/file.txt")` | Reads file content with syntax highlighting |
| "Write content to file" | `write_file("/path/to/output.txt", "content")` | Writes text content to specified file |
| "Create directory" | `create_directory("/path/to/new/folder")` | Creates directory structure recursively |
| "Delete file" | `delete_file("/path/to/unwanted.txt")` | Removes specified file |
| "Delete directory" | `delete_directory("/path/to/folder")` | Removes directory and all contents |
| "Move file" | `move_file("/old/path/file.txt", "/new/path/file.txt")` | Moves file to new location |
| "Show file tree" | `list_file_tree("/project")` | Displays hierarchical directory structure |
| "View images" | `observe_images("/path/to/images")` | Displays image files in directory |
| "Read PDF file" | `read_pdf("/path/to/document.pdf")` | Extracts and displays PDF content |

---

## Text Editing & File Modification

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Edit specific line in file" | `edit_file("/path/file.py", line_number=10, new_content="updated line")` | Modifies specific lines in text files |
| "Search for text in file" | `search_in_file("/path/file.txt", "search_term")` | Searches for specific text within a file |
| "Insert text at line" | `insert_at_line("/path/file.py", line_number=5, content="import pandas")` | Inserts new content at specified line number |
| "Delete lines from file" | `delete_lines("/path/file.txt", start_line=10, end_line=15)` | Removes specified range of lines from file |
| "Create new file" | `create_file("/path/to/new_file.txt", "initial content")` | Creates new file with specified content |

---

## Code Search & Analysis

### Basic Search Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Search files by pattern" | `glob("/project", "*.py")` | Finds files matching glob pattern |
| "Search text in files" | `grep("/project", "TODO", "*.py")` | Searches for text patterns in files |
| "List directory contents" | `ls("/path/to/directory")` | Lists files and directories with details |

---

## Code Quality & Validation

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Check Python syntax" | `validate_python_code("/path/script.py")` | Validates Python code syntax and structure |
| "Validate shell command" | `validate_command("complex_command")` | Validates shell command syntax |
| "Find common code errors" | `detect_common_errors("/path/script.py")` | Detects common programming errors |
| "Suggest function alternatives" | `suggest_function_alternatives("old_function")` | Recommends modern alternatives to functions |
| "Validate function call" | `validate_function_call("function_name", ["arg1", "arg2"])` | Checks function call syntax and arguments |
| "Check import statements" | `validate_imports("/path/script.py")` | Validates import statements in Python files |
| "Check code style" | `check_code_style("/path/script.py")` | Analyzes code style and formatting |

---

## Jupyter Notebook Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Show notebook contents" | `read_notebook("/path/analysis.ipynb")` | Displays notebook with beautiful formatting |
| "Edit notebook cell" | `edit_notebook_cell("/path/notebook.ipynb", cell_index=2, new_content="code")` | Modifies specific notebook cells |
| "Add new cell" | `add_notebook_cell("/path/notebook.ipynb", content="# Analysis", cell_type="markdown")` | Inserts new cells at specified positions |
| "Delete notebook cell" | `delete_notebook_cell("/path/notebook.ipynb", cell_index=3)` | Removes specific cell from notebook |
| "Create new notebook" | `create_notebook("/path/new_analysis.ipynb")` | Creates new empty Jupyter notebook |
| "Copy notebook cell" | `copy_notebook_cell("/path/notebook.ipynb", source_index=1, target_index=5)` | Copies cell to different position |
| "Move notebook cell" | `move_notebook_cell("/path/notebook.ipynb", from_index=2, to_index=8)` | Moves cell to different position |
| "Add notebook template" | `add_notebook_template("/path/notebook.ipynb", template_type="data_analysis")` | Adds predefined template cells |

---

## Web Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Fetch web page content" | `web_fetch("https://example.com")` | Retrieves and displays web page content |
| "Search the internet" | `web_search("Python data analysis tutorials")` | Searches web using search engine |
| "Google search for information" | `google_search("ATAC-seq analysis methods")` | Performs Google search and returns results |
| "Scrape web page data" | `fetch_web_page("https://site.com")` | Extracts structured data from web pages |

---

## System & Shell Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Run shell command directly" | `!ls -la` | **NEW**: Direct shell execution with `!` prefix |
| "Execute system command" | `!pwd` | **NEW**: Instant system command execution |
| "Create new shell session" | `new_shell("my_shell_session")` | Creates new isolated shell environment |
| "Close shell session" | `close_shell("shell_session_id")` | Terminates specific shell session |
| "Run command in shell" | `run_command_in_shell("shell_id", "ls -la")` | Executes command in specific shell session |
| "Get shell output" | `get_shell_output("shell_session_id")` | Retrieves output from shell session |
| "Run single command" | `run_command("ps aux")` | Executes single shell command |

---

## File Transfer & Synchronization

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Open file for writing" | `open_file_for_write("/path/to/file.txt")` | Opens file handle for streaming write operations |
| "Write data chunk" | `write_chunk("file_handle", "data_chunk")` | Writes data chunk to open file handle |
| "Close file handle" | `close_file("file_handle")` | Closes file handle and finalizes write |

---

## Vector Database & RAG Operations

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Search vector database" | `query_vector_db("database_name", "search query")` | Performs semantic search in vector database |
| "Get database information" | `get_vector_db_info("database_name")` | Retrieves metadata and statistics about vector database |

---

## AI-Powered External Toolset Generation

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "Generate web scraping toolset" | `generate_toolset("web_scraper", "web_scraping", "scrape e-commerce sites")` | AI creates custom web scraping tools |
| "Create blockchain analysis tools" | `generate_toolset("crypto_analyzer", "cryptocurrency", "analyze DeFi protocols")` | AI generates blockchain analysis toolset |
| "Build ML pipeline toolset" | `generate_toolset("ml_pipeline", "machine_learning", "automate model training")` | AI creates machine learning workflow tools |
| "Generate custom domain toolset" | `generate_toolset("custom_tools", "domain", "description")` | AI adapts to ANY domain automatically |
| "List existing toolsets" | `list_existing_toolsets()` | Shows all generated external toolsets |
| "Remove specific toolset" | `remove_toolset("toolset_name")` | Removes specific external toolset |
| "Clear all external toolsets" | `clear_all_toolsets()` | Removes all external toolsets (use with caution) |
| "Get generation help" | `get_generation_help()` | Get help and examples for toolset generation |

**Key Features:**
- **AI-Powered**: Each toolset includes intelligent workflow guidance
- **Domain Agnostic**: Works for any field - web scraping, blockchain, ML, bioinformatics, etc.
- **Auto-Integration**: Generated toolsets automatically work with TodoList management
- **Smart Tools**: AI determines appropriate tools based on domain and description

---

## Service & Endpoint Management

| **Natural Language Request** | **Tool Call** | **Function Description** |
|------------------------------|---------------|--------------------------|
| "List all services" | `list_services()` | Shows all registered service endpoints |
| "Get image as base64" | `fetch_image_base64("image_url")` | Retrieves image and converts to base64 format |
| "Add new service" | `add_service("service_name", "http://localhost:8080")` | Registers new service endpoint |
| "Get service details" | `get_service("service_name")` | Retrieves information about specific service |
| "Check services status" | `services_ready()` | Checks if all registered services are available |

---

## Common Workflow Patterns

### Pattern 1: Complete Data Analysis Workflow
```
User: "I want to do a complete RNA-seq analysis"
CLI: add_todo("Complete RNA-seq analysis from FASTQ to results")
CLI: execute_current_task()  # Get guidance for first step
CLI: run_python("# preprocessing script")  # Execute analysis
CLI: mark_task_done("Preprocessing completed")  # Mark progress
```

### Pattern 2: ATAC-seq Analysis from Start
```
User: "Analyze ATAC-seq data in my folder"
CLI: /bio atac init  # Enter ATAC mode
CLI: /bio atac upstream ./data_folder  # Start upstream analysis
CLI: bio_atac_scan_folder("./data_folder")  # Scan contents
CLI: bio_atac_auto_detect_species("./data_folder")  # Detect species
CLI: bio_atac_setup_genome_resources("human", "hg38")  # Setup genome
CLI: bio_atac_auto_align_fastq("./data_folder")  # Align reads
```

### Pattern 3: Code Quality & Validation
```
User: "Check and fix my Python code"
CLI: validate_python_code("/path/script.py")  # Check syntax
CLI: check_code_style("/path/script.py")  # Quality check
CLI: edit_file("/path/script.py", line_number=10, new_content="fixed_line")  # Fix issues
CLI: validate_python_code("/path/script.py")  # Re-validate
```

### Pattern 4: Research & Documentation
```
User: "Research best practices for my analysis"
CLI: web_search("ATAC-seq analysis best practices 2024")  # Search web
CLI: google_search("quality control methods genomics")  # Google search
CLI: query_vector_db("research_docs", "QC metrics importance")  # Search docs
CLI: write_file("/notes/research.md", "findings")  # Save notes
```

### Pattern 5: Single-Cell RNA-seq Analysis
```
User: "Analyze single-cell RNA data with cell type annotation"
CLI: /bio scrna init  # Enter scRNA mode
CLI: /bio scrna analysis ./data.h5ad  # Load and analyze data
CLI: bio_scrna_run_quality_control("./data.h5ad")  # QC analysis
CLI: bio_scrna_run_preprocessing("./data.h5ad")  # Normalization
CLI: bio_scrna_run_cell_type_annotation("./data.h5ad")  # Cell annotation
```

### Pattern 6: Julia High-Performance Computing
```
User: "Solve differential equations with Julia"
CLI: ]using DifferentialEquations  # Load packages directly
CLI: ]prob = ODEProblem(f, u0, tspan)  # Define problem
CLI: ]sol = solve(prob)  # Solve equations
CLI: ]using Plots; plot(sol)  # Plot solution
```

### Pattern 7: Multi-Language REPL Development
```
User: "Quick analysis across languages"
CLI: !ls data/  # Check available data files
CLI: %import pandas as pd; df = pd.read_csv("data/data.csv")  # Load in Python
CLI: %df.describe()  # Python descriptive stats
CLI: >library(ggplot2)  # Switch to R for plotting
CLI: >data <- read.csv("data/data.csv")  # Load in R
CLI: >ggplot(data, aes(x=var1, y=var2)) + geom_point()  # R visualization
CLI: ]using DataFrames, CSV  # Switch to Julia for performance
CLI: ]df = CSV.read("data/data.csv", DataFrame)  # Load in Julia
CLI: ]mean(df.var1)  # Julia statistical computation
```

### Pattern 8: System Administration Workflow
```
User: "Check system status and manage processes"
CLI: !ps aux | grep python  # Check Python processes
CLI: !df -h  # Check disk usage
CLI: !free -m  # Check memory usage
CLI: !top -n 1 | head -20  # Check CPU usage
CLI: %import psutil; psutil.cpu_percent()  # Python system monitoring
CLI: !systemctl status nginx  # Check service status
```

## Pro Tips

1. **Start with TODO** - Use `add_todo()` for complex multi-step tasks
2. **Get Guidance** - Use `execute_current_task()` when unsure what to do next  
3. **Use Direct Execution** - **NEW**: Use `!`, `%`, `>`, `]` for instant code execution
4. **Quick Commands** - Use `!ls`, `%import pandas`, `>summary()`, `]println()` for fast operations
5. **Validate First** - Check syntax and files before running expensive operations
6. **Track Progress** - Always use `mark_task_done()` to track completion
7. **Use Sessions** - Create interpreter/shell sessions for complex workflows
8. **Combine Tools** - Chain multiple tools together for powerful workflows
9. **Search Smart** - Use `glob()` for files, `grep()` for content, `web_search()` for research
10. **Manage Quality** - Validate code with quality tools before execution
11. **Choose Right Language** - Python for general analysis, R for statistics/Seurat, Julia for HPC
12. **Switch Languages Seamlessly** - **NEW**: Mix `%python`, `>R`, `]julia`, `!shell` in same session
13. **Bio Mode Entry** - Use `/bio <tool> init` to enter focused analysis modes
14. **Generate Custom Tools** - Use `generate_toolset()` for domain-specific workflows
15. **Leverage AI Guidance** - Generated toolsets include intelligent workflow suggestions
16. **REPL-Style Development** - **NEW**: Use direct execution for interactive development and testing

This guide maps your natural language intentions to specific tool calls, making it easy to accomplish any task in the Pantheon CLI efficiently.