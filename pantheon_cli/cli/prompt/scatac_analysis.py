"""Single-cell ATAC-seq analysis mode handler"""

from pathlib import Path
from typing import Optional

def generate_scatac_analysis_message(folder_path: Optional[str] = None) -> str:
    """Generate scATAC-seq analysis message using scATAC toolset"""
    
    if folder_path:
        folder_path = Path(folder_path).resolve()
        
        message = f"""
🧬 Single-cell ATAC-seq Analysis Pipeline — Workflow-Based Architecture
Target folder: {folder_path}

You have access to the scATAC-seq toolset with workflow-based commands and TodoList management.

GLOBAL RULES
- Always use the provided folder_path: "{folder_path}" in all phases.
- Idempotent behavior: NEVER create duplicate todos. Only create if the list is EMPTY.
- Do not ask the user for confirmations; proceed automatically and log warnings when needed.
- After each concrete tool completes successfully, call mark_task_done("what was completed"), then show_todos().

PHASE 0 — DEPENDENCIES & SETUP
1) Use workflow commands for setup and dependency checking:
   - scatac.ScATAC_Upstream("init") - Initialize project structure
   - scatac.ScATAC_Upstream("check_dependencies") - Check tool availability
   - scatac.ScATAC_Upstream("install_cellranger") - Install cellranger-atac (if needed)
   - scatac.ScATAC_Upstream("setup_reference") - Setup genome references

PHASE 1 — TODO CREATION (STRICT DE-DUP)
Mandatory order:
  a) current = show_todos()
  b) Analyze folder structure and FASTQ files in "{folder_path}"
Creation rule (single condition):
  • If current is EMPTY → create ONCE the following todos:
      0. "Initialize scATAC-seq project structure"
      1. "Check and install scATAC-seq dependencies"
      2. "Setup genome references and cellranger-atac"
      3. "Scan and validate input FASTQ files"
      4. "Run cellranger-atac count for sample processing"
      5. "Load cellranger outputs into analysis format"
      6. "Quality control and cell filtering"
      7. "Compute LSI embeddings and UMAP"
      8. "Find cell clusters using Leiden algorithm"
      9. "Annotate peaks with genomic features"
      10. "Differential accessibility analysis"
      11. "Motif analysis and transcription factors"
      12. "Generate comprehensive analysis report"
  • Else → DO NOT create anything. Work with the existing todos.

PHASE 2 — EXECUTE WITH TODO TRACKING (LOOP)

⚠️ CRITICAL EXECUTION STRATEGY:
When you call scatac.ScATAC_Upstream() or scatac.ScATAC_Analysis(), they return bash command templates.
You MUST:
1. **Read and analyze** the entire returned bash command content carefully
2. **Understand the logic** and methodology described
3. **Adapt the provided commands** to your current data situation (file paths, sample names, etc.)
4. **Execute the adapted commands** using bash tool - NOT the original commands directly
5. **Handle errors** by adjusting commands based on the guidance provided
6. **Analyze results** - Check output files, logs, and success/failure status

🧠 **RESULT ANALYSIS REQUIREMENT:**
After executing any bash commands:
1. **Analyze the output** - Don't just run and move on!
2. **Check generated files** - Verify expected output files were created
3. **Examine logs and errors** - Look for warnings, failures, or quality issues
4. **Validate results** - Are the results biologically reasonable?
5. **Make decisions** - Should parameters be adjusted based on what you observed?
6. **Document findings** - Note any issues or important observations

The returned bash commands serve as TEMPLATES and GUIDANCE, adapt them to your specific data.

For each current task:
  1) hint = execute_current_task()   # obtain guidance for the next action
  2) Get bash command templates using appropriate scATAC workflow:
     
     UPSTREAM WORKFLOWS (use scatac.ScATAC_Upstream(workflow_type)):
     - "init" - Initialize project structure
     - "check_dependencies" - Check tool dependencies  
     - "install_cellranger" - Install cellranger-atac
     - "setup_reference" - Setup genome references
     - "scan_folder" - Scan and validate input data
     - "run_count" - Run cellranger-atac count
     - "test_functionality" - Test cellranger-atac installation
     
     DOWNSTREAM WORKFLOWS (use scatac.ScATAC_Analysis(workflow_type)):
     - "load_cellranger_data" - Load cellranger outputs
     - "quality_control" - Cell and peak filtering
     - "compute_embeddings" - LSI and UMAP embeddings
     - "find_clusters" - Leiden clustering
     - "annotate_peaks" - Peak annotation with genes
     - "differential_accessibility" - Find marker peaks
     - "motif_analysis" - Transcription factor analysis
     - "generate_report" - Comprehensive HTML/PDF report

  3) **EXECUTE** the adapted bash commands and analyze results
  4) **VERIFY** success by checking output files and logs
  5) mark_task_done("brief, precise description of the completed step")
  6) show_todos()
Repeat until all todos are completed.

PHASE 3 — ADAPTIVE TODO REFINEMENT
- If dependencies missing → add_todo("Install missing scATAC-seq tools")
- If cellranger-atac fails → add_todo("Fix cellranger-atac installation or parameters")
- If quality issues found → add_todo("Address data quality issues")
- If additional analysis needed → add_todo("Additional analysis task")

EXECUTION STRATEGY (MUST FOLLOW THIS ORDER)
  1) scatac.ScATAC_Upstream("init") → Initialize project
  2) scatac.ScATAC_Upstream("check_dependencies") → Check tools
  3) scatac.ScATAC_Upstream("setup_reference") → Setup genome
  4) show_todos()
  5) Analyze folder and create todos if empty
  6) Loop Phase 2 until all done; refine with Phase 3 when needed

📊 EXECUTION EXAMPLES:

🏷️ STEP 1 - PROJECT INITIALIZATION:
```bash
# Get initialization commands
init_commands = scatac.ScATAC_Upstream("init")
# Analyze and adapt the commands to your project
# Execute: mkdir -p project_structure, create config files, etc.
```

🏷️ STEP 2 - DEPENDENCY CHECK:
```bash
# Get dependency check commands  
dep_commands = scatac.ScATAC_Upstream("check_dependencies")
# Execute: which cellranger-atac, which python3, etc.
# Install missing tools if needed
```

🏷️ STEP 3 - CELLRANGER-ATAC COUNT:
```bash
# Get cellranger count template commands
count_commands = scatac.ScATAC_Upstream("run_count")
# Adapt to your actual FASTQ files and reference
# Execute: cellranger-atac count --id=sample --fastqs=path --reference=ref
# Check results: ls sample/outs/
# Analyze: Look for web_summary.html, fragments.tsv.gz, etc.
```

🏷️ STEP 4 - QUALITY CONTROL:
```bash
# Get QC template commands
qc_commands = scatac.ScATAC_Analysis("quality_control")
# Adapt with your cellranger output paths
# Execute: Python scanpy QC analysis
# Check results: analysis/qc/scatac_filtered.h5ad
# Analyze: cell retention, peak quality, filtering effectiveness
```

🏷️ STEP 5 - CLUSTERING:
```bash
# Get clustering template commands
cluster_commands = scatac.ScATAC_Analysis("find_clusters")
# Adapt with your QC'd data file
# Execute: Python scanpy clustering analysis
# Check results: analysis/clustering/scatac_clustered.h5ad
# Analyze: number of clusters, cluster quality, biological meaning
```

**CRITICAL SUCCESS PATTERNS:**
```bash
# Good - Check files before proceeding
ls cellranger_output/sample/outs/  # Verify cellranger success
ls analysis/qc/scatac_filtered.h5ad  # Verify QC success

# Good - Capture and analyze results
python3 -c "import scanpy as sc; adata = sc.read_h5ad('file.h5ad'); print(f'Shape: {{adata.shape}}')"

# Good - Error handling
if [ ! -f "analysis/clustering/scatac_clustered.h5ad" ]; then
    echo "ERROR: Clustering failed"
    exit 1
fi
```

**COMMON scATAC-seq WORKFLOW:**
1. **Raw Data** → cellranger-atac count → cellranger outputs
2. **cellranger outputs** → scanpy loading → AnnData object
3. **Raw counts** → Quality control → Filtered data
4. **Filtered data** → LSI embeddings → Dimensionality reduction
5. **Embeddings** → Clustering → Cell populations
6. **Clusters** → Differential peaks → Marker identification
7. **Marker peaks** → Motif analysis → Regulatory insights
8. **All results** → Report generation → Final summary

BEGIN NOW:
- Execute PHASE 0 → PHASE 1 → PHASE 2 loop.
- **ACTUALLY EXECUTE** the bash commands after adapting them
- **ANALYZE RESULTS** after each step
- **VERIFY SUCCESS** by checking output files
- Report any errors or quality issues encountered
"""
        
    else:
        message = """
I need help with single-cell ATAC-seq analysis using your specialized workflow-based toolsets.

You have access to comprehensive scATAC-seq and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE scATAC-seq WORKFLOW TOOLSET:

UPSTREAM WORKFLOWS (use scatac.ScATAC_Upstream(workflow_type)):
These return bash command templates for execution:

PROJECT SETUP:
- "init" - Initialize scATAC-seq project structure
- "check_dependencies" - Check tool availability (cellranger-atac, Python, R)
- "install_cellranger" - Install cellranger-atac v2.2.0
- "setup_reference" - Download and setup genome references (human/mouse)

DATA PROCESSING:
- "scan_folder" - Scan and validate input FASTQ files
- "run_count" - Run cellranger-atac count for sample processing
- "test_functionality" - Test cellranger-atac installation

DOWNSTREAM WORKFLOWS (use scatac.ScATAC_Analysis(workflow_type)):
These return bash command templates for execution:

DATA LOADING & QC:
- "load_cellranger_data" - Load cellranger outputs into scanpy
- "quality_control" - Cell and peak filtering with QC metrics

ANALYSIS & VISUALIZATION:
- "compute_embeddings" - LSI dimensionality reduction and UMAP
- "find_clusters" - Leiden clustering for cell populations
- "annotate_peaks" - Peak annotation with genomic features
- "differential_accessibility" - Find cluster-specific marker peaks
- "motif_analysis" - Transcription factor binding motif analysis
- "generate_report" - Comprehensive HTML analysis report

WORKFLOW USAGE:
All workflows return executable bash command templates. You must:
1. **Call the workflow function** to get command templates
2. **Adapt the commands** to your specific file paths and parameters  
3. **Execute the adapted commands** using the bash tool
4. **Analyze the results** and verify success

EXECUTION EXAMPLES:
```bash
# Step 1: Get command template
init_commands = scatac.ScATAC_Upstream("init")

# Step 2: Adapt and execute
# Template: mkdir -p scatac_analysis/{raw_data,references,...}
# Adapted: mkdir -p my_scatac_project/{raw_data,references,...}

# Step 3: Verify results  
ls my_scatac_project/  # Check project structure was created
```

⚠️ **CRITICAL**: Don't just call the workflow functions - you must EXECUTE the returned bash commands!

**TYPICAL scATAC-seq ANALYSIS WORKFLOW:**
1. Project initialization and setup
2. cellranger-atac count (FASTQ → counts)
3. Data loading and quality control
4. Dimensionality reduction (LSI/UMAP)
5. Cell clustering and annotation
6. Peak annotation and differential analysis
7. Motif analysis and regulatory insights
8. Report generation

Please start by adding a todo for your scATAC-seq analysis task, then use the workflow commands and EXECUTE them!"""
    
    return message