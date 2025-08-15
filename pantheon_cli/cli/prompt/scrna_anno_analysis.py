"""Single-cell RNA-seq analysis mode handler with omicverse integration"""

from pathlib import Path
from typing import Optional

def generate_scrna_analysis_message(folder_path: Optional[str] = None) -> str:
    """Generate scRNA-seq analysis message using scrna toolset with omicverse"""
    
    if folder_path:
        data_path = Path(folder_path).resolve()
        
        # Determine if it's a file or folder
        if data_path.is_file():
            target_description = f"Target data file: {data_path}"
            path_instruction = f'Always use the provided data_path: "{data_path}" for data loading and analysis.'
        else:
            target_description = f"Target folder: {data_path}"
            path_instruction = f'Always use the provided folder_path: "{data_path}" to scan for scRNA-seq data files.'
        
        message = f"""
🧬 Single-cell RNA-seq Analysis Pipeline with omicverse Integration
{target_description}
⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
- **PERSISTENT STATE**: Python interpreter maintains ALL variables across calls! 
- **MEMORY OPTIMIZATION**: Variables persist! NEVER re-read or re-import data that already exists in memory!
- **SMART VARIABLE CHECKING**: Use `try/except` or `'var' in globals()` to check existence - NO redundant file I/O!
- **EFFICIENCY FIRST**: 
  - Check if adata exists before loading: `if 'adata' not in globals()`
  - Use existing results: `if 'pca_result' in adata.obsm`
  - Reuse computed values: `if 'marker_genes' in locals()`
- **ERROR RECOVERY**: If code fails, analyze error and fix - don't reload everything!
- **NO REPETITION**: Each import/load/compute happens ONCE per session unless explicitly needed
- **After each step**: mark_task_done("description"), then show_todos()
- **AUTOMATIC EXECUTION**: Proceed automatically without confirmations; log warnings when needed.

{path_instruction}

PHASE 0 — SETUP & VALIDATION
1) Environment check: scrna.check_dependencies()
2) Data discovery: scrna.scan_folder() if folder, or proceed with file

PHASE 1 — TODO CREATION (ONCE ONLY)
Execute: current = show_todos()
IF current is EMPTY, create these todos ONCE:
1. "Check Python environment and load initial data"
2. "Inspect data structure and determine processing pipeline"  
3. "Apply quality control with omicverse.pp.qc"
4. "Perform preprocessing with omicverse.pp.preprocess"
5. "Compute PCA with omicverse.pp.pca"
6. "Apply batch correction if needed"
7. "Run clustering analysis"
8. "Ask user for data context (tissue/condition)"
9. "Generate context-specific cell types and markers from description"
10. "Find cluster-specific marker genes from data"
11. "Calculate AUCell scores for cell type markers"
12. "Annotate cell type with LLM"
13. "Conduct downstream analysis"
14. "Generate analysis report"

⚡ AUTOMATIC WORKFLOW MODE:
- Execute each todo task automatically without asking for confirmation
- After successful completion of any step, immediately call mark_task_done("description") and proceed to next
- Continue the workflow seamlessly until all tasks complete or user intervenes

PHASE 2 — ADAPTIVE EXECUTION WORKFLOW

⚠️ CRITICAL EXECUTION STRATEGY:
When you call scrna.run_workflow(), it returns guidance, explanations, and example Python code using toolset function run_python_code.
You MUST:
1. **Read and analyze** the entire returned content carefully
2. **Understand the logic** and methodology described
3. **Adapt the provided code** to your current data situation (adata shape, available columns, etc.)
4. **Modify parameters** based on your actual data characteristics
5. **Execute the adapted code** - NOT the original code directly
6. **Handle errors** by adjusting code based on the guidance provided

🧠 **RESULT ANALYSIS REQUIREMENT:**
After executing any code:
1. **Analyze the output** - Don't just print and move on!
2. **Interpret the results** - What do the numbers, plots, and warnings mean?
3. **Check for issues** - Are there data quality problems or unexpected patterns?
4. **Make decisions** - Should parameters be adjusted based on what you observed?
5. **Document findings** - Save key insights to results directory
6. **Proceed intelligently** - Use results to inform next steps

The returned content serves as GUIDANCE and TEMPLATES, not direct execution scripts.

📊 STEP 1 - DATA LOADING, INSPECTION & PROJECT SETUP:
```python
# EFFICIENT DATA LOADING - Check memory first!
try:
    # Check if adata exists and is valid
    print(f"Using existing adata: {{adata.shape}} (n_obs, n_var)")
    data_already_loaded = True
except NameError:
    # Only load if not in memory
    print("Loading data for the first time...")
    import scanpy as sc
    import omicverse as ov
    import pandas as pd
    import numpy as np
    
    # Load data based on detected format
    adata = sc.read_xxx("path")  # .h5ad, .h5, .mtx, etc.
    print(f"Loaded: {{adata.shape}} (n_obs, n_var)")
    data_already_loaded = False

# Only import libraries once
if 'sc' not in globals():
    import scanpy as sc
    import omicverse as ov
    import pandas as pd
    import numpy as np
    print("Libraries imported")

# Create structured output directory (only if not exists)
print("\\n📁 Setting up project structure...")
try:
    # Check if we already have a results directory
    if 'results_dir' not in globals():
        import os
        from datetime import datetime
        
        # Create main results directory with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = f"scrna_analysis_results_{{timestamp}}"
        os.makedirs(results_dir, exist_ok=True)
    else:
        print(f"Using existing results directory: {{results_dir}}")
    
    # Create subdirectories for different analysis components
    subdirs = [
        "01_data_loading",
        "02_quality_control",
        "03_preprocessing", 
        "04_dimensionality_reduction",
        "05_batch_correction",
        "06_clustering",
        "07_cell_type_annotation",
        "08_visualization",
        "09_downstream_analysis",
        "10_reports",
        "logs"
    ]
    
    for subdir in subdirs:
        os.makedirs(os.path.join(results_dir, subdir), exist_ok=True)
    
    # Store results directory in adata for later use (if adata exists)
    if 'adata' in globals():
        adata.uns['results_directory'] = results_dir
        print(f"✅ Project structure ready: {{results_dir}}")
    else:
        print(f"✅ Project structure created: {{results_dir}} (will link to adata after loading)")
    
except Exception as e:
    print(f"❌ Failed to create project structure: {{e}}")

# Initial data inspection using unified toolset
print("\\n🔍 Running initial data inspection...")
data_inspection = scrna.load_and_inspect_data(data_path="{data_path}", output_dir=results_dir)
print("✅ Data inspection complete")
```

🏷️ STEP 2 - QUALITY CONTROL:
Get QC guidance and adapt the code to your data:
scrna.run_workflow(workflow_type="qc")
Then analyze the returned guidance and implement adapted QC code based on your adata structure.
**CRITICAL**: Analyze QC results - cell counts, gene expression distributions, mitochondrial percentages. Interpret plots and decide on filtering thresholds.

🏷️ STEP 3 - PREPROCESSING:
Get preprocessing guidance and adapt the code to your data:
scrna.run_workflow(workflow_type="preprocessing")
Then analyze the returned guidance and implement adapted preprocessing code based on your adata characteristics.
**CRITICAL**: Examine normalization results, highly variable genes selection. Check if the data distribution looks appropriate.

🏷️ STEP 4 - PCA:
Get PCA guidance and adapt the code to your data:
scrna.run_workflow(workflow_type="pca")
Then analyze the returned guidance and implement adapted PCA code based on your adata dimensions.
**CRITICAL**: Analyze PCA results - variance explained, elbow plots. Determine optimal number of components to use.

🏷️ STEP 5 - BATCH CORRECTION (if needed):
```python
# EFFICIENT BATCH CHECK - Only compute if not done before
if 'real_batch_keys' not in globals():
    # Check if batch correction is needed - look for REAL batch keys
    # Real batch keys are typically: 'batch', 'sample', 'donor', 'experiment', 'plate', 'condition'
    # NOT QC metrics like 'passing_mt', 'passing_ngenes', 'n_genes', 'total_counts', etc.
    
    real_batch_keys = []
    potential_batch_keys = ['batch', 'sample', 'donor', 'experiment', 'plate', 'condition', 'library_id']
    
    for key in potential_batch_keys:
        if key in adata.obs.columns:
            # Check if it has multiple unique values and is categorical
            unique_vals = adata.obs[key].nunique()
            if unique_vals > 1 and unique_vals < adata.n_obs * 0.5:  # Not too many unique values
                real_batch_keys.append(key)
                print(f"Found real batch key: {{key}} with {{unique_vals}} unique values")
else:
    print(f"Using previously identified batch keys: {{real_batch_keys}}")

if real_batch_keys:
    print(f"\\n🔧 Real batch keys detected: {{real_batch_keys}}")
    print("Proceeding with batch correction...")
    # Only proceed if real batch keys exist
else:
    print("\\n✅ No real batch keys found - skipping batch correction")
    print("Note: QC metrics like 'passing_mt', 'passing_ngenes' are NOT batch keys")
```
Only if real_batch_keys were found, get guidance and adapt the code:
If real_batch_keys: scrna.run_workflow(workflow_type="batch_correction")
Then implement adapted batch correction code based on your specific batch keys.
**CRITICAL**: Only apply batch correction if there are REAL batch effects, not QC filtering metrics.

🏷️ STEP 6 - CLUSTERING:
Get clustering guidance and adapt the code to your data:
scrna.run_workflow(workflow_type="clustering")
Then analyze the returned guidance and implement adapted clustering code based on your adata.

🏷️ STEP 7 - VISUALIZATION:
Get UMAP guidance and adapt the code to your data:
scrna.run_workflow(workflow_type="umap")
Then analyze the returned guidance and implement adapted visualization code.

🏷️ STEP 8 - DATA CONTEXT COLLECTION:
```python
# EFFICIENT CONTEXT COLLECTION - Check if already collected
if 'user_data_context' not in globals():
    print("\\n📝 **DATA CONTEXT COLLECTION**")
    user_data_context = input("Please briefly describe your data (tissue, condition, experiment): ").strip()
    print(f"Data context recorded: {{user_data_context}}")
    
    # Store context in adata
    adata.uns['user_data_context'] = user_data_context
else:
    print(f"Using existing data context: {{user_data_context}}")
```

🏷️ STEP 9 - Marker from description:
Get marker generation guidance based on data context:
scrna.run_workflow(workflow_type="marker_from_desc", description=user_data_context)
Then adapt and implement marker generation code based on the returned guidance and your tissue context.

🏷️ STEP 10 - Marker from data:
Get data-driven marker analysis guidance:
scrna.run_workflow(workflow_type="marker_from_data")
Then adapt and implement marker analysis code based on your actual cluster structure.
**CRITICAL**: Evaluate marker genes - fold changes, p-values, specificity. Identify the most discriminative markers per cluster.

🏷️ STEP 11 - AUCELL CELL TYPE SCORING:
Get AUCell scoring guidance and methodology:
scrna.run_workflow(workflow_type="aucell")
Then adapt and implement AUCell scoring code based on your marker genes and cell clusters.
**CRITICAL**: Examine AUCell scores distribution, thresholds, and how well they separate cell types. Validate scoring results.

🏷️ STEP 12 - ANNOTATION:
Get LLM-powered annotation guidance and workflow:
scrna.run_workflow(workflow_type="llm_anno", description=user_data_context)
Then adapt the annotation workflow based on your specific clustering results and evidence.
**CRITICAL**: Carefully review LLM annotations against marker evidence. Verify biological plausibility of assigned cell types.

🏷️ STEP 13 - DOWNSTREAM ANALYSIS:
```python
print("\\n🧬 Conducting downstream analysis...")

# Generate comprehensive analysis report
report_generation = scrna.generate_report(
    data_path="{data_path}",
    output_dir=results_dir,
    include_qc=True,
    include_clustering=True,
    include_annotation=True
)

print("✅ Downstream analysis and reporting complete")
```


🔧 **AVAILABLE TOOLSET FUNCTIONS:**

**UNIFIED WORKFLOW ENGINE:**
- `scrna.run_workflow(workflow_type="qc")` - Quality control with omicverse
- `scrna.run_workflow(workflow_type="preprocessing")` - Preprocessing with omicverse
- `scrna.run_workflow(workflow_type="pca")` - PCA with omicverse
- `scrna.run_workflow(workflow_type="clustering")` - Clustering analysis
- `scrna.run_workflow(workflow_type="umap")` - Calculate UMAP
- `scrna.run_workflow(workflow_type="aucell")` - AUCell scoring

**EXECUTION STRATEGY:**
1. Load data and create project structure
2. Execute todos in sequence using appropriate workflow functions
3. Use modular functions for specialized analysis steps
4. Leverage omicverse integration with scanpy fallbacks
5. Interactive LLM annotation for expert cell type assignment
6. Comprehensive result saving and reporting

**EFFICIENCY PRINCIPLES:**
1. **CHECK BEFORE COMPUTE**: Always check if variables/results exist before recomputing
2. **USE TRY/EXCEPT**: Gracefully handle missing variables without re-reading files
3. **MEMORY-FIRST**: Trust the persistent Python interpreter - no redundant I/O
4. **SMART RECOVERY**: Fix errors in-place, don't restart entire analysis
5. **INCREMENTAL PROGRESS**: Each step builds on previous results

Example patterns:
```python
# Good - Check memory first
try:
    print(f"Using existing adata: {{adata.shape}}")
except NameError:
    adata = sc.read_h5ad(path)

# Good - Check computed results
if 'X_pca' not in adata.obsm:
    sc.tl.pca(adata)
else:
    print("PCA already computed")

# Bad - Redundant file I/O
adata = sc.read_h5ad(path)  # Don't do this if adata exists!
```

**Remember:** Maintain persistent state, avoid redundant operations, and mark tasks complete with mark_task_done()!
"""
        
    else:
        message = """
I need help with single-cell RNA-seq analysis using your specialized toolsets.

You have access to comprehensive scRNA-seq and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE scRNA-seq TOOLSET:

**UNIFIED WORKFLOW ENGINE:**
- `scrna.run_workflow(workflow_type="qc")` - Quality control with omicverse
- `scrna.run_workflow(workflow_type="preprocessing")` - Preprocessing with omicverse
- `scrna.run_workflow(workflow_type="pca")` - PCA with omicverse
- `scrna.run_workflow(workflow_type="clustering")` - Clustering analysis
- `scrna.run_workflow(workflow_type="umap")` - Calculate UMAP
- `scrna.run_workflow(workflow_type="aucell")` - AUCell scoring

**EXECUTION STRATEGY:**
1. Load data and create project structure
2. Execute todos in sequence using appropriate workflow functions
3. Use modular functions for specialized analysis steps
4. Leverage omicverse integration with scanpy fallbacks
5. Interactive LLM annotation for expert cell type assignment
6. Comprehensive result saving and reporting

**GUIDANCE:**
- scrna.suggest_next_step() - Smart recommendations

Please start by adding a todo for your scRNA-seq analysis task, then use the appropriate scRNA tools!"""
    
    return message