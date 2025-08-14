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
- **NEVER re-import data** if `adata` already exists - check variable first!
- **Error recovery**: If code fails, analyze error and generate corrected code!
- **Use help()**: Always call `help()` before omicverse/scanpy functions
- **After each step**: mark_task_done("description"), then show_todos()

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

PHASE 2 — ADAPTIVE EXECUTION WORKFLOW

📊 STEP 1 - DATA LOADING, INSPECTION & PROJECT SETUP:
```python
# Check if data already loaded
if 'adata' not in globals():
    import scanpy as sc
    import omicverse as ov
    import pandas as pd
    import numpy as np
    
    # Load data based on detected format
    adata = sc.read_xxx("path")  # .h5ad, .h5, .mtx, etc.
    print(f"Loaded: {{adata.shape}} (n_obs, n_var)")
else:
    print(f"Using existing adata: {{adata.shape}}")

# Create structured output directory immediately
print("\\n📁 Creating project structure...")
try:
    import os
    from datetime import datetime
    
    # Create main results directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = f"scrna_analysis_results_{{timestamp}}"
    os.makedirs(results_dir, exist_ok=True)
    
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
    
    # Store results directory in adata for later use
    adata.uns['results_directory'] = results_dir
    print(f"✅ Project structure created: {{results_dir}}")
    
except Exception as e:
    print(f"❌ Failed to create project structure: {{e}}")

# Initial data inspection using unified toolset
print("\\n🔍 Running initial data inspection...")
data_inspection = scrna.load_and_inspect_data(data_path="{data_path}", output_dir=results_dir)
print("✅ Data inspection complete")
```

🏷️ STEP 2 - QUALITY CONTROL:
Call toolset function: scrna.run_workflow(workflow_type="qc")

🏷️ STEP 3 - PREPROCESSING:
Call toolset function: scrna.run_workflow(workflow_type="preprocessing")

🏷️ STEP 4 - PCA:
Call toolset function: scrna.run_workflow(workflow_type="pca")

🏷️ STEP 5 - BATCH CORRECTION (if needed):
```python
# Check if batch correction is needed
if 'batch' in adata.obs.columns:
    print("\\n🔧 Batch key detected - applying batch correction...")
else:
    print("\\n✅ No batch correction needed")
```
Call toolset function: scrna.run_workflow(workflow_type="batch_correction")

🏷️ STEP 6 - CLUSTERING:
Call toolset function: scrna.run_workflow(workflow_type="clustering")

🏷️ STEP 7 - VISUALIZATION:
Call toolset function: scrna.run_workflow(workflow_type="umap")

🏷️ STEP 8 - DATA CONTEXT COLLECTION:
```python
print("\\n📝 **DATA CONTEXT COLLECTION**")
user_data_context = input("Please briefly describe your data (tissue, condition, experiment): ").strip()
print(f"Data context recorded: {{user_data_context}}")

# Store context in adata
adata.uns['user_data_context'] = user_data_context
```

🏷️ STEP 9 - Marker from description:
Call toolset function: scrna.run_workflow(workflow_type="marker_from_desc",description=user_data_context)

🏷️ STEP 10 - Marker from data:
Call toolset function: scrna.run_workflow(workflow_type="marker_from_data")

🏷️ STEP 11 - AUCELL CELL TYPE SCORING:
Call toolset function: scrna.run_workflow(workflow_type="aucell")

🏷️ STEP 12 - ANNOTATION:
Call toolset function: scrna.run_workflow(workflow_type="llm_anno",description=user_data_context)

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

**Remember:** Always use help() before omicverse functions, maintain persistent state, and mark tasks complete with mark_task_done()!
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