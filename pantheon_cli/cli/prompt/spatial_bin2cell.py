


from pathlib import Path
from typing import Optional

def generate_spatial_workflow_message(workflow_type: str) -> str:
    """Generate spatial workflow message using spatial toolset with omicverse"""
    bin2cell_message = f"""
🧬 Spatial Analysis Pipeline with omicverse Integration 

🔧 AVAILABLE SPATIAL ANALYSIS TOOL:
You have access to the Spatial_Bin2Cell_Analysis tool which provides workflow guidance.
Use it by calling: Spatial_Bin2Cell_Analysis(workflow_type="<type>")
Available workflow types: cellpose_he, expand_labels, cellpose_gex, salvage_secondary_labels, bin2cell

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

PHASE 0 — SETUP & VALIDATION
1) Data discovery: Use ls command to check folder contents, then proceed with file loading
2) Environment check will be done automatically within data loading step

PHASE 1 - TODO CREATION (ONCE ONLY)
Execute: current = show_todos()
IF current is EMPTY, create these todos ONCE:
1. "Check Python environment and load initial data"
2. "Run read_visium workflow"
3. "Run cellpose_he workflow"
4. "Run expand_labels workflow"
5. "Run cellpose_gex workflow"
6. "Run salvage_secondary_labels workflow"
7. "Run bin2cell workflow"

⚡ AUTOMATIC WORKFLOW MODE:
- Execute each todo task automatically without asking for confirmation
- After successful completion of any step, immediately call mark_task_done("description") and proceed to next
- Continue the workflow seamlessly until all tasks complete or user intervenes


PHASE 2 — ADAPTIVE EXECUTION WORKFLOW
⚠️ CRITICAL EXECUTION STRATEGY:
The Spatial_Bin2Cell_Analysis tool is a registered toolset function that returns guidance and example code.
To call it, you MUST use the tool calling mechanism:
- Tool name: Spatial_Bin2Cell_Analysis
- Required parameter: workflow_type (string) - one of: "cellpose_he", "expand_labels", "cellpose_gex", "salvage_secondary_labels", "bin2cell"

When you call Spatial_Bin2Cell_Analysis(), it returns guidance, explanations, and example Python code.
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
The path is the path of the data,you can use ls command to check the folder contents,then proceed with file loading.
if you don't know the path, you need to ask the user for the path and use ls command to explore the folder contents and structure.
after you have the path, you need to find some look like `"binned_outputs/square_002um/"` in the path,and then load the data.
the img should be stored in any location in the path, you need to find the correct img and ask the user to confirm the img path.
the img should end with btf or tiff, you need to find the correct img and ask the user to confirm the img path.


```python
# EFFICIENT DATA LOADING - Check memory first!
try:
    # Check if adata exists and is valid
    print(f"Using existing adata: {{adata.shape}} (n_obs, n_var)")
    data_already_loaded = True
except NameError:
    # Only load if not in memory
    print("Loading data for the first time...")
    
    # Check and install required packages
    import subprocess
    import sys
    print("Checking required packages...")
    required = ['scanpy', 'omicverse', 'pandas', 'numpy']
    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            print(f"Installing {{pkg}}...")
            subprocess.run([sys.executable, '-m', 'pip', 'install', pkg])
    
    import scanpy as sc
    import omicverse as ov
    import pandas as pd
    import numpy as np
    
    # Load data based on detected format
    adata = ov.space.read_visium_10x("path")  # .h5ad, .h5, .mtx, etc.
    print(f"Loaded adata: {{adata.shape}} (n_obs, n_var)")

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
        results_dir = f"spatial_analysis_results_{{timestamp}}"
        os.makedirs(results_dir, exist_ok=True)
    else:
        print(f"Using existing results directory: {{results_dir}}")
    
    # Create subdirectories for different analysis components
    subdirs = [
        "01_data_loading",
        "02_stardist",
        "03_bin2cell", 
    ]
    for subdir in subdirs:
        os.makedirs(os.path.join(results_dir, subdir), exist_ok=True)
    print(f"📂 Project structure created in: {{results_dir}}")
except Exception as e:
    print(f"❌ Error creating project structure: {{e}}")

```

🏷️ STEP 2 - cellpose_he:
Get cellpose_he guidance and adapt the code to your data:
All function using in omicverse, should use help() to check the function parameters!!!

First, call the spatial analysis tool to get guidance:
**Execute this tool call:** Spatial_Bin2Cell_Analysis workflow_type="cellpose_he"

Then analyze the returned guidance and implement adapted cellpose_he code based on your adata structure.
**CRITICAL**: Analyze cellpose_he results

🏷️ STEP 3 - expand_labels:
Get expand_labels guidance and adapt the code to your data:
All function using in omicverse, should use help() to check the function parameters!!!

First, call the spatial analysis tool to get guidance:
**Execute this tool call:** Spatial_Bin2Cell_Analysis workflow_type="expand_labels"

Then analyze the returned guidance and implement adapted expand_labels code based on your adata structure.
**CRITICAL**: Analyze expand_labels results

🏷️ STEP 4 - cellpose_gex:
Get cellpose_gex guidance and adapt the code to your data:
All function using in omicverse, should use help() to check the function parameters!!!

First, call the spatial analysis tool to get guidance:
**Execute this tool call:** Spatial_Bin2Cell_Analysis workflow_type="cellpose_gex"

Then analyze the returned guidance and implement adapted cellpose_gex code based on your adata structure.
**CRITICAL**: Analyze cellpose_gex results

🏷️ STEP 5 - salvage_secondary_labels:
Get salvage_secondary_labels guidance and adapt the code to your data:
All function using in omicverse, should use help() to check the function parameters!!!

First, call the spatial analysis tool to get guidance:
**Execute this tool call:** Spatial_Bin2Cell_Analysis workflow_type="salvage_secondary_labels"

Then analyze the returned guidance and implement adapted salvage_secondary_labels code based on your adata structure.
**CRITICAL**: Analyze salvage_secondary_labels results

🏷️ STEP 6 - bin2cell:
Get bin2cell guidance and adapt the code to your data:
All function using in omicverse, should use help() to check the function parameters!!!

First, call the spatial analysis tool to get guidance:
**Execute this tool call:** Spatial_Bin2Cell_Analysis workflow_type="bin2cell"

Then analyze the returned guidance and implement adapted bin2cell code based on your adata structure.
**CRITICAL**: Analyze bin2cell results

🔧 **AVAILABLE TOOLSET FUNCTIONS:**

**UNIFIED WORKFLOW ENGINE:**
- `Spatial_Bin2Cell_Analysis(workflow_type="cellpose_he")` - Cellpose HE with omicverse
- `Spatial_Bin2Cell_Analysis(workflow_type="expand_labels")` - Expand labels with omicverse
- `Spatial_Bin2Cell_Analysis(workflow_type="cellpose_gex")` - Cellpose GEX with omicverse
- `Spatial_Bin2Cell_Analysis(workflow_type="salvage_secondary_labels")` - Salvage secondary labels with omicverse
- `Spatial_Bin2Cell_Analysis(workflow_type="bin2cell")` - Bin2cell with omicverse

**EXECUTION STRATEGY:**
1. Load data and create project structure
2. Execute todos in sequence using appropriate workflow functions
3. Use modular functions for specialized analysis steps
4. Leverage omicverse integration with scanpy fallbacks
5. Comprehensive result saving and reporting

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
    if workflow_type == "bin2cell":
        return bin2cell_message
    else:
        return None