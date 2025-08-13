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
🧬 Single-cell RNA-seq Analysis Pipeline — Persistent Python Environment with omicverse Integration
{target_description}

⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
- **PERSISTENT STATE**: Python interpreter maintains ALL variables across calls! 
- **NEVER re-import data** if `adata` already exists - check variable first!
- **ALWAYS use help()** before calling ANY omicverse/scanpy function!
- **Error recovery**: If code fails, analyze error and generate corrected code!

GLOBAL EXECUTION RULES:
- {path_instruction}
- Use persistent Python variables - check existence before recreating
- After each step: mark_task_done("specific completion description"), then show_todos()
- MANDATORY: Use help(function_name) before every omicverse/scanpy call ONLY

🔧 PERSISTENT PYTHON WORKFLOW:

**STATE MANAGEMENT PATTERN (USE FOR EVERY STEP):**
```python
# 1. Check if variable exists
if 'adata' in locals() or 'adata' in globals():
    print(f"Using existing adata: {{adata.shape}}")
else:
    # Only load if not exists
    print("Loading data...")
    adata = sc.read_xxx(...)

# 2. Check current processing state
print("Current adata state:")
print(f"- Layers: {{list(adata.layers.keys())}}")
print(f"- .obsm keys: {{list(adata.obsm.keys())}}")
print(f"- .var columns: {{list(adata.var.columns)}}")

# 3. Use help() ONLY before omicverse/scanpy function calls
help(omicverse.pp.function_name)  # ONLY for omicverse functions
# help(scanpy.tl.function_name)   # ONLY for scanpy functions

# 4. Apply function with error handling
try:
    result = omicverse.pp.function_name(adata, param1=value1)
    print("✅ Step completed successfully")
except Exception as e:
    print(f"❌ Error occurred: {{e}}")
    # Analyze error and retry with corrected parameters
```

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
8. "Perform cell type annotation"
9. "Conduct downstream analysis"
10. "Generate analysis report"

PHASE 2 — ADAPTIVE EXECUTION WORKFLOW

📊 STEP 1 - DATA LOADING & INSPECTION:
```python
# Check if data already loaded
if 'adata' not in locals() and 'adata' not in globals():
    import scanpy as sc
    import omicverse as ov
    import pandas as pd
    import numpy as np
    
    # Load data based on detected format
    # [load appropriate format: .h5ad, .h5, .mtx, etc.]
    adata = sc.read_xxx("path")
    print(f"Loaded: {{adata.shape}} (n_obs, n_var)")
else:
    print(f"Using existing adata: {{adata.shape}}")

# Always inspect current state
print("\\n🔍 Current Data State:")
print(f"- Observations: {{adata.n_obs}} cells")
print(f"- Variables: {{adata.n_vars}} genes")
print(f"- Layers: {{list(adata.layers.keys())}}")
print(f"- Embeddings: {{list(adata.obsm.keys())}}")
print(f"- Metadata columns: {{list(adata.obs.columns)}}")
print(f"- Data type: {{type(adata.X)}}, Max value: {{adata.X.max()}}")
```

🔬 STEP 2 - QUALITY CONTROL (CONDITIONAL):
```python
# Check if QC already done
if 'pct_counts_mt' not in adata.obs.columns:
    print("\\n📊 Running Quality Control...")
    
    # MANDATORY: Check help first
    help(ov.pp.qc)
    
    try:
        # Apply QC with parameters from help()
        ov.pp.qc(adata, 
                target_sum=1e4,
                mt_gene_regex="^MT-|^mt-", 
                rp_gene_regex="^RP[SL]|^Rp[sl]")
        print("✅ QC completed successfully")
    except Exception as e:
        print(f"❌ QC failed: {{e}}")
        # Retry with different parameters based on error
        
else:
    print("✅ QC already completed - skipping")
```

🧬 STEP 3 - PREPROCESSING (CONDITIONAL):
```python
# Check if preprocessing needed
needs_preprocessing = (
    'highly_variable' not in adata.var.columns or 
    'scaled' not in adata.layers or
    adata.X.max() > 50  # Raw counts detected
)

if needs_preprocessing:
    print("\\n🧬 Running Preprocessing...")
    
    # MANDATORY: Check help first  
    help(ov.pp.preprocess)
    
    try:
        adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', n_top_genes=2000)
        print("✅ Preprocessing completed")
    except Exception as e:
        print(f"❌ Preprocessing failed: {{e}}")
        # Analyze error and retry with corrected approach
else:
    print("✅ Data already preprocessed - skipping")
```

🔢 STEP 4 - PCA (CONDITIONAL):
```python
# Check if PCA needed
if 'X_pca' not in adata.obsm.keys():
    print("\\n🔢 Computing PCA...")
    
    # MANDATORY: Check help first
    help(ov.pp.pca)
    
    try:
        adata = ov.pp.pca(adata, n_comps=50)
        print("✅ PCA completed")
    except Exception as e:
        print(f"❌ PCA failed: {{e}}")
        # Error recovery logic
else:
    print("✅ PCA already computed - skipping")
```

🔗 STEP 5 - BATCH CORRECTION (CONDITIONAL):
```python
# Check if batch correction needed and possible
has_batch = 'batch' or 'sample' in adata.obs.columns
has_corrected = any('scVI' in k or 'harmony' in k for k in adata.obsm.keys())

if has_batch and not has_corrected:
    print("\\n🔗 Applying Batch Correction...")
    
    help(ov.single.batch_correction)
    
    try:
        ov.single.batch_correction(adata, batch_key='batch', methods=['scVI'])
        print("✅ Batch correction completed")
    except Exception as e:
        print(f"❌ Batch correction failed: {{e}}")
else:
    print("✅ Batch correction not needed or already done")
```

🎯 STEP 6 - CLUSTERING (CONDITIONAL):
```python
# Check if clustering needed
if 'leiden' not in adata.obs.columns:
    print("\\n🎯 Running Clustering...")
    
    help(ov.utils.clusters)
    
    try:
        ov.utils.clusters(adata, method='leiden', resolution=0.5)
        print("✅ Clustering completed")
    except Exception as e:
        print(f"❌ Clustering failed: {{e}}")
else:
    print("✅ Clustering already done - skipping")
```

🏷️ STEP 7 - CELL TYPE ANNOTATION:
```python
print("\\n🏷️ Cell Type Annotation...")

# ONLY use help() for omicverse/scanpy functions:
help(ov.single.CellOntologyMapper)  # omicverse function
help(sc.tl.rank_genes_groups)       # scanpy function
help(ov.single.get_celltype_marker) # omicverse function

# Regular Python functions don't need help():
# print(), len(), list(), dict(), etc. - no help() needed

# [Implement annotation workflow based on help() output]
```

📈 STEP 8 - DOWNSTREAM ANALYSIS:
```python
print("\\n📈 Downstream Analysis with pertpy...")

# Use help() for pertpy functions when needed:
# help(pertpy.tl.some_function)  # Only for pertpy-specific functions
# help(scanpy.tl.some_function)  # Only for scanpy functions

# Regular Python/pandas operations don't need help():
# adata.obs['new_column'] = values  # No help() needed
# df.groupby('column').mean()       # No help() needed

# [Implement pertpy-based analysis]
```

💾 CRITICAL ERROR RECOVERY PATTERN:
```python
# For ANY failed operation, use this pattern:
try:
    # Original function call
    result = some_function(adata, param=value)
except Exception as e:
    print(f"❌ Error: {{e}}")
    print("🔧 Analyzing error and retrying...")
    
    # Analyze the error message
    if "parameter" in str(e):
        # Adjust parameters and retry
        result = some_function(adata, different_param=new_value)
    elif "data format" in str(e):
        # Convert data format and retry
        adata.X = adata.X.astype('float32')
        result = some_function(adata, param=value)
    # Add more error handling as needed
```

🚀 EXECUTION ORDER:
1. Setup environment
2. show_todos() and create if empty
3. For each todo: execute_current_task() → run step → mark_task_done() → show_todos()
4. Use persistent variables throughout - NEVER reload data unnecessarily
5. ALWAYS use help() before omicverse/scanpy/pertpy function calls ONLY
6. Implement error recovery for failed operations

BEGIN EXECUTION NOW:
Start with environment setup and persistent data loading workflow.
"""
        
    else:
        message = """
I need help with single-cell RNA-seq analysis using your specialized toolsets with omicverse integration.

You have access to comprehensive scRNA-seq and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE scRNA-seq TOOLSET (OMICVERSE INTEGRATION):

ENVIRONMENT & SETUP:
- scrna.check_dependencies() - Verify Python environment and packages
- scrna.install_missing_packages() - Install missing omicverse, scanpy, pertpy packages
- scrna.scan_folder() - Comprehensive scRNA-seq data analysis
- scrna.init() - Create scRNA project structure

DATA LOADING & INSPECTION:
- scrna.load_and_inspect_data() - Load and comprehensively analyze data structure
  * Examine adata.obs, adata.var, adata.obsm for existing analysis
  * Check QC metrics, cell type annotations, batch information
  * Detect data type and preprocessing state

QUALITY CONTROL (OMICVERSE INTEGRATION):
- scrna.run_quality_control() - QC analysis with omicverse.pp.qc
  * Calculate mitochondrial/ribosomal gene percentages
  * Filter low-quality cells and genes
  * Generate QC visualizations

PREPROCESSING (OMICVERSE INTEGRATION):
- scrna.run_preprocessing() - Normalization with omicverse.pp.preprocess
  * Target sum normalization (default: 1e4)
  * Log1p transformation
  * Highly variable gene selection
  * Data scaling for PCA

DIMENSIONALITY REDUCTION (OMICVERSE INTEGRATION):
- scrna.run_pca() - Principal component analysis with omicverse.pp.pca
  * Compute PCA on scaled/normalized data
  * Variance explained analysis
  * Prepare for downstream analysis

BATCH CORRECTION (OMICVERSE INTEGRATION):
- scrna.run_batch_correction() - Integration with omicverse.single.batch_correction
  * scVI-based batch correction
  * Harmony integration
  * Corrected UMAP computation

CLUSTERING & ANNOTATION:
- scrna.run_clustering() - Graph-based clustering with omicverse.utils.clusters
  * Leiden/Louvain clustering algorithms
  * Resolution optimization
  * Cluster validation metrics

- scrna.run_cell_type_annotation() - Comprehensive annotation workflow
  * CellOntologyMapper integration (ov.single.CellOntologyMapper)
  * Marker gene analysis (omicverse.single.get_celltype_marker)
  * Differential expression (scanpy.tl.rank_genes_groups)
  * COSG marker identification (omicverse.single.cosg)
  * Broad cell type assignment

DOWNSTREAM ANALYSIS (PERTPY INTEGRATION):
- Differential cell type analysis with pertpy
- Pathway enrichment analysis
- Trajectory inference capabilities
- Comparative analysis between conditions

REPORTING:
- scrna.generate_report() - Comprehensive analysis report
  * HTML/PDF/Markdown formats
  * Integrated visualizations
  * Methods documentation

🚀 WORKFLOW: 
1. Start by checking your Python environment and installing any missing packages
2. Scan your data folder to identify scRNA-seq files
3. Add todos for your analysis pipeline
4. Follow the adaptive workflow that checks data state at each step:
   - Load & inspect → QC (if needed) → Preprocess (if needed) → PCA (if needed)
   - Batch correction (if needed) → Clustering (if needed) → Annotation → Downstream analysis
5. Each step adapts based on what's already present in your data!

The toolset uses the latest omicverse, scanpy, and pertpy integration for state-of-the-art scRNA-seq analysis!"""
    
    return message