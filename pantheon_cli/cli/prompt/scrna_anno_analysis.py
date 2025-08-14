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
9. "Generate context-specific cell types and markers"
10. "Calculate AUCell scores for cell type markers"
11. "Analyze cluster-celltype associations"
12. "Find cluster-specific marker genes"
13. "Integrate AUCell + markers for final cell type annotation"
14. "Conduct downstream analysis"
15. "Generate analysis report"

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
Use unified workflow: `scrna.run_workflow(workflow_type="qc", output_dir=results_dir)`

🏷️ STEP 3 - PREPROCESSING:
Use unified workflow: `scrna.run_workflow(workflow_type="preprocessing", output_dir=results_dir)`

🏷️ STEP 4 - PCA:
Use unified workflow: `scrna.run_workflow(workflow_type="pca", output_dir=results_dir)`

🏷️ STEP 5 - BATCH CORRECTION (if needed):
```python
# Check if batch correction is needed
if 'batch' in adata.obs.columns:
    print("\\n🔧 Batch key detected - applying batch correction...")
    # Batch correction would go here - currently handled by preprocessing workflow
else:
    print("\\n✅ No batch correction needed")
```

🏷️ STEP 6 - CLUSTERING:
Use unified workflow: `scrna.run_workflow(workflow_type="clustering", output_dir=results_dir)`

🏷️ STEP 7 - VISUALIZATION:
Use unified workflow: `scrna.run_workflow(workflow_type="visualization", output_dir=results_dir)`

🏷️ STEP 8 - DATA CONTEXT COLLECTION:
```python
print("\\n📝 **DATA CONTEXT COLLECTION**")
user_data_context = input("Please briefly describe your data (tissue, condition, experiment): ").strip()
print(f"Data context recorded: {{user_data_context}}")

# Store context in adata
adata.uns['user_data_context'] = user_data_context
```

🏷️ STEP 9 - AUCELL CELL TYPE SCORING:
Use unified workflow: `scrna.run_workflow(workflow_type="aucell", output_dir=results_dir)`

🏷️ STEP 10 - CLUSTER-CELLTYPE ASSOCIATIONS:
```python
print("\\n🔍 Analyzing cluster-celltype associations...")

# Find AUCell columns
celltype_columns = [col for col in adata.obs.columns if 'AUCell' in col]

if celltype_columns:
    # Calculate mean AUCell scores per cluster for each cell type
    cluster_celltype_scores = dict()
    
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_scores = dict()
        for celltype_col in celltype_columns:
            mean_score = adata.obs[adata.obs['leiden'] == cluster][celltype_col].mean()
            cluster_scores[celltype_col] = mean_score
        cluster_celltype_scores[cluster] = cluster_scores
    
    # Find best cell type for each cluster
    cluster_annotations = dict()
    for cluster, scores in cluster_celltype_scores.items():
        if scores:
            best_celltype_col = max(scores.items(), key=lambda x: x[1])[0]
            best_score = scores[best_celltype_col]
            
            # Extract cell type name from column name
            best_celltype = best_celltype_col.replace('_AUCell', '').replace('AUCell_', '')
            
            cluster_annotations[cluster] = best_celltype
            print(f"Cluster {{cluster}} -> {{best_celltype}} (score: {{round(best_score, 3)}})")
    
    # Store preliminary annotations
    adata.uns['preliminary_cluster_annotations'] = cluster_annotations
    print("\\n✅ Preliminary cluster annotations completed")
else:
    print("⚠️ No AUCell columns found - run AUCell scoring first")
```

🏷️ STEP 11 - EXTRACT CLUSTER MARKERS:
Use modular approach: `scrna.extract_markers(method="omicverse", log2fc_min=1.0, pval_cutoff=0.05, top_n=10)`

🏷️ STEP 12 - INTERACTIVE LLM-POWERED CLUSTER ANNOTATION:
```python
print("\\n🤖 Starting interactive LLM-powered cluster annotation workflow...")

# Get user data context if not already set
user_data_context = adata.uns.get('user_data_context', '')
if not user_data_context:
    user_data_context = input("\\n📝 **DATA CONTEXT:** Please briefly describe your data (tissue, condition, experiment): ").strip()
    adata.uns['user_data_context'] = user_data_context

# Launch comprehensive annotation workflow using modular toolset functions
annotation_workflow = scrna.llm_anno(
    cluster_id="all_clusters", 
    user_context=user_data_context,
    confidence_threshold=0.7
)

print("\\n🎯 Comprehensive cluster annotation workflow initiated")
print("💬 Follow the interactive prompts for each cluster annotation")
```

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

🏷️ STEP 14 - FINAL DATA EXPORT:
```python
print("\\n💾 Final data export...")

# Save final annotated dataset
final_export = scrna.sc_save(save_type="adata", output_dir=results_dir)

print("\\n🎉 **ANALYSIS PIPELINE COMPLETE**")
print(f"📁 **Results saved to:** {{results_dir}}")
print("📊 **Pipeline summary:**")
print("1. ✅ Data loading and inspection")
print("2. ✅ Quality control with omicverse")
print("3. ✅ Preprocessing and normalization")
print("4. ✅ PCA and dimensionality reduction")
print("5. ✅ Clustering analysis")
print("6. ✅ Visualization generation")
print("7. ✅ AUCell cell type scoring")
print("8. ✅ Interactive LLM-powered annotation")
print("9. ✅ Comprehensive reporting")
print("10. ✅ Final data export")
```

🔧 **AVAILABLE TOOLSET FUNCTIONS:**

**UNIFIED WORKFLOW ENGINE:**
- `scrna.run_workflow(workflow_type="qc")` - Quality control with omicverse
- `scrna.run_workflow(workflow_type="preprocessing")` - Preprocessing with omicverse
- `scrna.run_workflow(workflow_type="pca")` - PCA with omicverse
- `scrna.run_workflow(workflow_type="clustering")` - Clustering analysis
- `scrna.run_workflow(workflow_type="visualization")` - Generate plots
- `scrna.run_workflow(workflow_type="aucell")` - AUCell scoring

**SPECIALIZED FUNCTIONS:**
- `scrna.load_and_inspect_data()` - Data loading and inspection
- `scrna.extract_markers()` - Extract cluster markers
- `scrna.analyze_aucell()` - Analyze AUCell for specific clusters
- `scrna.prepare_evidence()` - Prepare annotation evidence
- `scrna.llm_anno()` - Interactive LLM annotation workflow
- `scrna.sc_save()` - Save analysis results
- `scrna.generate_report()` - Generate comprehensive reports

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
- scrna.run_workflow(workflow_type="qc") - Complete quality control workflow
- scrna.run_workflow(workflow_type="preprocessing") - Normalization and preprocessing
- scrna.run_workflow(workflow_type="pca") - Principal component analysis
- scrna.run_workflow(workflow_type="clustering") - Clustering with UMAP and Leiden
- scrna.run_workflow(workflow_type="visualization") - Generate visualizations
- scrna.run_workflow(workflow_type="aucell") - AUCell cell type scoring

**SPECIALIZED FUNCTIONS:**
- scrna.load_and_inspect_data() - Load and inspect scRNA-seq data
- scrna.extract_markers() - Cluster-specific marker gene extraction  
- scrna.analyze_aucell() - AUCell analysis for specific clusters
- scrna.prepare_evidence() - Prepare annotation evidence summaries
- scrna.llm_anno() - Interactive LLM annotation workflow
- scrna.sc_save() - Save analysis results
- scrna.generate_report() - Generate comprehensive reports

**GUIDANCE:**
- scrna.suggest_next_step() - Smart recommendations

Please start by adding a todo for your scRNA-seq analysis task, then use the appropriate scRNA tools!"""
    
    return message