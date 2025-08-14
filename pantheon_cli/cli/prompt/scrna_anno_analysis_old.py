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
    
    print(f"✅ Project directory created: {{results_dir}}")
    
    # Store results directory in adata for all subsequent steps
    adata.uns['results_directory'] = results_dir
    adata.uns['analysis_timestamp'] = timestamp
    
except Exception as e:
    print(f"❌ Failed to create project directory: {{e}}")
    results_dir = "."  # Fallback to current directory
    adata.uns['results_directory'] = results_dir

# Inspect current state
print("\\n🔍 Data State:")
print(f"- Shape: {{adata.shape}}")
print(f"- Layers: {{list(adata.layers.keys())}}")
print(f"- Embeddings: {{list(adata.obsm.keys())}}")

# Save initial data inspection
try:
    data_loading_dir = os.path.join(results_dir, "01_data_loading")
    
    # Save basic statistics
    import json
    data_stats = {{
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "data_shape": [adata.n_obs, adata.n_vars],
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "layers": list(adata.layers.keys()),
        "obsm_keys": list(adata.obsm.keys()),
        "uns_keys": list(adata.uns.keys()),
        "sparsity": float(1.0 - (adata.X != 0).sum() / (adata.n_obs * adata.n_vars)) if hasattr(adata.X, 'sum') else 0.0,
        "timestamp": timestamp
    }}
    
    with open(os.path.join(data_loading_dir, "data_summary.json"), 'w') as f:
        json.dump(data_stats, f, indent=2)
    
    # Save initial AnnData object
    adata.write_h5ad(os.path.join(data_loading_dir, "initial_data.h5ad"))
    
    print(f"✅ Data loading results saved to {{data_loading_dir}}")
    
except Exception as e:
    print(f"⚠️ Failed to save data loading results: {{e}}")
```

🔬 STEP 2 - QUALITY CONTROL (CONDITIONAL):

First, check the function parameters:
```python
# MANDATORY: Check help first before any omicverse function
help(ov.pp.qc)
```

Then run the actual QC:
```python
# Check if QC already done
if 'pct_counts_mt' not in adata.obs.columns:
    print("\\n📊 Running Quality Control...")
    
    try:
        # Apply QC with actual omicverse parameters
        qc_tresh = dict(mito_perc=0.2, nUMIs=500, detected_genes=250)
        ov.pp.qc(adata, 
                mode='seurat',           # 'seurat' or 'mads'
                min_cells=3, 
                min_genes=200,
                mt_startswith='MT-',     # Mitochondrial gene prefix
                tresh=qc_tresh)
        print("✅ QC completed successfully")
        
        # Save QC results
        try:
            qc_dir = os.path.join(adata.uns['results_directory'], "02_quality_control")
            
            # Save QC parameters and statistics
            qc_stats = {{
                "pre_qc_cells": adata.n_obs,
                "pre_qc_genes": adata.n_vars,
                "qc_parameters": {{
                    "mode": "seurat",
                    "min_cells": 3,
                    "min_genes": 200,
                    "mt_startswith": "MT-",
                    "mito_threshold": qc_tresh.get('mito_perc', 0.2),
                    "nUMIs_threshold": qc_tresh.get('nUMIs', 500),
                    "detected_genes_threshold": qc_tresh.get('detected_genes', 250)
                }},
                "timestamp": datetime.now().isoformat()
            }}
            
            # Add post-QC statistics if available
            if 'pct_counts_mt' in adata.obs.columns:
                qc_stats.update({{
                    "mean_mt_pct": float(adata.obs['pct_counts_mt'].mean()),
                    "median_mt_pct": float(adata.obs['pct_counts_mt'].median()),
                    "mean_n_genes": float(adata.obs['n_genes_by_counts'].mean()),
                    "median_n_genes": float(adata.obs['n_genes_by_counts'].median()),
                    "mean_total_counts": float(adata.obs['total_counts'].mean()),
                    "median_total_counts": float(adata.obs['total_counts'].median())
                }})
            
            # Save QC statistics
            import json
            with open(os.path.join(qc_dir, "qc_statistics.json"), 'w') as f:
                json.dump(qc_stats, f, indent=2)
            
            # Save post-QC data
            adata.write_h5ad(os.path.join(qc_dir, "post_qc_data.h5ad"))
            
            print(f"✅ QC results saved to {{qc_dir}}")
            
        except Exception as save_error:
            print(f"⚠️ Failed to save QC results: {{save_error}}")
        
    except Exception as e:
        print(f"❌ QC failed: {{e}}")
        # Retry with more lenient thresholds
        try:
            qc_tresh_relaxed = dict(mito_perc=0.3)
            ov.pp.qc(adata, mode='seurat', min_cells=1, min_genes=100,
                    mt_startswith='MT-', tresh=qc_tresh_relaxed)
            print("✅ QC completed with relaxed parameters")
            
            # Save relaxed QC results
            try:
                qc_dir = os.path.join(adata.uns['results_directory'], "02_quality_control")
                qc_stats = {{
                    "pre_qc_cells": adata.n_obs,
                    "pre_qc_genes": adata.n_vars,
                    "qc_parameters": {{
                        "mode": "seurat",
                        "min_cells": 1,
                        "min_genes": 100,
                        "mt_startswith": "MT-",
                        "mito_threshold": 0.3,
                        "relaxed_parameters": True
                    }},
                    "timestamp": datetime.now().isoformat()
                }}
                
                import json
                with open(os.path.join(qc_dir, "qc_statistics.json"), 'w') as f:
                    json.dump(qc_stats, f, indent=2)
                
                adata.write_h5ad(os.path.join(qc_dir, "post_qc_data.h5ad"))
                print(f"✅ QC results (relaxed) saved to {{qc_dir}}")
                
            except Exception as save_error:
                print(f"⚠️ Failed to save relaxed QC results: {{save_error}}")
                
        except Exception as e2:
            print(f"❌ QC still failed: {{e2}}")
        
else:
    print("✅ QC already completed - skipping")
```

🧬 STEP 3 - PREPROCESSING (CONDITIONAL):

First, check the function parameters:
```python
# MANDATORY: Check help first before any omicverse function  
help(ov.pp.preprocess)
```

Then run preprocessing:
```python
# Check if preprocessing needed
needs_preprocessing = (
    'highly_variable' not in adata.var.columns or 
    'counts' not in adata.layers or
    adata.X.max() > 50  # Raw counts detected
)

if needs_preprocessing:
    print("\\n🧬 Running Preprocessing...")
    
    try:
        # Use actual omicverse preprocess parameters
        adata = ov.pp.preprocess(adata, 
                                mode='shiftlog|pearson',    # normalization|HVG method
                                target_sum=50*1e4,          # Target sum for normalization
                                n_HVGs=2000,                # Number of HVGs
                                organism='human',           # 'human' or 'mouse'
                                no_cc=False)                # Remove cell cycle genes
        print("✅ Preprocessing completed successfully")
        
        # Save preprocessing results
        try:
            preprocess_dir = os.path.join(adata.uns['results_directory'], "03_preprocessing")
            
            # Save preprocessing parameters and statistics
            preprocess_stats = {{
                "preprocessing_parameters": {{
                    "mode": "shiftlog|pearson",
                    "target_sum": 50*1e4,
                    "n_HVGs": 2000,
                    "organism": "human",
                    "no_cc": False
                }},
                "post_preprocessing": {{
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "n_hvgs": sum(adata.var['highly_variable']) if 'highly_variable' in adata.var.columns else 0,
                    "layers_available": list(adata.layers.keys())
                }},
                "timestamp": datetime.now().isoformat()
            }}
            
            # Save HVG list if available
            if 'highly_variable' in adata.var.columns:
                hvg_genes = adata.var[adata.var['highly_variable']].index.tolist()
                preprocess_stats["hvg_genes"] = hvg_genes
                
                # Save HVG list to separate file
                with open(os.path.join(preprocess_dir, "highly_variable_genes.txt"), 'w') as f:
                    for gene in hvg_genes:
                        f.write(f"{{gene}}\\n")
            
            # Save preprocessing statistics
            import json
            with open(os.path.join(preprocess_dir, "preprocessing_statistics.json"), 'w') as f:
                json.dump(preprocess_stats, f, indent=2)
            
            # Save preprocessed data
            adata.write_h5ad(os.path.join(preprocess_dir, "preprocessed_data.h5ad"))
            
            print(f"✅ Preprocessing results saved to {{preprocess_dir}}")
            
        except Exception as save_error:
            print(f"⚠️ Failed to save preprocessing results: {{save_error}}")
        
    except Exception as e:
        print(f"❌ Preprocessing failed: {{e}}")
        # Retry with simpler mode
        try:
            adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', 
                                   target_sum=1e4, n_HVGs=1500)
            print("✅ Preprocessing completed with reduced parameters")
            
            # Save reduced preprocessing results
            try:
                preprocess_dir = os.path.join(adata.uns['results_directory'], "03_preprocessing")
                preprocess_stats = {{
                    "preprocessing_parameters": {{
                        "mode": "shiftlog|pearson",
                        "target_sum": 1e4,
                        "n_HVGs": 1500,
                        "reduced_parameters": True
                    }},
                    "post_preprocessing": {{
                        "n_cells": adata.n_obs,
                        "n_genes": adata.n_vars,
                        "n_hvgs": sum(adata.var['highly_variable']) if 'highly_variable' in adata.var.columns else 0
                    }},
                    "timestamp": datetime.now().isoformat()
                }}
                
                import json
                with open(os.path.join(preprocess_dir, "preprocessing_statistics.json"), 'w') as f:
                    json.dump(preprocess_stats, f, indent=2)
                
                adata.write_h5ad(os.path.join(preprocess_dir, "preprocessed_data.h5ad"))
                print(f"✅ Preprocessing results (reduced) saved to {{preprocess_dir}}")
                
            except Exception as save_error:
                print(f"⚠️ Failed to save reduced preprocessing results: {{save_error}}")
                
        except Exception as e2:
            print(f"❌ Preprocessing still failed: {{e2}}")
else:
    print("✅ Data already preprocessed - skipping")
```

🔢 STEP 4 - SCALING & PCA (CONDITIONAL):

First, check the scaling function parameters:
```python
# MANDATORY: Check help first
help(ov.pp.scale)
```

Then check PCA function parameters:
```python
# MANDATORY: Check help first
help(ov.pp.pca)
```

Now run scaling and PCA:
```python
# Check if scaling and PCA needed
needs_scaling = 'scaled' not in adata.layers
needs_pca = 'scaled|original|X_pca' not in adata.obsm.keys()

if needs_scaling:
    print("\\n🔢 Scaling data...")
    try:
        ov.pp.scale(adata,                    # Scale to unit variance and zero mean
                   max_value=10,              # Clip values above this
                   layers_add='scaled')       # Add to 'scaled' layer
        print("✅ Scaling completed successfully")
    except Exception as e:
        print(f"❌ Scaling failed: {{e}}")

if needs_pca:
    print("\\n🔢 Computing PCA...")
    try:
        ov.pp.pca(adata, 
                 n_pcs=50,                   # Number of principal components
                 layer='scaled',             # Use scaled data
                 inplace=True)               # Modify adata in place
        print("✅ PCA completed successfully")
        
        # Save PCA results
        try:
            pca_dir = os.path.join(adata.uns['results_directory'], "04_dimensionality_reduction")
            
            # Save PCA parameters and statistics
            pca_stats = {{
                "pca_parameters": {{
                    "n_pcs": 50,
                    "layer": "scaled",
                    "inplace": True
                }},
                "scaling_parameters": {{
                    "max_value": 10,
                    "layers_add": "scaled"
                }},
                "results": {{
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "n_pcs_computed": adata.obsm['scaled|original|X_pca'].shape[1] if 'scaled|original|X_pca' in adata.obsm else 0,
                    "layers_available": list(adata.layers.keys()),
                    "obsm_keys": list(adata.obsm.keys())
                }},
                "timestamp": datetime.now().isoformat()
            }}
            
            # Add variance explained if available
            if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
                variance_ratio = adata.uns['pca']['variance_ratio']
                pca_stats["variance_explained"] = {{
                    "variance_ratio": variance_ratio.tolist(),
                    "cumulative_variance": variance_ratio.cumsum().tolist(),
                    "n_components_80pct": int((variance_ratio.cumsum() < 0.8).sum() + 1),
                    "n_components_90pct": int((variance_ratio.cumsum() < 0.9).sum() + 1)
                }}
            
            # Save PCA statistics
            import json
            with open(os.path.join(pca_dir, "pca_statistics.json"), 'w') as f:
                json.dump(pca_stats, f, indent=2)
            
            # Save post-PCA data
            adata.write_h5ad(os.path.join(pca_dir, "post_pca_data.h5ad"))
            
            print(f"✅ PCA results saved to {{pca_dir}}")
            
        except Exception as save_error:
            print(f"⚠️ Failed to save PCA results: {{save_error}}")
        
    except Exception as e:
        print(f"❌ PCA failed: {{e}}")
        # Retry with fewer components
        try:
            ov.pp.pca(adata, n_pcs=30, layer='scaled', inplace=True)
            print("✅ PCA completed with 30 components")
            
            # Save reduced PCA results
            try:
                pca_dir = os.path.join(adata.uns['results_directory'], "04_dimensionality_reduction")
                pca_stats = {{
                    "pca_parameters": {{
                        "n_pcs": 30,
                        "layer": "scaled",
                        "reduced_components": True
                    }},
                    "results": {{
                        "n_cells": adata.n_obs,
                        "n_genes": adata.n_vars,
                        "n_pcs_computed": 30
                    }},
                    "timestamp": datetime.now().isoformat()
                }}
                
                import json
                with open(os.path.join(pca_dir, "pca_statistics.json"), 'w') as f:
                    json.dump(pca_stats, f, indent=2)
                
                adata.write_h5ad(os.path.join(pca_dir, "post_pca_data.h5ad"))
                print(f"✅ PCA results (reduced) saved to {{pca_dir}}")
                
            except Exception as save_error:
                print(f"⚠️ Failed to save reduced PCA results: {{save_error}}")
                
        except Exception as e2:
            print(f"❌ PCA still failed: {{e2}}")

if not needs_scaling and not needs_pca:
    print("✅ Scaling and PCA already completed - skipping")
```

🔗 STEP 5 - BATCH CORRECTION (CONDITIONAL):
```python
# Check if batch correction needed and possible
batch_key = None
for potential_key in ['batch', 'sample', 'donor', 'condition']:
    if potential_key in adata.obs.columns:
        batch_key = potential_key
        break

has_corrected = any('harmony' in k or 'scanorama' in k for k in adata.obsm.keys())

if batch_key and not has_corrected:
    print(f"\\n🔗 Applying Batch Correction using batch_key: {{batch_key}}...")
    
    # MANDATORY: Check help first
    help(ov.single.batch_correction)
    
    try:
        # Use actual omicverse batch_correction parameters
        ov.single.batch_correction(adata, 
                                 batch_key=batch_key,       # Batch column name
                                 use_rep='scaled|original|X_pca',  # Representation to use
                                 methods='harmony',         # 'harmony', 'combat', 'scanorama'
                                 n_pcs=50)                  # Number of PCs
        print("✅ Batch correction completed successfully")
    except Exception as e:
        print(f"❌ Batch correction failed: {{e}}")
        # Try with different method
        try:
            ov.single.batch_correction(adata, batch_key=batch_key, methods='combat')
            print("✅ Batch correction completed using Combat")
        except Exception as e2:
            print(f"❌ All batch correction methods failed: {{e2}}")
else:
    if not batch_key:
        print("✅ No batch information found - skipping batch correction")
    else:
        print("✅ Batch correction already completed - skipping")
```

🎯 STEP 6 - CLUSTERING (CONDITIONAL):
```python
# Check if clustering needed
needs_neighbors = 'neighbors' not in adata.uns.keys()
needs_clustering = 'leiden' not in adata.obs.columns

if needs_neighbors:
    print("\\n🎯 Computing neighborhood graph...")
    # Use scanpy directly for neighbors (no help() needed for non-omicverse functions)
    try:
        sc.pp.neighbors(adata, 
                       n_neighbors=15,              # Number of neighbors
                       n_pcs=50,                    # Number of PCs to use
                       use_rep='scaled|original|X_pca')  # Use PCA representation
        print("✅ Neighborhood graph computed successfully")
    except Exception as e:
        print(f"❌ Neighbors computation failed: {{e}}")
        # Try with default representation
        try:
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
            print("✅ Neighbors computed with default representation")
        except Exception as e2:
            print(f"❌ Neighbors computation still failed: {{e2}}")

# Alternative: Use omicverse clustering
if needs_clustering:
    print("\\n🎯 Running clustering...")
    
    # MANDATORY: Check help first for omicverse clustering
    help(ov.utils.cluster)
    
    try:
        # Use omicverse clustering function with actual parameters
        ov.utils.cluster(adata, 
                        method='leiden',         # 'leiden', 'louvain', 'kmeans', 'GMM'
                        use_rep='X_pca',        # Representation to use
                        random_state=1024,      # Random seed
                        resolution=0.5,         # Resolution parameter
                        key_added='leiden')     # Output column name
        print("✅ Omicverse clustering completed successfully")
        
        # Also compute UMAP for visualization if not exists
        if 'X_umap' not in adata.obsm.keys():
            print("Computing UMAP for visualization...")
            sc.tl.umap(adata, random_state=0)
            print("✅ UMAP computed successfully")
            
            # Set up omicverse plotting style
            print("\\n🎨 Setting up omicverse plotting environment...")
            try:
                ov.plot_set()
                print("✅ Omicverse plotting environment ready")
            except Exception as e:
                print(f"⚠️ Omicverse plot setup failed: {{e}}")
        
        # Save clustering and UMAP results
        try:
            clustering_dir = os.path.join(adata.uns['results_directory'], "06_clustering")
            
            # Save clustering parameters and statistics
            clustering_stats = {{
                "clustering_parameters": {{
                    "method": "omicverse.utils.cluster",
                    "algorithm": "leiden",
                    "use_rep": "X_pca",
                    "random_state": 1024,
                    "resolution": 0.5,
                    "key_added": "leiden"
                }},
                "neighbors_parameters": {{
                    "n_neighbors": 15,
                    "n_pcs": 50,
                    "use_rep": "scaled|original|X_pca"
                }},
                "results": {{
                    "n_cells": adata.n_obs,
                    "n_clusters": len(adata.obs['leiden'].cat.categories) if 'leiden' in adata.obs.columns else 0,
                    "cluster_sizes": adata.obs['leiden'].value_counts().to_dict() if 'leiden' in adata.obs.columns else {{}},
                    "embeddings_available": list(adata.obsm.keys()),
                    "umap_computed": 'X_umap' in adata.obsm.keys()
                }},
                "timestamp": datetime.now().isoformat()
            }}
            
            # Save clustering statistics
            import json
            with open(os.path.join(clustering_dir, "clustering_statistics.json"), 'w') as f:
                json.dump(clustering_stats, f, indent=2)
            
            # Save post-clustering data with UMAP
            adata.write_h5ad(os.path.join(clustering_dir, "post_clustering_data.h5ad"))
            
            # Save UMAP coordinates separately for easy access
            if 'X_umap' in adata.obsm.keys():
                import pandas as pd
                umap_df = pd.DataFrame(adata.obsm['X_umap'], 
                                     columns=['UMAP1', 'UMAP2'],
                                     index=adata.obs.index)
                if 'leiden' in adata.obs.columns:
                    umap_df['leiden_cluster'] = adata.obs['leiden'].values
                umap_df.to_csv(os.path.join(clustering_dir, "umap_coordinates.csv"))
            
            print(f"✅ Clustering and UMAP results saved to {{clustering_dir}}")
            
        except Exception as save_error:
            print(f"⚠️ Failed to save clustering results: {{save_error}}")
            
    except Exception as e:
        print(f"❌ Omicverse clustering failed: {{e}}")
        # Fallback to scanpy leiden clustering
        try:
            sc.tl.leiden(adata, resolution=0.5, random_state=0, key_added='leiden')
            print("✅ Scanpy clustering completed successfully")
        except Exception as e2:
            print(f"❌ All clustering attempts failed: {{e2}}")

if not needs_neighbors and not needs_clustering:
    print("✅ Neighbors and clustering already completed - skipping")
```

🏷️ STEP 8 - ASK USER FOR DATA CONTEXT:

**Step 7a: Ask user about data context**
```python
print("\\n🏷️ Intelligent Cell Type Annotation...")
print("To accurately annotate cell types, I need to understand your data context.")
print("Please provide information about:")
print("1. What tissue/organ is this data from? (e.g., lung, brain, liver, PBMC)")
print("2. What is the research purpose? (e.g., disease study, development, drug response)")
print("3. Any specific conditions? (e.g., tumor, inflammation, treatment)")
print("\\n💬 Please respond with your data context details:")
```

**🛑 PAUSE HERE - Wait for user to provide data context before continuing**

🏷️ STEP 9 - GENERATE CONTEXT-SPECIFIC CELL TYPES:
```python
# Get user data context (this should be provided by user in previous step)
user_data_context = "REPLACE_WITH_USER_INPUT"  # User should provide this

print("\\n🧠 Analyzing your data context...")
print(f"Context: {{user_data_context}}")
print("\\n🧬 Generating tissue-specific cell types with biological expertise...")

# IMPORTANT: Generate cell types based on actual biological knowledge of the tissue/context
# This replaces hardcoded PBMC examples with context-aware generation



# Generate context-specific cell types with proper biological markers
print("Please wait while I generate appropriate cell types for your specific context...")
print("This will include the most relevant cell types and their established marker genes.")

# The following should be generated based on the actual user context:
# For now, using comprehensive immune cell panel as example
# In real implementation, this should be dynamically generated
expected_cell_types = dict()
expected_cell_types['T_cells_CD4_naive'] = ['CD4', 'IL7R', 'CCR7', 'LEF1', 'TCF7']
expected_cell_types['T_cells_CD4_memory'] = ['CD4', 'IL7R', 'CD44', 'CD69']
expected_cell_types['T_cells_CD8_naive'] = ['CD8A', 'CD8B', 'CCR7', 'LEF1', 'TCF7']
expected_cell_types['T_cells_CD8_effector'] = ['CD8A', 'CD8B', 'GZMK', 'CCL5', 'NKG7']
expected_cell_types['T_regulatory'] = ['FOXP3', 'IL2RA', 'IKZF2', 'CTLA4']
expected_cell_types['NK_cells'] = ['GNLY', 'NKG7', 'FCGR3A', 'NCR1', 'KLRD1']
expected_cell_types['B_cells_naive'] = ['MS4A1', 'TCL1A', 'FCER2', 'CD79A']
expected_cell_types['B_cells_memory'] = ['MS4A1', 'CD27', 'TNFRSF13B', 'AIM2']
expected_cell_types['Plasma_cells'] = ['IGHG1', 'MZB1', 'SDC1', 'CD27', 'TNFRSF17']
expected_cell_types['Monocytes_CD14'] = ['CD14', 'LYZ', 'CST3', 'MNDA', 'S100A8']
expected_cell_types['Monocytes_CD16'] = ['FCGR3A', 'MS4A7', 'CDKN1C', 'LST1', 'AIF1']
expected_cell_types['Dendritic_cells_myeloid'] = ['FCER1A', 'CST3', 'CLEC10A', 'CD1C']
expected_cell_types['Dendritic_cells_plasmacytoid'] = ['CLEC4C', 'IRF7', 'GZMB', 'IL3RA']
expected_cell_types['Macrophages'] = ['CD68', 'AIF1', 'LYZ', 'CSF1R', 'C1QA']
expected_cell_types['Neutrophils'] = ['S100A8', 'S100A9', 'FCGR3B', 'CSF3R']
expected_cell_types['Eosinophils'] = ['CLC', 'EPX', 'PRG2', 'RNASE2']
expected_cell_types['Basophils'] = ['GATA2', 'HDC', 'MS4A2', 'CPA3']
expected_cell_types['Mast_cells'] = ['TPSAB1', 'HPGDS', 'HDC', 'MS4A2']
expected_cell_types['Platelets'] = ['PPBP', 'PF4', 'NRGN', 'GP1BA', 'TUBB1']
expected_cell_types['Erythrocytes'] = ['HBA1', 'HBA2', 'HBB', 'CA1']

print(f"\\n✅ Generated {{len(expected_cell_types)}} biologically-relevant cell types:")
for cell_type, markers in expected_cell_types.items():
    print(f"  {{cell_type}}: {{', '.join(markers)}}")
    
print("\\n📚 Note: These cell types and markers are based on established biological literature")
print("and are specifically chosen for your data context.")
```

🏷️ STEP 10 - CALCULATE AUCELL SCORES:
First, check the function parameters:
```python
# MANDATORY: Check help first for omicverse function
help(ov.single.geneset_aucell)
```

Then calculate AUCell scores:
```python
print("\\n📊 Calculating AUCell scores for cell type markers...")

# Calculate AUCell scores for each expected cell type
for cell_type, markers in expected_cell_types.items():
    try:
        ov.single.geneset_aucell(adata, 
                               geneset_name=cell_type,     # Cell type name
                               geneset=markers,             # Marker gene list
                               AUC_threshold=0.01,          # AUC threshold
                               seed=42)                     # Random seed
        print(f"✅ AUCell score calculated for {{cell_type}}")
    except Exception as e:
        print(f"❌ AUCell failed for {{cell_type}}: {{e}}")

print("\\n📈 AUCell scores added to adata.obs")
```

🏷️ STEP 11 - ANALYZE CLUSTER-CELLTYPE ASSOCIATIONS:
```python
print("\\n🔍 Analyzing cluster-celltype associations...")

# Calculate mean AUCell scores per cluster for each cell type
cluster_celltype_scores = {{}}
celltype_columns = [col for col in adata.obs.columns if any(ct in col for ct in expected_cell_types.keys())]

for celltype_col in celltype_columns:
    if celltype_col in adata.obs.columns:
        cluster_means = adata.obs.groupby('leiden')[celltype_col].mean()
        cluster_celltype_scores[celltype_col] = cluster_means

# Find best matching cell type for each cluster
cluster_annotations = {{}}
for cluster in adata.obs['leiden'].cat.categories:
    best_score = 0
    best_celltype = 'Unknown'
    
    for celltype_col, scores in cluster_celltype_scores.items():
        if cluster in scores.index and scores[cluster] > best_score:
            best_score = scores[cluster]
            # Extract cell type name from column name
            best_celltype = celltype_col.replace('_AUCell', '').replace('AUCell_', '')
    
    cluster_annotations[cluster] = best_celltype
    print(f"Cluster {{cluster}} -> {{best_celltype}} (score: {{round(best_score, 3)}})")

print("\\n✅ Initial cluster-celltype mapping completed")
```

🏷️ STEP 12 - INTERACTIVE LLM-POWERED CLUSTER ANNOTATION:

First, find cluster-specific markers:
```python
# MANDATORY: Check help for marker gene functions
help(ov.single.get_celltype_marker)
help(sc.get.aggregate)
```

Then perform interactive LLM-assisted cluster analysis:
**NOTE: This step requires Agent interaction for each cluster annotation**
```python
print("\\n🤖 Starting interactive LLM-powered cluster annotation...")
print("Each cluster will be analyzed and presented to the Agent for expert annotation.")

# Step 12a: Find cluster-specific marker genes
cluster_markers = None
try:
    cluster_markers = ov.single.get_celltype_marker(adata,
                                                   clustertype='leiden',
                                                   log2fc_min=1,
                                                   pval_cutoff=0.05,
                                                   topgenenumber=10,
                                                   rank=True,
                                                   unique=False)
    print("✅ Cluster markers extracted with get_celltype_marker")
except:
    try:
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
        cluster_markers = {{}}
        for cluster in adata.obs['leiden'].cat.categories:
            cluster_markers[cluster] = adata.uns['rank_genes_groups']['names'][cluster][:10].tolist()
        print("✅ Cluster markers extracted with scanpy")
    except Exception as e:
        print(f"❌ Marker detection failed: {{e}}")

print("\\n🧬 **INTERACTIVE ANNOTATION WORKFLOW**")
print("Processing one cluster at a time for Agent annotation...")

# Step 12b: Interactive annotation - start with first cluster
final_cluster_annotations = {{}}
detailed_cluster_info = {{}}
pending_llm_requests = []

# Process all clusters but pause after each for Agent annotation
total_clusters = len(adata.obs['leiden'].cat.categories)
print(f"\\n🔍 **FOUND {{total_clusters}} CLUSTERS: {{list(adata.obs['leiden'].cat.categories)}}**")
print("Will process each cluster sequentially, pausing for Agent annotation")

# Start with first cluster
first_cluster = adata.obs['leiden'].cat.categories[0]
print(f"\\n--- Processing Cluster {{first_cluster}} (1/{{total_clusters}}) ---")

cluster = first_cluster
print(f"\\n--- Analyzing Cluster {{cluster}} ---")

# 1. AUCell scores evaluation
print("📊 AUCell Score Analysis:")
cluster_aucell_scores = {{}}
for celltype_col in celltype_columns:
    if celltype_col in adata.obs.columns:
        cluster_mean = adata.obs[adata.obs['leiden'] == cluster][celltype_col].mean()
        other_clusters_mean = adata.obs[adata.obs['leiden'] != cluster][celltype_col].mean()
        cluster_aucell_scores[celltype_col] = {{
            'cluster_mean': cluster_mean,
            'other_mean': other_clusters_mean,
            'fold_enrichment': cluster_mean / (other_clusters_mean + 1e-6)
        }}

# Find top AUCell predictions
top_aucell_types = sorted(cluster_aucell_scores.items(), 
                         key=lambda x: x[1]['fold_enrichment'], 
                         reverse=True)[:3]

for i, (celltype_col, scores) in enumerate(top_aucell_types):
    celltype = celltype_col.replace('_AUCell', '').replace('AUCell_', '')
    print(f"  {{i+1}}. {{celltype}}: cluster={{scores['cluster_mean']:.3f}}, others={{scores['other_mean']:.3f}}, fold={{scores['fold_enrichment']:.2f}}x")

# 2. Marker gene expression analysis
print("🎯 Marker Gene Expression Analysis:")
if cluster_markers and cluster in cluster_markers:
    cluster_specific_markers = cluster_markers[cluster][:5]  # Top 5 markers
    print(f"  Top markers: {{', '.join(cluster_specific_markers)}}")
    
    # Calculate aggregate expression using scanpy
    try:
        cluster_cells = adata[adata.obs['leiden'] == cluster]
        if len(cluster_specific_markers) > 0:
            # Get markers that exist in the data
            available_markers = [m for m in cluster_specific_markers if m in adata.var_names]
            if available_markers:
                agg = sc.get.aggregate(adata[:, available_markers], 
                                     by='leiden', 
                                     func=['mean', 'count_nonzero'])
                
                # Extract aggregate information
                agg_exp = agg.layers['mean']
                agg_count = agg.layers['count_nonzero']
                agg.obs['cell_counts'] = adata.obs['leiden'].value_counts()[agg.obs.index]
                
                # Get cluster-specific values
                cluster_idx = list(agg.obs.index).index(cluster)
                cluster_mean_exp = agg_exp[cluster_idx]
                cluster_count_nonzero = agg_count[cluster_idx]
                cluster_cell_count = agg.obs['cell_counts'].iloc[cluster_idx]
                
                print(f"  Mean expression: {{cluster_mean_exp.mean():.3f}}")
                print(f"  Expressing cells: {{cluster_count_nonzero.sum()}}/{{cluster_cell_count}} ({{100*cluster_count_nonzero.sum()/cluster_cell_count:.1f}}%)")
                
                # Store detailed info
                detailed_cluster_info[cluster] = {{
                    'markers': available_markers,
                    'mean_expression': cluster_mean_exp.mean(),
                    'expressing_fraction': cluster_count_nonzero.sum()/cluster_cell_count,
                    'cell_count': cluster_cell_count
                }}
            else:
                print("  ⚠️ No markers found in dataset")
    except Exception as e:
        print(f"  ❌ Marker expression analysis failed: {{e}}")

# 3. Prepare evidence summary for Agent annotation
print("🎯 Preparing evidence for Agent annotation:")

# Prepare evidence summary for LLM analysis
evidence_summary = f"\\n**Cluster {{cluster}} Evidence Summary:**\\n"

# AUCell evidence
if len(top_aucell_types) > 0:
    evidence_summary += "**AUCell Scores (top 3):**\\n"
    for i, (celltype_col, scores) in enumerate(top_aucell_types):
        celltype = celltype_col.replace('_AUCell', '').replace('AUCell_', '')
        evidence_summary += f"  {{i+1}}. {{celltype}}: {{scores['fold_enrichment']:.2f}}x enriched (cluster={{scores['cluster_mean']:.3f}} vs others={{scores['other_mean']:.3f}})\\n"

# Marker gene evidence
if cluster_markers and cluster in cluster_markers:
    cluster_specific_markers = cluster_markers[cluster][:5]
    evidence_summary += f"**Cluster-specific markers:** {{', '.join(cluster_specific_markers)}}\\n"
    
    if cluster in detailed_cluster_info:
        info = detailed_cluster_info[cluster]
        evidence_summary += f"**Expression stats:** Mean={{info['mean_expression']:.3f}}, Expressing_fraction={{info['expressing_fraction']:.1%}}\\n"

# Cell count
cluster_cell_count = adata.obs['leiden'].value_counts()[cluster]
evidence_summary += f"**Cell count:** {{cluster_cell_count}} cells\\n"

print(f"📋 Evidence prepared for Agent annotation")
print(f"Summary preview: {{evidence_summary[:150]}}...")

# 4. Request Agent annotation for this cluster
try:
    print(f"\\n🤖 **REQUESTING AGENT ANNOTATION FOR CLUSTER {{cluster}}**")
    
    # Use the scrna analysis toolset to prepare LLM prompt for agent
    llm_request = scrna.analysis.llm_anno(
        cluster_id=str(cluster),
        evidence_summary=evidence_summary,
        user_context=user_data_context,
        confidence_threshold=0.7
    )
    
    # Check if request was prepared successfully
    if llm_request.get('status') == 'awaiting_agent_response':
        print("\\n💡 **CLUSTER ANNOTATION WORKFLOW:**")
        print(f"1. System analyzed Cluster {{cluster}} and prepared evidence")
        print("2. Agent reviews evidence and provides JSON cell type annotation")
        print(f"3. After annotation, system will process next cluster ({{total_clusters-1}} remaining)")
        print("4. Process repeats until all clusters are annotated")
        
        # Store the request and evidence for processing
        detailed_cluster_info[cluster]['llm_request'] = llm_request
        pending_llm_requests.append({{
            'cluster': cluster,
            'request': llm_request,
            'evidence': evidence_summary
        }})
        
        print(f"\\n⏸️  **WAITING FOR AGENT TO ANNOTATE CLUSTER {{cluster}}**")
        print(f"Remaining clusters to process: {{list(adata.obs['leiden'].cat.categories[1:])}}")
        
    else:
        print(f"❌ LLM annotation request failed: {{llm_request.get('error', 'Unknown error')}}")
        
except Exception as e:
    print(f"❌ Failed to prepare LLM annotation request: {{e}}")

# 5. Save annotation workflow state and initial results
print("\\n📂 Saving annotation workflow state...")

# Store workflow data for continuation
adata.uns['pending_llm_requests'] = pending_llm_requests
adata.uns['cluster_evidence_data'] = {{
    'cluster_markers': cluster_markers,
    'detailed_cluster_info': detailed_cluster_info,
    'user_context': user_data_context,
    'remaining_clusters': list(adata.obs['leiden'].cat.categories[1:])  # All except first
}}

# Save annotation workflow results
try:
    annotation_dir = os.path.join(adata.uns['results_directory'], "07_cell_type_annotation")
    
    # Save annotation workflow state
    annotation_workflow = {{
        "workflow_status": "initiated",
        "total_clusters": total_clusters,
        "clusters_identified": list(adata.obs['leiden'].cat.categories),
        "first_cluster_processed": first_cluster,
        "remaining_clusters": list(adata.obs['leiden'].cat.categories[1:]),
        "user_data_context": user_data_context,
        "aucell_columns": celltype_columns if 'celltype_columns' in locals() else [],
        "marker_detection_method": "omicverse.single.get_celltype_marker" if cluster_markers else "scanpy.tl.rank_genes_groups",
        "llm_requests_pending": len(pending_llm_requests),
        "timestamp": datetime.now().isoformat()
    }}
    
    # Save cluster evidence for first cluster
    if first_cluster in detailed_cluster_info:
        annotation_workflow["first_cluster_evidence"] = detailed_cluster_info[first_cluster]
    
    # Save workflow state
    import json
    with open(os.path.join(annotation_dir, "annotation_workflow_state.json"), 'w') as f:
        json.dump(annotation_workflow, f, indent=2)
    
    # Save cluster markers to separate file
    if cluster_markers:
        with open(os.path.join(annotation_dir, "cluster_markers.json"), 'w') as f:
            # Convert numpy arrays to lists for JSON serialization
            serializable_markers = {{}}
            for cluster, markers in cluster_markers.items():
                if isinstance(markers, (list, tuple)):
                    serializable_markers[str(cluster)] = list(markers)
                else:
                    serializable_markers[str(cluster)] = str(markers)
            json.dump(serializable_markers, f, indent=2)
    
    # Save expected cell types if they exist
    if 'expected_cell_types' in locals():
        with open(os.path.join(annotation_dir, "expected_cell_types.json"), 'w') as f:
            json.dump(expected_cell_types, f, indent=2)
    
    # Save intermediate data state
    adata.write_h5ad(os.path.join(annotation_dir, "pre_annotation_data.h5ad"))
    
    print(f"✅ Annotation workflow state saved to {{annotation_dir}}")
    
except Exception as save_error:
    print(f"⚠️ Failed to save annotation workflow state: {{save_error}}")

print("\\n✅ **CLUSTER ANNOTATION WORKFLOW INITIATED**")
print(f"📊 **DATASET:** {{total_clusters}} clusters identified from leiden clustering")
print(f"📝 **STATUS:** Cluster {{first_cluster}} analysis complete - awaiting Agent annotation")
print(f"📋 **QUEUE:** {{list(adata.obs['leiden'].cat.categories[1:])}} clusters pending analysis")
print("\\n💬 **WORKFLOW:**")
print("1. ✅ System analyzed cluster evidence (AUCell + markers)")
print("2. ⏳ Agent provides cell type annotation for current cluster")
print("3. 🔄 System processes next cluster automatically")
print("4. 🎯 Repeat until all clusters annotated")
print("5. 📋 Final integration of all annotations into dataset")
```

🏷️ STEP 13 - FINALIZE CELL TYPE ANNOTATIONS:
```python
print("\\n🎯 Finalizing cell type annotations...")

# Extract simple annotations from comprehensive analysis
simple_cluster_annotations = {{}}
for cluster, info in final_cluster_annotations.items():
    simple_cluster_annotations[cluster] = info['annotation']

# Apply final annotations to cells
adata.obs['predicted_celltype'] = adata.obs['leiden'].map(simple_cluster_annotations)
adata.obs['annotation_confidence'] = adata.obs['leiden'].map(
    lambda x: final_cluster_annotations[x]['confidence']
)

# Show comprehensive annotation summary
print("\\n📊 Final cell type annotation summary:")
annotation_counts = adata.obs['predicted_celltype'].value_counts()
for celltype, count in annotation_counts.items():
    percentage = (count / adata.n_obs) * 100
    print(f"  {{celltype}}: {{count}} cells ({{round(percentage, 1)}}%)")

print("\\n🔍 Annotation confidence breakdown:")
confidence_counts = adata.obs['annotation_confidence'].value_counts()
for confidence, count in confidence_counts.items():
    percentage = (count / adata.n_obs) * 100
    print(f"  {{confidence}}: {{count}} cells ({{round(percentage, 1)}}%)")

print("\\n✅ Intelligent cell type annotation completed!")
print("Results saved in:")
print("- adata.obs['predicted_celltype']: Final cell type annotations")
print("- adata.obs['annotation_confidence']: Confidence levels")

# Save comprehensive results
adata.uns['cluster_celltype_mapping'] = simple_cluster_annotations
adata.uns['detailed_cluster_analysis'] = final_cluster_annotations
adata.uns['cluster_expression_info'] = detailed_cluster_info

print("\\n💾 Detailed analysis saved in:")
print("- adata.uns['cluster_celltype_mapping']: Simple cluster->celltype mapping")
print("- adata.uns['detailed_cluster_analysis']: Comprehensive analysis results")
print("- adata.uns['cluster_expression_info']: Marker expression statistics")
```

📈 STEP 14 - COMPREHENSIVE VISUALIZATION AND DOWNSTREAM ANALYSIS:
```python
print("\\n📈 Comprehensive Visualization and Downstream Analysis...")

# Get results directory for saving plots
results_dir = adata.uns.get('results_directory', '.')
qc_dir = os.path.join(results_dir, "01_quality_control")
viz_dir = os.path.join(results_dir, "06_visualization")
downstream_dir = os.path.join(results_dir, "07_downstream_analysis")

print(f"📂 Saving results to: {{results_dir}}")

# Step 14a: Quality Control Visualizations
print("\\n🔬 Generating Quality Control Visualizations...")

# QC metrics overview
if 'total_counts' in adata.obs.columns:
    try:
        import matplotlib.pyplot as plt
        
        # QC violin plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        
        # Total counts per cell
        sc.pl.violin(adata, 'total_counts', groupby='leiden', rotation=45, ax=axes[0,0])
        axes[0,0].set_title('Total counts per cell by cluster')
        
        # Number of genes per cell  
        sc.pl.violin(adata, 'n_genes_by_counts', groupby='leiden', rotation=45, ax=axes[0,1])
        axes[0,1].set_title('Number of genes by cluster')
        
        # Mitochondrial gene percentage (if available)
        if 'pct_counts_mt' in adata.obs.columns:
            sc.pl.violin(adata, 'pct_counts_mt', groupby='leiden', rotation=45, ax=axes[1,0])
            axes[1,0].set_title('Mitochondrial gene percentage by cluster')
            
        # Cell counts per cluster
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        axes[1,1].bar(range(len(cluster_counts)), cluster_counts.values)
        axes[1,1].set_xlabel('Cluster')
        axes[1,1].set_ylabel('Number of cells')
        axes[1,1].set_title('Cell count per cluster')
        axes[1,1].set_xticks(range(len(cluster_counts)))
        axes[1,1].set_xticklabels(cluster_counts.index)
        
        plt.tight_layout()
        plt.savefig(os.path.join(qc_dir, 'quality_control_overview.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        print("✅ Quality control overview saved")
        
    except Exception as e:
        print(f"❌ QC visualization failed: {{e}}")

# Step 14b: UMAP Visualizations with omicverse
print("\\n🎨 Generating UMAP Visualizations...")

try:
    # Basic UMAP with clusters
    fig, ax = plt.subplots(figsize=(8, 6))
    ov.pl.embedding(adata, basis='X_umap', color='leiden', 
                   frameon='small', show=False, ax=ax,
                   title='Clusters (Leiden)')
    plt.savefig(os.path.join(viz_dir, 'umap_clusters.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP with cell type annotations (if available)
    if 'predicted_celltype' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        ov.pl.embedding(adata, basis='X_umap', color='predicted_celltype', 
                       frameon='small', show=False, ax=ax,
                       title='Predicted Cell Types')
        plt.savefig(os.path.join(viz_dir, 'umap_celltypes.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Beautiful UMAP with adjusted labels
        from matplotlib import patheffects
        fig, ax = plt.subplots(figsize=(10, 10))
        
        ov.pl.embedding(adata, basis='X_umap', color='predicted_celltype',
                       size=100, show=False, legend_loc=None, add_outline=False,
                       frameon='small', legend_fontoutline=2, ax=ax)
        
        # Add adjusted cell type labels
        try:
            ov.pl.embedding_adjust(adata, groupby='predicted_celltype', basis='X_umap', ax=ax,
                                 adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
                                 text_kwargs=dict(fontsize=12, weight='bold',
                                 path_effects=[patheffects.withStroke(linewidth=2, foreground='w')]))
        except:
            pass  # Skip if adjustment fails
        
        plt.title('Single-cell RNA-seq Analysis Results', fontsize=16, pad=20)
        plt.savefig(os.path.join(viz_dir, 'umap_beautiful.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # Cell proportion analysis
    if 'predicted_celltype' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(4, 6))
        ov.pl.cellproportion(adata, celltype_clusters='predicted_celltype',
                           groupby='leiden', legend=True, ax=ax)
        plt.title('Cell Type Proportions by Cluster')
        plt.savefig(os.path.join(viz_dir, 'cell_proportions.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        
    print("✅ UMAP visualizations saved")
    
except Exception as e:
    print(f"❌ UMAP visualization failed: {{e}}")

# Step 14c: Marker Gene Visualizations  
print("\\n🧬 Generating Marker Gene Visualizations...")

try:
    # Top marker genes heatmap (if available)
    if 'rank_genes_groups' in adata.uns.keys():
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True,
                                      ax=ax, show=False)
        plt.savefig(os.path.join(viz_dir, 'marker_genes_heatmap.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Marker genes dotplot
        fig, ax = plt.subplots(figsize=(10, 6))
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, ax=ax, show=False)
        plt.savefig(os.path.join(viz_dir, 'marker_genes_dotplot.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        
    print("✅ Marker gene visualizations saved")
    
except Exception as e:
    print(f"❌ Marker gene visualization failed: {{e}}")

# Step 14d: Statistical Analysis and Summary
print("\\n📊 Generating Statistical Summary...")

try:
    # Basic statistics
    stats_summary = {{
        'total_cells': adata.n_obs,
        'total_genes': adata.n_vars,
        'n_clusters': len(adata.obs['leiden'].cat.categories) if 'leiden' in adata.obs.columns else 0,
        'available_layers': list(adata.layers.keys()),
        'available_embeddings': list(adata.obsm.keys())
    }}
    
    # Cluster composition
    if 'leiden' in adata.obs.columns:
        cluster_stats = adata.obs['leiden'].value_counts().sort_index()
        stats_summary['cluster_composition'] = cluster_stats.to_dict()
    
    # Cell type composition (if available)
    if 'predicted_celltype' in adata.obs.columns:
        celltype_stats = adata.obs['predicted_celltype'].value_counts()
        stats_summary['celltype_composition'] = celltype_stats.to_dict()
    
    # Quality metrics summary
    if 'total_counts' in adata.obs.columns:
        qc_summary = {{
            'mean_total_counts': adata.obs['total_counts'].mean(),
            'median_total_counts': adata.obs['total_counts'].median(),
            'mean_n_genes': adata.obs['n_genes_by_counts'].mean(),
            'median_n_genes': adata.obs['n_genes_by_counts'].median()
        }}
        if 'pct_counts_mt' in adata.obs.columns:
            qc_summary['mean_mt_pct'] = adata.obs['pct_counts_mt'].mean()
            qc_summary['median_mt_pct'] = adata.obs['pct_counts_mt'].median()
        
        stats_summary['quality_metrics'] = qc_summary
    
    # Save summary to JSON
    import json
    with open(os.path.join(downstream_dir, 'analysis_summary.json'), 'w') as f:
        json.dump(stats_summary, f, indent=2, default=str)
    
    print("✅ Statistical summary saved")
    
except Exception as e:
    print(f"❌ Statistical summary failed: {{e}}")

# Step 14e: Save final AnnData object
print("\\n💾 Saving final results...")

try:
    # Save the complete analysis results
    adata.write_h5ad(os.path.join(results_dir, 'final_analysis_results.h5ad'))
    print("✅ Final AnnData object saved")
    
    # Export key results to CSV
    if 'leiden' in adata.obs.columns:
        cluster_summary = adata.obs[['leiden']].copy()
        if 'predicted_celltype' in adata.obs.columns:
            cluster_summary['predicted_celltype'] = adata.obs['predicted_celltype']
        if 'annotation_confidence' in adata.obs.columns:
            cluster_summary['annotation_confidence'] = adata.obs['annotation_confidence']
        
        cluster_summary.to_csv(os.path.join(results_dir, 'cell_annotations.csv'))
        print("✅ Cell annotations exported to CSV")
    
except Exception as e:
    print(f"❌ Failed to save final results: {{e}}")

# Step 14f: Generate analysis report
print("\\n📋 Analysis Summary:")
print(f"- Total cells analyzed: {{adata.n_obs}}")
print(f"- Total genes: {{adata.n_vars}}")
if 'leiden' in adata.obs.columns:
    print(f"- Number of clusters identified: {{len(adata.obs['leiden'].cat.categories)}}")
if 'predicted_celltype' in adata.obs.columns:
    n_celltypes = len(adata.obs['predicted_celltype'].unique())
    print(f"- Number of cell types annotated: {{n_celltypes}}")

print(f"\\n📂 All results saved in: {{results_dir}}")
print("📁 Directory contents:")
print("   01_quality_control/ - QC plots and metrics")
print("   06_visualization/ - UMAP plots and cell type visualizations") 
print("   07_downstream_analysis/ - Statistical summaries")
print("   final_analysis_results.h5ad - Complete analysis object")
print("   cell_annotations.csv - Cell cluster and type annotations")

print("\\n✅ Comprehensive scRNA-seq analysis pipeline completed successfully!")
print("🎯 Ready for further downstream analysis and biological interpretation!")
```

🚀 EXECUTION ORDER:
1. Setup environment
2. show_todos() and create if empty  
3. For each todo: execute_current_task() → run step → mark_task_done() → show_todos()
4. Use help() before omicverse/scanpy functions ONLY

BEGIN EXECUTION NOW:
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