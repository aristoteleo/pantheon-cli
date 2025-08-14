"""Single-cell RNA-seq analysis mode handler with omicverse integration"""
#TODO: This is a placeholder for the subtype analysis.
#TODO: It is not implemented yet.
#TODO: It is a placeholder for the subtype analysis.
#TODO: It is a placeholder for the subtype analysis.

from pathlib import Path
from typing import Optional

def generate_scrna_subtype_analysis_message(folder_path: Optional[str] = None) -> str:
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
8. "Perform cell type annotation"
9. "Conduct downstream analysis"
10. "Generate analysis report"

PHASE 2 — ADAPTIVE EXECUTION WORKFLOW

📊 STEP 1 - DATA LOADING & INSPECTION:
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

# Inspect current state
print("\\n🔍 Data State:")
print(f"- Shape: {{adata.shape}}")
print(f"- Layers: {{list(adata.layers.keys())}}")
print(f"- Embeddings: {{list(adata.obsm.keys())}}")
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
        ov.pp.qc(adata, 
                mode='seurat',           # 'seurat' or 'mads'
                min_cells=3, 
                min_genes=200,
                mt_startswith='MT-',     # Mitochondrial gene prefix
                tresh={{'mito_perc': 20, 'nUMIs': 1000, 'detected_genes': 2500}})
        print("✅ QC completed successfully")
    except Exception as e:
        print(f"❌ QC failed: {{e}}")
        # Retry with more lenient thresholds
        try:
            ov.pp.qc(adata, mode='seurat', min_cells=1, min_genes=100,
                    mt_startswith='MT-', tresh={{'mito_perc': 30}})
            print("✅ QC completed with relaxed parameters")
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
    except Exception as e:
        print(f"❌ Preprocessing failed: {{e}}")
        # Retry with simpler mode
        try:
            adata = ov.pp.preprocess(adata, mode='shiftlog|pearson', 
                                   target_sum=1e4, n_HVGs=1500)
            print("✅ Preprocessing completed with reduced parameters")
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
    except Exception as e:
        print(f"❌ PCA failed: {{e}}")
        # Retry with fewer components
        try:
            ov.pp.pca(adata, n_pcs=30, layer='scaled', inplace=True)
            print("✅ PCA completed with 30 components")
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

🏷️ STEP 7 - CELL TYPE ANNOTATION:
```python
print("\\n🏷️ Cell Type Annotation...")

# Step 7a: Find marker genes for each cluster
if 'rank_genes_groups' not in adata.uns.keys():
    print("Finding marker genes...")
    # Use scanpy rank_genes_groups (no help() needed - scanpy function)
    sc.tl.rank_genes_groups(adata, 
                           groupby='leiden',        # Cluster column
                           method='wilcoxon',       # Statistical test
                           key_added='rank_genes_groups',
                           n_genes=25)              # Number of top genes
    print("✅ Marker genes computed")

# Step 7b: Extract marker genes with omicverse
print("\\nExtracting top marker genes...")

# MANDATORY: Check help first for omicverse function
help(ov.single.get_celltype_marker)

try:
    # Extract marker genes using omicverse with actual parameters
    marker_dict = ov.single.get_celltype_marker(adata,
                                               clustertype='leiden',        # Cluster column
                                               log2fc_min=2,              # Min log2 fold change
                                               scores_type='scores',       # 'scores' or 'logfoldchanges'
                                               pval_cutoff=0.05,           # P-value threshold
                                               rank=False,                 # Whether to rank
                                               key='rank_genes_groups',    # Rank genes key
                                               method='wilcoxon',          # Test method
                                               foldchange=None,            # Fold change threshold
                                               topgenenumber=10,           # Top genes per cluster
                                               unique=True,                # Remove duplicates within cluster
                                               global_unique=False)        # Remove duplicates globally
    print("✅ Marker gene extraction completed")
    
    # Display marker genes for each cluster
    for cluster, genes in marker_dict.items():
        print(f"Cluster {{cluster}}: {{', '.join(genes[:5])}}...")  # Show top 5 genes
        
except Exception as e:
    print(f"❌ Marker gene extraction failed: {{e}}")
    # Fallback: manually extract from rank_genes_groups
    marker_dict = {{}}
    for cluster in adata.obs['leiden'].cat.categories:
        marker_dict[cluster] = adata.uns['rank_genes_groups']['names'][cluster][:10].tolist()
    print("✅ Fallback marker gene extraction completed")

# Step 7c: Cell ontology mapping (optional)
print("\\nInitializing Cell Ontology Mapper...")

try:
    # MANDATORY: Check help first for omicverse function
    help(ov.single.CellOntologyMapper)
    
    # Initialize Cell Ontology Mapper with actual parameters
    cell_mapper = ov.single.CellOntologyMapper(
        cl_obo_file=None,                    # Path to CL ontology file
        embeddings_path=None,               # Path to pre-computed embeddings
        model_name="all-mpnet-base-v2",    # Sentence transformer model
        local_model_dir=None                # Local model directory
    )
    
    # Setup for cell type mapping if needed
    # cell_mapper.setup_llm_expansion(
    #     api_type="openai",
    #     tissue_context="your_tissue",  # e.g., "lung", "brain"
    #     species="human"                 # or "mouse"
    # )
    
    print("✅ Cell Ontology Mapper initialized")
    
except Exception as e:
    print(f"❌ Cell Ontology Mapper initialization failed: {{e}}")
    print("Continuing with manual annotation based on marker genes...")

# Step 7d: Add cluster annotations to adata (placeholder)
adata.obs['cell_type_annotation'] = 'Unknown'  # Initialize
print("✅ Cell type annotation column initialized")
print("📝 Manual annotation required based on marker genes and biological knowledge")
```

📈 STEP 8 - DOWNSTREAM ANALYSIS:
```python
print("\\n📈 Downstream Analysis...")

# Step 8a: Basic statistics and summary
print("\\n📊 Analysis Summary:")
print(f"- Total cells: {{adata.n_obs}}")
print(f"- Total genes: {{adata.n_vars}}")
print(f"- Number of clusters: {{len(adata.obs['leiden'].cat.categories)}}")
print(f"- Available layers: {{list(adata.layers.keys())}}")
print(f"- Available embeddings: {{list(adata.obsm.keys())}}")

# Step 8b: Cluster composition analysis
if 'leiden' in adata.obs.columns:
    print("\\n🔢 Cluster composition:")
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    for cluster, count in cluster_counts.items():
        percentage = (count / adata.n_obs) * 100
        print(f"  Cluster {{cluster}}: {{count}} cells ({{percentage:.1f}}%)")

# Step 8c: Quality metrics per cluster
if 'leiden' in adata.obs.columns and 'total_counts' in adata.obs.columns:
    print("\\n📈 Quality metrics by cluster:")
    qc_metrics = ['total_counts', 'n_genes_by_counts']
    if 'pct_counts_mt' in adata.obs.columns:
        qc_metrics.append('pct_counts_mt')
    
    cluster_qc = adata.obs.groupby('leiden')[qc_metrics].mean()
    print(cluster_qc)

# Step 8d: Save intermediate results
print("\\n💾 Analysis state saved to adata object")
print("Available for further analysis:")
print("- adata.obs['leiden']: Cluster assignments")
print("- adata.obsm['X_pca']: PCA coordinates")
if 'X_umap' in adata.obsm.keys():
    print("- adata.obsm['X_umap']: UMAP coordinates")
if 'rank_genes_groups' in adata.uns.keys():
    print("- adata.uns['rank_genes_groups']: Differential expression results")

# Step 8e: Pertpy analysis (optional if pertpy available)
try:
    import pertpy as pt
    print("\\n🧪 Pertpy available for perturbation analysis")
    # Add pertpy-specific analyses here if needed
    # help(pt.tl.some_function)  # Use help() for pertpy functions
except ImportError:
    print("\\n📝 Pertpy not available - basic analysis completed")

print("\\n✅ scRNA-seq analysis pipeline completed successfully!")
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