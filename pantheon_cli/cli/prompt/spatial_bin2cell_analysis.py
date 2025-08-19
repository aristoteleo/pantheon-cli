"""Spatial transcriptomics analysis mode handler with bin2cell integration"""

from pathlib import Path
from typing import Optional

def generate_spatial_analysis_message(folder_path: Optional[str] = None, source_image_path: Optional[str] = None) -> str:
    """Generate spatial transcriptomics analysis message using spatial toolset with bin2cell"""
    
    if folder_path:
        data_path = Path(folder_path).resolve()
        
        # Determine if it's a file or folder
        if data_path.is_file():
            target_description = f"Target data file: {data_path}"
            path_instruction = f'Always use the provided data_path: "{data_path}" for data loading and analysis.'
        else:
            target_description = f"Target folder: {data_path}"
            path_instruction = f'Always use the provided folder_path: "{data_path}" to scan for spatial transcriptomics data files.'
        
        image_instruction = ""
        if source_image_path:
            image_instruction = f'\nSource H&E image: "{source_image_path}" - Required for Visium HD bin2cell analysis.'
        
        message = f"""
🔬 Spatial Transcriptomics Analysis Pipeline with bin2cell Integration
{target_description}{image_instruction}

⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
- **PERSISTENT STATE**: Python interpreter maintains ALL variables across calls! 
- **MEMORY OPTIMIZATION**: Variables persist! NEVER re-read or re-import data that already exists in memory!
- **SMART VARIABLE CHECKING**: Use `try/except` or `'var' in globals()` to check existence - NO redundant file I/O!
- **EFFICIENCY FIRST**: 
  - Check if adata exists before loading: `if 'adata' not in globals()`
  - Use existing results: `if 'labels_he' in adata.obs`
  - Reuse computed values: `if 'stardist_done' in locals()`
- **ERROR RECOVERY**: If code fails, analyze error and fix - don't reload everything!
- **NO REPETITION**: Each import/load/compute happens ONCE per session unless explicitly needed
- **After each step**: mark_task_done("description"), then show_todos()
- **AUTOMATIC EXECUTION**: Proceed automatically without confirmations; log warnings when needed.

{path_instruction}

PHASE 0 — SETUP & VALIDATION
1) Environment check: Use pip list to check scanpy and bin2cell packages
2) Data discovery: Use ls command to check folder contents, then proceed with file loading

PHASE 1 — TODO CREATION (ONCE ONLY)
Execute: current = show_todos()
IF current is EMPTY, create these todos ONCE:
1. "Check Python environment and load spatial data"
2. "Inspect spatial data structure and determine bin2cell compatibility"
3. "Apply basic filtering to remove low-quality spots"
4. "Run destriping to correct variable bin dimensions (Visium HD)"
5. "Check and adjust spatial alignment between coordinates and H&E image"
6. "Perform H&E image segmentation using StarDist"
7. "Perform gene expression image segmentation using StarDist"
8. "Expand nuclear labels to capture full cell bodies"
9. "Convert bins to cells using joint segmentation results"

⚡ AUTOMATIC WORKFLOW MODE:
- Execute each todo task automatically without asking for confirmation
- After successful completion of any step, immediately call mark_task_done("description") and proceed to next
- Continue the workflow seamlessly until all tasks complete or user intervenes

PHASE 2 — ADAPTIVE EXECUTION WORKFLOW

⚠️ CRITICAL EXECUTION STRATEGY:
When you call spatial.run_spatial_workflow(), it returns guidance, explanations, and example Python code using toolset function run_python_code.
You MUST:
1. **Read and analyze** the entire returned content carefully
2. **Understand the bin2cell methodology** and spatial analysis logic described
3. **Adapt the provided code** to your current spatial data situation (adata shape, coordinate systems, etc.)
4. **Modify parameters** based on your actual tissue type and data characteristics
5. **Execute the adapted code** - NOT the original code directly
6. **Handle segmentation errors** by adjusting parameters based on the guidance provided

🧠 **SPATIAL RESULT ANALYSIS REQUIREMENT:**
After executing any spatial code:
1. **Analyze the output** - Don't just print and move on!
2. **Interpret spatial patterns** - What do the segmentation results, spatial distributions mean?
3. **Check alignment** - Are spatial coordinates properly aligned with H&E images?
4. **Validate segmentation** - Are cells properly detected and boundaries reasonable?
5. **Assess bin-to-cell conversion** - Does the aggregation make biological sense?
6. **Document spatial findings** - Save key insights and visualizations
7. **Proceed with spatial context** - Use spatial information to inform next steps

The returned content serves as GUIDANCE and TEMPLATES for spatial analysis, not direct execution scripts.

📊 STEP 1 - SPATIAL DATA LOADING, INSPECTION & PROJECT SETUP:
```python
# STEP 1A: Check Python environment with pip list
print("\\n🔍 Checking Python environment...")
import subprocess
import sys

# Check for required packages
required_packages = ['scanpy', 'bin2cell', 'pandas', 'numpy']
result = subprocess.run([sys.executable, '-m', 'pip', 'list'], capture_output=True, text=True)
installed_packages = result.stdout.lower()

missing_packages = []
for pkg in required_packages:
    if pkg.lower() not in installed_packages:
        missing_packages.append(pkg)

if missing_packages:
    print(f"❌ Missing packages: {{', '.join(missing_packages)}}")
    print("Installing missing packages...")
    for pkg in missing_packages:
        subprocess.run([sys.executable, '-m', 'pip', 'install', pkg])
    print("✅ Package installation complete")
else:
    print("✅ All required packages are installed")

# STEP 1B: EFFICIENT SPATIAL DATA LOADING - Check memory first!
try:
    # Check if adata exists and is valid
    print(f"Using existing spatial adata: {{adata.shape}} (n_obs, n_var)")
    data_already_loaded = True
except NameError:
    # Only load if not in memory
    print("Loading spatial data for the first time...")
    import bin2cell as b2c
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import os
    
    # Create directory for stardist files
    os.makedirs("stardist", exist_ok=True)
    
    # Load spatial data based on detected format
    # For Visium HD with bin2cell (requires both path and source_image_path)
    {"if source_image_path:" if source_image_path else ""}
    path = "{data_path}"  # Path to square_002um folder
    source_image_path = "{source_image_path}"  # BTF or TIF H&E image
    adata = b2c.read_visium(path, source_image_path=source_image_path)
    {"else:" if source_image_path else ""}
    {"# For non-Visium HD data" if source_image_path else ""}
    {"adata = sc.read_h5ad('{data_path}')  # .h5ad files" if source_image_path else ""}
    adata.var_names_make_unique()
    print(f"Loaded spatial data: {{adata.shape}} (n_obs, n_var)")
    print(f"Available obs columns: {{list(adata.obs.columns)}}")
    print(f"Spatial keys: {{list(adata.obsm.keys())}}")
    data_already_loaded = False

# Only import libraries once
if 'b2c' not in globals():
    import bin2cell as b2c
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import os
    print("Spatial analysis libraries imported")

# Create structured output directory for spatial analysis (only if not exists)
print("\\n📁 Setting up spatial analysis project structure...")
try:
    # Check if we already have a results directory
    if 'results_dir' not in globals():
        from datetime import datetime
        
        # Create main results directory with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = f"spatial_analysis_results_{{timestamp}}"
        os.makedirs(results_dir, exist_ok=True)
    else:
        print(f"Using existing results directory: {{results_dir}}")
    
    # Create subdirectories for spatial analysis components
    subdirs = [
        "01_data_loading",
        "02_basic_filtering", 
        "03_destriping",
        "04_alignment_check",
        "05_he_segmentation",
        "06_gex_segmentation", 
        "07_label_expansion",
        "08_bin_to_cell",
        "stardist",  # StarDist segmentation files
        "images",    # H&E and visualization images
    ]
    
    for subdir in subdirs:
        os.makedirs(os.path.join(results_dir, subdir), exist_ok=True)
    
    # Store results directory in adata for later use (if adata exists)
    if 'adata' in globals():
        adata.uns['results_directory'] = results_dir
        print(f"✅ Spatial project structure ready: {{results_dir}}")
    else:
        print(f"✅ Spatial project structure created: {{results_dir}} (will link to adata after loading)")
    
except Exception as e:
    print(f"❌ Failed to create spatial project structure: {{e}}")

# Initial directory inspection with ls command
print("\\n🔍 Checking directory contents...")
import subprocess
import os

# Check if path is file or directory
if os.path.isfile("{data_path}"):
    # If it's a file, show the file and its parent directory
    parent_dir = os.path.dirname("{data_path}")
    print(f"Target file: {data_path}")
    if parent_dir:
        result = subprocess.run(['ls', '-la', parent_dir], capture_output=True, text=True)
        print(f"Parent directory contents:\\n{{result.stdout}}")
else:
    # If it's a directory, show its contents
    result = subprocess.run(['ls', '-la', "{data_path}"], capture_output=True, text=True)
    print(f"Directory contents:\\n{{result.stdout}}")
    
    # Check for common spatial data files
    spatial_files = []
    if os.path.exists(os.path.join("{data_path}", "filtered_feature_bc_matrix.h5")):
        spatial_files.append("filtered_feature_bc_matrix.h5")
    if os.path.exists(os.path.join("{data_path}", "spatial")):
        spatial_files.append("spatial/ directory")
    if os.path.exists(os.path.join("{data_path}", "tissue_positions.parquet")):
        spatial_files.append("tissue_positions.parquet")
    
    if spatial_files:
        print(f"✅ Detected spatial data files: {{', '.join(spatial_files)}}")
    else:
        print("⚠️ No standard spatial data files detected")

print("✅ Directory inspection complete")
```

🧹 STEP 2 - BASIC FILTERING:
Get basic filtering guidance and adapt the code to your spatial data:
spatial.run_spatial_workflow(workflow_type="basic_filter")
Then analyze the returned guidance and implement adapted filtering code for spatial spots.
**CRITICAL**: Filter genes in at least 3 spots and spots with at least 1 count. VisiumHD data is extremely sparse initially.
**TUTORIAL EXAMPLE**: sc.pp.filter_genes(adata, min_cells=3); sc.pp.filter_cells(adata, min_counts=1)

🎯 STEP 3 - DESTRIPING (Visium HD specific):
Get destriping guidance to correct variable bin dimensions:
spatial.run_spatial_workflow(workflow_type="destripe")
Then analyze the returned guidance and implement destriping correction for Visium HD data.
**CRITICAL**: Visium HD has ~10% variability in bin width/height causing striped patterns. Destriping identifies 0.99 quantile for each row/column and rescales.
**TUTORIAL FUNCTION**: b2c.destripe(adata) creates 'n_counts_adjusted' and 'destripe_factor' columns.

🔧 STEP 4 - ALIGNMENT CHECK:
Get alignment checking guidance and methodology:
spatial.run_spatial_workflow(workflow_type="alignment_check")
Then adapt and implement alignment verification between spatial coordinates and H&E image.
**CRITICAL**: Check fine-grain alignment between spatial coordinates and H&E image. Use test region (e.g., rows 2050-2250, cols 1350-1550).
**TUTORIAL ADJUSTMENT**: Apply coordinate shifts like up_down=-5, right_left=+0 to correct misalignment before segmentation.

🎨 STEP 5 - H&E SEGMENTATION:
Get H&E segmentation guidance using StarDist:
spatial.run_spatial_workflow(workflow_type="he_segmentation")
Then adapt and implement H&E image segmentation code based on your tissue type and image quality.
**CRITICAL**: Use mpp=0.3 for most tissues, prob_thresh=0.1 for more lenient nuclear detection. Model: '2D_versatile_he'.
**TUTORIAL WORKFLOW**: 1) b2c.scaled_he_image(), 2) b2c.stardist(), 3) b2c.insert_labels() with basis='spatial'

🧬 STEP 6 - GENE EXPRESSION SEGMENTATION:
Get gene expression segmentation guidance:
spatial.run_spatial_workflow(workflow_type="gex_segmentation")
Then adapt and implement GEX-based segmentation code for your spatial expression patterns.
**CRITICAL**: Use 'n_counts_adjusted' with sigma=5 smoothing. Model: '2D_versatile_fluo', prob_thresh=0.01, nms_thresh=0.1.
**TUTORIAL WORKFLOW**: 1) b2c.grid_image() with array coordinates, 2) b2c.stardist(), 3) b2c.insert_labels() with basis='array'

📏 STEP 7 - EXPAND LABELS:
Get label expansion guidance and methodology:
spatial.run_spatial_workflow(workflow_type="expand_labels")
Then adapt and implement label expansion to capture full cell bodies beyond just nuclei.
**CRITICAL**: Use max_bin_distance=4 for most tissues. Algorithm uses PCA space for equidistant bin assignment.
**TUTORIAL WORKFLOW**: 1) b2c.expand_labels() for H&E, 2) b2c.salvage_secondary_labels() to combine H&E+GEX segmentation

🔄 STEP 8 - BIN TO CELL CONVERSION:
Get bin-to-cell conversion guidance:
spatial.run_spatial_workflow(workflow_type="bin_to_cell")
Then adapt and implement the conversion from spatial bins to reconstructed cells.
**CRITICAL**: Use 'labels_joint' for best results. Preserves spatial coordinates as centroids, aggregates expression.
**TUTORIAL RESULTS**: Expect major reduction in observations (e.g., 8M bins → 257K cells) with preserved spatial structure.

✅ BIN2CELL WORKFLOW COMPLETE:
```python
print("\\n🎉 Bin2cell workflow completed successfully!")
print(f"- Original bins: {{adata.n_obs}}")
print(f"- Reconstructed cells: {{cdata.n_obs}}")
print(f"- Data saved: spatial_bins.h5ad and spatial_cells.h5ad")
print("\\n📊 Results available for downstream analysis:")
print("- adata: Original bin-level data with segmentation labels")
print("- cdata: Reconstructed cell-level data ready for scRNA-seq analysis")
```

🔧 **AVAILABLE SPATIAL TOOLSET FUNCTIONS:**

**BIN2CELL WORKFLOW ENGINE:**
- `spatial.run_spatial_workflow(workflow_type="basic_filter")` - Basic spot filtering
- `spatial.run_spatial_workflow(workflow_type="destripe")` - Visium HD destriping correction
- `spatial.run_spatial_workflow(workflow_type="alignment_check")` - Spatial alignment verification
- `spatial.run_spatial_workflow(workflow_type="he_segmentation")` - H&E segmentation with StarDist
- `spatial.run_spatial_workflow(workflow_type="gex_segmentation")` - Gene expression segmentation
- `spatial.run_spatial_workflow(workflow_type="expand_labels")` - Label expansion for full cell capture
- `spatial.run_spatial_workflow(workflow_type="bin_to_cell")` - Bin to cell conversion

**BIN2CELL EXECUTION STRATEGY:**
1. Load Visium HD spatial data with bin2cell integration
2. Apply destriping correction for variable bin dimensions
3. Verify and adjust spatial-image alignment
4. Perform dual segmentation (H&E + gene expression) with StarDist
5. Expand nuclear labels to capture full cell bodies
6. Convert high-resolution bins to reconstructed cells
7. Save both bin-level and cell-level data for downstream analysis

**BIN2CELL EFFICIENCY PRINCIPLES:**
1. **CHECK SEGMENTATION STATE**: Verify if segmentation already done before re-running StarDist
2. **PRESERVE SPATIAL INFO**: Maintain spatial coordinates through all transformations
3. **VALIDATE ALIGNMENT**: Always check spatial-image alignment before segmentation
4. **BIN2CELL OPTIMIZATION**: Use appropriate parameters for tissue type and resolution
5. **DESTRIPING FIRST**: Apply destriping before segmentation for Visium HD data

Example spatial patterns:
```python
# Good - Check segmentation status
if 'labels_he' not in adata.obs:
    # Run H&E segmentation
    spatial.run_spatial_workflow(workflow_type="he_segmentation")
else:
    print("H&E segmentation already completed")

# Good - Verify spatial coordinates
if 'spatial' in adata.obsm:
    print(f"Spatial coordinates available: {{adata.obsm['spatial'].shape}}")
else:
    print("Warning: No spatial coordinates found")

# Good - Check bin2cell compatibility
if 'array_row' in adata.obs and 'array_col' in adata.obs:
    print("Data is bin2cell compatible (Visium HD format)")
else:
    print("Standard spatial format detected")
```

**Remember:** bin2cell converts Visium HD 2μm bins to reconstructed cells through image segmentation. The final cdata object is ready for standard scRNA-seq analysis!
"""
        
    else:
        message = """
I need help with spatial transcriptomics analysis using your specialized bin2cell toolsets.

You have access to comprehensive spatial transcriptomics and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🔬 BIN2CELL TOOLSET:

**BIN2CELL WORKFLOW ENGINE:**
- `spatial.run_spatial_workflow(workflow_type="basic_filter")` - Basic spot filtering
- `spatial.run_spatial_workflow(workflow_type="destripe")` - Visium HD destriping correction
- `spatial.run_spatial_workflow(workflow_type="alignment_check")` - Spatial alignment verification
- `spatial.run_spatial_workflow(workflow_type="he_segmentation")` - H&E image segmentation with StarDist
- `spatial.run_spatial_workflow(workflow_type="gex_segmentation")` - Gene expression segmentation
- `spatial.run_spatial_workflow(workflow_type="expand_labels")` - Label expansion for full cell capture
- `spatial.run_spatial_workflow(workflow_type="bin_to_cell")` - Bin to cell conversion

**BIN2CELL INTEGRATION:**
- Full Visium HD support with 2μm resolution
- Destriping correction for variable bin dimensions  
- H&E and gene expression dual segmentation with StarDist
- Nuclear label expansion using PCA-based similarity
- Joint segmentation combining H&E (primary) + GEX (secondary)
- Automated bin-to-cell aggregation with centroid coordinates
- Quality metrics: bin_count, segmentation source tracking

**BIN2CELL EXECUTION STRATEGY:**
1. Load Visium HD data with bin2cell compatibility check
2. Apply destriping correction for variable bin dimensions
3. Verify and adjust spatial-image alignment
4. Perform dual segmentation (H&E + gene expression) for cell detection
5. Expand nuclear labels to capture full cell bodies
6. Convert high-resolution bins to reconstructed cells
7. Save bin-level and cell-level data for downstream analysis

**GUIDANCE:**
- spatial.suggest_next_step() - Smart spatial analysis recommendations

Please start by adding a todo for your bin2cell analysis task, then use pip list to check packages and ls to examine directory contents before proceeding!

**Key Features:**
• Visium HD bin2cell integration
• H&E and gene expression dual segmentation  
• StarDist-based nuclear detection
• Destriping for variable bin correction
• Label expansion for full cell capture
• Seamless bin-to-cell conversion
• Output ready for scRNA-seq analysis"""
    
    return message