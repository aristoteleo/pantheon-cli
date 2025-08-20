"""Interactive molecular docking workflow handler - no path input required"""

def generate_interactive_dock_message() -> str:
    """Generate interactive molecular docking workflow message"""
    
    message = """
🧬 Interactive Molecular Docking Pipeline — AutoDock Vina Workflow

🔧 AVAILABLE WORKFLOW TOOLS:
You have access to the Dock_Workflow tool for getting official templates.
💡 Use the tool when you need guidance, official examples, or best practices for a specific workflow step.
You need to ask the user for the path of the data, and then use ls command to check the folder contents,and then proceed with file loading.


⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
- **PERSISTENT STATE**: Python interpreter maintains ALL variables across calls! 
- **MEMORY OPTIMIZATION**: Variables persist! NEVER re-read or re-import data that already exists in memory!
- **SMART VARIABLE CHECKING**: Use `try/except` or `'var' in globals()` to check existence - NO redundant file I/O!
- **EFFICIENCY FIRST**: 
  - Check if data exists before loading: `if 'receptor_data' not in globals()`
  - Use existing results: `if 'docking_results' in locals()`
  - Reuse computed values: `if 'binding_energies' in locals()`
- **ERROR RECOVERY**: If code fails, analyze error and fix - don't reload everything!
- **NO REPETITION**: Each import/load/compute happens ONCE per session unless explicitly needed
- **After each step**: mark_task_done("description"), then show_todos()
- **AUTOMATIC EXECUTION**: Proceed automatically without confirmations; log warnings when needed.

PHASE 0 — SETUP & VALIDATION
1) Data discovery: Use ls command to check current folder contents, then proceed with file loading
2) Environment check will be done automatically within data loading step

PHASE 1 - TODO CREATION (ONCE ONLY)
Execute: current = show_todos()
IF current is EMPTY, create these todos ONCE:
1. "Check Python environment and molecular docking dependencies"
2. "Initialize docking project structure"
3. "Discover and prepare receptor PDB files"
4. "Discover and prepare ligand SDF/MOL2 files"
5. "Calculate binding site centers and configure docking parameters"
6. "Run AutoDock Vina molecular docking"
7. "Analyze protein-ligand interactions and binding results"
8. "Generate docking results summary and visualization"

⚡ AUTOMATIC WORKFLOW MODE:
- Execute each todo task automatically without asking for confirmation
- After successful completion of any step, immediately call mark_task_done("description") and proceed to next
- Continue the workflow seamlessly until all tasks complete or user intervenes

PHASE 2 — INTELLIGENT EXECUTION STRATEGY
🧠 SMART DECISION MAKING:

**ASSESS CURRENT SITUATION FIRST:**
- What data files do you have in the current directory?
- What analysis steps are completed?
- What specific guidance do you need?

**CHOOSE YOUR APPROACH:**

**Option A - Use Template Tool (Recommended for new users or complex steps):**
- Call Dock_Workflow(workflow_type="<type>") to get official template
- Study the returned guidance and code patterns
- Adapt the template to your specific data
- Execute the adapted code

**Option B - Direct Implementation (For experienced users with clear requirements):**
- Directly write and execute code based on vina documentation
- Use help() functions to check parameters
- Follow established best practices

**WHEN TO USE TEMPLATE TOOL:**
✅ When you need official guidance for a workflow step
✅ When you want to see best practices and parameter examples  
✅ When you're unsure about the correct approach
✅ When you want standardized, tested code patterns

**WHEN DIRECT IMPLEMENTATION IS OK:**
✅ When you have clear requirements and know the approach
✅ When adapting previous successful code
✅ When making minor parameter adjustments

🧠 **RESULT ANALYSIS REQUIREMENT:**
After executing any code:
1. **Analyze the output** - Don't just print and move on!
2. **Interpret the results** - What do the binding energies, poses, and warnings mean?
3. **Check for issues** - Are there data quality problems or unexpected patterns?
4. **Make decisions** - Should parameters be adjusted based on what you observed?
5. **Document findings** - Save key insights to results directory
6. **Proceed intelligently** - Use results to inform next steps

🏷️ STEP 1 - DEPENDENCY CHECK:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on dependency installation?
- Are you familiar with meeko and vina packages?
- Do you have the required molecular docking tools?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="check_dependencies")
```
Then adapt the returned template to your environment.

**Path B - Direct Implementation:**
Direct implementation using pip if you know the approach:
```python
import subprocess
import sys
# Check and install required packages
result = subprocess.run([sys.executable, '-m', 'pip', 'install', 'meeko', 'vina'], capture_output=True, text=True)
```

🎯 **GOAL:** Ensure all molecular docking dependencies are installed

🏷️ STEP 2 - PROJECT INITIALIZATION:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on project structure setup?
- Are you familiar with docking project organization?
- Do you have specific directory requirements?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="init")
```
Then adapt the returned template to your workspace.

**Path B - Direct Implementation:**
Direct implementation using os.makedirs if you know the structure:
```python
import os
directories = ['receptors', 'ligands', 'prepare', 'output', 'scripts', 'analysis']
for dir_name in directories:
    os.makedirs(dir_name, exist_ok=True)
```

🎯 **GOAL:** Create organized project structure for docking workflow

🏷️ STEP 3 - RECEPTOR PREPARATION:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on receptor preparation workflow?
- Are you familiar with meeko and PDBQT conversion?
- Do you have the required PDB files?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="prepare_receptor")
```
Then adapt the returned template to your receptor files.

**Path B - Direct Implementation:**
Direct implementation using meeko if you know the approach:
```python
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
# Process receptor files
```

🎯 **GOAL:** Convert PDB receptor files to PDBQT format for docking

🏷️ STEP 4 - LIGAND PREPARATION:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on ligand preparation workflow?
- Are you familiar with SDF/MOL2 to PDBQT conversion?
- Do you have the required ligand files?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="prepare_ligand")
```
Then adapt the returned template to your ligand files.

**Path B - Direct Implementation:**
Direct implementation using meeko if you know the approach:
```python
from meeko import MoleculePreparation
from rdkit import Chem
# Process ligand files
```

🎯 **GOAL:** Convert SDF/MOL2 ligand files to PDBQT format for docking

🏷️ STEP 5 - DOCKING EXECUTION:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on AutoDock Vina docking?
- Are you familiar with binding site definition and parameters?
- Do you have prepared PDBQT files?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="docking_vina")
```
Then adapt the returned template to your specific receptor-ligand pairs.

**Path B - Direct Implementation:**
Direct implementation using vina if you know the approach:
```python
from vina import Vina
# Setup docking parameters and execute
```

🎯 **GOAL:** Perform molecular docking and generate binding poses

🏷️ STEP 6 - BATCH DOCKING (Optional):
🤔 **ASSESS YOUR NEEDS:**
- Do you have multiple receptor-ligand pairs?
- Do you need guidance on batch processing?
- Are you familiar with automated docking workflows?

# Get docking template commands (generates Python scripts)
dock_commands = dock.Dock_Workflow("batch_docking")
# Execute: creates scripts/batch_docking.py and runs it
# Execute: creates scripts/batch_docking.py and runs it
# Execute: creates scripts/batch_docking.py and runs it
Don't run code using run_python_code, use bash command to run the script. This is important.
You don't need to write any subprocess.run or os.system in the script. follow the dock.Dock_Workflow("batch_docking")'s example to write the script.
# Python script uses: from vina import Vina (NOT command line vina)
# Check results: grep "REMARK VINA RESULT" output/*.pdbqt

🎯 **GOAL:** Efficiently dock multiple protein-ligand combinations

🏷️ STEP 7 - INTERACTION ANALYSIS:
🤔 **ASSESS YOUR NEEDS:**
- Do you need guidance on interaction analysis?
- Are you familiar with PDBQT result interpretation?
- Do you want to analyze binding energies and poses?

🛠️ **CHOOSE YOUR PATH:**

**Path A - Get Official Template:**
```
Dock_Workflow(workflow_type="analyze_interactions")
```
Then adapt the returned template to your results.

**Path B - Direct Implementation:**
Direct implementation using file parsing if you know the approach:
```python
import glob
import pandas as pd
# Parse PDBQT results and analyze binding energies
```

🎯 **GOAL:** Analyze protein-ligand interactions and binding affinities

🔧 **AVAILABLE GUIDANCE TOOLS:**

**MOLECULAR DOCKING TEMPLATE ENGINE (Optional):**
- `Dock_Workflow(workflow_type="check_dependencies")` - Get dependency installation template
- `Dock_Workflow(workflow_type="init")` - Get project initialization template
- `Dock_Workflow(workflow_type="prepare_receptor")` - Get receptor preparation template
- `Dock_Workflow(workflow_type="prepare_ligand")` - Get ligand preparation template
- `Dock_Workflow(workflow_type="docking_vina")` - Get single docking template
- `Dock_Workflow(workflow_type="batch_docking")` - Get batch docking template
- `Dock_Workflow(workflow_type="analyze_interactions")` - Get interaction analysis template

**FLEXIBLE EXECUTION STRATEGY:**
1. **ASSESS**: Evaluate your current needs and knowledge level
2. **CHOOSE**: Use templates for guidance OR implement directly
3. **ADAPT**: Modify any template code to fit your specific data
4. **EXECUTE**: Run your adapted code with run_python_code
5. **PROGRESS**: Mark tasks complete and continue workflow

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
    print(f"Using existing receptor_data: {len(receptor_files)} files")
except NameError:
    receptor_files = glob.glob("*.pdb")

# Good - Check computed results
if 'docking_results' not in locals():
    # Run docking
    pass
else:
    print("Docking already completed")

# Bad - Redundant file I/O
receptor_files = glob.glob("*.pdb")  # Don't do this if receptor_files exists!
```

**Remember:** Maintain persistent state, avoid redundant operations, and mark tasks complete with mark_task_done()!
"""
    
    return message