"""Simplified molecular docking mode handler"""

from pathlib import Path
from typing import Optional

def generate_dock_analysis_message(folder_path: Optional[str] = None) -> str:
    """Generate molecular docking analysis message using dock toolset"""
    
    if folder_path:
        folder_path = Path(folder_path).resolve()
        
        message = f"""
🧬 Molecular Docking Pipeline — AutoDock Vina Workflow (NEW VERSION)

Target folder: {folder_path}

⚡ **IMPORTANT**: This is the NEW UPDATED dock toolset that uses Python vina API (not command line).
The docking workflows generate Python scripts and execute them with bash, NOT direct vina commands.

You have access to the molecular docking toolset with workflow-based commands and TodoList management.

GLOBAL RULES
- Always use the provided folder_path: "{folder_path}" in all phases.
- Idempotent behavior: NEVER create duplicate todos. Only create if the list is EMPTY.
- Do not ask the user for confirmations; proceed automatically and log warnings when needed.
- After each concrete tool completes successfully, call mark_task_done("what was completed"), then show_todos().

PHASE 0 — DEPENDENCY CHECK & PROJECT SETUP
1) Use workflow commands for setup:
   - dock.Dock_Workflow("check_dependencies") - Install meeko, vina
   - dock.Dock_Workflow("init") - Create docking project structure

PHASE 1 — TODO CREATION (STRICT DE-DUP)
Mandatory order:
  a) current = show_todos()
  b) Analyze folder structure and files in "{folder_path}"
Creation rule (single condition):
  • If current is EMPTY → create ONCE the following todos:
      0. "Check and install molecular docking dependencies"
      1. "Initialize docking project structure"
      2. "Prepare receptor PDB files to PDBQT format"
      3. "Prepare ligand SDF files to PDBQT format"
      4. "Calculate binding site centers and grid parameters"
      5. "Run AutoDock Vina molecular docking"
      6. "Analyze protein-ligand interactions from PDBQT files"
      7. "Generate docking results summary and visualization"
  • Else → DO NOT create anything. Work with the existing todos.

PHASE 2 — EXECUTE WITH TODO TRACKING (LOOP)

⚠️ CRITICAL EXECUTION STRATEGY - NEW VERSION:
When you call dock.Dock_Workflow(), they return bash command templates that GENERATE AND EXECUTE Python scripts.
You MUST:
1. **Read and analyze** the entire returned bash command content carefully  
2. **Understand that the docking uses Python vina API**, not command line vina
3. **Adapt the provided commands** to your current data situation (file paths, receptor names, ligand names, etc.)
4. **Execute the bash commands** which will create Python scripts in scripts/ folder and run them
5. **Handle errors** by checking Python vina package installation and script execution
6. **Analyze results** - Check output files, binding energies, and PDBQT poses

🧠 **RESULT ANALYSIS REQUIREMENT:**
After executing any bash commands:
1. **Analyze binding energies** - Check for reasonable docking scores
2. **Verify output files** - Ensure PDBQT, PDB files were created
3. **Examine interaction maps** - Review PLIP analysis results
4. **Validate poses** - Are the binding poses chemically reasonable?
5. **Check convergence** - Did the docking converge properly?
6. **Document findings** - Note binding affinities and key interactions

For each current task:
  1) hint = execute_current_task()   # obtain guidance for the next action
  2) Get bash command templates using appropriate dock workflow:
     
     DOCKING WORKFLOWS (use dock.Dock_Workflow(workflow_type)):
     - "check_dependencies" - Check and install required packages
     - "init" - Initialize docking project structure
     - "prepare_receptor" - Convert PDB to PDBQT (receptor preparation)
     - "prepare_ligand" - Convert SDF/MOL2 to PDBQT (ligand preparation)
     - "docking_vina" - Single protein-ligand docking with Vina
     - "batch_docking" - Multiple protein-ligand pairs docking
     - "analyze_interactions" - Basic interaction analysis from PDBQT files
     
  3) **EXECUTE** the adapted bash commands and analyze results
  4) **VERIFY** success by checking binding energies and output files
  5) mark_task_done("brief, precise description of the completed step")
  6) show_todos()
Repeat until all todos are completed.

PHASE 3 — ADAPTIVE TODO REFINEMENT
- If dependencies missing → add_todo("Install missing docking tools")
- If binding sites unclear → add_todo("Define binding site coordinates")
- If poor docking scores → add_todo("Optimize docking parameters")

EXECUTION STRATEGY - NEW VERSION (MUST FOLLOW THIS ORDER)
  1) dock.Dock_Workflow("check_dependencies") → Install Python vina and meeko packages
  2) dock.Dock_Workflow("init") → Create project structure with scripts/ folder
  3) show_todos()
  4) Analyze folder and create todos if empty
  5) Loop Phase 2 until all done; refine with Phase 3 when needed
  
🔧 **KEY DIFFERENCES IN NEW VERSION:**
- Uses Python vina API instead of vina command line tool
- Generates Python scripts in scripts/ folder before execution  
- No pymol or plip dependencies required
- All docking operations use "from vina import Vina"

📊 EXECUTION EXAMPLES:

🏷️ STEP 1 - DEPENDENCY CHECK:
```bash
# Get dependency check commands
dep_commands = dock.Dock_Workflow("check_dependencies")
# Execute: pip install meeko vina
```

🏷️ STEP 2 - PROJECT SETUP:
```bash
# Get initialization commands
init_commands = dock.Dock_Workflow("init")
# Execute: mkdir -p docking/receptors, ligands, prepare, output, etc.
```

🏷️ STEP 3 - RECEPTOR PREPARATION:
```bash
# Get receptor preparation commands
prep_commands = dock.Dock_Workflow("prepare_receptor")
# Execute: mk_prepare_receptor.py -i protein.pdb -o prepare/protein.pdbqt
```

🏷️ STEP 4 - LIGAND PREPARATION:
```bash
# Get ligand preparation commands
ligand_commands = dock.Dock_Workflow("prepare_ligand")
# Execute: mk_prepare_ligand.py -i ligand.sdf -o prepare/ligand.pdbqt
```

🏷️ STEP 5 - DOCKING EXECUTION (NEW VERSION):
```bash
# Get docking template commands (generates Python scripts)
dock_commands = dock.Dock_Workflow("batch_docking")
# Execute: creates scripts/batch_docking.py and runs it
# Python script uses: from vina import Vina (NOT command line vina)
# Check results: grep "REMARK VINA RESULT" output/*.pdbqt
```

🏷️ STEP 6 - INTERACTION ANALYSIS:
```bash
# Get interaction analysis commands
analysis_commands = dock.Dock_Workflow("analyze_interactions")
# Execute basic interaction analysis
# Check: binding site residues and interaction distances
```

**CRITICAL SUCCESS PATTERNS (NEW VERSION):**
```bash
# Good - Check Python packages (NEW - no plip/pymol needed)
pip list | grep -E "(meeko|vina)"

# Good - Verify file formats
file receptors/*.pdb ligands/*.sdf

# Good - Check scripts generation (NEW)
ls scripts/*.py

# Good - Check docking results (same as before)
grep "REMARK VINA RESULT" output/docked_poses.pdbqt | head -5

# Good - Analyze binding energies
awk '/REMARK VINA RESULT/ {{print $4}}' output/*.pdbqt | sort -n
```

**GRID CENTER DETECTION:**
The toolset automatically detects binding site centers using:
1. pythonsh prepare_gpf.py (if available)
2. Geometric center calculation from PDBQT coordinates
3. Geometric center from receptor coordinates
4. User-specified active site residues

BEGIN NOW:
- Execute PHASE 0 → PHASE 1 → PHASE 2 loop.
- **ACTUALLY EXECUTE** the bash commands after adapting them
- **ANALYZE BINDING RESULTS** after each step
- **VERIFY SUCCESS** by checking binding energies and poses
- Report any convergence issues or unreasonable binding scores
"""
        
    else:
        message = """
I need help with molecular docking analysis using your specialized workflow-based toolsets.

⚡ **NEW VERSION**: This dock toolset uses Python vina API (not command line) and generates Python scripts.

You have access to comprehensive molecular docking and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE MOLECULAR DOCKING WORKFLOW TOOLSET:

DOCKING WORKFLOWS (use dock.Dock_Workflow(workflow_type)):
These return bash command templates for execution:

PROJECT SETUP:
- "check_dependencies" - Check and install docking tools
- "init" - Initialize docking project structure

STRUCTURE PREPARATION:
- "prepare_receptor" - Convert PDB to PDBQT format
- "prepare_ligand" - Convert SDF/MOL2 to PDBQT format

MOLECULAR DOCKING:
- "docking_vina" - Single receptor-ligand docking
- "batch_docking" - Multiple protein-ligand pairs docking

ANALYSIS:
- "analyze_interactions" - Basic interaction analysis from PDBQT files

WORKFLOW USAGE (NEW VERSION):
All workflows return executable bash command templates that generate and run Python scripts. You must:
1. **Call the workflow function** to get bash command templates
2. **Execute the bash commands** which will create Python scripts in scripts/ folder and run them
3. **The Python scripts use "from vina import Vina"** (not command line vina)  
4. **Analyze the binding results** from PDBQT output files

EXECUTION EXAMPLES:
```bash
# Step 1: Get command template
dep_commands = dock.Dock_Workflow("check_dependencies")

# Step 2: Adapt and execute
# Template: pip install meeko vina
# Execute the installation commands

# Step 3: Verify installation
pip list | grep -E "(meeko|vina)"
```

⚠️ **CRITICAL (NEW VERSION)**: Don't just call the workflow functions - you must EXECUTE the returned bash commands which will generate and run Python scripts!

**SUPPORTED FILE FORMATS:**
- Receptors: PDB files
- Ligands: SDF, MOL2 files
- Output: PDBQT files for docking, PDB files for visualization

**DOCKING PARAMETERS:**
- AutoDock Vina scoring function
- Exhaustiveness: 32 (thorough search)
- Number of poses: 20 per ligand
- Box size: 20x20x20 Å (adjustable)

Please start by adding a todo for your molecular docking task, then use the workflow commands and EXECUTE them!"""
    
    return message