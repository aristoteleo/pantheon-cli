from pathlib import Path


def generate_single_cell_workflow_message(workflow_type: str) -> str:
    """Generate SingleCellAgent workflow message with persistent-state rules and templates"""
    message = f"""
🧬 Single-Cell Analysis Pipeline

🔧 AVAILABLE WORKFLOW TOOLS:
You have access to the SingleCellAgent_Analysis tool for getting official templates.
💡 Use the tool when you need guidance, official examples, or best practices for a specific workflow step.
Ask for dataset path if unknown, use ls to inspect folders, then proceed with file loading.

⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
- PERSISTENT STATE: Python interpreter maintains ALL variables across calls!
- MEMORY OPTIMIZATION: Variables persist! NEVER re-read or re-import data that already exists in memory!
- SMART VARIABLE CHECKING: Use try/except or 'var' in globals() to check existence — NO redundant I/O!
- EFFICIENCY FIRST:
  - Check if adata exists before loading: if 'adata' not in globals()
  - Use existing results: if 'X_pca' in adata.obsm
  - Reuse computed values: if 'marker_genes' in locals()
- ERROR RECOVERY: If code fails, analyze error and fix — don't reload everything!
- NO REPETITION: Each import/load/compute happens ONCE per session unless needed
- After each step: mark_task_done("description"), then show_todos()
- AUTOMATIC EXECUTION: Proceed automatically without confirmations; log warnings when needed.

PHASE 0 — SETUP & VALIDATION
1) Data discovery: Ask for path or use ls to inspect, then load an .h5ad dataset
2) Environment checks happen during data loading

PHASE 1 — TODO CREATION (ONCE ONLY)
Execute: current = show_todos()
IF current is EMPTY, create these todos ONCE:
1. "Check Python environment and load initial data"
2. "Run qc workflow"
3. "Run annotation workflow"
4. "Run trajectory workflow"
5. "Run differential workflow"
6. "Run visualization workflow"
7. "Run clustering workflow"
8. "Run batch_integration workflow"

⚡ AUTOMATIC WORKFLOW MODE:
- Execute each todo automatically; after success, mark_task_done() and proceed to next

PHASE 2 — INTELLIGENT EXECUTION STRATEGY
🧠 SMART DECISION MAKING:

ASSESS:
- What data is loaded? What results are available?
- What guidance is needed?

CHOOSE YOUR APPROACH:

Option A — Use Template Tool (Recommended):
- Call SingleCellAgent_Analysis(workflow_type="<type>") to get official templates
- Study returned guidance and patterns
- Adapt template to your data
- Execute adapted code

Option B — Direct Implementation (Experienced users):
- Write and run code based on your preferred toolkit's docs
- Use help() to verify parameters
- Follow best practices

WHEN TO USE TEMPLATES:
✅ Need official guidance or best practices
✅ Unsure about parameters
✅ Want standardized patterns

RESULT ANALYSIS REQUIREMENT:
1. Analyze outputs (numbers, plots, warnings)
2. Interpret results
3. Check for issues
4. Adjust parameters if needed
5. Document key insights to results dir
6. Proceed based on findings

The returned content serves as GUIDANCE and TEMPLATES, not execution scripts.

🏷️ STEP EXAMPLES:

— annotation —
Use: SingleCellAgent_Analysis(workflow_type="annotation")
Then adapt to your data and execute.

— trajectory —
Use: SingleCellAgent_Analysis(workflow_type="trajectory")

— differential —
Use: SingleCellAgent_Analysis(workflow_type="differential")

— visualization —
Use: SingleCellAgent_Analysis(workflow_type="visualization")

— qc —
Use: SingleCellAgent_Analysis(workflow_type="qc")

— clustering —
Use: SingleCellAgent_Analysis(workflow_type="clustering")

— batch_integration —
Use: SingleCellAgent_Analysis(workflow_type="batch_integration")

— communication | grn | drug | metacell —
Use the corresponding workflow_type with SingleCellAgent_Analysis

"""

    return message
