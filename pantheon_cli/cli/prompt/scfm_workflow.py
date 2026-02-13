"""Prompt generator for SCFM (Single Cell Foundation Model) workflows.

Supports foundation model-based single-cell analysis using models such as
scGPT, scBERT, Geneformer, scFoundation, and UCE.
"""

from textwrap import dedent


SUPPORTED_MODELS = {
    "scgpt": "scGPT - Generative pre-trained transformer for single-cell",
    "scbert": "scBERT - BERT-based cell type annotation",
    "geneformer": "Geneformer - Transfer learning for single-cell genomics",
    "scfoundation": "scFoundation - Large-scale foundation model",
    "uce": "UCE - Universal Cell Embedding",
}


def generate_scfm_workflow_message(
    dataset_path: str,
    model_name: str = "auto",
) -> str:
    """Create an instruction message for the LLM to run an SCFM workflow.

    Args:
        dataset_path: Path to .h5ad dataset.
        model_name: Foundation model to use (auto|scgpt|scbert|geneformer|scfoundation|uce).

    Returns:
        A multi-section prompt for the LLM to execute autonomously.
    """

    model_line = (
        f"Model: {model_name}"
        if model_name != "auto"
        else "Model: auto (LLM will select the best model based on data characteristics)"
    )

    return dedent(
        f"""\
        🧬 Single Cell Foundation Model (SCFM) Analysis Pipeline

        Dataset: {dataset_path}
        {model_line}

        ⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
        - PERSISTENT STATE: Python interpreter keeps variables across calls.
        - MEMORY OPTIMIZATION: NEVER re-read or re-import if already in memory.
        - SMART CHECKS: Use try/except or `'var' in globals()` to avoid redundant I/O.
        - ERROR RECOVERY: Fix errors in-place; do not restart from scratch.
        - INCREMENTAL PROGRESS: Reuse computed embeddings/annotations if present.

        PHASE 0 — SETUP & VALIDATION
        1) Verify dataset exists and is readable.
        2) Inspect dataset: number of cells, genes, obs keys, batch information.
        3) Probe environment: check if foundation model libraries are available
           (scgpt, geneformer, transformers, torch, scanpy, anndata).

        PHASE 1 — TODO CREATION (ONCE ONLY)
        Execute: current = show_todos()
        IF empty, create todos ONCE in this order:
        1. "Load and inspect dataset"
        2. "Preprocess data for foundation model input"
        3. "Run foundation model embedding/annotation"
        4. "Evaluate and visualize results"
        5. "Save results and generate report"

        ⚡ AUTOMATIC WORKFLOW MODE:
        - Execute each todo task automatically without asking for confirmation
        - After successful completion, call mark_task_done("description") and proceed
        - Continue until all tasks complete or user intervenes

        PHASE 2 — INTELLIGENT EXECUTION STRATEGY
        🧠 SMART DECISION MAKING:

        **ASSESS CURRENT SITUATION FIRST:**
        - What data do you have loaded?
        - What preprocessing steps are completed?
        - What model is selected and available?

        **MODEL SELECTION (if auto):**
        - For cell type annotation → prefer scGPT or scBERT
        - For general embeddings → prefer scGPT or Geneformer
        - For large datasets (>100k cells) → prefer scFoundation or UCE
        - Fall back to available model if preferred is not installed

        **EXECUTION:**
        - Load data with scanpy/anndata
        - Preprocess following model-specific requirements
        - Generate embeddings or annotations using the foundation model
        - Visualize with UMAP colored by predictions
        - Save annotated h5ad and summary figures
        """
    )
