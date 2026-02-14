"""Prompt generator for SCFM (Single Cell Foundation Model) workflows.

Supports foundation model-based single-cell analysis using models such as
scGPT, scBERT, Geneformer, scFoundation, and UCE — integrated with the
Pantheon Agent's SingleCellAgent toolset and OmicVerse APIs.
"""

from textwrap import dedent


SUPPORTED_MODELS = {
    "scgpt": "scGPT - Generative pre-trained transformer for single-cell",
    "scbert": "scBERT - BERT-based cell type annotation",
    "geneformer": "Geneformer - Transfer learning for single-cell genomics",
    "scfoundation": "scFoundation - Large-scale foundation model",
    "uce": "UCE - Universal Cell Embedding",
}

# Analysis types provided by the SingleCellAgent toolset in pantheon-agents
ANALYSIS_TYPES = {
    "comprehensive": "Full analysis: annotation, trajectory, DE, pathways, visualization",
    "annotation": "Cell type identification (pySCSA / gptcelltype / CellVote via OmicVerse)",
    "trajectory": "Pseudotime and developmental trajectory (TrajInfer via OmicVerse)",
    "differential": "Differential expression analysis (DCT via OmicVerse)",
    "visualization": "UMAP/t-SNE embedding and gene expression plots",
    "qc": "Quality control and preprocessing summary",
    "clustering": "Leiden/Louvain clustering",
    "batch_integration": "Batch effect correction (Harmony/Combat via OmicVerse)",
    "communication": "Cell-cell communication (CellPhoneDB via OmicVerse)",
    "grn": "Gene regulatory network / TF activity (SCENIC/AUCell)",
    "drug": "Drug response prediction (scDrug)",
    "metacell": "Metacell construction and summary",
    "custom": "Free-form analysis guided by user research question",
}


def generate_scfm_workflow_message(
    dataset_path: str,
    model_name: str = "auto",
    question: str | None = None,
    analysis_type: str | None = None,
) -> str:
    """Create an instruction message for the LLM to run an SCFM workflow.

    The generated prompt instructs the Agent to prefer the registered
    ``SingleCellAgent`` tool (backed by OmicVerse) for standard analysis
    types, and fall back to ``run_python_code`` for direct foundation-model
    operations (scGPT / Geneformer / etc.) when needed.

    Args:
        dataset_path: Path to .h5ad dataset.
        model_name: Foundation model to use (auto|scgpt|scbert|geneformer|scfoundation|uce).
        question: Optional user-specified analysis question or goal.
        analysis_type: Optional analysis type matching SingleCellAgent
            (e.g. annotation, trajectory, comprehensive).

    Returns:
        A multi-section prompt for the LLM to execute autonomously.
    """

    model_line = (
        f"Model: {model_name}"
        if model_name != "auto"
        else "Model: auto (LLM will select the best model based on data characteristics)"
    )

    question_block = ""
    if question:
        question_block = f"""
        USER ANALYSIS GOAL:
        {question}
        (Tailor the entire workflow — model selection, preprocessing, and evaluation — to address this goal.)
        """

    analysis_type_block = ""
    if analysis_type and analysis_type in ANALYSIS_TYPES:
        analysis_type_block = f"""
        REQUESTED ANALYSIS TYPE: {analysis_type}
        Description: {ANALYSIS_TYPES[analysis_type]}
        """

    return dedent(
        f"""\
        🧬 Single Cell Foundation Model (SCFM) Analysis Pipeline

        Dataset: {dataset_path}
        {model_line}
        {question_block}
        {analysis_type_block}

        ⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
        - PERSISTENT STATE: Python interpreter keeps variables across calls.
        - MEMORY OPTIMIZATION: NEVER re-read or re-import if already in memory.
        - SMART CHECKS: Use try/except or `'var' in globals()` to avoid redundant I/O.
        - ERROR RECOVERY: Fix errors in-place; do not restart from scratch.
        - INCREMENTAL PROGRESS: Reuse computed embeddings/annotations if present.

        PHASE 0 — SETUP & VALIDATION
        1) Verify dataset exists and is readable.
        2) Inspect dataset: number of cells, genes, obs keys, batch information.
        3) Probe environment — check BOTH paths:
           a) SingleCellAgent tool availability (registered via BioToolsetManager).
              If the ``SingleCellAgent`` tool is available, it provides OmicVerse-powered
              analysis: annotation (pySCSA/gptcelltype/CellVote), trajectory (TrajInfer),
              DE (DCT), clustering, batch integration, communication, GRN, drug, metacell.
           b) Foundation model libraries: scgpt, geneformer, transformers, torch,
              scanpy, anndata, omicverse.

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

        **TOOL SELECTION — TWO-TIER STRATEGY:**

        TIER 1 (PREFERRED) — Use the ``SingleCellAgent`` tool when:
        - The analysis type matches one of: comprehensive, annotation, trajectory,
          differential, visualization, qc, clustering, batch_integration,
          communication, grn, drug, metacell, custom.
        - Call it via: SingleCellAgent(dataset="{dataset_path}",
          analysis_type="<type>", research_question="<question>")
        - This uses OmicVerse internally and handles the entire workflow.

        TIER 2 (FALLBACK) — Use ``run_python_code`` when:
        - A specific foundation model (scGPT, scBERT, Geneformer, scFoundation, UCE)
          is requested via --model and is NOT covered by SingleCellAgent.
        - The SingleCellAgent tool is not available or returned an error.
        - Custom low-level FM operations are needed (e.g. extracting raw embeddings,
          fine-tuning, or model-specific preprocessing).

        **ASSESS CURRENT SITUATION FIRST:**
        - What data do you have loaded?
        - What preprocessing steps are completed?
        - What model is selected and available?

        **MODEL SELECTION (if auto):**
        - For cell type annotation → prefer SingleCellAgent(analysis_type="annotation")
          which uses OmicVerse pySCSA; OR scGPT / scBERT via run_python_code
        - For general embeddings → prefer scGPT or Geneformer via run_python_code
        - For large datasets (>100k cells) → prefer scFoundation or UCE
        - For trajectory / DE / clustering / other standard analyses →
          prefer SingleCellAgent with the matching analysis_type
        - Fall back to available model if preferred is not installed

        **EXECUTION:**
        - Load data with scanpy/anndata
        - Preprocess following model-specific requirements
        - Generate embeddings or annotations using the foundation model
        - Visualize with UMAP colored by predictions
        - Save annotated h5ad and summary figures
        """
    )
