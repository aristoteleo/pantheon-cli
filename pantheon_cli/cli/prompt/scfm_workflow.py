"""Prompt generator for SCFM (Single Cell Foundation Model) workflows.

Generates concise user-intent messages that let the Agent's LLM
autonomously discover and call the registered SingleCellAgent tool
(backed by OmicVerse) or use run_python_code for direct foundation-model
operations (scGPT, Geneformer, etc.).
"""


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
    """Create a concise user-intent message for scFM analysis.

    The message expresses what the user wants and lets the Agent's LLM
    autonomously select the appropriate tools (SingleCellAgent,
    run_python_code, etc.) based on its registered tool descriptions.

    Args:
        dataset_path: Path to .h5ad dataset.
        model_name: Foundation model preference
            (auto|scgpt|scbert|geneformer|scfoundation|uce).
        question: Optional user-specified analysis question or goal.
        analysis_type: Optional analysis type
            (e.g. annotation, trajectory, comprehensive).

    Returns:
        A concise intent message for the LLM.
    """
    parts = [
        f"Run a Single Cell Foundation Model (SCFM) analysis on dataset: {dataset_path}"
    ]

    if model_name != "auto":
        parts.append(f"Model: {model_name}")
    else:
        parts.append(
            "Model: auto — select the best foundation model based on data characteristics."
        )

    if analysis_type and analysis_type in ANALYSIS_TYPES:
        parts.append(
            f"Analysis type: {analysis_type} — {ANALYSIS_TYPES[analysis_type]}"
        )

    if question:
        parts.append(f"User analysis goal: {question}")

    return "\n".join(parts)
