"""Prompt generator for SCFM (Single Cell Foundation Model) workflows.

The CLI accepts free-form natural language from the user and passes it
directly to the Pantheon Agent.  The Agent's LLM router inspects its
registered tools (SingleCellAgent, run_python_code, etc.) and
autonomously selects the right scFM and analysis parameters.

This module only adds minimal scFM context so the router knows the
request is in the single-cell foundation model domain.
"""


# Reference catalogues — kept for /bio scfm list_models and list_analysis_types
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


def generate_scfm_workflow_message(user_query: str) -> str:
    """Wrap the user's natural-language request with scFM context.

    The Pantheon Agent's LLM router will read this message, recognise
    the scFM intent, and autonomously select the appropriate tool
    (SingleCellAgent, run_python_code, etc.) and parameters.

    Args:
        user_query: The user's free-form natural-language request,
            e.g. "annotate cell types in pbmc3k.h5ad using scGPT".

    Returns:
        A message for the Agent that preserves the user's intent and
        adds minimal scFM context for the router.
    """
    return (
        f"[scFM request] {user_query}\n"
        "\n"
        "Use the registered SingleCellAgent tool or run_python_code to "
        "handle this single-cell foundation model request. Select the "
        "appropriate scFM model and analysis type based on the user's intent."
    )
