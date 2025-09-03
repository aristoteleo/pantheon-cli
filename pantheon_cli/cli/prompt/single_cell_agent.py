"""Prompt generator for SingleCellAgent workflows (LLM-first orchestration).

This module follows the Biology_toolsets_build_guide.md standards:
- Single dispatcher entry under /bio singlecellagent
- Phase-based guidance (setup, todo, execution)
- Persistent interpreter rules and automatic workflow mode
"""

from textwrap import dedent


def generate_singlecell_workflow_message(
    dataset_path: str,
    profile: str | None = "auto",
    question: str | None = None,
    save_results: bool = True,
    visualize: bool = True,
) -> str:
    """Create an instruction message for the LLM to run SingleCellAgent.

    Args:
        dataset_path: Path to .h5ad dataset
        profile: auto|quick|standard|comprehensive (hint only; LLM decides)
        question: Optional research question to bias interpretation
        save_results: Whether to save JSON/figures
        visualize: Whether to generate plots

    Returns:
        A multi-section prompt for the LLM to execute autonomously.
    """

    question_line = f"Research Question: {question}" if question else "Research Question: General exploration"

    return dedent(
        f"""
        🧬 Single-Cell RNA-seq Analysis Pipeline (OmicVerse + Pantheon)

        {question_line}
        Dataset: {dataset_path}
        Profile hint: {profile}
        Save results: {str(save_results).lower()}
        Visualize: {str(visualize).lower()}

        ⚠️ CRITICAL PYTHON ENVIRONMENT RULES:
        - PERSISTENT STATE: Python interpreter keeps variables across calls.
        - MEMORY OPTIMIZATION: NEVER re-read or re-import if already in memory.
        - SMART CHECKS: Use try/except or `'var' in globals()` to avoid redundant I/O.
        - ERROR RECOVERY: Fix errors in-place; do not restart from scratch.
        - INCREMENTAL PROGRESS: Reuse computed PCA/UMAP/clusters if present.

        PHASE 0 — SETUP & VALIDATION
        1) Verify dataset exists and is readable (no long prints).
        2) Infer dataset profile (cells, genes, obs keys, batch) minimally.
        3) Probe capabilities: check if optional libs (velocity/communication/GRN/drug/metacell) are available.

        PHASE 1 — TODO CREATION (ONCE ONLY)
        Execute: current = show_todos()
        IF empty, create todos ONCE in this order (auto mode):
        1. "QC & Preprocess"
        2. "Clustering"
        3. "Annotation"
        4. "Differential Expression"
        5. "Pathway Enrichment"
        6. "Trajectory"
        7. "Optional: Integration"
        8. "Optional: Velocity"
        9. "Optional: Communication"
        10. "Optional: GRN"
        11. "Optional: Drug Response"
        12. "Optional: Metacell"
        13. "Report Summary"

        ⚡ AUTOMATIC WORKFLOW MODE:
        - Execute each todo without asking for confirmation.
        - After success, call mark_task_done("<task>") and continue immediately.
        - If a step is unsupported, SKIP with a short explanation and proceed.

        PHASE 2 — INTELLIGENT EXECUTION STRATEGY
        🧠 TOOL SELECTION (LLM decides; no user assignment):
        - Use the SingleCellAgent tool with appropriate analysis_type in this progression:
          qc → clustering → annotation → differential → enrichment → trajectory
        - Run optional steps only if capabilities are detected (integration/velocity/communication/GRN/drug/metacell).
        - Prefer fast, robust defaults; escalate to advanced methods if feasible.
        - Persist outputs to singlecell_results/ (JSON + figures). Avoid verbose prints.

        ✅ REQUIRED TOOL CALLS (examples; adapt parameters as needed):
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="qc", save_results={str(save_results).lower()}, visualize={str(visualize).lower()})
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="clustering", save_results={str(save_results).lower()}, visualize={str(visualize).lower()})
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="annotation", save_results={str(save_results).lower()}, visualize={str(visualize).lower()})
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="differential", save_results={str(save_results).lower()}, visualize=false)
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="visualization", save_results={str(save_results).lower()}, visualize={str(visualize).lower()})
        - SingleCellAgent(dataset="{dataset_path}", analysis_type="trajectory", save_results={str(save_results).lower()}, visualize={str(visualize).lower()})

        📊 REPORTING
        - Summarize key findings (cell types, DE genes, top pathways, trajectory branches).
        - Save a concise Markdown report under singlecell_results/report.md.
        - Include figure paths; avoid embedding large binary content in messages.

        🔒 SAFETY & FALLBACKS
        - If a dependency is missing, log a one-line suggestion, SKIP the step, and continue.
        - If memory pressure is detected, sub-sample for visualization while preserving core stats.
        - Always continue to the next task unless a critical failure occurs (e.g., dataset unreadable).

        Return only the necessary tool calls and minimal narrative to proceed efficiently.
        """
    ).strip()


__all__ = ["generate_singlecell_workflow_message"]

