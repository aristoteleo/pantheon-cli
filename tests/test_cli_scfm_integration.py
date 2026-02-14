"""Integration tests for CLI -> Pantheon Agent -> scFM call chain.

Validates the full path: user command -> REPL -> BioCommandHandler ->
natural language pass-through -> agent message dispatch.

The scFM design follows the pantheon-agents architecture: the user provides
natural language after `/bio scfm`, and the Pantheon Agent's LLM router
reads its registered tool descriptions (SingleCellAgent, run_python_code,
etc.) and autonomously selects the right scFM model and analysis parameters.

The CLI does NOT parse --model / --analysis_type flags.
"""

import importlib
import importlib.util
import sys
import asyncio
import pytest
from unittest.mock import AsyncMock, MagicMock, patch, PropertyMock
from rich.console import Console

from pantheon_cli.repl.bio_handler import (
    BioCommandHandler,
    BIO_COMMAND_MAP,
    get_bio_command_suggestions,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _import_scfm_workflow():
    """Import scfm_workflow module directly, bypassing pantheon_cli.cli.__init__."""
    spec = importlib.util.spec_from_file_location(
        "scfm_workflow",
        "/home/user/pantheon-cli/pantheon_cli/cli/prompt/scfm_workflow.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _patch_scfm_modules():
    """Return a patch context that makes the relative scfm_workflow import work."""
    scfm_mod = _import_scfm_workflow()
    cli_mock = MagicMock()
    prompt_mock = MagicMock()
    prompt_mock.scfm_workflow = scfm_mod
    cli_mock.prompt = prompt_mock
    return patch.dict(
        sys.modules,
        {
            "pantheon_cli.cli": cli_mock,
            "pantheon_cli.cli.prompt": prompt_mock,
            "pantheon_cli.cli.prompt.scfm_workflow": scfm_mod,
        },
    )


@pytest.fixture
def handler():
    """Create a BioCommandHandler with a mocked console."""
    console = MagicMock(spec=Console)
    return BioCommandHandler(console)


# ===========================================================================
# 1. scFM Workflow Prompt – Natural Language Pass-Through
# ===========================================================================


class TestScfmPromptTemplate:
    """Verify generate_scfm_workflow_message wraps natural language correctly."""

    def test_user_query_preserved_in_output(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            "annotate cell types in pbmc3k.h5ad"
        )
        assert "annotate cell types in pbmc3k.h5ad" in msg

    def test_scfm_context_tag_present(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("analyze my data")
        assert "scFM" in msg

    def test_routing_hint_mentions_tools(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("find DEGs between clusters")
        assert "SingleCellAgent" in msg

    def test_prompt_is_concise(self):
        """The generated prompt should be concise, not a verbose instruction manual."""
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("analyze x.h5ad")
        assert "PHASE 0" not in msg
        assert "PHASE 1" not in msg
        assert "PHASE 2" not in msg
        assert "TIER 1" not in msg
        assert "PERSISTENT STATE" not in msg
        assert len(msg) < 500

    def test_supported_models_completeness(self):
        mod = _import_scfm_workflow()
        expected = {"scgpt", "scbert", "geneformer", "scfoundation", "uce"}
        assert set(mod.SUPPORTED_MODELS.keys()) == expected

    def test_analysis_types_dict_exists(self):
        mod = _import_scfm_workflow()
        assert hasattr(mod, "ANALYSIS_TYPES")
        expected_types = {
            "comprehensive", "annotation", "trajectory", "differential",
            "visualization", "qc", "clustering", "batch_integration",
            "communication", "grn", "drug", "metacell", "custom",
        }
        assert set(mod.ANALYSIS_TYPES.keys()) == expected_types

    def test_various_natural_language_queries(self):
        """Various natural-language queries should all be preserved."""
        mod = _import_scfm_workflow()
        queries = [
            "annotate cell types in pbmc3k.h5ad using scGPT",
            "run trajectory analysis on my_data.h5ad",
            "find DEGs between clusters in ./data.h5ad",
            "comprehensive analysis of single_cell.h5ad with geneformer",
            "cluster and visualize neurons.h5ad",
        ]
        for q in queries:
            msg = mod.generate_scfm_workflow_message(q)
            assert q in msg, f"Query not preserved: {q}"
            assert "scFM" in msg


# ===========================================================================
# 2. BioCommandHandler – SCFM Routing
# ===========================================================================


class TestScfmCommandRouting:
    """Verify the bio_handler routes /bio scfm commands correctly."""

    @pytest.mark.asyncio
    async def test_scfm_help_shows_natural_language_usage(self, handler):
        """'/bio scfm' should show help mentioning natural language usage."""
        result = await handler.handle_bio_command("/bio scfm")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "init" in printed
        assert "list_models" in printed

    @pytest.mark.asyncio
    async def test_scfm_init_returns_strict_clear(self, handler):
        result = await handler.handle_bio_command("/bio scfm init")
        assert result is not None
        assert "SCFM INIT MODE" in result
        assert "clear_all_todos" in result
        assert "show_todos" in result

    @pytest.mark.asyncio
    async def test_scfm_list_models_prints_all_models(self, handler):
        result = await handler.handle_bio_command("/bio scfm list_models")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        for model in ("scgpt", "scbert", "geneformer", "scfoundation", "uce"):
            assert model in printed

    @pytest.mark.asyncio
    async def test_scfm_natural_language_returns_prompt(self, handler):
        """Natural language after /bio scfm should return a prompt for agent."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm annotate cell types in ./data.h5ad"
            )
        assert result is not None
        assert "annotate cell types in ./data.h5ad" in result

    @pytest.mark.asyncio
    async def test_scfm_natural_language_model_mention_preserved(self, handler):
        """Mentioning a model in natural language should be preserved."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm use geneformer to analyze ./d.h5ad"
            )
        assert "geneformer" in result

    @pytest.mark.asyncio
    async def test_scfm_natural_language_analysis_mention_preserved(self, handler):
        """Mentioning an analysis type in natural language should be preserved."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm trajectory analysis on ./d.h5ad"
            )
        assert "trajectory" in result

    @pytest.mark.asyncio
    async def test_scfm_list_analysis_types_returns_none(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command("/bio scfm list_analysis_types")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "comprehensive" in printed
        assert "annotation" in printed
        assert "trajectory" in printed

    @pytest.mark.asyncio
    async def test_scfm_list_models_shows_omicverse_methods(self, handler):
        result = await handler.handle_bio_command("/bio scfm list_models")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "pySCSA" in printed
        assert "gptcelltype" in printed
        assert "CellVote" in printed
        assert "OmicVerse" in printed


# ===========================================================================
# 3. CLI -> REPL -> Agent Integration (mocked Agent)
# ===========================================================================


class TestCliToAgentIntegration:
    """Simulate the full REPL dispatch: /bio scfm -> handler -> agent.run()."""

    @pytest.mark.asyncio
    async def test_natural_language_prompt_reaches_agent(self, handler):
        """The natural-language prompt should be sent to agent.run()."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm annotate cell types in ./pbmc3k.h5ad using scGPT"
            )

        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="SCFM analysis complete.")

        assert prompt is not None
        assert len(prompt) > 20

        response = await mock_agent.run(prompt)
        mock_agent.run.assert_called_once_with(prompt)
        assert response == "SCFM analysis complete."

    @pytest.mark.asyncio
    async def test_prompt_contains_user_intent(self, handler):
        """The prompt sent to agent should contain the user's original intent."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm find DEGs in ./big.h5ad"
            )

        assert prompt is not None
        assert "find DEGs in ./big.h5ad" in prompt
        assert "scFM" in prompt

    @pytest.mark.asyncio
    async def test_scfm_init_prompt_reaches_agent(self, handler):
        """The init mode strict prompt should also be sent to agent."""
        prompt = await handler.handle_bio_command("/bio scfm init")
        assert prompt is not None
        assert "SCFM INIT MODE" in prompt

        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="SCFM init ready")
        await mock_agent.run(prompt)
        mock_agent.run.assert_called_once_with(prompt)


# ===========================================================================
# 4. Command Map & Suggestions Completeness
# ===========================================================================


class TestScfmCommandMapCompleteness:
    """Ensure BIO_COMMAND_MAP and suggestions include SCFM entries."""

    def test_command_map_has_scfm_keys(self):
        expected_keys = {"scfm_init", "scfm_list_models", "scfm_list_analysis_types"}
        assert expected_keys.issubset(set(BIO_COMMAND_MAP.keys()))

    def test_suggestions_contain_scfm_commands(self):
        suggestions = get_bio_command_suggestions()
        assert "/bio scfm init" in suggestions
        assert "/bio scfm list_models" in suggestions
        assert "/bio scfm list_analysis_types" in suggestions

    def test_suggestions_list_is_not_empty(self):
        suggestions = get_bio_command_suggestions()
        assert len(suggestions) > 10


# ===========================================================================
# 5. SCFM Does Not Break Other Bio Tools
# ===========================================================================


class TestScfmNoSideEffects:
    """Verify that SCFM routing does not interfere with other bio commands."""

    @pytest.mark.asyncio
    async def test_scrna_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio scrna init")
        assert result is not None
        assert "scRNA INIT MODE" in result

    @pytest.mark.asyncio
    async def test_atac_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio atac init")
        assert result is not None
        assert "ATAC INIT MODE" in result

    @pytest.mark.asyncio
    async def test_spatial_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio spatial init")
        assert result is not None
        assert "spatial INIT MODE" in result

    @pytest.mark.asyncio
    async def test_rna_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio rna init")
        assert result is not None
        assert "RNA INIT MODE" in result

    @pytest.mark.asyncio
    async def test_hic_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio hic init")
        assert result is not None
        assert "HiC INIT MODE" in result

    @pytest.mark.asyncio
    async def test_dock_init_still_works(self, handler):
        result = await handler.handle_bio_command("/bio dock init")
        assert result is not None
        assert "DOCK INIT MODE" in result

    @pytest.mark.asyncio
    async def test_bio_list_still_works(self, handler):
        result = await handler.handle_bio_command("/bio list")
        assert result == "bio list"

    @pytest.mark.asyncio
    async def test_bio_help_mentions_scfm(self, handler):
        result = await handler.handle_bio_command("/bio")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "scfm" in printed


# ===========================================================================
# 6. Edge Cases & Error Handling
# ===========================================================================


class TestScfmEdgeCases:
    """Edge-case and error-handling tests for SCFM commands."""

    @pytest.mark.asyncio
    async def test_scfm_import_error_handled(self, handler):
        """If scfm_workflow import fails, handler should return None gracefully."""
        with patch.dict(
            sys.modules,
            {
                "pantheon_cli.cli": None,
                "pantheon_cli.cli.prompt": None,
                "pantheon_cli.cli.prompt.scfm_workflow": None,
            },
        ):
            result = await handler.handle_bio_command(
                "/bio scfm analyze ./d.h5ad"
            )
        assert result is None

    @pytest.mark.asyncio
    async def test_various_h5ad_paths_accepted(self, handler):
        """Various .h5ad file paths in natural language should be accepted."""
        with _patch_scfm_modules():
            for path in (
                "./data.h5ad",
                "/absolute/path/to/cells.h5ad",
                "relative/pbmc3k.h5ad",
                "../parent/data.h5ad",
            ):
                result = await handler.handle_bio_command(
                    f"/bio scfm analyze {path}"
                )
                assert result is not None, f"Failed for path: {path}"
                assert path in result

    @pytest.mark.asyncio
    async def test_multiple_consecutive_commands(self, handler):
        """Handler should work correctly for multiple sequential calls."""
        # init
        r1 = await handler.handle_bio_command("/bio scfm init")
        assert "SCFM INIT MODE" in r1

        # list_models
        r2 = await handler.handle_bio_command("/bio scfm list_models")
        assert r2 is None

        # natural language query
        with _patch_scfm_modules():
            r3 = await handler.handle_bio_command(
                "/bio scfm annotate cell types in ./d.h5ad"
            )
        assert "annotate cell types in ./d.h5ad" in r3

    @pytest.mark.asyncio
    async def test_query_with_special_characters(self, handler):
        """Queries with special characters should still work."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm find T-cell subtypes (CD4+/CD8+) in data.h5ad"
            )
        assert result is not None
        assert "T-cell" in result


# ===========================================================================
# 7. End-to-End Workflow Validation
# ===========================================================================


class TestScfmEndToEnd:
    """Validate the complete expected workflow: init -> query -> agent receives prompt."""

    @pytest.mark.asyncio
    async def test_full_workflow_init_then_query(self, handler):
        """Simulate user doing /bio scfm init, then a natural language query."""
        # Step 1: init
        init_prompt = await handler.handle_bio_command("/bio scfm init")
        assert init_prompt is not None
        assert "clear_all_todos" in init_prompt

        # Step 2: Mock agent processes init
        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="SCFM init ready")
        init_response = await mock_agent.run(init_prompt)

        # Step 3: natural language query
        with _patch_scfm_modules():
            run_prompt = await handler.handle_bio_command(
                "/bio scfm annotate T cells in ./pbmc3k.h5ad using scGPT"
            )
        assert run_prompt is not None
        assert "annotate T cells in ./pbmc3k.h5ad using scGPT" in run_prompt
        assert "scFM" in run_prompt

        # Step 4: Mock agent processes query
        mock_agent.run = AsyncMock(
            return_value="SCFM analysis complete. Results saved to ./results/"
        )
        run_response = await mock_agent.run(run_prompt)
        assert "complete" in run_response.lower()

    @pytest.mark.asyncio
    async def test_prompt_structure_for_agent_consumption(self, handler):
        """The generated prompt should contain user intent and scFM context."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm comprehensive analysis of ./data.h5ad"
            )

        assert "comprehensive analysis of ./data.h5ad" in prompt
        assert "scFM" in prompt
        assert "SingleCellAgent" in prompt
        # Should NOT contain verbose instruction phases
        assert "PHASE 0" not in prompt
        assert "PHASE 1" not in prompt
        assert "PHASE 2" not in prompt

    @pytest.mark.asyncio
    async def test_agent_receives_natural_language_not_structured_flags(self, handler):
        """The agent should receive natural language, not parsed --model/--analysis_type."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm use scBERT to annotate cell types in my_data.h5ad"
            )
        # Natural language preserved
        assert "use scBERT to annotate cell types in my_data.h5ad" in prompt
        # No structured flag artifacts
        assert "Model: auto" not in prompt
        assert "Analysis type:" not in prompt
