"""Integration tests for CLI -> Pantheon Agent -> scFM call chain.

Validates the full path: user command -> REPL -> BioCommandHandler -> scfm_workflow
prompt generation -> agent message dispatch.  All external dependencies (pantheon
Agent, heavy toolsets) are mocked so the tests run without installing the full
pantheon stack.

The prompt is now a concise user-intent message (not a multi-phase instruction
manual).  The Agent's LLM is expected to recognise the request and call the
registered SingleCellAgent tool autonomously.
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
# 1. scFM Workflow Prompt – Template Correctness
# ===========================================================================


class TestScfmPromptTemplate:
    """Verify generate_scfm_workflow_message produces correct, concise prompts."""

    def test_auto_model_contains_dataset_and_auto(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="pbmc3k.h5ad")
        assert "pbmc3k.h5ad" in msg
        assert "auto" in msg.lower()

    def test_specific_model_embedded(self):
        mod = _import_scfm_workflow()
        for model in ("scgpt", "scbert", "geneformer", "scfoundation", "uce"):
            msg = mod.generate_scfm_workflow_message(
                dataset_path="cells.h5ad", model_name=model
            )
            assert model in msg

    def test_question_block_present_when_provided(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="x.h5ad", question="annotate T cells"
        )
        assert "User analysis goal" in msg
        assert "annotate T cells" in msg

    def test_question_block_absent_when_none(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="x.h5ad")
        assert "User analysis goal" not in msg

    def test_prompt_is_concise(self):
        """The generated prompt should be concise, not a verbose instruction manual."""
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="x.h5ad")
        # Should NOT contain verbose instruction phases
        assert "PHASE 0" not in msg
        assert "PHASE 1" not in msg
        assert "PHASE 2" not in msg
        assert "TIER 1" not in msg
        assert "TIER 2" not in msg
        assert "PERSISTENT STATE" not in msg
        # Should be short
        assert len(msg) < 500

    def test_prompt_mentions_scfm(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="x.h5ad")
        assert "SCFM" in msg or "Foundation Model" in msg

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

    def test_analysis_type_embedded_in_prompt(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="x.h5ad", analysis_type="annotation"
        )
        assert "Analysis type" in msg
        assert "annotation" in msg

    def test_analysis_type_absent_when_none(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="x.h5ad")
        assert "Analysis type" not in msg

    def test_analysis_type_invalid_not_embedded(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="x.h5ad", analysis_type="nonexistent_type"
        )
        assert "Analysis type" not in msg

    def test_all_params_combined(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="data.h5ad",
            model_name="scgpt",
            question="find DEGs",
            analysis_type="differential",
        )
        assert "data.h5ad" in msg
        assert "scgpt" in msg
        assert "find DEGs" in msg
        assert "differential" in msg


# ===========================================================================
# 2. BioCommandHandler – SCFM Routing
# ===========================================================================


class TestScfmCommandRouting:
    """Verify the bio_handler routes /bio scfm commands to the correct handlers."""

    @pytest.mark.asyncio
    async def test_scfm_help_shows_all_subcommands(self, handler):
        """'/bio scfm' should show help listing init, run, list_models."""
        result = await handler.handle_bio_command("/bio scfm")
        assert result is None  # help only prints, returns None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "init" in printed
        assert "run" in printed
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
    async def test_scfm_run_no_dataset_errors(self, handler):
        result = await handler.handle_bio_command("/bio scfm run")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "Error" in printed

    @pytest.mark.asyncio
    async def test_scfm_run_with_dataset_returns_prompt(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command("/bio scfm run ./data.h5ad")
        assert result is not None
        assert "./data.h5ad" in result

    @pytest.mark.asyncio
    async def test_scfm_run_model_flag_forwarded(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --model geneformer"
            )
        assert "geneformer" in result

    @pytest.mark.asyncio
    async def test_scfm_run_question_flag_forwarded(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --question find_markers"
            )
        assert "find_markers" in result
        assert "User analysis goal" in result

    @pytest.mark.asyncio
    async def test_scfm_run_model_and_question_combined(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --model scgpt --question batch_correction"
            )
        assert "scgpt" in result
        assert "batch_correction" in result
        assert "User analysis goal" in result

    @pytest.mark.asyncio
    async def test_scfm_run_default_model_is_auto(self, handler):
        """Without --model, the prompt should use auto model selection."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command("/bio scfm run ./d.h5ad")
        assert "auto" in result.lower()

    @pytest.mark.asyncio
    async def test_scfm_run_analysis_type_flag_forwarded(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --analysis_type annotation"
            )
        assert "Analysis type" in result
        assert "annotation" in result

    @pytest.mark.asyncio
    async def test_scfm_run_all_three_flags(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --model scgpt --question find_DEGs --analysis_type differential"
            )
        assert "scgpt" in result
        assert "find_DEGs" in result
        assert "Analysis type" in result
        assert "differential" in result

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
    async def test_scfm_help_mentions_analysis_type(self, handler):
        result = await handler.handle_bio_command("/bio scfm")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "analysis_type" in printed

    @pytest.mark.asyncio
    async def test_scfm_list_models_shows_omicverse_methods(self, handler):
        result = await handler.handle_bio_command("/bio scfm list_models")
        assert result is None
        printed = " ".join(str(c) for c in handler.console.print.call_args_list)
        assert "pySCSA" in printed
        assert "gptcelltype" in printed
        assert "CellVote" in printed
        assert "OmicVerse" in printed

    @pytest.mark.asyncio
    async def test_scfm_unknown_subcommand_falls_through(self, handler):
        result = await handler.handle_bio_command("/bio scfm foobar")
        assert result is not None
        assert "bio_scfm_foobar" in result

    @pytest.mark.asyncio
    async def test_scfm_unknown_subcommand_with_params(self, handler):
        result = await handler.handle_bio_command("/bio scfm custom x y z")
        assert result == "bio_scfm_custom x y z"


# ===========================================================================
# 3. CLI -> REPL -> Agent Integration (mocked Agent)
# ===========================================================================


class TestCliToAgentIntegration:
    """Simulate the full REPL dispatch: /bio scfm -> handler -> agent.run()."""

    @pytest.mark.asyncio
    async def test_scfm_prompt_reaches_agent(self, handler):
        """The prompt returned by the handler would be sent to agent.run().

        We verify the handler returns a non-empty string that an agent could
        process, containing the essential dataset and model information.
        """
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm run ./pbmc3k.h5ad --model scgpt --question annotate_cells"
            )

        # Simulate what REPL does: send prompt to agent
        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="SCFM analysis complete.")

        assert prompt is not None
        assert len(prompt) > 20  # concise but non-trivial

        # Fire mock agent — verifies prompt is agent-compatible
        response = await mock_agent.run(prompt)
        mock_agent.run.assert_called_once_with(prompt)
        assert response == "SCFM analysis complete."

    @pytest.mark.asyncio
    async def test_scfm_auto_model_prompt_reaches_agent(self, handler):
        """Auto model selection prompt should also be agent-compatible."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command("/bio scfm run ./big.h5ad")

        assert prompt is not None
        assert "auto" in prompt.lower()
        assert "./big.h5ad" in prompt

        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="Done.")
        await mock_agent.run(prompt)
        mock_agent.run.assert_called_once()

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
    """Ensure BIO_COMMAND_MAP and suggestions include all SCFM entries."""

    def test_command_map_has_all_scfm_keys(self):
        expected_keys = {"scfm_init", "scfm_run", "scfm_list_models", "scfm_list_analysis_types"}
        assert expected_keys.issubset(set(BIO_COMMAND_MAP.keys()))

    def test_suggestions_contain_scfm_commands(self):
        suggestions = get_bio_command_suggestions()
        assert "/bio scfm init" in suggestions
        assert "/bio scfm run" in suggestions
        assert "/bio scfm list_models" in suggestions
        assert "/bio scfm list_analysis_types" in suggestions

    def test_suggestions_list_is_not_empty(self):
        suggestions = get_bio_command_suggestions()
        assert len(suggestions) > 10  # reasonable minimum


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
    async def test_scfm_run_import_error_handled(self, handler):
        """If scfm_workflow import fails, handler should return None gracefully."""
        with patch.dict(sys.modules, {"pantheon_cli.cli": None}):
            result = await handler.handle_bio_command("/bio scfm run ./d.h5ad")
        assert result is None

    @pytest.mark.asyncio
    async def test_scfm_run_with_h5ad_extension(self, handler):
        """Various .h5ad file paths should be accepted."""
        with _patch_scfm_modules():
            for path in (
                "./data.h5ad",
                "/absolute/path/to/cells.h5ad",
                "relative/pbmc3k.h5ad",
                "../parent/data.h5ad",
            ):
                result = await handler.handle_bio_command(f"/bio scfm run {path}")
                assert result is not None, f"Failed for path: {path}"
                assert path in result

    @pytest.mark.asyncio
    async def test_scfm_run_model_only_no_question(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --model uce"
            )
        assert "uce" in result
        assert "User analysis goal" not in result

    @pytest.mark.asyncio
    async def test_scfm_run_question_only_no_model(self, handler):
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --question cluster_analysis"
            )
        assert "cluster_analysis" in result
        assert "User analysis goal" in result
        assert "auto" in result.lower()

    @pytest.mark.asyncio
    async def test_scfm_run_flags_order_reversed(self, handler):
        """--question before --model should still work."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --question find_DEGs --model scbert"
            )
        assert "scbert" in result
        assert "find_DEGs" in result

    @pytest.mark.asyncio
    async def test_scfm_multiple_consecutive_commands(self, handler):
        """Handler should work correctly for multiple sequential calls."""
        # init
        r1 = await handler.handle_bio_command("/bio scfm init")
        assert "SCFM INIT MODE" in r1

        # list_models
        r2 = await handler.handle_bio_command("/bio scfm list_models")
        assert r2 is None

        # run
        with _patch_scfm_modules():
            r3 = await handler.handle_bio_command("/bio scfm run ./d.h5ad")
        assert "./d.h5ad" in r3


# ===========================================================================
# 7. End-to-End Workflow Validation
# ===========================================================================


class TestScfmEndToEnd:
    """Validate the complete expected workflow: init -> run -> agent receives prompt."""

    @pytest.mark.asyncio
    async def test_full_workflow_init_then_run(self, handler):
        """Simulate user doing /bio scfm init, then /bio scfm run."""
        # Step 1: init
        init_prompt = await handler.handle_bio_command("/bio scfm init")
        assert init_prompt is not None
        assert "clear_all_todos" in init_prompt

        # Step 2: Mock agent processes init
        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="SCFM init ready • todos=0")
        init_response = await mock_agent.run(init_prompt)
        assert "todos=0" in init_response

        # Step 3: run
        with _patch_scfm_modules():
            run_prompt = await handler.handle_bio_command(
                "/bio scfm run ./pbmc3k.h5ad --model scgpt --question annotate_T_cells"
            )
        assert run_prompt is not None
        assert "scgpt" in run_prompt
        assert "annotate_T_cells" in run_prompt
        assert "./pbmc3k.h5ad" in run_prompt

        # Step 4: Mock agent processes run
        mock_agent.run = AsyncMock(
            return_value="SCFM analysis complete. Results saved to ./results/"
        )
        run_response = await mock_agent.run(run_prompt)
        assert "complete" in run_response.lower()

    @pytest.mark.asyncio
    async def test_prompt_structure_for_agent_consumption(self, handler):
        """The generated prompt should be concise — dataset, model, and intent."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm run ./data.h5ad --model scgpt"
            )

        # Agent needs dataset path and model
        assert "./data.h5ad" in prompt
        assert "scgpt" in prompt
        # Should mention SCFM or Foundation Model
        assert "SCFM" in prompt or "Foundation Model" in prompt
        # Should NOT contain verbose instruction phases
        assert "PHASE 0" not in prompt
        assert "PHASE 1" not in prompt
        assert "PHASE 2" not in prompt

    @pytest.mark.asyncio
    async def test_full_workflow_with_analysis_type(self, handler):
        """Simulate user using --analysis_type for OmicVerse-based analysis."""
        with _patch_scfm_modules():
            prompt = await handler.handle_bio_command(
                "/bio scfm run ./pbmc3k.h5ad --analysis_type trajectory --question pseudotime"
            )
        assert prompt is not None
        assert "Analysis type" in prompt
        assert "trajectory" in prompt
        assert "pseudotime" in prompt

        mock_agent = MagicMock()
        mock_agent.run = AsyncMock(return_value="Trajectory analysis complete.")
        response = await mock_agent.run(prompt)
        assert "complete" in response.lower()
