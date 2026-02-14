"""Tests for SCFM (Single Cell Foundation Model) router in bio_handler.py"""

import importlib
import sys
import pytest
from unittest.mock import MagicMock, patch
from rich.console import Console

from pantheon_cli.repl.bio_handler import BioCommandHandler


@pytest.fixture
def handler():
    """Create a BioCommandHandler with a mocked console."""
    console = MagicMock(spec=Console)
    return BioCommandHandler(console)


def _import_scfm_workflow():
    """Import scfm_workflow module directly, bypassing pantheon_cli.cli.__init__."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "scfm_workflow",
        "/home/user/pantheon-cli/pantheon_cli/cli/prompt/scfm_workflow.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Routing tests ────────────────────────────────────────────────────────


class TestScfmRouting:
    """Verify that /bio scfm commands are routed to the SCFM handler."""

    @pytest.mark.asyncio
    async def test_bio_scfm_shows_help(self, handler):
        """'/bio scfm' with no subcommand should show help and return None."""
        result = await handler.handle_bio_command("/bio scfm")
        assert result is None
        # Verify help was printed
        handler.console.print.assert_called()

    @pytest.mark.asyncio
    async def test_bio_scfm_init_returns_init_message(self, handler):
        """'/bio scfm init' should return the SCFM init mode strict message."""
        result = await handler.handle_bio_command("/bio scfm init")
        assert result is not None
        assert "SCFM INIT MODE" in result
        assert "clear_all_todos" in result

    @pytest.mark.asyncio
    async def test_bio_scfm_list_models_returns_none(self, handler):
        """'/bio scfm list_models' prints model list and returns None."""
        result = await handler.handle_bio_command("/bio scfm list_models")
        assert result is None
        # Should have printed model information
        handler.console.print.assert_called()

    @pytest.mark.asyncio
    async def test_bio_scfm_run_without_dataset_returns_none(self, handler):
        """'/bio scfm run' with no dataset should print error and return None."""
        result = await handler.handle_bio_command("/bio scfm run")
        assert result is None
        # Should have printed an error
        calls = [str(c) for c in handler.console.print.call_args_list]
        assert any("Error" in c for c in calls)

    @pytest.mark.asyncio
    async def test_bio_scfm_run_with_dataset(self, handler):
        """'/bio scfm run ./data.h5ad' should return a concise intent message."""
        scfm_mod = _import_scfm_workflow()
        cli_mock = MagicMock()
        prompt_mock = MagicMock()
        prompt_mock.scfm_workflow = scfm_mod
        cli_mock.prompt = prompt_mock
        with patch.dict(
            sys.modules,
            {
                "pantheon_cli.cli": cli_mock,
                "pantheon_cli.cli.prompt": prompt_mock,
                "pantheon_cli.cli.prompt.scfm_workflow": scfm_mod,
            },
        ):
            result = await handler.handle_bio_command("/bio scfm run ./data.h5ad")
        assert result is not None
        assert "SCFM" in result or "Foundation Model" in result
        assert "./data.h5ad" in result

    @pytest.mark.asyncio
    async def test_bio_scfm_run_with_model_flag(self, handler):
        """'/bio scfm run ./data.h5ad --model scgpt' should pass model to prompt."""
        scfm_mod = _import_scfm_workflow()
        cli_mock = MagicMock()
        prompt_mock = MagicMock()
        prompt_mock.scfm_workflow = scfm_mod
        cli_mock.prompt = prompt_mock
        with patch.dict(
            sys.modules,
            {
                "pantheon_cli.cli": cli_mock,
                "pantheon_cli.cli.prompt": prompt_mock,
                "pantheon_cli.cli.prompt.scfm_workflow": scfm_mod,
            },
        ):
            result = await handler.handle_bio_command(
                "/bio scfm run ./data.h5ad --model scgpt"
            )
        assert result is not None
        assert "scgpt" in result

    @pytest.mark.asyncio
    async def test_bio_scfm_unknown_subcommand(self, handler):
        """'/bio scfm unknown' should fall through to generic handler."""
        result = await handler.handle_bio_command("/bio scfm unknown")
        assert result is not None
        assert "bio_scfm_unknown" in result

    @pytest.mark.asyncio
    async def test_bio_scfm_unknown_with_params(self, handler):
        """'/bio scfm unknown foo bar' should pass params to generic handler."""
        result = await handler.handle_bio_command("/bio scfm unknown foo bar")
        assert result == "bio_scfm_unknown foo bar"

    @pytest.mark.asyncio
    async def test_bio_scfm_run_with_question_flag(self, handler):
        """'/bio scfm run ./data.h5ad --question ...' should pass question to prompt."""
        scfm_mod = _import_scfm_workflow()
        cli_mock = MagicMock()
        prompt_mock = MagicMock()
        prompt_mock.scfm_workflow = scfm_mod
        cli_mock.prompt = prompt_mock
        with patch.dict(
            sys.modules,
            {
                "pantheon_cli.cli": cli_mock,
                "pantheon_cli.cli.prompt": prompt_mock,
                "pantheon_cli.cli.prompt.scfm_workflow": scfm_mod,
            },
        ):
            result = await handler.handle_bio_command(
                "/bio scfm run ./data.h5ad --question identify_T_cells"
            )
        assert result is not None
        assert "identify_T_cells" in result
        assert "User analysis goal" in result

    @pytest.mark.asyncio
    async def test_bio_scfm_run_with_model_and_question(self, handler):
        """'/bio scfm run ./d.h5ad --model scgpt --question ...' passes both."""
        scfm_mod = _import_scfm_workflow()
        cli_mock = MagicMock()
        prompt_mock = MagicMock()
        prompt_mock.scfm_workflow = scfm_mod
        cli_mock.prompt = prompt_mock
        with patch.dict(
            sys.modules,
            {
                "pantheon_cli.cli": cli_mock,
                "pantheon_cli.cli.prompt": prompt_mock,
                "pantheon_cli.cli.prompt.scfm_workflow": scfm_mod,
            },
        ):
            result = await handler.handle_bio_command(
                "/bio scfm run ./d.h5ad --model scgpt --question find_markers"
            )
        assert result is not None
        assert "scgpt" in result
        assert "find_markers" in result
        assert "User analysis goal" in result


# ── Prompt template tests ────────────────────────────────────────────────


class TestScfmWorkflowPrompt:
    """Verify the SCFM workflow prompt generator."""

    def test_generate_auto_model(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="test.h5ad")
        assert "test.h5ad" in msg
        assert "auto" in msg.lower()

    def test_generate_specific_model(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="my_data.h5ad",
            model_name="geneformer",
        )
        assert "my_data.h5ad" in msg
        assert "geneformer" in msg

    def test_generate_with_question(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="cells.h5ad",
            question="annotate T cell subtypes",
        )
        assert "cells.h5ad" in msg
        assert "annotate T cell subtypes" in msg
        assert "User analysis goal" in msg

    def test_generate_without_question(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="cells.h5ad")
        assert "User analysis goal" not in msg

    def test_supported_models_dict(self):
        mod = _import_scfm_workflow()
        assert "scgpt" in mod.SUPPORTED_MODELS
        assert "scbert" in mod.SUPPORTED_MODELS
        assert "geneformer" in mod.SUPPORTED_MODELS
        assert "scfoundation" in mod.SUPPORTED_MODELS
        assert "uce" in mod.SUPPORTED_MODELS

    def test_analysis_types_dict(self):
        mod = _import_scfm_workflow()
        assert hasattr(mod, "ANALYSIS_TYPES")
        assert "comprehensive" in mod.ANALYSIS_TYPES
        assert "annotation" in mod.ANALYSIS_TYPES
        assert "trajectory" in mod.ANALYSIS_TYPES

    def test_generate_with_analysis_type(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(
            dataset_path="cells.h5ad",
            analysis_type="annotation",
        )
        assert "Analysis type" in msg
        assert "annotation" in msg

    def test_prompt_is_concise(self):
        """The prompt should be concise, not a multi-phase instruction manual."""
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message(dataset_path="cells.h5ad")
        assert "SCFM" in msg or "Foundation Model" in msg
        # Should NOT contain verbose instruction phases
        assert "PHASE 0" not in msg
        assert "PHASE 1" not in msg
        assert "PHASE 2" not in msg
        # Should be short (under 500 chars for a simple request)
        assert len(msg) < 500


# ── Integration with existing routing ────────────────────────────────────


class TestScfmIntegration:
    """Verify SCFM integrates correctly with the existing routing system."""

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_scrna(self, handler):
        """SCFM routing should not break existing scrna routing."""
        result = await handler.handle_bio_command("/bio scrna init")
        assert result is not None
        assert "scRNA INIT MODE" in result

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_spatial(self, handler):
        """SCFM routing should not break existing spatial routing."""
        result = await handler.handle_bio_command("/bio spatial init")
        assert result is not None
        assert "spatial INIT MODE" in result

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_atac(self, handler):
        """SCFM routing should not break existing atac routing."""
        result = await handler.handle_bio_command("/bio atac init")
        assert result is not None
        assert "ATAC INIT MODE" in result

    @pytest.mark.asyncio
    async def test_bio_list_still_works(self, handler):
        """'/bio list' should still work."""
        result = await handler.handle_bio_command("/bio list")
        assert result == "bio list"

    @pytest.mark.asyncio
    async def test_bio_help_shows_scfm(self, handler):
        """'/bio' help should mention scfm."""
        result = await handler.handle_bio_command("/bio")
        assert result is None
        # Check that scfm was mentioned in help output
        calls = [str(c) for c in handler.console.print.call_args_list]
        assert any("scfm" in c for c in calls)


# ── Command map tests ────────────────────────────────────────────────────


class TestScfmCommandMap:
    """Verify SCFM entries in BIO_COMMAND_MAP and suggestions."""

    def test_scfm_in_command_map(self):
        from pantheon_cli.repl.bio_handler import BIO_COMMAND_MAP

        assert "scfm_init" in BIO_COMMAND_MAP
        assert "scfm_run" in BIO_COMMAND_MAP
        assert "scfm_list_models" in BIO_COMMAND_MAP
        assert "scfm_list_analysis_types" in BIO_COMMAND_MAP

    def test_scfm_in_suggestions(self):
        from pantheon_cli.repl.bio_handler import get_bio_command_suggestions

        suggestions = get_bio_command_suggestions()
        assert "/bio scfm init" in suggestions
        assert "/bio scfm run" in suggestions
        assert "/bio scfm list_models" in suggestions
        assert "/bio scfm list_analysis_types" in suggestions
