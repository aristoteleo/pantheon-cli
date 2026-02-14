"""Tests for SCFM (Single Cell Foundation Model) router in bio_handler.py

The scFM design follows the pantheon-agents architecture: the user provides
natural language, and the Pantheon Agent's LLM router selects the right scFM
model and analysis parameters.  The CLI does NOT parse --model/--analysis_type
flags; instead it wraps the user's natural language with minimal scFM context.
"""

import importlib
import importlib.util
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


# ── Routing tests ────────────────────────────────────────────────────────


class TestScfmRouting:
    """Verify that /bio scfm commands are routed correctly."""

    @pytest.mark.asyncio
    async def test_bio_scfm_shows_help(self, handler):
        """'/bio scfm' with no subcommand should show help and return None."""
        result = await handler.handle_bio_command("/bio scfm")
        assert result is None
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
        handler.console.print.assert_called()

    @pytest.mark.asyncio
    async def test_bio_scfm_list_analysis_types_returns_none(self, handler):
        """'/bio scfm list_analysis_types' prints types and returns None."""
        result = await handler.handle_bio_command("/bio scfm list_analysis_types")
        assert result is None
        handler.console.print.assert_called()

    @pytest.mark.asyncio
    async def test_natural_language_query_forwarded(self, handler):
        """'/bio scfm annotate cells in pbmc3k.h5ad' should forward to agent."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm annotate cell types in pbmc3k.h5ad using scGPT"
            )
        assert result is not None
        assert "annotate cell types in pbmc3k.h5ad using scGPT" in result
        assert "scFM" in result

    @pytest.mark.asyncio
    async def test_natural_language_with_dataset_path(self, handler):
        """Natural language mentioning a dataset path should be preserved."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm run trajectory analysis on ./my_data.h5ad"
            )
        assert result is not None
        assert "./my_data.h5ad" in result

    @pytest.mark.asyncio
    async def test_natural_language_with_model_name(self, handler):
        """Natural language mentioning a model name should be preserved."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm use geneformer to embed cells.h5ad"
            )
        assert result is not None
        assert "geneformer" in result

    @pytest.mark.asyncio
    async def test_single_word_query(self, handler):
        """Even a single-word query after /bio scfm should work."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command("/bio scfm analyze")
        assert result is not None
        assert "analyze" in result

    @pytest.mark.asyncio
    async def test_agent_routing_hint_present(self, handler):
        """The forwarded message should hint at SingleCellAgent tool usage."""
        with _patch_scfm_modules():
            result = await handler.handle_bio_command(
                "/bio scfm cluster my dataset"
            )
        assert result is not None
        assert "SingleCellAgent" in result


# ── Prompt template tests ────────────────────────────────────────────────


class TestScfmWorkflowPrompt:
    """Verify the SCFM workflow prompt generator."""

    def test_generate_preserves_user_query(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("annotate cell types in test.h5ad")
        assert "annotate cell types in test.h5ad" in msg

    def test_generate_adds_scfm_context(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("analyze my data")
        assert "scFM" in msg

    def test_generate_mentions_tool_routing(self):
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("find DEGs")
        assert "SingleCellAgent" in msg

    def test_prompt_is_concise(self):
        """The prompt should be concise, not a multi-phase instruction manual."""
        mod = _import_scfm_workflow()
        msg = mod.generate_scfm_workflow_message("analyze cells.h5ad")
        assert "PHASE 0" not in msg
        assert "PHASE 1" not in msg
        assert "PHASE 2" not in msg
        assert len(msg) < 500

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


# ── Integration with existing routing ────────────────────────────────────


class TestScfmIntegration:
    """Verify SCFM integrates correctly with the existing routing system."""

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_scrna(self, handler):
        result = await handler.handle_bio_command("/bio scrna init")
        assert result is not None
        assert "scRNA INIT MODE" in result

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_spatial(self, handler):
        result = await handler.handle_bio_command("/bio spatial init")
        assert result is not None
        assert "spatial INIT MODE" in result

    @pytest.mark.asyncio
    async def test_scfm_does_not_interfere_with_atac(self, handler):
        result = await handler.handle_bio_command("/bio atac init")
        assert result is not None
        assert "ATAC INIT MODE" in result

    @pytest.mark.asyncio
    async def test_bio_list_still_works(self, handler):
        result = await handler.handle_bio_command("/bio list")
        assert result == "bio list"

    @pytest.mark.asyncio
    async def test_bio_help_shows_scfm(self, handler):
        result = await handler.handle_bio_command("/bio")
        assert result is None
        calls = [str(c) for c in handler.console.print.call_args_list]
        assert any("scfm" in c for c in calls)


# ── Command map tests ────────────────────────────────────────────────────


class TestScfmCommandMap:
    """Verify SCFM entries in BIO_COMMAND_MAP and suggestions."""

    def test_scfm_in_command_map(self):
        from pantheon_cli.repl.bio_handler import BIO_COMMAND_MAP

        assert "scfm_init" in BIO_COMMAND_MAP
        assert "scfm_list_models" in BIO_COMMAND_MAP
        assert "scfm_list_analysis_types" in BIO_COMMAND_MAP

    def test_scfm_in_suggestions(self):
        from pantheon_cli.repl.bio_handler import get_bio_command_suggestions

        suggestions = get_bio_command_suggestions()
        assert "/bio scfm init" in suggestions
        assert "/bio scfm list_models" in suggestions
        assert "/bio scfm list_analysis_types" in suggestions
