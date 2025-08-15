"""Pantheon-CLI: Command Line Interface for Pantheon Agent Framework"""

__version__ = "0.1.1"

# Entry point function for CLI command
def cli_main():
    """Entry point for pantheon-cli command"""
    from .cli.core import cli
    cli()