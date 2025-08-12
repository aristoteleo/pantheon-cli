"""Pantheon REPL - Entry point for python -m pantheon.repl"""

import asyncio
from pantheon.agent import Agent
from .core import Repl

def main():
    """Main entry point for pantheon-repl command"""
    agent = Agent(
        "agent",
        "You are a helpful assistant."
    )
    repl = Repl(agent)
    asyncio.run(repl.run())

if __name__ == "__main__":
    main()