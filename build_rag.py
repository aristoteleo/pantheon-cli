#!/usr/bin/env python
"""Build RAG database with Pantheon API keys automatically loaded"""

import sys
import os
from pathlib import Path

# Load Pantheon API keys first
sys.path.insert(0, str(Path(__file__).parent))
from pantheon_cli.utils.load_global_keys import load_pantheon_api_keys

print("Loading Pantheon API keys...")
load_pantheon_api_keys()

# Now run the RAG build command
if __name__ == "__main__":
    import subprocess
    
    # Default paths
    config_file = Path(__file__).parent / "pantheon_cli/cli/rag_system_config.yaml"
    output_dir = "tmp/pantheon_cli_tools_rag"
    
    # Allow custom paths from command line
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    
    print(f"\nBuilding RAG database...")
    print(f"Config: {config_file}")
    print(f"Output: {output_dir}\n")
    
    # Run the build command
    cmd = [sys.executable, "-m", "pantheon.toolsets.utils.rag", "build", 
           str(config_file), str(output_dir)]
    
    result = subprocess.run(cmd)
    sys.exit(result.returncode)