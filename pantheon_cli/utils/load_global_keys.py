#!/usr/bin/env python
"""Load API keys from Pantheon global config and set them as environment variables"""

import os
import json
from pathlib import Path
import sys

# Add parent directory to path to import api_key_manager
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from pantheon_cli.cli.manager.api_key_manager import APIKeyManager


def load_pantheon_api_keys():
    """Load API keys from Pantheon global/local config and set them as environment variables"""
    
    # Initialize API key manager with dummy local path
    local_config = Path.cwd() / ".pantheon_config.json"
    api_manager = APIKeyManager(local_config)
    
    # The manager already loads and decrypts keys in __init__
    # Just sync them to environment
    api_manager.sync_environment_variables()
    
    # Print loaded keys for confirmation
    loaded_keys = []
    for key in ["OPENAI_API_KEY", "ANTHROPIC_API_KEY", "GOOGLE_API_KEY", 
                "DEEPSEEK_API_KEY", "QWEN_API_KEY", "MOONSHOT_API_KEY", "GROK_API_KEY"]:
        if os.environ.get(key):
            loaded_keys.append(key)
    
    if loaded_keys:
        print(f"✅ Loaded API keys from Pantheon config: {', '.join(loaded_keys)}")
        return True
    else:
        print("⚠️  No API keys found in Pantheon config")
        return False


if __name__ == "__main__":
    load_pantheon_api_keys()