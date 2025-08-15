"""Configuration file loader for Pantheon CLI.

Loads user-defined prompts and instructions from local PANTHEON.md and PANTHEON.local.md files,
similar to Claude Code's CLAUDE.md functionality.
"""

import os
from pathlib import Path
from typing import Optional, Tuple


class ConfigLoader:
    """Loads configuration from PANTHEON.md and PANTHEON.local.md files."""
    
    CONFIG_FILES = ["PANTHEON.md", "PANTHEON.local.md"]
    
    def __init__(self, search_paths: Optional[list] = None):
        """Initialize the config loader.
        
        Args:
            search_paths: List of paths to search for config files. 
                         If None, uses current directory and parent directories.
        """
        self.search_paths = search_paths or self._get_default_search_paths()
    
    def _get_default_search_paths(self) -> list:
        """Get default search paths: current directory and all parent directories."""
        paths = []
        current = Path.cwd()
        
        # Add current directory and all parent directories
        for path in [current] + list(current.parents):
            paths.append(path)
        
        return paths
    
    def _find_config_files(self) -> dict:
        """Find all available config files in search paths.
        
        Returns:
            dict: Mapping of config filename to full path
        """
        found_files = {}
        
        for search_path in self.search_paths:
            for config_file in self.CONFIG_FILES:
                config_path = Path(search_path) / config_file
                if config_path.exists() and config_path.is_file():
                    # Use the first occurrence found (closest to current directory)
                    if config_file not in found_files:
                        found_files[config_file] = config_path
        
        return found_files
    
    def _load_file_content(self, file_path: Path) -> str:
        """Load content from a config file.
        
        Args:
            file_path: Path to the config file
            
        Returns:
            str: File content, or empty string if loading fails
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read().strip()
            return content
        except Exception as e:
            print(f"[Warning] Failed to load config file {file_path}: {e}")
            return ""
    
    def load_config(self) -> Tuple[str, list]:
        """Load configuration from PANTHEON.md and PANTHEON.local.md files.
        
        Returns:
            Tuple[str, list]: (combined_config_content, found_config_files)
        """
        found_files = self._find_config_files()
        
        if not found_files:
            return "", []
        
        config_sections = []
        found_file_paths = []
        
        # Load PANTHEON.md first (general config)
        if "PANTHEON.md" in found_files:
            file_path = found_files["PANTHEON.md"]
            content = self._load_file_content(file_path)
            if content:
                config_sections.append(f"# Configuration from {file_path.name}")
                config_sections.append(content)
                found_file_paths.append(str(file_path))
        
        # Load PANTHEON.local.md second (local overrides)
        if "PANTHEON.local.md" in found_files:
            file_path = found_files["PANTHEON.local.md"]
            content = self._load_file_content(file_path)
            if content:
                config_sections.append(f"# Local Configuration from {file_path.name}")
                config_sections.append(content)
                found_file_paths.append(str(file_path))
        
        # Combine all sections
        combined_config = "\n\n".join(config_sections) if config_sections else ""
        
        return combined_config, found_file_paths
    
    def get_config_for_agent(self) -> str:
        """Get formatted configuration content for agent instructions.
        
        Returns:
            str: Formatted configuration content ready to append to agent instructions
        """
        config_content, found_files = self.load_config()
        
        if not config_content:
            return ""
        
        # Format for agent instructions
        formatted_config = f"""

# User-Defined Configuration

{config_content}

# End User-Defined Configuration

"""
        return formatted_config
    
    def print_loaded_config_info(self) -> None:
        """Print information about loaded configuration files."""
        config_content, found_files = self.load_config()
        
        if found_files:
            print(f"📋 Loaded configuration from: {', '.join(found_files)}")
        else:
            # Silent when no config files are found - this is normal
            pass


def load_user_config(search_paths: Optional[list] = None) -> str:
    """Convenience function to load user configuration.
    
    Args:
        search_paths: Optional list of paths to search for config files
        
    Returns:
        str: Formatted configuration content for agent instructions
    """
    loader = ConfigLoader(search_paths)
    return loader.get_config_for_agent()


def print_config_info(search_paths: Optional[list] = None) -> None:
    """Convenience function to print config file info.
    
    Args:
        search_paths: Optional list of paths to search for config files
    """
    loader = ConfigLoader(search_paths)
    loader.print_loaded_config_info()