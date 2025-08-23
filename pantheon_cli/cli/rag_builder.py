"""Simple RAG database builder for Pantheon CLI"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Optional
from .manager.api_key_manager import APIKeyManager


class RAGBuilder:
    """Build RAG database with automatic API key loading from Pantheon config"""
    
    def __init__(self):
        # Initialize API key manager to load keys
        local_config = Path.cwd() / ".pantheon_config.json"
        self.api_manager = APIKeyManager(local_config)
        self.api_manager.sync_environment_variables()
    
    def build(self, config_path: Optional[str] = None, output_dir: Optional[str] = None) -> bool:
        """
        Build RAG database
        
        Args:
            config_path: Path to RAG config YAML file (default: use built-in config)
            output_dir: Output directory for RAG database (default: tmp/pantheon_cli_tools_rag)
        
        Returns:
            bool: True if build successful, False otherwise
        """
        # Use defaults if not provided
        if config_path is None:
            config_path = Path(__file__).parent / "rag_system_config.yaml"
        if output_dir is None:
            output_dir = "tmp/pantheon_cli_tools_rag"
        
        # Check if OpenAI key is available
        if not os.environ.get("OPENAI_API_KEY"):
            print("❌ OpenAI API key is required for building RAG database")
            print("💡 Please set it using: /api-key openai sk-your-key-here")
            return False
        
        print(f"🔨 Building RAG database...")
        print(f"   Config: {config_path}")
        print(f"   Output: {output_dir}")
        print(f"   API Key: ✅ OpenAI key loaded from Pantheon config\n")
        
        # Build the command
        cmd = [
            sys.executable, "-m", 
            "pantheon.toolsets.utils.rag", "build",
            str(config_path), str(output_dir)
        ]
        
        try:
            # Run the build command
            result = subprocess.run(cmd, capture_output=False, text=True)
            
            if result.returncode == 0:
                print("\n✅ RAG database built successfully!")
                print(f"📁 Database location: {output_dir}")
                return True
            else:
                print("\n❌ Failed to build RAG database")
                return False
                
        except Exception as e:
            print(f"\n❌ Error building RAG database: {e}")
            return False
    
    def check_status(self) -> dict:
        """Check RAG build status and requirements"""
        status = {
            "openai_key": bool(os.environ.get("OPENAI_API_KEY")),
            "config_exists": (Path(__file__).parent / "rag_system_config.yaml").exists(),
            "default_db_exists": Path("tmp/pantheon_cli_tools_rag").exists()
        }
        
        return status
    
    def handle_command(self, command: str) -> str:
        """Handle /rag-build commands"""
        parts = command.strip().split()
        
        if len(parts) == 1 or parts[1] == "status":
            # Show status
            status = self.check_status()
            result = "📊 RAG Database Status:\n\n"
            result += f"{'✅' if status['openai_key'] else '❌'} OpenAI API Key: {'Available' if status['openai_key'] else 'Not configured'}\n"
            result += f"{'✅' if status['config_exists'] else '❌'} Config File: {'Found' if status['config_exists'] else 'Missing'}\n"
            result += f"{'✅' if status['default_db_exists'] else '❌'} Default Database: {'Built' if status['default_db_exists'] else 'Not built'}\n"
            
            if not status['openai_key']:
                result += "\n💡 To build RAG database, first set OpenAI key:\n"
                result += "   /api-key openai sk-your-key-here\n"
            elif not status['default_db_exists']:
                result += "\n💡 To build RAG database, run:\n"
                result += "   /rag-build\n"
            
            return result
        
        elif parts[1] == "build" or len(parts) == 1:
            # Build with optional custom paths
            config_path = None
            output_dir = None
            
            if len(parts) > 2:
                config_path = parts[2]
            if len(parts) > 3:
                output_dir = parts[3]
            
            success = self.build(config_path, output_dir)
            
            if success:
                return "✅ RAG database built successfully!"
            else:
                return "❌ Failed to build RAG database. Check the error messages above."
        
        elif parts[1] == "help":
            return """🔧 RAG Database Builder Commands:

/rag-build              - Build RAG database with default config
/rag-build status       - Check RAG database status
/rag-build help         - Show this help message
/rag-build <config> <output>  - Build with custom config and output dir

Examples:
  /rag-build                     # Build with defaults
  /rag-build my_config.yaml      # Use custom config
  /rag-build my_config.yaml ./rag_db  # Custom config and output
"""
        
        else:
            return f"❌ Unknown subcommand '{parts[1]}'. Try: /rag-build help"


# Convenience function for CLI integration
def build_rag_database(config_path: Optional[str] = None, output_dir: Optional[str] = None) -> bool:
    """Build RAG database with automatic API key loading"""
    builder = RAGBuilder()
    return builder.build(config_path, output_dir)