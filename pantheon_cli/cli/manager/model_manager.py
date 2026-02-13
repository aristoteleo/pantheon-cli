"""Model Management Module for Pantheon CLI"""

import json
from pathlib import Path
from typing import Optional, Tuple

from .api_key_manager import APIKeyManager

# Available models configuration
AVAILABLE_MODELS = {
    # OpenAI Models - GPT-5 Series (Latest)
    "gpt-5": "OpenAI GPT-5 (Latest)",
    "gpt-5-mini": "OpenAI GPT-5 Mini",
    "gpt-5-nano": "OpenAI GPT-5 Nano",
    "gpt-5-chat-latest": "OpenAI GPT-5 Chat Latest",
    # OpenAI Models - GPT-4 Series
    "gpt-4.1": "OpenAI GPT-4.1", 
    "gpt-4.1-mini": "OpenAI GPT-4.1 Mini",
    "gpt-4.1-nano": "OpenAI GPT-4.1 Nano",
    "gpt-4o": "OpenAI GPT-4o",
    "gpt-4o-2024-05-13": "OpenAI GPT-4o (2024-05-13)",
    "gpt-4o-audio-preview": "OpenAI GPT-4o Audio Preview",
    "gpt-4o-realtime-preview": "OpenAI GPT-4o Realtime Preview",
    "gpt-4o-mini": "OpenAI GPT-4o Mini",
    "gpt-4o-mini-audio-preview": "OpenAI GPT-4o Mini Audio Preview",
    "gpt-4o-mini-realtime-preview": "OpenAI GPT-4o Mini Realtime Preview",
    # OpenAI Models - o-Series (Reasoning)
    "o1": "OpenAI o1 (Reasoning)",
    "o1-pro": "OpenAI o1 Pro (Reasoning)",
    "o3-pro": "OpenAI o3 Pro (Reasoning)",
    "o3": "OpenAI o3 (Reasoning)",
    "o3-deep-research": "OpenAI o3 Deep Research",
    "o4-mini": "OpenAI o4 Mini (Reasoning)",
    "o4-mini-deep-research": "OpenAI o4 Mini Deep Research",
    "o3-mini": "OpenAI o3 Mini (Reasoning)",
    "o1-mini": "OpenAI o1 Mini (Reasoning)",
    # OpenAI Models - Codex Series
    "codex-mini-latest": "OpenAI Codex Mini Latest",
    # Anthropic Models - Claude 4 Series (Latest)
    "anthropic/claude-opus-4-1-20250805": "Claude Opus 4.1 (Latest)",
    "anthropic/claude-opus-4-20250514": "Claude Opus 4",
    "anthropic/claude-sonnet-4-20250514": "Claude Sonnet 4",
    "anthropic/claude-3-7-sonnet-20250219": "Claude Sonnet 3.7",
    "anthropic/claude-3-5-haiku-20241022": "Claude Haiku 3.5",
    # Anthropic Models - Claude 3 Series (Legacy)
    "anthropic/claude-3-opus-20240229": "Claude 3 Opus (Legacy)",
    "anthropic/claude-3-sonnet-20240229": "Claude 3 Sonnet (Legacy)", 
    "anthropic/claude-3-haiku-20240307": "Claude 3 Haiku (Legacy)",
    # Google Models
    "gemini/gemini-2.5-pro": "Gemini 2.5 Pro",
    "gemini/gemini-2.5-flash": "Gemini 2.5 Flash",
    "gemini/gemini-2.0-pro": "Gemini 2.0 Pro",
    "gemini/gemini-2.0-flash": "Gemini 2.0 Flash",
    "gemini/gemini-pro": "Gemini Pro",
    # DeepSeek Models
    "deepseek/deepseek-chat": "DeepSeek Chat",
    "deepseek/deepseek-reasoner": "DeepSeek Reasoner",
    # Qwen/Alibaba Models - Latest 2025 Series
    "qwq-plus": "QwQ Plus (Reasoning)",
    "qwen-max": "Qwen Max (Latest)",
    "qwen-max-latest": "Qwen Max Latest",
    "qwen-max-2025-01-25": "Qwen Max 2025-01-25",
    "qwen-plus": "Qwen Plus (Latest)",
    "qwen-plus-latest": "Qwen Plus Latest", 
    "qwen-plus-2025-04-28": "Qwen Plus 2025-04-28",
    "qwen-plus-2025-01-25": "Qwen Plus 2025-01-25",
    "qwen-turbo": "Qwen Turbo (Latest)",
    "qwen-turbo-latest": "Qwen Turbo Latest",
    "qwen-turbo-2025-04-28": "Qwen Turbo 2025-04-28", 
    "qwen-turbo-2024-11-01": "Qwen Turbo 2024-11-01",
    "qvq-max": "QVQ Max (Visual Reasoning)",
    "qvq-max-latest": "QVQ Max Latest",
    "qvq-max-2025-03-25": "QVQ Max 2025-03-25",
    # Qwen/Alibaba Models - Legacy
    "qwen/qwen-2.5-72b-instruct": "Qwen 2.5 72B (Legacy)",
    # Kimi/Moonshot Models - Latest K2 Series
    "moonshot/kimi-k2-0711-preview": "Kimi K2 (Preview)",
    "moonshot/kimi-k2-turbo-preview": "Kimi K2 Turbo (Preview)",
    # Kimi/Moonshot Models - Latest Series
    "moonshot/kimi-latest": "Kimi Latest (Auto Context)",
    "moonshot/kimi-latest-8k": "Kimi Latest 8K",
    "moonshot/kimi-latest-32k": "Kimi Latest 32K",
    "moonshot/kimi-latest-128k": "Kimi Latest 128K",
    # Kimi/Moonshot Models - Moonshot V1 Series
    "moonshot/moonshot-v1-8k": "Moonshot V1 8K",
    "moonshot/moonshot-v1-32k": "Moonshot V1 32K",
    "moonshot/moonshot-v1-128k": "Moonshot V1 128K",
    "moonshot/moonshot-v1-8k-vision-preview": "Moonshot V1 8K Vision",
    "moonshot/moonshot-v1-32k-vision-preview": "Moonshot V1 32K Vision",
    "moonshot/moonshot-v1-128k-vision-preview": "Moonshot V1 128K Vision",
    # Kimi/Moonshot Models - Thinking Series
    "moonshot/kimi-thinking-preview": "Kimi Thinking (Preview)",
    # Kimi/Moonshot Models - Legacy
    "moonshot/moonshot-v1-8k": "Kimi 8K (Legacy)",
    "moonshot/moonshot-v1-32k": "Kimi 32K (Legacy)", 
    "moonshot/moonshot-v1-128k": "Kimi 128K (Legacy)",
    # Grok/xAI Models
    "grok/grok-beta": "Grok Beta",
    "grok/grok-2": "Grok 2",
    # Zhipu AI/GLM Models (User-friendly names, internally converted to openai/ format)
    "zhipu/glm-4.5": "GLM-4.5 (Zhipu AI - Latest)",
    "zhipu/glm-4.5-air": "GLM-4.5 Air (Zhipu AI - Latest)",
    "zhipu/glm-4.5-flash": "GLM-4.5 Flash (Zhipu AI - Latest)",
    "zhipu/glm-4": "GLM-4 (Zhipu AI)",
    "zhipu/glm-4-plus": "GLM-4 Plus (Zhipu AI)", 
    "zhipu/glm-4-air": "GLM-4 Air (Zhipu AI)",
    "zhipu/glm-4-flash": "GLM-4 Flash (Zhipu AI - Free)",
    "zhipu/glm-4-long": "GLM-4 Long (Zhipu AI)",
    
    # Local/Other Models
    "ollama/llama3.2": "Llama 3.2 (Local)",
}

REASONING_EFFORT_LEVELS = {"minimal", "low", "medium", "high"}


class ModelManager:
    """Manages model selection and switching for Pantheon CLI"""
    
    def __init__(self, config_file_path: Path, api_key_manager: APIKeyManager):
        self.config_file_path = config_file_path
        self.api_key_manager = api_key_manager
        self.current_model = "gpt-5"
        self.current_agent = None
        self.reasoning_effort: Optional[str] = None
        self._load_model_config()

    def _load_model_config(self) -> str:
        """Load saved model configuration"""
        if self.config_file_path and self.config_file_path.exists():
            try:
                with open(self.config_file_path, 'r') as f:
                    config = json.load(f)
                    self.current_model = config.get('model', 'gpt-4.1')
                    effort = config.get('reasoning_effort')
                    if isinstance(effort, str) and effort.lower() in REASONING_EFFORT_LEVELS:
                        self.reasoning_effort = effort.lower()
            except Exception:
                pass
        return self.current_model

    def save_model_config(self, model: str):
        """Save current model configuration"""
        if not self.config_file_path:
            return
        
        # Load existing config to preserve API keys
        config = {'model': model}
        if self.config_file_path.exists():
            try:
                with open(self.config_file_path, 'r') as f:
                    config = json.load(f)
                    config['model'] = model
            except Exception:
                pass

        if self.reasoning_effort:
            config['reasoning_effort'] = self.reasoning_effort
        else:
            config.pop('reasoning_effort', None)

        try:
            with open(self.config_file_path, 'w') as f:
                json.dump(config, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save model config: {e}")

    def set_agent(self, agent):
        """Set the current agent reference for model updates"""
        self.current_agent = agent
        if hasattr(self.current_agent, "set_reasoning_effort"):
            applied_effort, _ = self._resolve_reasoning_effort_for_model(self.current_model)
            self.current_agent.set_reasoning_effort(applied_effort)


    def switch_model(self, new_model: str) -> str:
        """Switch to a new model"""
        if new_model not in AVAILABLE_MODELS:
            available = "\n".join([f"  {k}: {v}" for k, v in AVAILABLE_MODELS.items()])
            return f"❌ Model '{new_model}' not available. Available models:\n{available}"
        
        # Check API key availability
        key_available, key_message = self.api_key_manager.check_api_key_for_model(new_model)
        if not key_available:
            return f"❌ Cannot switch to {new_model}: {key_message}"
        
        old_model = self.current_model
        self.current_model = new_model

        # Update agent's model
        if self.current_agent:
            if isinstance(new_model, str):
                self.current_agent.models = [new_model]
                if new_model != "gpt-5-mini":
                    self.current_agent.models.append("gpt-5-mini")
            else:
                self.current_agent.models = new_model

            if hasattr(self.current_agent, "set_reasoning_effort"):
                applied_effort, inactive_preference = self._resolve_reasoning_effort_for_model(new_model)
                self.current_agent.set_reasoning_effort(applied_effort)

        # Save configuration
        self.save_model_config(new_model)

        extra_hint = ""
        allowed_levels = self._get_allowed_reasoning_levels(new_model)
        if allowed_levels:
            levels_str = "/".join(sorted(allowed_levels))
            applied_effort, inactive_preference = self._resolve_reasoning_effort_for_model(new_model)
            if inactive_preference and applied_effort:
                extra_hint = (
                    "\n🧠 This model supports reasoning effort"
                    f" (levels: {levels_str})."
                    f"\n   Requested '{inactive_preference}' isn't supported; currently using {applied_effort}."
                    "\n   Adjust with `/reasoning effort <minimal|low|medium|high>` before your next request."
                )
            else:
                display_effort = self.reasoning_effort or applied_effort or "medium"
                descriptor = "default" if self.reasoning_effort is None else "active"
                extra_hint = (
                    "\n🧠 This model supports reasoning effort"
                    f" (levels: {levels_str})."
                    f"\n   Current {descriptor} setting: {display_effort}."
                    "\n   Adjust with `/reasoning effort <minimal|low|medium|high>` before your next request."
                )
        elif self.reasoning_effort:
            extra_hint = (
                f"\n⚠️ Reasoning effort '{self.reasoning_effort}' is stored but this model does not"
                " support custom effort levels."
            )

        return (
            f"✅ Switched from {AVAILABLE_MODELS.get(old_model, old_model)} to {AVAILABLE_MODELS[new_model]} ({new_model})"
            f"\nℹ️ {key_message}{extra_hint}"
        )

    def list_models(self) -> str:
        """List all available models with API key status"""
        result = "🤖 Available Models (Top 3 per provider):\n\n"

        # Group models by provider with correct categorization
        providers = {}
        for model_id, description in AVAILABLE_MODELS.items():
            # Determine provider based on model naming patterns
            if model_id.startswith("anthropic/"):
                provider = "Anthropic"
            elif model_id.startswith(("qwq-", "qwen-", "qvq-")) or model_id.startswith("qwen/"):
                provider = "Qwen"
            elif model_id.startswith(("kimi-", "moonshot-")) or model_id.startswith("moonshot/"):
                provider = "Kimi"
            elif model_id.startswith("grok/"):
                provider = "Grok"
            elif model_id.startswith("gemini/"):
                provider = "Google"
            elif model_id.startswith("deepseek/"):
                provider = "DeepSeek"
            elif model_id.startswith("ollama/"):
                provider = "Local"
            elif model_id.startswith("zhipu/"):
                provider = "Zhipu AI"
            else:
                provider = "OpenAI"
            
            if provider not in providers:
                providers[provider] = []
            providers[provider].append((model_id, description))
        
        for provider, models in providers.items():
            result += f"{provider}:\n"
            # Show only first 3 models per provider
            top_models = models[:3]
            for model_id, description in top_models:
                current_indicator = " ← Current" if model_id == self.current_model else ""
                
                # Check API key status
                key_available, _ = self.api_key_manager.check_api_key_for_model(model_id)
                from .api_key_manager import PROVIDER_API_KEYS
                if PROVIDER_API_KEYS.get(model_id) is None:
                    key_status = " 🟢"  # Green circle for no key needed
                elif key_available:
                    key_status = " ✅"  # Checkmark for available key
                else:
                    key_status = " ❌"  # X for missing key

                reasoning_hint = ""
                allowed_levels = self._get_allowed_reasoning_levels(model_id)
                if allowed_levels:
                    levels_str = "/".join(sorted(allowed_levels))
                    if model_id == self.current_model:
                        applied_effort, inactive_preference = self._resolve_reasoning_effort_for_model(model_id)
                        if inactive_preference and applied_effort:
                            reasoning_hint = (
                                f" (effort: {self.reasoning_effort} inactive; using {applied_effort}; levels: {levels_str})"
                            )
                        else:
                            display_effort = self.reasoning_effort or applied_effort or "medium"
                            reasoning_hint = f" (effort: {display_effort}; levels: {levels_str})"
                    else:
                        reasoning_hint = f" (supports effort levels: {levels_str})"

                result += (
                    f"  • {model_id}: {description}{key_status}{current_indicator}{reasoning_hint}\n"
                )
            
            # Show count if there are more models
            if len(models) > 3:
                result += f"  ... and {len(models) - 3} more models\n"
            result += "\n"
        
        result += "Legend: 🟢 No API key needed | ✅ API key available | ❌ API key missing\n\n"
        result += f"💡 Usage: /model <model_id> | /api-key <provider> <key>\n"
        result += f"📝 Current: {AVAILABLE_MODELS.get(self.current_model, self.current_model)} ({self.current_model})"
        
        return result

    def get_current_model_status(self) -> str:
        """Get current model with API key status"""
        key_available, key_message = self.api_key_manager.check_api_key_for_model(self.current_model)
        key_status = "✅" if key_available else "❌"
        allowed_levels = self._get_allowed_reasoning_levels()
        reasoning_info = ""
        if allowed_levels:
            applied_effort, inactive_preference = self._resolve_reasoning_effort_for_model()
            levels_str = "/".join(sorted(allowed_levels))
            if inactive_preference and applied_effort:
                reasoning_line = (
                    f"🧠 Reasoning effort: {inactive_preference} (inactive; using {applied_effort})."
                )
            else:
                if self.reasoning_effort:
                    reasoning_line = f"🧠 Reasoning effort: {self.reasoning_effort}."
                else:
                    reasoning_line = f"🧠 Reasoning effort: {applied_effort} (default)."
            reasoning_info = (
                f"\n{reasoning_line}\n   Supported levels: {levels_str}."
                "\n   Adjust with `/reasoning effort <minimal|low|medium|high>`."
            )
        elif self.reasoning_effort:
            reasoning_info = (
                f"\n⚠️ Stored reasoning effort '{self.reasoning_effort}' is inactive"
                " for this model (no custom effort support)."
            )

        return (
            f"📱 Current Model: {AVAILABLE_MODELS.get(self.current_model, self.current_model)} ({self.current_model})"
            f"\n{key_status} {key_message}{reasoning_info}"
        )

    def _flatten_model_ids(self, model: Optional[str]) -> list[str]:
        if model is None:
            return []
        if isinstance(model, str):
            return [model]
        if isinstance(model, (list, tuple)):
            result: list[str] = []
            for item in model:
                result.extend(self._flatten_model_ids(item))
            return result
        return [str(model)]

    @staticmethod
    def _normalize_model_name(model_id: str) -> str:
        normalized = str(model_id)
        if "/" in normalized:
            parts = normalized.split("/", 1)
            if len(parts) == 2:
                normalized = parts[1]
        return normalized

    def _allowed_levels_for_model_id(self, model_id: str) -> set[str]:
        normalized = self._normalize_model_name(model_id)
        if normalized.startswith("gpt-5"):
            return {"minimal", "low", "medium", "high"}
        if normalized.startswith(("o1", "o3", "o4")):
            return {"low", "medium", "high"}
        return set()

    def _get_allowed_reasoning_levels(self, model: Optional[str] = None) -> set[str]:
        allowed: set[str] = set()
        for model_id in self._flatten_model_ids(model if model is not None else self.current_model):
            allowed |= self._allowed_levels_for_model_id(model_id)
        return allowed

    def _resolve_reasoning_effort_for_model(self, model: Optional[str] = None) -> Tuple[Optional[str], Optional[str]]:
        allowed_levels = self._get_allowed_reasoning_levels(model)
        if not allowed_levels:
            return None, None
        preferred = self.reasoning_effort.lower() if isinstance(self.reasoning_effort, str) else None
        if preferred and preferred in allowed_levels:
            return preferred, None
        fallback = "medium" if "medium" in allowed_levels else sorted(allowed_levels)[0]
        return fallback, preferred

    def supports_reasoning_effort(self, model: Optional[str] = None) -> bool:
        """Check if reasoning effort is supported for the given model"""
        return bool(self._get_allowed_reasoning_levels(model))

    def set_reasoning_effort(self, effort: Optional[str]) -> str:
        """Update reasoning effort and sync to agent/config"""
        allowed_levels = self._get_allowed_reasoning_levels()

        if effort is None:
            self.reasoning_effort = None
        else:
            effort_lower = effort.lower()
            if effort_lower not in REASONING_EFFORT_LEVELS:
                allowed = ", ".join(sorted(REASONING_EFFORT_LEVELS))
                return f"❌ Invalid effort '{effort}'. Choose from: {allowed}."
            if not allowed_levels:
                return "⚠️ Current model does not support reasoning effort adjustments."
            if effort_lower not in allowed_levels:
                allowed = "/".join(sorted(allowed_levels))
                return (
                    f"⚠️ Effort '{effort_lower}' not supported by {self.current_model}."
                    f" Available levels: {allowed}."
                )
            self.reasoning_effort = effort_lower

        if self.current_agent and hasattr(self.current_agent, "set_reasoning_effort"):
            applied_effort, _ = self._resolve_reasoning_effort_for_model()
            self.current_agent.set_reasoning_effort(applied_effort)

        # Persist new value alongside model configuration
        self.save_model_config(self.current_model)

        allowed_levels = self._get_allowed_reasoning_levels()
        if self.reasoning_effort:
            return (
                f"✅ Reasoning effort set to {self.reasoning_effort}."
                f" Supported levels: {'/'.join(sorted(allowed_levels))}."
            )
        if allowed_levels:
            return (
                "✅ Reasoning effort cleared. Using default level "
                f"{self._resolve_reasoning_effort_for_model()[0]} for current model."
            )
        return "✅ Reasoning effort cleared (current model does not support custom effort)."

    def handle_model_command(self, command: str) -> str:
        """Handle /model commands"""
        parts = command.strip().split()
        
        if len(parts) == 1:  # Just "/model"
            return self.list_models()
        
        subcommand = parts[1].lower()
        
        if subcommand == "list":
            return self.list_models()
        elif subcommand == "current":
            return self.get_current_model_status()
        elif subcommand in AVAILABLE_MODELS:
            return self.switch_model(subcommand)
        else:
            # Try to match partial model names
            matches = [m for m in AVAILABLE_MODELS.keys() if subcommand in m.lower()]
            if len(matches) == 1:
                return self.switch_model(matches[0])
            elif len(matches) > 1:
                match_list = "\n".join([f"  • {m}: {AVAILABLE_MODELS[m]}" for m in matches])
                return f"🔍 **Multiple matches found:**\n{match_list}\n\n💡 Use the full model ID: `/model <model_id>`"
            else:
                return f"❌ Model '{subcommand}' not found. Use `/model list` to see available models."

    def handle_reasoning_command(self, command: str) -> str:
        """Handle /reasoning commands"""
        parts = command.strip().split()

        if len(parts) == 1 or parts[1].lower() in {"help", "?"}:
            return (
                "🧠 Reasoning Effort Control:\n"
                "  /reasoning status       - Show current effort and levels\n"
                "  /reasoning effort high  - Set effort to high\n"
                "  /reasoning minimal      - Shortcut to minimal (GPT-5 only)\n"
                "  /reasoning clear        - Remove custom setting\n"
                "\n💡 Supported models: GPT-5 reasoning models (minimal/low/medium/high)"
                " and OpenAI o-series (low/medium/high)."
            )

        subcommand = parts[1].lower()

        if subcommand == "status" or subcommand == "current":
            allowed_levels = self._get_allowed_reasoning_levels()
            applied_effort, inactive_preference = self._resolve_reasoning_effort_for_model()
            if allowed_levels:
                levels_str = "/".join(sorted(allowed_levels))
                if inactive_preference and applied_effort:
                    current_desc = f"{inactive_preference} (inactive; using {applied_effort})"
                elif self.reasoning_effort:
                    current_desc = self.reasoning_effort
                elif applied_effort:
                    current_desc = f"{applied_effort} (default)"
                else:
                    current_desc = "default"
                return (
                    f"🧠 Reasoning effort: {current_desc}."
                    f" Supported levels: {levels_str}."
                )
            return "🧠 Reasoning effort: not supported by current model."

        if subcommand == "clear":
            return self.set_reasoning_effort(None)

        if subcommand in REASONING_EFFORT_LEVELS:
            if not self.supports_reasoning_effort():
                return "⚠️ Current model does not support reasoning effort adjustments."
            return self.set_reasoning_effort(subcommand)

        if subcommand in {"effort", "set"} and len(parts) > 2:
            level = parts[2].lower()
            if level not in REASONING_EFFORT_LEVELS:
                allowed = ", ".join(sorted(REASONING_EFFORT_LEVELS))
                return f"❌ Invalid effort '{level}'. Choose from: {allowed}."
            if not self.supports_reasoning_effort():
                return "⚠️ Current model does not support reasoning effort adjustments."
            return self.set_reasoning_effort(level)

        return "❌ Unknown /reasoning command. Use `/reasoning help` for options."
