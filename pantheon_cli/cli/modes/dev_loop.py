"""Developer loop mode: plan → code → review."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any
import json
import re

from rich.console import Console

from pantheon.agent import Agent
from pantheon.memory import Memory
from pantheon.utils.log import logger


@dataclass
class DevLoopResult:
    success: bool
    iterations: int
    final_message: str


def _safe_json_extract(text: str) -> Optional[dict]:
    """Try to parse JSON, falling back to the first balanced {...} block."""

    text = (text or "").strip()
    if not text:
        return None
    try:
        return json.loads(text)
    except Exception:
        pass

    stack = 0
    start = -1
    for idx, ch in enumerate(text):
        if ch == "{":
            if stack == 0:
                start = idx
            stack += 1
        elif ch == "}":
            if stack > 0:
                stack -= 1
                if stack == 0 and start != -1:
                    candidate = text[start : idx + 1]
                    try:
                        return json.loads(candidate)
                    except Exception:
                        start = -1
                        continue
    return None


async def run_devloop(
    agent: Agent,
    console: Console,
    workspace_path: Path,
    goal: str,
    max_iters: int = 10,
) -> DevLoopResult:
    """Run the developer loop until completion or iteration limit."""

    memory = Memory(name="dev-loop")

    PLAN_MODEL, PLAN_EFFORT = "gpt-5", "high"
    CODE_MODEL = "gpt-4.1"
    REVIEW_MODEL, REVIEW_EFFORT = "gpt-5", "minimal"

    console.print("\n[bold cyan]🔁 Developer Loop[/bold cyan] plan → code → review")
    console.print(f"Workspace: {workspace_path}")
    console.print(f"Goal: {goal}")

    next_task: Optional[str] = None
    acceptance: Optional[Any] = None
    iters = 0

    while iters < max_iters:
        iters += 1
        console.print(f"\n[bold]Iteration {iters}/{max_iters}[/bold]")

        # ----------------- PLAN -----------------
        agent.set_reasoning_effort(PLAN_EFFORT)
        plan_prompt = (
            "You are the Planner. Break the overall goal into the smallest next actionable subtask.\n"
            "Respond STRICTLY in JSON with keys: next_task (str), acceptance (list[str]), rationale (str), done (bool).\n"
            f"Goal: {goal}\n"
            "Use existing context in memory to avoid repetition."
        )
        plan_resp = await agent.run(
            plan_prompt,
            model=PLAN_MODEL,
            tool_use=False,
            memory=memory,
            update_memory=True,
        )
        plan_text = str(plan_resp.content or "").strip()
        plan_json = _safe_json_extract(plan_text) or {}

        done = bool(plan_json.get("done")) if isinstance(plan_json, dict) else False
        if done:
            console.print("[green]Planner marked goal complete.[/green]")
            return DevLoopResult(True, iters, "Planner concluded task is complete.")

        next_task = plan_json.get("next_task") if isinstance(plan_json, dict) else None
        acceptance = plan_json.get("acceptance") if isinstance(plan_json, dict) else None

        if not next_task:
            match = re.search(r"next_task\s*[:=]\s*\"([^\"]+)\"", plan_text, re.IGNORECASE)
            if match:
                next_task = match.group(1)

        if not next_task:
            if plan_text:
                console.print(f"[yellow]Planner raw output:[/yellow] {plan_text}")
            else:
                console.print("[yellow]Planner returned empty output; check API access or model configuration.[/yellow]")
            logger.warning("Planner did not provide a valid next_task. Stopping.")
            return DevLoopResult(False, iters, "Planner did not produce a next_task.")

        # ----------------- CODE -----------------
        agent.set_reasoning_effort(None)
        code_prompt = (
            "You are the Coder. Implement the subtask using available tools (file_editor, file_manager, code_search, python_interpreter).\n"
            f"Workspace root: {workspace_path}.\n"
            "- Make precise, minimal changes.\n"
            "- Summarize modifications at the end in JSON: {\"summary\": str}.\n"
            f"Subtask: {next_task}\n"
            f"Acceptance criteria: {acceptance or '[]'}\n"
        )
        await agent.run(
            code_prompt,
            model=CODE_MODEL,
            tool_use=True,
            memory=memory,
            update_memory=True,
        )

        # ----------------- REVIEW -----------------
        agent.set_reasoning_effort(REVIEW_EFFORT)
        review_prompt = (
            "You are the Reviewer. Assess whether the latest changes meet the acceptance criteria.\n"
            "Respond STRICTLY in JSON with keys: passed (bool), issues (list[str]), next_hint (str|null).\n"
            f"Subtask: {next_task}\n"
            f"Acceptance criteria: {acceptance or '[]'}"
        )
        review_resp = await agent.run(
            review_prompt,
            model=REVIEW_MODEL,
            tool_use=False,
            memory=memory,
            update_memory=True,
        )
        review_text = str(review_resp.content or "").strip()
        review_json = _safe_json_extract(review_text) or {}
        passed = bool(review_json.get("passed")) if isinstance(review_json, dict) else False

        if passed:
            console.print("[green]✅ Subtask accepted by Reviewer.[/green]")
            continue

        hint = None
        if isinstance(review_json, dict):
            hint = review_json.get("next_hint")
        if hint:
            memory.add_messages([
                {"role": "user", "content": f"Reviewer hint for next subtask: {hint}"}
            ])
        console.print("[yellow]⚠️ Subtask not accepted; iterating based on reviewer feedback.[/yellow]")

    return DevLoopResult(False, iters, "Reached iteration limit without completion.")
