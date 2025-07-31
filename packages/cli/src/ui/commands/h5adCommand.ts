/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import {
  SlashCommand,
  CommandKind,
  CommandContext,
  ToolActionReturn,
} from './types.js';
import { MessageType } from '../types.js';

export const h5adCommand: SlashCommand = {
  name: 'h5ad',
  description:
    'Analyze h5ad single-cell files with RAG-enhanced code generation',
  kind: CommandKind.BUILT_IN,
  action: async (
    context: CommandContext,
    args: string,
  ): Promise<ToolActionReturn | void> => {
    const argList = args.trim().split(/\s+/);

    if (argList.length < 2) {
      context.ui.addItem(
        {
          type: MessageType.INFO,
          text: `**Usage:** \`/h5ad <file_path> "<analysis_task>" [--params parameters] [--output directory]\`

**Examples:**
- \`/h5ad data.h5ad "perform basic preprocessing and clustering"\`
- \`/h5ad pbmc3k.h5ad "find marker genes and create UMAP visualization"\`
- \`/h5ad sample.h5ad "trajectory analysis with RNA velocity" --output results/\`

**Description:**
This command analyzes h5ad single-cell files using RAG-enhanced dynamic code generation. It combines knowledge from the RAG system with the Gemini API to generate and execute appropriate analysis code.`,
        },
        Date.now(),
      );
      return;
    }

    const filePath = argList[0];
    let analysisTask = argList[1];

    // Task mapping for short commands
    const taskMappings: Record<string, string> = {
      qc: '"perform quality control and generate QC plots"',
    };

    if (analysisTask in taskMappings) {
      analysisTask = taskMappings[analysisTask];
    }

    // Parse optional parameters
    let parameters: string | undefined;
    let outputDir: string | undefined;

    for (let i = 2; i < argList.length; i++) {
      if (argList[i] === '--output' && i + 1 < argList.length) {
        outputDir = argList[i + 1];
        i++; // Skip next arg
      } else if (argList[i] === '--params' && i + 1 < argList.length) {
        parameters = argList[i + 1];
        i++; // Skip next arg
      }
    }

    // Return tool action to execute the H5AD analyzer
    const toolArgs: Record<string, unknown> = {
      h5adPath: filePath,
      analysisTask,
    };

    if (parameters) toolArgs.parameters = parameters;
    if (outputDir) toolArgs.outputDir = outputDir;

    return {
      type: 'tool',
      toolName: 'h5ad_analyzer',
      toolArgs,
    };
  },
};
