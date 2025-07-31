/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import { agentRegistry, AgentDefinition } from './agent-loader.js';
import { existsSync } from 'fs';
import { join } from 'path';

export class AgentService {
  private initialized = false;

  async initialize() {
    if (this.initialized) return;

    // Define reference folders to scan
    const referenceFolders = [
      join(process.cwd(), 'pantheon-single-cell'),
      join(process.cwd(), 'pantheon-agents'),
      // Also scan any existing agent directories
      join(process.cwd(), 'packages/core/src/agents'),
    ];

    // Filter to only existing folders
    const existingFolders = referenceFolders.filter((folder) =>
      existsSync(folder),
    );

    if (existingFolders.length > 0) {
      console.log(
        `Initializing agent registry with folders: ${existingFolders.join(', ')}`,
      );
      await agentRegistry.initialize(existingFolders);
      console.log(`Registered ${agentRegistry.getAllAgents().length} agents`);
    }

    this.initialized = true;
  }

  async processQuery(query: string): Promise<{
    agent?: AgentDefinition;
    handler?: string;
    description?: string;
    fallback?: boolean;
  }> {
    await this.initialize();

    const match = agentRegistry.matchIntent(query);

    if (match.agent) {
      return { agent: match.agent };
    } else if (match.handler) {
      return {
        handler: match.handler,
        description: match.description,
        fallback: true,
      };
    }

    return { fallback: true };
  }

  async getAvailableAgents(): Promise<AgentDefinition[]> {
    await this.initialize();
    return agentRegistry.getAllAgents();
  }

  addIntentMapping(pattern: RegExp, handler: string, description: string) {
    agentRegistry.addIntentMapping(pattern, handler, description);
  }
}

export const agentService = new AgentService();
