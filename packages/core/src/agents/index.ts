/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

// Agent system exports
export {
  AgentRegistry,
  AgentDefinition,
  IntentMapping,
} from './agent-loader.js';
export { AgentService, agentService } from './agent-service.js';
export {
  ChatAgentIntegration,
  chatAgentIntegration,
  ChatIntegrationResult,
} from './chat-agent-integration.js';
export {
  testReferenceFolderDeletion,
  runDeletionTest,
} from './deletion-test.js';
