import { agentService } from './agent-service.js';

export interface ChatIntegrationResult {
  shouldContinueWithNormalChat: boolean;
  agentResponse?: string;
  functionCall?: string;
  metadata?: Record<string, any>;
}

export class ChatAgentIntegration {
  async processUserInput(input: string): Promise<ChatIntegrationResult> {
    try {
      const result = await agentService.processQuery(input);
      
      if (result.agent) {
        // Found a matching agent
        return {
          shouldContinueWithNormalChat: true,
          agentResponse: `🤖 Matched agent "${result.agent.name}": ${result.agent.instructions}`,
          metadata: {
            agent_name: result.agent.name,
            agent_model: result.agent.model,
            source_path: result.agent.source_path
          }
        };
      } else if (result.handler && result.fallback) {
        // Found a function mapping
        return {
          shouldContinueWithNormalChat: true,
          agentResponse: `🔧 Intent matched: ${result.description}`,
          functionCall: result.handler,
          metadata: {
            intent_handler: result.handler,
            description: result.description
          }
        };
      }
      
      // No match found, continue with normal chat
      return {
        shouldContinueWithNormalChat: true
      };
    } catch (error) {
      console.error('Error in agent integration:', error);
      return {
        shouldContinueWithNormalChat: true
      };
    }
  }

  async getAgentSuggestions(): Promise<string[]> {
    const agents = await agentService.getAvailableAgents();
    return agents.map((agent: any) => `${agent.name}: ${agent.instructions.substring(0, 100)}...`);
  }
}

export const chatAgentIntegration = new ChatAgentIntegration();