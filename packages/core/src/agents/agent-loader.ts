import { existsSync, readdirSync, statSync } from 'fs';
import { join, extname } from 'path';
import { spawn } from 'child_process';

export interface AgentDefinition {
  name: string;
  instructions: string;
  model?: string;
  icon?: string;
  tools?: string[];
  source_path: string;
}

export interface IntentMapping {
  pattern: RegExp;
  handler: string;
  description: string;
}

export class AgentRegistry {
  private agents: Map<string, AgentDefinition> = new Map();
  private intentMappings: IntentMapping[] = [];
  private initialized = false;

  constructor() {
    this.setupDefaultIntentMappings();
  }

  private setupDefaultIntentMappings() {
    this.intentMappings = [
      {
        pattern: /normalize.*counts?\.h5ad/i,
        handler: 'scanpy.pp.normalize_total',
        description: 'Normalize single-cell RNA-seq counts'
      },
      {
        pattern: /log.*transform/i,
        handler: 'scanpy.pp.log1p',
        description: 'Log transform expression data'
      },
      {
        pattern: /map.*reads?.*to.*(hg38|mm10|genome)/i,
        handler: 'star_aligner',
        description: 'Map sequencing reads to reference genome'
      },
      {
        pattern: /cluster.*cells?/i,
        handler: 'scanpy.tl.leiden',
        description: 'Cluster single cells'
      },
      {
        pattern: /find.*marker.*genes?/i,
        handler: 'scanpy.tl.rank_genes_groups',
        description: 'Find marker genes for clusters'
      },
      {
        pattern: /umap.*visualization/i,
        handler: 'scanpy.tl.umap',
        description: 'Generate UMAP visualization'
      }
    ];
  }

  async initialize(referenceFolders: string[] = []) {
    if (this.initialized) return;

    for (const folder of referenceFolders) {
      if (existsSync(folder)) {
        await this.scanAgentFolder(folder);
      }
    }

    this.initialized = true;
  }

  private async scanAgentFolder(folderPath: string) {
    const files = readdirSync(folderPath);
    
    for (const file of files) {
      const filePath = join(folderPath, file);
      const stat = statSync(filePath);
      
      if (stat.isDirectory()) {
        await this.scanAgentFolder(filePath);
      } else if (extname(file) === '.py') {
        await this.extractPythonAgents(filePath);
      }
    }
  }

  private async extractPythonAgents(filePath: string): Promise<void> {
    return new Promise((resolve, reject) => {
      const process = spawn('python3', ['-c', `
import ast
import sys
import os

def extract_agent_info(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        tree = ast.parse(content)
        agents = []
        
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if hasattr(node.func, 'id') and node.func.id == 'Agent':
                    agent_info = {'source_path': file_path}
                    
                    for keyword in node.keywords:
                        if keyword.arg in ['name', 'instructions', 'model', 'icon']:
                            if isinstance(keyword.value, ast.Constant):
                                agent_info[keyword.arg] = keyword.value.value
                    
                    if 'name' in agent_info:
                        agents.append(agent_info)
        
        return agents
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return []

if __name__ == "__main__":
    file_path = sys.argv[1]
    agents = extract_agent_info(file_path)
    import json
    print(json.dumps(agents))
      `, filePath]);

      let stdout = '';
      let stderr = '';

      process.stdout.on('data', (data) => {
        stdout += data.toString();
      });

      process.stderr.on('data', (data) => {
        stderr += data.toString();
      });

      process.on('close', (code) => {
        if (code === 0) {
          try {
            const agents = JSON.parse(stdout);
            agents.forEach((agent: AgentDefinition) => {
              this.registerAgent(agent);
            });
            resolve();
          } catch (e) {
            console.warn(`Failed to parse agent info from ${filePath}:`, e);
            resolve();
          }
        } else {
          console.warn(`Failed to extract agents from ${filePath}:`, stderr);
          resolve();
        }
      });

      process.on('error', (err) => {
        console.warn(`Error running Python script for ${filePath}:`, err);
        resolve();
      });
    });
  }

  registerAgent(agent: AgentDefinition) {
    this.agents.set(agent.name, agent);
  }

  getAgent(name: string): AgentDefinition | undefined {
    return this.agents.get(name);
  }

  getAllAgents(): AgentDefinition[] {
    return Array.from(this.agents.values());
  }

  matchIntent(query: string): { agent?: AgentDefinition; handler?: string; description?: string } {
    // First try to match registered agents by name or instructions
    for (const agent of this.agents.values()) {
      if (query.toLowerCase().includes(agent.name.toLowerCase()) ||
          (agent.instructions && query.toLowerCase().includes(agent.instructions.toLowerCase().substring(0, 50)))) {
        return { agent };
      }
    }

    // Then try intent mappings
    for (const mapping of this.intentMappings) {
      if (mapping.pattern.test(query)) {
        return { 
          handler: mapping.handler, 
          description: mapping.description 
        };
      }
    }

    return {};
  }

  addIntentMapping(pattern: RegExp, handler: string, description: string) {
    this.intentMappings.push({ pattern, handler, description });
  }
}

export const agentRegistry = new AgentRegistry();