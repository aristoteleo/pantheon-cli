import { readdirSync, statSync, readFileSync } from 'fs';
import { join, extname } from 'path';

export interface NotebookCell {
  cell_type: 'code' | 'markdown';
  source: string[];
  metadata?: Record<string, any>;
  outputs?: any[];
}

export interface NotebookData {
  cells: NotebookCell[];
  metadata?: Record<string, any>;
}

export interface RagSnippet {
  id: string;
  content: string;
  metadata: {
    filename: string;
    cell_index: number;
    cell_type: 'code' | 'markdown';
    notebook_path: string;
  };
}

export interface RagClient {
  indexSnippet(id: string, content: string, metadata: Record<string, any>): Promise<void>;
}

export class NotebookRagUpdater {
  constructor(private ragClient: RagClient) {}

  async updateFromDirectories(directories: string[]): Promise<void> {
    const notebooks = this.findNotebooks(directories);
    
    for (const notebookPath of notebooks) {
      await this.processNotebook(notebookPath);
    }
  }

  private findNotebooks(directories: string[]): string[] {
    const notebooks: string[] = [];
    
    for (const dir of directories) {
      this.collectNotebooks(dir, notebooks);
    }
    
    return notebooks;
  }

  private collectNotebooks(directory: string, notebooks: string[]) {
    try {
      const files = readdirSync(directory);
      
      for (const file of files) {
        const filePath = join(directory, file);
        const stat = statSync(filePath);
        
        if (stat.isDirectory()) {
          this.collectNotebooks(filePath, notebooks);
        } else if (extname(file) === '.ipynb') {
          notebooks.push(filePath);
        }
      }
    } catch (error) {
      console.warn(`Warning: Could not read directory ${directory}:`, error);
    }
  }

  private async processNotebook(notebookPath: string): Promise<void> {
    try {
      const content = readFileSync(notebookPath, 'utf8');
      const notebook: NotebookData = JSON.parse(content);
      
      const filename = notebookPath.split('/').pop() || 'unknown.ipynb';
      
      for (let i = 0; i < notebook.cells.length; i++) {
        const cell = notebook.cells[i];
        const cellContent = Array.isArray(cell.source) 
          ? cell.source.join('') 
          : cell.source;
        
        if (cellContent.trim()) {
          const snippet: RagSnippet = {
            id: `${notebookPath}:${i}`,
            content: cellContent,
            metadata: {
              filename,
              cell_index: i,
              cell_type: cell.cell_type,
              notebook_path: notebookPath
            }
          };
          
          await this.ragClient.indexSnippet(
            snippet.id,
            snippet.content,
            snippet.metadata
          );
        }
      }
      
      console.log(`Processed notebook: ${filename} (${notebook.cells.length} cells)`);
    } catch (error) {
      console.error(`Error processing notebook ${notebookPath}:`, error);
    }
  }
}

export class MockRagClient implements RagClient {
  private snippets: Map<string, { content: string; metadata: Record<string, any> }> = new Map();

  async indexSnippet(id: string, content: string, metadata: Record<string, any>): Promise<void> {
    this.snippets.set(id, { content, metadata });
    console.log(`Indexed snippet: ${id} (${content.length} chars)`);
  }

  getSnippets(): Map<string, { content: string; metadata: Record<string, any> }> {
    return this.snippets;
  }

  clear(): void {
    this.snippets.clear();
  }
}

export function createRagClient(): RagClient {
  // In a real implementation, this would connect to the actual RAG system
  // For now, return a mock client
  return new MockRagClient();
}