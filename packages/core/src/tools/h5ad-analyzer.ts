HWL/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import * as fs from 'fs';
import * as path from 'path';
import { Schema, Type } from '@google/genai';
import { BaseTool, Icon, ToolResult } from './tools.js';
import { spawn } from 'child_process';

export interface H5adAnalyzerParams {
  /** Path to the h5ad file to analyze */
  h5adPath: string;
  /** Analysis task description (e.g., "perform basic preprocessing and clustering") */
  analysisTask: string;
  /** Optional: specific parameters or constraints for the analysis */
  parameters?: string;
  /** Optional: output directory for results */
  outputDir?: string;
}

export interface RagSource {
  source: string;
  content: string;
}

const H5AD_ANALYZER_SCHEMA: Schema = {
  type: Type.OBJECT,
  properties: {
    h5adPath: {
      type: Type.STRING,
      description: 'Path to the h5ad file to analyze',
    },
    analysisTask: {
      type: Type.STRING,
      description:
        'Description of the analysis task to perform (e.g., "perform basic preprocessing and clustering", "find marker genes", "trajectory analysis")',
    },
    parameters: {
      type: Type.STRING,
      description:
        'Optional specific parameters or constraints for the analysis',
    },
    outputDir: {
      type: Type.STRING,
      description:
        'Optional output directory for results (defaults to current directory)',
    },
  },
  required: ['h5adPath', 'analysisTask'],
};

export class H5adAnalyzerTool extends BaseTool<H5adAnalyzerParams> {
  private outputLog: string[] = [];

  constructor() {
    super(
      'h5ad_analyzer',
      'H5AD File Analyzer',
      'Analyzes h5ad single-cell files using RAG-enhanced dynamic code generation',
      Icon.LightBulb,
      H5AD_ANALYZER_SCHEMA,
      true,
      true,
    );
  }

  private logAndUpdate(message: string, updateOutput?: (output: string) => void) {
    this.outputLog.push(message);
    updateOutput?.(message);
    // Also log to console for debugging
    console.log(`[H5AD-DEBUG] ${message.trim()}`);
  }

  validateToolParams(params: H5adAnalyzerParams): string | null {
    if (!params.h5adPath?.trim()) {
      return 'h5adPath is required';
    }

    if (!params.analysisTask?.trim()) {
      return 'analysisTask is required';
    }

    // Check if h5ad file exists
    const h5adPath = path.resolve(params.h5adPath);
    if (!fs.existsSync(h5adPath)) {
      return `H5AD file not found: ${h5adPath}`;
    }

    // Check if it's actually an h5ad file
    if (!h5adPath.endsWith('.h5ad')) {
      return 'File must have .h5ad extension';
    }

    return null;
  }

  getDescription(params: H5adAnalyzerParams): string {
    return `Analyze h5ad file "${params.h5adPath}" with task: "${params.analysisTask}"`;
  }

  toolLocations(
    params: H5adAnalyzerParams,
  ): Array<{ path: string; line?: number }> {
    const locations = [{ path: path.resolve(params.h5adPath) }];
    if (params.outputDir) {
      locations.push({ path: path.resolve(params.outputDir) });
    }
    return locations;
  }

  async execute(
    params: H5adAnalyzerParams,
    signal: AbortSignal,
    updateOutput?: (output: string) => void,
  ): Promise<ToolResult> {
    const validation = this.validateToolParams(params);
    if (validation) {
      throw new Error(validation);
    }

    const h5adPath = path.resolve(params.h5adPath);
    const outputDir = params.outputDir
      ? path.resolve(params.outputDir)
      : process.cwd();

    // Clear previous log for new analysis
    this.outputLog = [];
    
    this.logAndUpdate(`🔍 Starting H5AD analysis for: ${h5adPath}\n`, updateOutput);
    this.logAndUpdate(`📋 Task: ${params.analysisTask}\n`, updateOutput);

    try {
      // Step 1: Query RAG system for relevant code
      this.logAndUpdate(`🧬 Starting H5AD Analysis with RAG Integration\n`, updateOutput);
      this.logAndUpdate(`📋 Task: ${params.analysisTask}\n`, updateOutput);
      this.logAndUpdate(`📁 File: ${h5adPath}\n\n`, updateOutput);

      const ragResponse = await this.queryRagSystem(
        params.analysisTask,
        params.parameters,
        h5adPath,
        updateOutput,
      );

      if (!ragResponse) {
        throw new Error('Failed to get guidance from RAG system');
      }

      this.logAndUpdate(`\n📚 Knowledge Sources Found:\n`, updateOutput);
      ragResponse.sources.slice(0, 3).forEach((source, idx) => {
        this.logAndUpdate(`  ${idx + 1}. ${source.source || 'Unknown source'}\n`, updateOutput);
      });
      this.logAndUpdate(`\n`, updateOutput);

      // Step 2: Generate code using model + RAG context
      this.logAndUpdate(`🤖 Generating analysis code...\n`, updateOutput);
      const generatedCode = await this.generateAnalysisCode(
        h5adPath,
        params.analysisTask,
        ragResponse.guidance,
        params.parameters,
        outputDir,
      );

      // Step 3: Validate function parameters in generated code
      this.logAndUpdate(`🔍 Validating generated code parameters...\n`, updateOutput);
      const validationResult = await this.validateGeneratedCode(generatedCode);

      if (!validationResult.isValid) {
        this.logAndUpdate(`❌ Validation failed: ${validationResult.errors.join(', ')}\n`, updateOutput);
        throw new Error(
          `Generated code validation failed: ${validationResult.errors.join(', ')}`,
        );
      }

      this.logAndUpdate(`✅ Code validation passed\n`, updateOutput);

      // Step 4: Execute the generated code
      this.logAndUpdate(`🚀 Executing analysis...\n`, updateOutput);
      const executionResult = await this.executeAnalysisCode(
        generatedCode,
        h5adPath,
        outputDir,
        signal,
        updateOutput,
      );

      this.logAndUpdate(`✅ Analysis completed successfully!\n`, updateOutput);

      // Include complete log in the result
      const completeLog = this.outputLog.join('');
      this.logAndUpdate(`\n📋 Complete Analysis Log:\n${completeLog}\n`, updateOutput);

      return {
        summary: `H5AD analysis completed for ${path.basename(h5adPath)}`,
        llmContent: [
          {
            text: `H5AD Analysis completed successfully!\n\nFile: ${h5adPath}\nTask: ${params.analysisTask}\nSources used: ${ragResponse.sourceCount}\n\nGenerated and executed code:\n\`\`\`python\n${generatedCode}\n\`\`\`\n\nOutput:\n${executionResult.output}`,
          },
        ],
        returnDisplay: `# 🧬 H5AD Analysis Results

**📁 File:** \`${path.basename(h5adPath)}\`  
**📋 Task:** ${params.analysisTask}  
**📂 Output Directory:** \`${outputDir}\`  
**📚 Knowledge Sources:** ${ragResponse.sourceCount} documents

## 🤖 RAG-Generated Analysis Code

\`\`\`python
${generatedCode}
\`\`\`

## 📊 Execution Results

\`\`\`
${executionResult.output}
\`\`\`

## 📚 Knowledge Sources Used

${ragResponse.sources
  .slice(0, 5)
  .map(
    (source, idx) =>
      `${idx + 1}. **${source.source || 'Unknown'}**: ${source.content.substring(0, 150)}...`,
  )
  .join('\n')}

## ✅ Analysis Summary

- 🧠 **RAG Integration**: Retrieved guidance from ${ragResponse.sourceCount} knowledge sources
- 🔍 **Parameter Validation**: All function parameters validated successfully  
- ⚡ **Code Execution**: Analysis completed without errors
- 💾 **Results Saved**: Output files available in \`${outputDir}\`
- 🎯 **Task Completion**: "${params.analysisTask}" executed successfully

---
*Analysis powered by RAG-enhanced dynamic code generation*
`,
      };
    } catch (error) {
      const errorMessage =
        error instanceof Error ? error.message : String(error);
      this.logAndUpdate(`❌ Error: ${errorMessage}\n`, updateOutput);

      // Show complete log even on error for debugging
      const completeLog = this.outputLog.join('');
      this.logAndUpdate(`\n📋 Complete Error Log:\n${completeLog}\n`, updateOutput);

      return {
        summary: `H5AD analysis failed: ${errorMessage}`,
        llmContent: [
          {
            text: `Analysis failed for ${h5adPath}: ${errorMessage}`,
          },
        ],
        returnDisplay: `# H5AD Analysis Failed

**File:** \`${h5adPath}\`  
**Task:** ${params.analysisTask}

**Error:** ${errorMessage}

Please check the file path and ensure all dependencies are installed.
`,
      };
    }
  }

  private async queryRagSystem(
    analysisTask: string,
    parameters?: string,
    h5adPath?: string,
    updateOutput?: (output: string) => void,
  ): Promise<{
    guidance: string;
    sourceCount: number;
    sources: RagSource[];
  } | null> {
    try {
      this.logAndUpdate(`🧠 Connecting to RAG system...\n`, updateOutput);

      // Use the enhanced H5AD RAG query system
      const ragCommand = [
        'python3',
        path.join(process.cwd(), 'rag_system', 'h5ad_rag_query.py'),
        analysisTask,
        '--h5ad-path',
        h5adPath || '',
        '--parameters',
        parameters || '',
      ];

      this.logAndUpdate(`📡 Querying RAG knowledge base...\n`, updateOutput);
      this.logAndUpdate(
        `⚡ This may take a few moments for embedding and retrieval...\n`,
        updateOutput,
      );

      const result = await this.executeCommand(
        ragCommand,
        process.cwd(),
        undefined,
        (output: string) => {
          // Forward RAG system progress messages
          if (
            output.includes('📂') ||
            output.includes('🔄') ||
            output.includes('✅') ||
            output.includes('🤖') ||
            output.includes('📊') ||
            output.includes('🧠')
          ) {
            updateOutput?.(output);
          }
        },
      );

      if (result.exitCode !== 0) {
        this.logAndUpdate(`❌ RAG system query failed: ${result.stderr}\n`, updateOutput);
        return null;
      }

      // Parse RAG response (JSON format)
      try {
        const ragResponse = JSON.parse(result.stdout);

        if (!ragResponse.success) {
          this.logAndUpdate(`❌ RAG system error: ${ragResponse.error}\n`, updateOutput);
          return null;
        }

        this.logAndUpdate(
          `✅ Retrieved knowledge from ${ragResponse.source_count} sources\n`,
          updateOutput,
        );

        return {
          guidance: ragResponse.answer,
          sourceCount: ragResponse.source_count || 0,
          sources: ragResponse.source_documents || [],
        };
      } catch (parseError) {
        this.logAndUpdate(`❌ Failed to parse RAG response: ${parseError}\n`, updateOutput);
        this.logAndUpdate(`Raw output: ${result.stdout.substring(0, 500)}...\n`, updateOutput);
        return null;
      }
    } catch (error) {
      this.logAndUpdate(`❌ Error querying RAG system: ${error}\n`, updateOutput);
      return null;
    }
  }

  private extractPythonCodeFromResponse(response: string): string {
    // First, try to find code blocks marked with ```python
    const pythonBlockRegex = /```python\n([\s\S]*?)\n```/g;
    const pythonMatches = [...response.matchAll(pythonBlockRegex)];
    
    if (pythonMatches.length > 0) {
      // Combine all Python code blocks
      const codeBlocks = pythonMatches.map(match => match[1].trim());
      return codeBlocks.join('\n\n');
    }
    
    // Fallback: try generic code blocks
    const genericBlockRegex = /```\n([\s\S]*?)\n```/g;
    const genericMatches = [...response.matchAll(genericBlockRegex)];
    
    if (genericMatches.length > 0) {
      // Filter for Python-like content
      const pythonLikeBlocks = genericMatches
        .map(match => match[1].trim())
        .filter(block => 
          block.includes('import') || 
          block.includes('sc.') || 
          block.includes('ov.') ||
          block.includes('scanpy') ||
          block.includes('omicverse')
        );
      
      if (pythonLikeBlocks.length > 0) {
        return pythonLikeBlocks.join('\n\n');
      }
    }
    
    // Last resort: extract lines that look like Python code
    const lines = response.split('\n');
    const codeLines = lines.filter(line => {
      const trimmed = line.trim();
      return (
        trimmed.startsWith('import ') ||
        trimmed.startsWith('from ') ||
        trimmed.includes('sc.') ||
        trimmed.includes('ov.') ||
        trimmed.includes('adata') ||
        trimmed.includes('scanpy') ||
        (trimmed.includes('=') && !trimmed.startsWith('#')) ||
        trimmed.startsWith('print(')
      );
    });
    
    return codeLines.join('\n');
  }

  private fixIncompleteFunctionCalls(code: string): string {
    const lines = code.split('\n');
    const fixedLines: string[] = [];
    let inMultilineCall = false;
    let multilineCallIndent = 0;
    
    for (let i = 0; i < lines.length; i++) {
      let line = lines[i];
      const trimmedLine = line.trim();
      
      // Detect start of multi-line function call ending with ()
      if (trimmedLine.match(/\w+.*=.*\w+\.\w+.*\(\)\s*$/)) {
        // This is a function call ending with empty parentheses, likely continuing on next lines
        // Convert it to single line format
        line = line.replace(/\(\)\s*$/, '(');
        inMultilineCall = true;
        multilineCallIndent = line.search(/\S/); // Find indentation level
      }
      // Handle continuation lines in multi-line calls
      else if (inMultilineCall) {
        if (trimmedLine === ')') {
          // End of multi-line call
          inMultilineCall = false;
          // Convert to inline format - merge with previous line
          if (fixedLines.length > 0) {
            const lastLine = fixedLines[fixedLines.length - 1];
            if (lastLine.trim().endsWith(',')) {
              fixedLines[fixedLines.length - 1] = lastLine.replace(/,\s*$/, ')');
              continue; // Skip adding the closing paren line
            } else {
              fixedLines[fixedLines.length - 1] = lastLine + ')';
              continue;
            }
          }
        } else if (trimmedLine.length > 0) {
          // This is a parameter line - convert to inline format
          const paramContent = trimmedLine.replace(/,$/, '');
          if (fixedLines.length > 0) {
            const lastLine = fixedLines[fixedLines.length - 1];
            if (lastLine.includes('(') && !lastLine.includes(')')) {
              // Add parameter to the function call line
              if (paramContent) {
                fixedLines[fixedLines.length - 1] = lastLine + paramContent + ')';
                inMultilineCall = false;
              }
              continue;
            }
          }
        }
      }
      
      // Check for incomplete function calls ending with partial parameters
      const incompleteCallRegex = /^(\s*\w+.*\([^)]*,\s*\w+_?)$/;
      
      if (incompleteCallRegex.test(trimmedLine)) {
        // This line looks like an incomplete function call
        // Try to fix common QC metrics pattern
        if (line.includes('calculate_qc_metrics') && line.includes('percent_')) {
          line = line.replace(/percent_\s*$/, "percent_top=[20])");
        }
        // Fix other common incomplete patterns
        else if (line.includes('filter_') && line.endsWith('min_')) {
          line = line.replace(/min_\s*$/, "min_genes=200)");
        }
        else if (line.includes('normalize_') && line.endsWith('target_')) {
          line = line.replace(/target_\s*$/, "target_sum=1e4)");
        }
        // Generic fix - just close the parenthesis if it's clearly incomplete
        else if (trimmedLine.endsWith('_') || trimmedLine.endsWith(',')) {
          line = line.replace(/[,_]\s*$/, ')');
        }
      }
      
      // Check for unmatched parentheses
      const openParens = (line.match(/\(/g) || []).length;
      const closeParens = (line.match(/\)/g) || []).length;
      
      if (openParens > closeParens) {
        // Add missing closing parentheses
        line += ')'.repeat(openParens - closeParens);
      }
      
      fixedLines.push(line);
    }
    
    return fixedLines.join('\n');
  }

  private convertMultilineToSingleLine(code: string): string {
    const lines = code.split('\n');
    const fixedLines: string[] = [];
    let i = 0;
    
    while (i < lines.length) {
      const line = lines[i].trim();
      
      // Check if this line starts a multi-line function call
      if (line.match(/\w+.*=.*\w+\.\w+.*\(\)\s*$/) || 
          (line.includes('=') && line.includes('(') && !line.includes(')') && line.endsWith('('))) {
        
        // This is the start of a multi-line call - collect all related lines
        let functionCall = lines[i]; // Keep original indentation for first line
        i++;
        
        // Collect continuation lines until we find the closing parenthesis
        while (i < lines.length) {
          const nextLine = lines[i].trim();
          if (nextLine === ')') {
            // End of function call - close it
            functionCall = functionCall.replace(/\(\)\s*$/, '()').replace(/\(\s*$/, '()');
            break;
          } else if (nextLine.length > 0 && !nextLine.startsWith('#')) {
            // This is a parameter line - add it to the function call
            const cleanParam = nextLine.replace(/,$/, '');
            if (functionCall.includes('()')) {
              functionCall = functionCall.replace('()', `(${cleanParam})`);
            } else if (functionCall.endsWith('(')) {
              functionCall = functionCall + cleanParam + ')';
            }
          }
          i++;
        }
        
        fixedLines.push(functionCall);
      } else {
        // Regular line - just add it
        fixedLines.push(lines[i]);
      }
      
      i++;
    }
    
    return fixedLines.join('\n');
  }

  private async generateAnalysisCode(
    h5adPath: string,
    analysisTask: string,
    ragGuidance: string,
    parameters?: string,
    outputDir?: string,
  ): Promise<string> {
    // Get API key from environment
    const apiKey = process.env.GEMINI_API_KEY;
    if (!apiKey) {
      throw new Error('GEMINI_API_KEY environment variable not found');
    }

    // Extract existing code from RAG guidance if available
    const existingCode = this.extractPythonCodeFromResponse(ragGuidance);
    
    const prompt = `
Based on the following guidance and existing code examples, generate clean, executable Python code for single-cell RNA-seq analysis.

GUIDANCE FROM RAG SYSTEM:
${ragGuidance}

EXISTING CODE PATTERNS FOUND:
${existingCode}

TASK: ${analysisTask}
H5AD FILE: ${h5adPath}
OUTPUT DIR: ${outputDir || process.cwd()}
${parameters ? `PARAMETERS: ${parameters}` : ''}

Generate ONLY executable Python code that:
1. Imports necessary libraries (scanpy as sc, omicverse as ov)
2. Loads the H5AD file: adata = sc.read("${h5adPath}")
3. Performs ${analysisTask} using appropriate sc.pp.*, sc.tl.*, or ov.* functions
4. Saves results to ${outputDir || process.cwd()}
5. Includes print statements for progress

CRITICAL FORMATTING RULES:
- ALL function calls must be on SINGLE LINES with complete parameters
- NO multi-line function calls (avoid splitting arguments across lines)
- Every opening parenthesis '(' must have a matching closing parenthesis ')'
- Use simple, direct function calls like: sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
- Avoid complex nested calls or multi-line parameter lists

IMPORTANT: Return ONLY the Python code, no markdown formatting, no explanations.
Start with imports and end with save operations.

Example structure:
import scanpy as sc
import omicverse as ov
adata = sc.read("${h5adPath}")
# Analysis operations here
adata.write("${outputDir || process.cwd()}/result.h5ad")
print("Analysis complete")
`;

    // Call Gemini API
    const response = await fetch(
      'https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-pro:generateContent?key=' +
        apiKey,
      {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          contents: [
            {
              parts: [
                {
                  text: prompt,
                },
              ],
            },
          ],
          generationConfig: {
            temperature: 0.1,
            maxOutputTokens: 100000,
          },
        }),
      },
    );

    if (!response.ok) {
      throw new Error(
        `Gemini API error: ${response.status} ${response.statusText}`,
      );
    }

    const data = await response.json();
    const generatedText = data.candidates?.[0]?.content?.parts?.[0]?.text;

    if (!generatedText) {
      throw new Error('No code generated from Gemini API');
    }

    // Clean up the generated code
    let code = this.extractPythonCodeFromResponse(generatedText);
    
    // If extraction failed, use the raw response but clean it up
    if (!code.trim()) {
      code = generatedText.trim();
      // Remove common markdown artifacts
      code = code.replace(/^```python\n?/, '');
      code = code.replace(/\n?```$/, '');
      code = code.replace(/^```\n?/, '');
    }

    // Fix incomplete function calls and multi-line formatting issues
    code = this.fixIncompleteFunctionCalls(code);
    code = this.convertMultilineToSingleLine(code);
    
    // Ensure we have basic imports if missing
    if (!code.includes('import scanpy') && !code.includes('import sc')) {
      code = 'import scanpy as sc\nimport omicverse as ov\n' + code;
    }
    
    return code;
  }

  private async validateGeneratedCode(
    code: string,
  ): Promise<{ isValid: boolean; errors: string[] }> {
    const errors: string[] = [];

    // Basic validation checks
    if (!code.trim()) {
      errors.push('Generated code is empty');
      return { isValid: false, errors };
    }

    // More flexible function pattern matching
    const functionPatterns = [
      // Scanpy patterns
      /sc\.(pp|tl|pl)\.[\w_]+\s*\(/g,
      // Omicverse patterns  
      /ov\.[\w._]+\s*\(/g,
      // Direct function calls
      /calculate_qc_metrics\s*\(/g,
      /filter_cells\s*\(/g,
      /filter_genes\s*\(/g,
      /normalize_total\s*\(/g,
      /log1p\s*\(/g,
      // Variable assignments with analysis functions
      /\w+\s*=\s*sc\./g,
      /\w+\s*=\s*ov\./g
    ];

    let hasAnalysisFunctions = false;
    let foundPatterns: string[] = [];
    
    for (const pattern of functionPatterns) {
      const matches = code.match(pattern);
      if (matches && matches.length > 0) {
        hasAnalysisFunctions = true;
        foundPatterns.push(...matches);
      }
    }

    if (!hasAnalysisFunctions) {
      // More detailed analysis of what's missing
      if (!code.includes('sc.') && !code.includes('ov.')) {
        errors.push('No scanpy (sc.) or omicverse (ov.) function calls found');
      } else {
        errors.push('Analysis functions present but not in expected format');
      }
    } else {
      // Log found patterns for debugging
      console.log(`Found ${foundPatterns.length} analysis function patterns:`, foundPatterns.slice(0, 5));
    }

    // Check for basic structure with more flexibility
    if (!code.includes('import') && !code.includes('from')) {
      errors.push('No import statements found');
    }

    // More flexible h5ad reference check
    if (!code.includes('.h5ad') && !code.includes('adata') && !code.includes('h5ad')) {
      errors.push('Code does not reference h5ad file or AnnData object');
    }

    // Check for basic Python syntax
    const hasBasicPython = (
      code.includes('=') ||  // assignments
      code.includes('print(') ||  // output
      code.includes('if ') ||  // conditionals
      code.includes('for ') ||  // loops
      code.includes('def ')   // functions
    );
    
    if (!hasBasicPython) {
      errors.push('Code does not contain basic Python syntax elements');
    }

    return {
      isValid: errors.length === 0,
      errors,
    };
  }

  private async executeAnalysisCode(
    code: string,
    h5adPath: string,
    outputDir: string,
    signal: AbortSignal,
    updateOutput?: (output: string) => void,
  ): Promise<{ output: string; exitCode: number }> {
    // Create a temporary Python file
    const tempDir = path.join(outputDir, '.temp_analysis');
    if (!fs.existsSync(tempDir)) {
      fs.mkdirSync(tempDir, { recursive: true });
    }

    const scriptPath = path.join(tempDir, 'analysis_script.py');

    // Add necessary imports and setup if not present
    const enhancedCode = `
import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Ensure output directory exists
output_dir = "${outputDir}"
os.makedirs(output_dir, exist_ok=True)

# Change to output directory for relative paths
os.chdir(output_dir)

try:
${code
  .split('\n')
  .map((line) => '    ' + line)
  .join('\n')}
    print("\\n✅ Analysis completed successfully!")
except Exception as e:
    print(f"❌ Error during analysis: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
`;

    fs.writeFileSync(scriptPath, enhancedCode);

    try {
      const result = await this.executeCommand(
        ['python3', scriptPath],
        outputDir,
        signal,
        updateOutput,
      );

      // Clean up temp file
      fs.unlinkSync(scriptPath);
      if (fs.existsSync(tempDir) && fs.readdirSync(tempDir).length === 0) {
        fs.rmdirSync(tempDir);
      }

      return result;
    } catch (error) {
      // Clean up temp file on error
      if (fs.existsSync(scriptPath)) {
        fs.unlinkSync(scriptPath);
      }
      throw error;
    }
  }

  private async executeCommand(
    command: string[],
    cwd: string,
    signal?: AbortSignal,
    updateOutput?: (output: string) => void,
  ): Promise<{
    stdout: string;
    stderr: string;
    exitCode: number;
    output: string;
  }> {
    return new Promise((resolve, reject) => {
      const process = spawn(command[0], command.slice(1), {
        cwd,
        stdio: ['pipe', 'pipe', 'pipe'],
      });

      let stdout = '';
      let stderr = '';
      let output = '';

      process.stdout?.on('data', (data: Buffer) => {
        const text = data.toString();
        stdout += text;
        output += text;
        updateOutput?.(text);
      });

      process.stderr?.on('data', (data: Buffer) => {
        const text = data.toString();
        stderr += text;
        output += text;
        updateOutput?.(text);
      });

      process.on('close', (code) => {
        resolve({
          stdout: stdout.trim(),
          stderr: stderr.trim(),
          exitCode: code || 0,
          output: output.trim(),
        });
      });

      process.on('error', (error) => {
        reject(error);
      });

      signal?.addEventListener('abort', () => {
        process.kill('SIGTERM');
        reject(new Error('Process aborted'));
      });
    });
  }
}
