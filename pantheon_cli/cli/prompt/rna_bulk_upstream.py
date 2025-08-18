"""Simplified RNA-seq mode handler"""

from pathlib import Path
from typing import Optional

def generate_rna_analysis_message(folder_path: Optional[str] = None) -> str:
    """Generate RNA-seq analysis message using RNA toolset"""
    
    if folder_path:
        folder_path = Path(folder_path).resolve()
        
        message = f"""
🧬 RNA-seq Analysis Pipeline — Workflow-Based Architecture
Target folder: {folder_path}

You have access to the RNA-seq toolset with workflow-based commands and TodoList management.

GLOBAL RULES
- Always use the provided folder_path: "{folder_path}" in all phases.
- Idempotent behavior: NEVER create duplicate todos. Only create if the list is EMPTY.
- Do not ask the user for confirmations; proceed automatically and log warnings when needed.
- After each concrete tool completes successfully, call mark_task_done("what was completed"), then show_todos().

PHASE 0 — SPECIES DETECTION & GENOME RESOURCES
1) Use workflow commands for species detection and setup:
   - rna.RNA_Upstream("init") - Initialize project structure
   - rna.RNA_Upstream("check_dependencies") - Check tool availability
   - rna.RNA_Upstream("setup_genome_resources") - Setup genome resources

PHASE 1 — TODO CREATION (STRICT DE-DUP)
Mandatory order:
  a) current = show_todos()
  b) Analyze folder structure and FASTQ files in "{folder_path}"
Creation rule (single condition):
  • If current is EMPTY → create ONCE the following todos:
      0. "Initialize RNA-seq project structure"
      1. "Check and install RNA-seq dependencies"
      2. "Setup genome resources and references"
      3. "RNA-seq Quality Control with FastQC"
      4. "RNA-seq Adapter Trimming"
      5. "RNA-seq Genome Alignment with STAR"
      6. "RNA-seq Expression Quantification"
      7. "RNA-seq BAM Processing and QC"
      8. "RNA-seq Coverage Track Generation"
      9. "RNA-seq QC Report Generation"
  • Else → DO NOT create anything. Work with the existing todos.

PHASE 2 — EXECUTE WITH TODO TRACKING (LOOP)

⚠️ CRITICAL EXECUTION STRATEGY:
When you call rna.RNA_Upstream() or rna.RNA_Analysis(), they return bash command templates.
You MUST:
1. **Read and analyze** the entire returned bash command content carefully
2. **Understand the logic** and methodology described
3. **Adapt the provided commands** to your current data situation (file paths, sample names, etc.)
4. **Execute the adapted commands** using bash tool - NOT the original commands directly
5. **Handle errors** by adjusting commands based on the guidance provided
6. **Analyze results** - Check output files, logs, and success/failure status

🧠 **RESULT ANALYSIS REQUIREMENT:**
After executing any bash commands:
1. **Analyze the output** - Don't just run and move on!
2. **Check generated files** - Verify expected output files were created
3. **Examine logs and errors** - Look for warnings, failures, or quality issues
4. **Validate results** - Are the results biologically reasonable?
5. **Make decisions** - Should parameters be adjusted based on what you observed?
6. **Document findings** - Note any issues or important observations

The returned bash commands serve as TEMPLATES and GUIDANCE, adapt them to your specific data.

For each current task:
  1) hint = execute_current_task()   # obtain guidance for the next action
  2) Get bash command templates using appropriate RNA workflow:
     
     UPSTREAM WORKFLOWS (use rna.RNA_Upstream(workflow_type)):
     - "init" - Initialize project structure
     - "check_dependencies" - Check tool dependencies  
     - "setup_genome_resources" - Setup genome resources
     - "run_fastqc" - Quality control analysis
     - "trim_adapters" - Adapter trimming
     - "align_star" - STAR alignment (recommended for RNA-seq)
     - "align_hisat2" - HISAT2 alignment (alternative)
     - "quantify_featurecounts" - featureCounts quantification (alignment-based)
     - "process_bam_smart" - Smart BAM processing pipeline
     - "rna_qc" - RNA-seq specific quality control
     
     DOWNSTREAM WORKFLOWS (use rna.RNA_Analysis(workflow_type)):
     - "differential_expression" - Differential expression analysis
     - "pathway_analysis" - Gene set enrichment analysis
     - "visualization" - Generate plots and visualizations

  3) **EXECUTE** the adapted bash commands and analyze results
  4) **VERIFY** success by checking output files and logs
  5) mark_task_done("brief, precise description of the completed step")
  6) show_todos()
Repeat until all todos are completed.

PHASE 3 — ADAPTIVE TODO REFINEMENT
- If dependencies missing → add_todo("Install missing RNA-seq tools")
- If quality issues found → add_todo("Address data quality issues")
- If additional analysis needed → add_todo("Additional analysis task")

EXECUTION STRATEGY (MUST FOLLOW THIS ORDER)
  1) rna.RNA_Upstream("init") → Initialize project
  2) rna.RNA_Upstream("check_dependencies") → Check tools
  3) rna.RNA_Upstream("setup_genome_resources") → Setup genome
  4) show_todos()
  5) Analyze folder and create todos if empty
  6) Loop Phase 2 until all done; refine with Phase 3 when needed

📊 EXECUTION EXAMPLES:

🏷️ STEP 1 - PROJECT INITIALIZATION:
```bash
# Get initialization commands
init_commands = rna.RNA_Upstream("init")
# Analyze and adapt the commands to your project
# Execute: mkdir -p project_structure, create config files, etc.
```

🏷️ STEP 2 - DEPENDENCY CHECK:
```bash
# Get dependency check commands  
dep_commands = rna.RNA_Upstream("check_dependencies")
# Execute: which fastqc, which STAR, which featureCounts, etc.
# Install missing tools if needed
```

🏷️ STEP 3 - FASTQC EXECUTION:
```bash
# Get FastQC template commands
fastqc_commands = rna.RNA_Upstream("run_fastqc")
# Adapt to your actual FASTQ files
# Execute: fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc/fastqc/
# Check results: ls qc/fastqc/*.html
# Analyze: Look for adapter contamination, quality issues
```

🏷️ STEP 4 - ALIGNMENT EXECUTION:
```bash
# Get STAR alignment template commands
align_commands = rna.RNA_Upstream("align_star") 
# Adapt with your actual file paths and genome index
# Execute: STAR --genomeDir genome/index/star/ --readFilesIn R1.fq.gz R2.fq.gz ...
# Check results: samtools flagstat sample_Aligned.sortedByCoord.out.bam
# Analyze: mapping rate, splice junction detection
```

🏷️ STEP 5 - QUANTIFICATION:
```bash
# Get featureCounts quantification template commands
quant_commands = rna.RNA_Upstream("quantify_featurecounts")
# Adapt with your aligned BAM files
# Execute: featureCounts -a annotation.gtf -o counts.txt sample.bam
# Check results: head counts.txt
# Analyze: number of quantified genes, count statistics
```

**CRITICAL SUCCESS PATTERNS:**
```bash
# Good - Check files before proceeding
ls *.fastq.gz  # Verify input files exist
fastqc *.fastq.gz -o qc/fastqc/
ls qc/fastqc/*.html  # Verify outputs created

# Good - Capture and analyze results
samtools flagstat sample_Aligned.sortedByCoord.out.bam > alignment_stats.txt
cat alignment_stats.txt  # Review mapping statistics

# Good - Error handling
if [ ! -f "sample_Aligned.sortedByCoord.out.bam" ]; then
    echo "ERROR: STAR alignment failed"
    exit 1
fi
```

BEGIN NOW:
- Execute PHASE 0 → PHASE 1 → PHASE 2 loop.
- **ACTUALLY EXECUTE** the bash commands after adapting them
- **ANALYZE RESULTS** after each step
- **VERIFY SUCCESS** by checking output files
- Report any errors or quality issues encountered
"""
        
    else:
        message = """
I need help with RNA-seq analysis using your specialized workflow-based toolsets.

You have access to comprehensive RNA-seq and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE RNA-seq WORKFLOW TOOLSET:

UPSTREAM WORKFLOWS (use rna.RNA_Upstream(workflow_type)):
These return bash command templates for execution:

PROJECT SETUP:
- "init" - Initialize project structure
- "check_dependencies" - Check tool availability
- "setup_genome_resources" - Setup genome resources

QUALITY CONTROL & PREPROCESSING:
- "run_fastqc" - Quality control analysis
- "trim_adapters" - Adapter trimming with Trim Galore

ALIGNMENT & QUANTIFICATION:
- "align_star" - STAR alignment (recommended for RNA-seq)
- "align_hisat2" - HISAT2 alignment (alternative)
- "quantify_featurecounts" - featureCounts quantification (alignment-based)
- "process_bam_smart" - Complete BAM processing pipeline
- "rna_qc" - RNA-seq specific quality control

DOWNSTREAM WORKFLOWS (use rna.RNA_Analysis(workflow_type)):
These return bash command templates for execution:

DIFFERENTIAL EXPRESSION:
- "differential_expression" - DESeq2 differential expression analysis
- "pathway_analysis" - Gene set enrichment analysis with clusterProfiler
- "visualization" - Generate PCA plots, volcano plots, heatmaps

WORKFLOW USAGE:
All workflows return executable bash command templates. You must:
1. **Call the workflow function** to get command templates
2. **Adapt the commands** to your specific file paths and parameters  
3. **Execute the adapted commands** using the bash tool
4. **Analyze the results** and verify success

EXECUTION EXAMPLES:
```bash
# Step 1: Get command template
fastqc_commands = rna.RNA_Upstream("run_fastqc")

# Step 2: Adapt and execute
# Template: fastqc *.fastq.gz -o qc/fastqc/
# Adapted: fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o qc/fastqc/

# Step 3: Verify results  
ls qc/fastqc/*.html  # Check HTML reports were generated
```

⚠️ **CRITICAL**: Don't just call the workflow functions - you must EXECUTE the returned bash commands!

Please start by adding a todo for your RNA-seq analysis task, then use the workflow commands and EXECUTE them!"""
    
    return message