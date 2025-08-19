"""Hi-C Analysis mode handler"""

from pathlib import Path
from typing import Optional

def generate_hic_analysis_message(folder_path: Optional[str] = None) -> str:
    """Generate Hi-C analysis message using Hi-C toolset"""
    
    if folder_path:
        folder_path = Path(folder_path).resolve()
        
        message = f"""
🧬 Hi-C Analysis Pipeline — Workflow-Based Architecture
Target folder: {folder_path}

You have access to the Hi-C toolset with comprehensive workflow-based commands and TodoList management.

GLOBAL RULES
- Always use the provided folder_path: "{folder_path}" in all phases.
- Idempotent behavior: NEVER create duplicate todos. Only create if the list is EMPTY.
- Do not ask the user for confirmations; proceed automatically and log warnings when needed.
- After each concrete tool completes successfully, call mark_task_done("what was completed"), then show_todos().

PHASE 0 — SPECIES DETECTION & GENOME RESOURCES
1) Use workflow commands for species detection and setup:
   - hic.init() - Initialize Hi-C project structure
   - hic.auto_detect_species("{folder_path}") - Auto-detect species and enzyme
   - hic.setup_genome_resources() - Setup genome and restriction sites
   (Tool availability will be checked automatically within each workflow)

PHASE 1 — TODO CREATION (STRICT DE-DUP)
Mandatory order:
  a) current = show_todos()
  b) Analyze folder structure and FASTQ files in "{folder_path}"
Creation rule (single condition):
  • If current is EMPTY → create ONCE the following todos:
      0. "Initialize Hi-C project structure"
      1. "Auto-detect species and restriction enzyme"
      2. "Check and install Hi-C dependencies (HiCExplorer, BWA, Cooler)"
      3. "Setup genome resources and restriction sites"
      4. "Hi-C Quality Control with FastQC"
      5. "Hi-C Adapter Trimming (minimal for Hi-C)"
      6. "Hi-C Reads Mapping with BWA (separate R1/R2)"
      7. "Hi-C BAM Processing and Filtering"
      8. "Hi-C Matrix Building at Multiple Resolutions"
      9. "Hi-C Matrix Correction (ICE normalization)"
      10. "Hi-C TAD Calling with Enhanced Methods"
      11. "Hi-C A/B Compartment Analysis"
      12. "Hi-C Chromatin Loop Detection"
      13. "Hi-C Visualization and Track Generation"
      14. "Hi-C Comprehensive QC Report"
  • Else → DO NOT create anything. Work with the existing todos.

PHASE 2 — EXECUTE WITH TODO TRACKING (LOOP)

⚠️ CRITICAL EXECUTION STRATEGY:
When you call hic.HiC_Upstream() or hic.HiC_Analysis(), they return bash command templates.
You MUST:
1. **Read and analyze** the entire returned bash command content carefully
2. **Understand the Hi-C methodology** and biology described
3. **Adapt the provided commands** to your current data situation (file paths, sample names, enzyme)
4. **Execute the adapted commands** using bash tool - NOT the original commands directly
5. **Handle errors** by adjusting commands based on Hi-C specific guidance
6. **Analyze results** - Check matrix quality, TAD boundaries, compartments

🧠 **Hi-C RESULT ANALYSIS REQUIREMENT:**
After executing any bash commands:
1. **Analyze the output** - Hi-C analysis has unique QC metrics!
2. **Check generated files** - Verify matrices, TADs, compartments were created
3. **Examine logs and errors** - Look for mapping rates, valid pairs, matrix sparsity
4. **Validate Hi-C results** - Are interaction patterns biologically reasonable?
5. **Check resolution quality** - Are different resolutions appropriate for analysis?
6. **Document findings** - Note restriction enzyme efficiency, library quality

The returned bash commands serve as TEMPLATES and GUIDANCE for Hi-C, adapt them to your specific data.

For each current task:
  1) hint = execute_current_task()   # obtain guidance for the next action
  2) Get bash command templates using appropriate Hi-C workflow:
     
     UPSTREAM WORKFLOWS (use hic.HiC_Upstream(workflow_type)):
     - "init" - Initialize Hi-C project structure
     - "setup_genome_resources" - Setup genome and restriction enzyme sites
     - "run_fastqc" - Quality control analysis for Hi-C
     - "trim_adapters" - Minimal adapter trimming
     - "align_reads" - BWA mapping (separate R1/R2 for Hi-C)
     - "process_bam" - BAM quality filtering
     - "build_matrix" - Build Hi-C contact matrix at multiple resolutions
     - "correct_matrix" - Apply ICE normalization
     - "generate_qc" - Generate Hi-C specific QC metrics
     
     DOWNSTREAM WORKFLOWS (use hic.HiC_Analysis(workflow_type)):
     - "call_tads" - Enhanced TAD calling with multiple methods
     - "find_compartments" - A/B compartment analysis with PCA
     - "call_loops" - Chromatin loop detection
     - "plot_matrix" - Visualize Hi-C matrices
     - "plot_tads" - Create TAD visualization with tracks
     - "differential_analysis" - Compare conditions
     - "integration_analysis" - Multi-omics integration
     - "generate_tracks" - Create genome browser tracks
     - "quality_control" - Comprehensive QC analysis
     - "convert_formats" - Convert between H5/Cool/HiC formats
     - "compare_matrices" - Statistical matrix comparisons

  3) **EXECUTE** the adapted bash commands and analyze Hi-C results
  4) **VERIFY** success by checking Hi-C matrices and quality metrics
  5) mark_task_done("brief, precise description of the completed Hi-C step")
  6) show_todos()
Repeat until all todos are completed.

PHASE 3 — ADAPTIVE TODO REFINEMENT
- If Hi-C tools missing → add_todo("Install missing Hi-C tools (HiCExplorer, Cooler)")
- If mapping rate low → add_todo("Investigate Hi-C library quality")
- If matrix sparsity high → add_todo("Adjust resolution or normalization")
- If TAD calling issues → add_todo("Optimize TAD calling parameters")
- If additional analysis needed → add_todo("Hi-C specific analysis task")

Hi-C EXECUTION STRATEGY (MUST FOLLOW THIS ORDER)
  1) hic.init() → Initialize Hi-C project (includes dependency check)
  2) hic.auto_detect_species("{folder_path}") → Detect species and enzyme
  3) hic.HiC_Upstream("setup_genome_resources") → Setup genome
  5) show_todos()
  6) Analyze folder and create Hi-C todos if empty
  7) Loop Phase 2 until all done; refine with Phase 3 when needed

📊 Hi-C EXECUTION EXAMPLES:

🏷️ STEP 1 - Hi-C PROJECT INITIALIZATION:
```bash
# Get Hi-C initialization commands
init_commands = hic.init(project_name="hic_analysis", enzyme="MboI")
# Analyze and adapt the commands to your Hi-C project
# Execute: mkdir -p hic_analysis/{{matrices,tads,compartments,loops}}, create config files
```

🏷️ STEP 2 - SPECIES & ENZYME DETECTION:
```bash
# Auto-detect species and restriction enzyme
species_info = hic.auto_detect_species("{folder_path}")
# Review detected species and enzyme recommendations
# Common enzymes: MboI (GATC), DpnII (GATC), HindIII (AAGCTT)
```

🏷️ STEP 3 - Hi-C DEPENDENCY CHECK:
```bash
# Get Hi-C dependency check commands  
dep_commands = hic.check_dependencies()
# Execute: which hicBuildMatrix, which cooler, which bwa, etc.
# Install missing Hi-C tools if needed
```

🏷️ STEP 4 - Hi-C MATRIX BUILDING:
```bash
# Get Hi-C matrix building template commands
matrix_commands = hic.HiC_Upstream("build_matrix")
# Adapt to your mapped BAM files and restriction enzyme
# Execute: hicBuildMatrix --samFiles R1.bam R2.bam --binSize 50000 --restrictionSequence GATC
# Check results: ls matrices/raw/*.h5, examine QC reports
# Analyze: valid pairs percentage, duplication rates
```

🏷️ STEP 5 - Hi-C MATRIX CORRECTION:
```bash
# Get matrix correction template commands
correct_commands = hic.HiC_Upstream("correct_matrix")
# Adapt with your raw matrix files
# Execute: hicCorrectMatrix correct --matrix raw_matrix.h5 --correctionMethod ICE
# Check results: hicPlotMatrix for visualization
# Analyze: convergence of normalization, matrix quality
```

🏷️ STEP 6 - TAD CALLING:
```bash
# Get TAD calling template commands
tad_commands = hic.HiC_Analysis("call_tads")
# Adapt with your corrected matrix
# Execute: hicFindTADs --matrix corrected_matrix.h5 --outPrefix tads/sample
# Check results: wc -l tads/sample_domains.bed
# Analyze: TAD boundary distribution, size distribution
```

**Hi-C CRITICAL SUCCESS PATTERNS:**
```bash
# Good - Check Hi-C file formats before proceeding
ls *.fastq.gz  # Verify paired-end Hi-C reads exist
file sample_R1.fastq.gz  # Confirm gzip format

# Good - Monitor Hi-C mapping rates
bwa mem genome.fa R1.fastq.gz | samtools flagstat  # Check mapping rate
# Hi-C mapping rates are typically 70-90%

# Good - Validate Hi-C matrix quality
hicInfo --matrix sample.h5  # Check matrix statistics
# Look for: total interactions, sparsity, resolution coverage

# Good - Verify Hi-C specific outputs
ls matrices/corrected/*.h5  # Corrected matrices
ls tads/*.bed  # TAD boundary files
ls compartments/*.bedgraph  # Compartment tracks
```

**Hi-C QUALITY INDICATORS:**
- Valid pairs: 25-40% of total reads (Hi-C specific)
- Mapping rate: >70% for each mate
- Matrix sparsity: <95% for good coverage
- TAD boundaries: 2000-5000 per chromosome
- Compartment signal: Clear A/B pattern in PC1

BEGIN NOW:
- Execute PHASE 0 → PHASE 1 → PHASE 2 loop for Hi-C analysis.
- **ACTUALLY EXECUTE** the bash commands after adapting them for Hi-C
- **ANALYZE Hi-C RESULTS** after each step (valid pairs, matrix quality, etc.)
- **VERIFY SUCCESS** by checking Hi-C specific output files
- Report any Hi-C library quality issues or analysis problems encountered
"""
        
    else:
        message = """
I need help with Hi-C (chromosome conformation capture) analysis using your specialized workflow-based toolsets.

You have access to comprehensive Hi-C and TODO management tools:

📋 TODO MANAGEMENT (use these for ALL tasks):
- add_todo() - Add tasks and auto-break them down
- show_todos() - Display current progress  
- execute_current_task() - Get smart guidance
- mark_task_done() - Mark tasks complete and progress

🧬 COMPLETE Hi-C WORKFLOW TOOLSET:

UPSTREAM WORKFLOWS (use hic.HiC_Upstream(workflow_type)):
These return bash command templates for execution:

PROJECT SETUP:
- "init" - Initialize Hi-C project structure
- "check_dependencies" - Check Hi-C tool availability (HiCExplorer, Cooler, BWA)
- "setup_genome_resources" - Setup genome and restriction enzyme sites

QUALITY CONTROL & PREPROCESSING:
- "run_fastqc" - Quality control analysis for Hi-C
- "trim_adapters" - Minimal adapter trimming (often not needed for Hi-C)

ALIGNMENT & MATRIX GENERATION:
- "align_reads" - BWA alignment (separate R1/R2 mapping for Hi-C)
- "process_bam" - BAM filtering and processing
- "build_matrix" - Build Hi-C contact matrix at multiple resolutions
- "correct_matrix" - Apply ICE normalization to matrices
- "generate_qc" - Generate Hi-C specific quality metrics

DOWNSTREAM WORKFLOWS (use hic.HiC_Analysis(workflow_type)):
These return bash command templates for execution:

STRUCTURAL ANALYSIS:
- "call_tads" - Enhanced TAD calling with multiple methods
- "find_compartments" - A/B compartment analysis using PCA
- "call_loops" - Chromatin loop detection

VISUALIZATION & INTEGRATION:
- "plot_matrix" - Visualize Hi-C contact matrices
- "plot_tads" - Create TAD visualization with tracks
- "differential_analysis" - Compare Hi-C between conditions
- "integration_analysis" - Multi-omics integration with ChIP/RNA-seq
- "generate_tracks" - Create genome browser tracks
- "quality_control" - Comprehensive QC analysis
- "convert_formats" - Convert between H5/Cool/HiC formats
- "compare_matrices" - Statistical matrix comparisons

WORKFLOW USAGE:
All workflows return executable bash command templates. You must:
1. **Call the workflow function** to get Hi-C command templates
2. **Adapt the commands** to your specific file paths, enzyme, and parameters  
3. **Execute the adapted commands** using the bash tool
4. **Analyze the Hi-C results** and verify biological relevance

EXECUTION EXAMPLES:
```bash
# Step 1: Get Hi-C command template
matrix_commands = hic.HiC_Upstream("build_matrix")

# Step 2: Adapt and execute for Hi-C
# Template: hicBuildMatrix --samFiles R1.bam R2.bam --binSize 50000 --restrictionSequence GATC
# Adapted: hicBuildMatrix --samFiles sample_R1.bam sample_R2.bam --binSize 50000 --restrictionSequence GATC --outFileName matrix.h5

# Step 3: Verify Hi-C results  
hicInfo --matrix matrix.h5  # Check Hi-C matrix statistics
ls qc/hicqc/*/  # Check QC reports for valid pairs, duplicates
```

⚠️ **CRITICAL**: Don't just call the workflow functions - you must EXECUTE the returned bash commands for Hi-C analysis!

Hi-C SPECIFIC CONSIDERATIONS:
- Restriction enzyme: Most common are MboI/DpnII (GATC), HindIII (AAGCTT)
- Valid pairs: 25-40% is typical for Hi-C (lower than other seq methods)
- Resolution: Start with 50kb-100kb for TADs, 10kb for loops
- Matrix correction: ICE normalization is essential for Hi-C
- Quality metrics: Focus on valid pairs, duplication rate, distance decay

Please start by adding a todo for your Hi-C analysis task, then use the workflow commands and EXECUTE them!"""
    
    return message