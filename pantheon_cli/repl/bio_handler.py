"""Bio Commands Handler for REPL"""

from rich.console import Console


class BioCommandHandler:
    """Handler for /bio commands in REPL"""
    
    def __init__(self, console: Console):
        self.console = console
    
    async def handle_bio_command(self, command: str) -> str:
        """
        Handle /bio commands for bioinformatics analysis
        
        Returns:
            str: Message to send to agent, or None if no message needed
        """
        # Parse command parts
        parts = command.split()
        
        if len(parts) == 1:
            # Just /bio - show help
            self._show_bio_help()
            return None
        
        # Route bio commands to the bio toolset
        if len(parts) >= 2:
            if parts[1] in ['list', 'info', 'help']:
                return self._handle_bio_manager_command(parts)
            else:
                return self._handle_tool_specific_command(parts)
        
        return None
    
    def _show_bio_help(self):
        """Show bio commands help"""
        self.console.print("\n[bold cyan]🧬 Bio Analysis Tools[/bold cyan]")
        self.console.print("[dim]/bio list[/dim] - List all available bio analysis tools")
        self.console.print("[dim]/bio info <tool>[/dim] - Get information about a specific tool")
        self.console.print("[dim]/bio help [tool][/dim] - Get help for bio tools")
        self.console.print("[dim]/bio <tool> <command>[/dim] - Run tool-specific commands")
        self.console.print("\n[dim]Examples:[/dim]")
        self.console.print("[dim]  /bio list                      # Show all available tools[/dim]")
        self.console.print("[dim]  /bio atac init                 # Initialize ATAC-seq project[/dim]")
        self.console.print("[dim]  /bio atac upstream ./data      # Run upstream ATAC analysis[/dim]")
        self.console.print("[dim]  /bio scatac init               # Initialize scATAC-seq project[/dim]")
        self.console.print("[dim]  /bio scatac upstream ./data    # Run cellranger-atac analysis[/dim]")
        self.console.print("[dim]  /bio scrna init                # Initialize scRNA-seq project[/dim]")
        self.console.print("[dim]  /bio scrna analysis ./data    # Load and analyze scRNA-seq data[/dim]")
        self.console.print("[dim]  /bio rna init                 # Initialize RNA-seq project[/dim]")
        self.console.print("[dim]  /bio rna upstream ./data      # Run RNA-seq upstream analysis[/dim]")
        self.console.print("[dim]  /bio spatial init              # Initialize spatial project[/dim]")
        self.console.print("[dim]  /bio spatial run_spatial_workflow <workflow_type> # Run spatial workflow[/dim]")
        self.console.print("[dim]  /bio dock init                # Initialize molecular docking project[/dim]")
        self.console.print("[dim]  /bio dock run_dock            # Interactive molecular docking workflow[/dim]")
        self.console.print("[dim]  /bio dock run ./data          # Run batch molecular docking on folder[/dim]")
        self.console.print("[dim]  /bio GeneAgent TP53,BRCA1,EGFR # Gene set analysis with AI[/dim]")
        self.console.print("")
    
    def _handle_bio_manager_command(self, parts) -> str:
        """Handle direct bio manager commands (list, info, help)"""
        method_name = parts[1]
        
        if len(parts) > 2 and parts[1] in ['info', 'help']:
            # Include tool name as parameter
            tool_name = parts[2]
            return f"bio {method_name} {tool_name}"
        else:
            return f"bio {method_name}"
    
    def _handle_tool_specific_command(self, parts) -> str:
        """Handle tool-specific commands"""
        tool_name = parts[1]
        
        # Handle ATAC commands with special logic
        if tool_name == "atac":
            return self._handle_atac_command(parts)
        
        # Handle scATAC commands with special logic
        if tool_name == "scatac":
            return self._handle_scatac_command(parts)
        
        # Handle scRNA commands with special logic
        if tool_name == "scrna":
            return self._handle_scrna_command(parts)
        
        
        # Handle GeneAgent commands with special logic
        if tool_name == "GeneAgent":
            return self._handle_gene_agent_command(parts)
        
        # Handle RNA commands with special logic
        if tool_name == "rna":
            return self._handle_rna_command(parts)
        
        # Handle dock commands with special logic
        if tool_name == "dock":
            return self._handle_dock_command(parts)
        
        # Handle hic commands with special logic
        if tool_name == "hic":
            return self._handle_hic_command(parts)
        
        if tool_name == "spatial":
            return self._handle_spatial_command(parts)
        
        # Generic handler for other tools
        if len(parts) > 2:
            method_name = parts[2]
            params = " ".join(parts[3:])  # Additional parameters
            if params:
                return f"bio_{tool_name}_{method_name} {params}"
            else:
                return f"bio_{tool_name}_{method_name}"
        else:
            # Just tool name, show tool help
            return f"bio info {tool_name}"
        
    def _handle_spatial_command(self, parts) -> str:
        """Handle spatial-specific commands with special logic"""
        if len(parts) == 2:
            # Just /bio spatial - show spatial help
            self.console.print("\n[bold]🧬 Spatial Analysis Helper[/bold]")
            self.console.print("[dim]/bio spatial init[/dim] - Initialize spatial analysis project")
            self.console.print("[dim]/bio spatial run_spatial_workflow <workflow_type>[/dim] - Run spatial workflow")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio spatial init                     # Initialize spatial analysis project[/dim]")
            self.console.print("[dim]  /bio spatial run_spatial_workflow bin2cell # Run bin2cell workflow[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter spatial mode - simple mode activation
            self.console.print("\n[bold cyan]🧬 Initializing spatial analysis project[/bold cyan]")
            
            
            # Clear all existing todos when entering spatial mode
            clear_message = """
spatial INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any spatial.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  spatial init ready • todos={len(todos)}
            """

            self.console.print("[dim]Clearing existing todos and preparing spatial environment...[/dim]")
            self.console.print("[dim]Ready for spatial analysis assistance...[/dim]")
            self.console.print("[dim]Spatial mode activated. You can now use spatial tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio spatial init - Enter spatial mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio spatial run_spatial_workflow <workflow_type> - Run spatial workflow[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "run_spatial_workflow":
            # Handle spatial workflow
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a workflow type[/red]")
                self.console.print("[dim]Usage: /bio spatial run_spatial_workflow <workflow_type>[/dim]")
                self.console.print("[dim]Example: /bio spatial run_spatial_workflow bin2cell[/dim]")
                return None
            
            workflow_type = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Starting spatial workflow: {workflow_type}[/bold cyan]")
            self.console.print("[dim]Preparing spatial analysis pipeline...[/dim]\n")
            
            
            # Generate the workflow message with workflow type
            from ..cli.prompt.spatial_bin2cell import generate_spatial_workflow_message
            spatial_message = generate_spatial_workflow_message(workflow_type=workflow_type)
            
            self.console.print("[dim]Sending spatial workflow request...[/dim]\n")
            
            return spatial_message
            
    
    def _handle_atac_command(self, parts) -> str:
        """Handle ATAC-specific commands with special logic like the original _handle_atac_command"""
        
        if len(parts) == 2:
            # Just /bio atac - show ATAC help
            self.console.print("\n[bold]🧬 ATAC-seq Analysis Helper[/bold]")
            self.console.print("[dim]/bio atac init[/dim] - Enter ATAC-seq analysis mode")
            self.console.print("[dim]/bio atac upstream <folder>[/dim] - Run upstream ATAC-seq analysis on folder")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio atac init                     # Enter ATAC mode[/dim]")
            self.console.print("[dim]  /bio atac upstream ./fastq_data   # Analyze FASTQ data[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter ATAC mode - simple mode activation without automation
            self.console.print("\n[bold cyan]🧬 Entering ATAC-seq Analysis Mode[/bold cyan]")
            
            # Clear all existing todos when entering ATAC mode
            clear_message = """
ATAC INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any atac.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  ATAC init ready • todos={len(todos)}
"""
            
            self.console.print("[dim]Clearing existing todos and preparing ATAC environment...[/dim]")
            self.console.print("[dim]Ready for ATAC-seq analysis assistance...[/dim]")
            self.console.print("[dim]ATAC-seq mode activated. You can now use ATAC tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio atac init - Enter ATAC mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio atac upstream <folder> - Run upstream analysis on specific folder[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "upstream":
            # Handle upstream analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a folder path[/red]")
                self.console.print("[dim]Usage: /bio atac upstream <folder_path>[/dim]")
                self.console.print("[dim]Example: /bio atac upstream ./fastq_data[/dim]")
                return None
                
            try:
                from ..cli.prompt.atac_bulk_upstream import generate_atac_analysis_message
                
                folder_path = parts[3]
                self.console.print(f"\n[bold cyan]🧬 Starting ATAC-seq Analysis[/bold cyan]")
                self.console.print(f"[dim]Target folder: {folder_path}[/dim]")
                self.console.print("[dim]Preparing analysis pipeline...[/dim]\n")
                
                # Generate the analysis message with folder
                atac_message = generate_atac_analysis_message(folder_path=folder_path)
                
                self.console.print("[dim]Sending ATAC-seq analysis request...[/dim]\n")
                
                return atac_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: ATAC module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing analysis: {str(e)}[/red]")
                return None
        
        else:
            # Handle other ATAC commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_atac_{command} {params}"
            else:
                return f"bio_atac_{command}"
    
    def _handle_scatac_command(self, parts) -> str:
        """Handle scATAC-specific commands with special logic"""
        
        if len(parts) == 2:
            # Just /bio scatac - show scATAC help
            self.console.print("\n[bold]🧬 Single-cell ATAC-seq Analysis Helper[/bold]")
            self.console.print("[dim]/bio scatac init[/dim] - Initialize scATAC-seq analysis project")
            self.console.print("[dim]/bio scatac install[/dim] - Download and install cellranger-atac")
            self.console.print("[dim]/bio scatac upstream <folder>[/dim] - Run cellranger-atac upstream analysis")
            self.console.print("[dim]/bio scatac count <sample>[/dim] - Run cellranger-atac count for single sample")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio scatac init                      # Initialize scATAC project[/dim]")
            self.console.print("[dim]  /bio scatac install                   # Download cellranger-atac v2.2.0[/dim]")
            self.console.print("[dim]  /bio scatac upstream ./fastq_data    # Analyze 10X Chromium data[/dim]")
            self.console.print("[dim]  /bio scatac count sample1             # Process single sample[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter scATAC mode - simple mode activation
            self.console.print("\n[bold cyan]🧬 Initializing scATAC-seq Project[/bold cyan]")
            
            # Clear all existing todos when entering scATAC mode
            clear_message = """
scATAC INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any scatac.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  scATAC init ready • todos={len(todos)}
"""
            
            self.console.print("[dim]Clearing existing todos and preparing scATAC environment...[/dim]")
            self.console.print("[dim]Ready for single-cell ATAC-seq analysis assistance...[/dim]")
            self.console.print("[dim]scATAC-seq mode activated. You can now use scATAC tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio scatac init - Enter scATAC mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio scatac upstream <folder> - Run cellranger-atac analysis on specific folder[/dim]")
            self.console.print("[dim]  - /bio scatac install - Download and install cellranger-atac[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "install":
            # Handle cellranger-atac installation
            self.console.print("\n[bold cyan]🔧 Installing cellranger-atac[/bold cyan]")
            self.console.print("[dim]Downloading and setting up cellranger-atac v2.2.0...[/dim]")
            self.console.print("[dim]This will download ~500MB and may take several minutes.[/dim]")
            self.console.print()
            return "scatac_install_cellranger_atac"
        
        elif command == "upstream":
            # Handle upstream analysis with cellranger-atac
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a folder path[/red]")
                self.console.print("[dim]Usage: /bio scatac upstream <folder_path>[/dim]")
                self.console.print("[dim]Example: /bio scatac upstream ./10x_data[/dim]")
                return None
            
            try:
                from ..cli.prompt.atac_sc_upstream import generate_scatac_analysis_message
                
                folder_path = parts[3]
                self.console.print(f"\n[bold cyan]🧬 Starting scATAC-seq Analysis[/bold cyan]")
                self.console.print(f"[dim]Target folder: {folder_path}[/dim]")
                self.console.print("[dim]Will scan for 10X Chromium ATAC data and run cellranger-atac pipeline...[/dim]")
                self.console.print("[dim]Preparing cellranger-atac analysis pipeline...[/dim]\n")
                
                # Generate the analysis message with folder
                scatac_message = generate_scatac_analysis_message(folder_path=folder_path)
                
                self.console.print("[dim]Sending scATAC-seq analysis request...[/dim]\n")
                
                return scatac_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: scATAC module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing scATAC analysis: {str(e)}[/red]")
                return None
        
        elif command == "count":
            # Handle cellranger-atac count for single sample
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify sample information[/red]")
                self.console.print("[dim]Usage: /bio scatac count <sample_id>[/dim]")
                self.console.print("[dim]Example: /bio scatac count Sample1[/dim]")
                return None
            
            sample_id = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Running cellranger-atac count[/bold cyan]")
            self.console.print(f"[dim]Sample ID: {sample_id}[/dim]")
            self.console.print("[dim]Processing single-cell ATAC-seq data...[/dim]")
            
            return f"scatac_count {sample_id}"
        
        else:
            # Handle other scATAC commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_scatac_{command} {params}"
            else:
                return f"bio_scatac_{command}"
    
    def _handle_scrna_command(self, parts) -> str:
        """Handle scRNA-specific commands with special logic"""
        
        if len(parts) == 2:
            # Just /bio scrna - show scRNA help
            self.console.print("\n[bold]🧬 Single-cell RNA-seq Analysis Helper[/bold]")
            self.console.print("[dim]/bio scrna init[/dim] - Initialize scRNA-seq analysis project")
            self.console.print("[dim]/bio scrna analysis <folder/file>[/dim] - Load and inspect scRNA-seq data")
            self.console.print("[dim]/bio scrna qc <folder/file>[/dim] - Run quality control analysis")
            self.console.print("[dim]/bio scrna preprocess <folder/file>[/dim] - Run preprocessing and normalization")
            self.console.print("[dim]/bio scrna annotate <folder/file>[/dim] - Perform cell type annotation")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio scrna init                        # Initialize scRNA project[/dim]")
            self.console.print("[dim]  /bio scrna load_data ./data.h5ad       # Load and inspect H5AD file[/dim]")
            self.console.print("[dim]  /bio scrna qc ./data.h5ad             # Run quality control[/dim]")
            self.console.print("[dim]  /bio scrna preprocess ./data.h5ad     # Normalize and preprocess[/dim]")
            self.console.print("[dim]  /bio scrna annotate ./data.h5ad       # Annotate cell types[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter scRNA mode - simple mode activation
            self.console.print("\n[bold cyan]🧬 Initializing scRNA-seq Project[/bold cyan]")
            
            # Clear all existing todos when entering scRNA mode
            clear_message = """
scRNA INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any scrna.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  scRNA init ready • todos={len(todos)}
"""
            
            self.console.print("[dim]Clearing existing todos and preparing scRNA environment...[/dim]")
            self.console.print("[dim]Ready for single-cell RNA-seq analysis assistance...[/dim]")
            self.console.print("[dim]scRNA-seq mode activated. You can now use scRNA tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio scrna init - Enter scRNA mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio scrna analysis <file/folder> - Load and analyze scRNA data[/dim]")
            self.console.print("[dim]  - /bio scrna qc <file/folder> - Run quality control analysis[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "analysis":
            # Handle data loading and inspection
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a data file path[/red]")
                self.console.print("[dim]Usage: /bio scrna analysis <file_path>[/dim]")
                self.console.print("[dim]Example: /bio scrna analysis ./data.h5ad[/dim]")
                return None
            
            try:
                from ..cli.prompt.scrna_anno_analysis import generate_scrna_analysis_message
                
                file_path = parts[3]
                self.console.print(f"\n[bold cyan]🧬 Starting scRNA-seq Data Analysis[/bold cyan]")
                self.console.print(f"[dim]Target file: {file_path}[/dim]")
                self.console.print("[dim]Will load and inspect scRNA-seq data with omicverse integration...[/dim]")
                self.console.print("[dim]Preparing comprehensive analysis pipeline...[/dim]\n")
                
                # Generate the analysis message with file path
                scrna_message = generate_scrna_analysis_message(folder_path=file_path)
                
                self.console.print("[dim]Sending scRNA-seq analysis request...[/dim]\n")
                
                return scrna_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: scRNA module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing scRNA analysis: {str(e)}[/red]")
                return None
        elif command == "subtype":
            # Handle subtype analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a data file path[/red]")
                self.console.print("[dim]Usage: /bio scrna subtype <file_path>[/dim]")
                self.console.print("[dim]Example: /bio scrna subtype ./data.h5ad[/dim]")
                return None
            try:
                from ..cli.prompt.scrna_subtype_analysis import generate_scrna_subtype_analysis_message
                
                file_path = parts[3]
                self.console.print(f"\n[bold cyan]🧬 Starting scRNA-seq Subtype Analysis[/bold cyan]")
                self.console.print(f"[dim]Target file: {file_path}[/dim]")
                self.console.print("[dim]Will load and inspect scRNA-seq data with omicverse integration...[/dim]")
                self.console.print("[dim]Preparing comprehensive analysis pipeline...[/dim]\n")
                
                # Generate the analysis message with file path
                scrna_message = generate_scrna_subtype_analysis_message(folder_path=file_path)

                self.console.print("[dim]Sending scRNA-seq subtype analysis request...[/dim]\n")
                
                return scrna_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: scRNA module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing scRNA subtype analysis: {str(e)}[/red]")
                return None
                
        elif command in ["qc", "quality_control"]:
            # Handle quality control
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a data file path[/red]")
                self.console.print("[dim]Usage: /bio scrna qc <file_path>[/dim]")
                self.console.print("[dim]Example: /bio scrna qc ./data.h5ad[/dim]")
                return None
            
            file_path = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Running scRNA-seq Quality Control[/bold cyan]")
            self.console.print(f"[dim]Data file: {file_path}[/dim]")
            self.console.print("[dim]Performing quality control with omicverse...[/dim]")
            
            return f"scrna_run_quality_control {file_path}"
        
        elif command in ["preprocess", "preprocessing"]:
            # Handle preprocessing
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a data file path[/red]")
                self.console.print("[dim]Usage: /bio scrna preprocess <file_path>[/dim]")
                self.console.print("[dim]Example: /bio scrna preprocess ./data.h5ad[/dim]")
                return None
            
            file_path = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Running scRNA-seq Preprocessing[/bold cyan]")
            self.console.print(f"[dim]Data file: {file_path}[/dim]")
            self.console.print("[dim]Normalizing and preprocessing with omicverse...[/dim]")
            
            return f"scrna_run_preprocessing {file_path}"
        
        elif command in ["annotate", "annotation", "cell_type"]:
            # Handle cell type annotation
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a data file path[/red]")
                self.console.print("[dim]Usage: /bio scrna annotate <file_path>[/dim]")
                self.console.print("[dim]Example: /bio scrna annotate ./data.h5ad[/dim]")
                return None
            
            file_path = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Running Cell Type Annotation[/bold cyan]")
            self.console.print(f"[dim]Data file: {file_path}[/dim]")
            self.console.print("[dim]Performing cell type annotation with omicverse + CellOntologyMapper...[/dim]")
            
            return f"scrna_run_cell_type_annotation {file_path}"
        
        else:
            # Handle other scRNA commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_scrna_{command} {params}"
            else:
                return f"bio_scrna_{command}"
    
    
    def _handle_gene_agent_command(self, parts) -> str:
        """Handle GeneAgent commands with special logic"""
        
        if len(parts) == 2:
            # Just /bio GeneAgent - show GeneAgent help
            self.console.print("\n[bold]🧬 GeneAgent - Gene Set Analysis[/bold]")
            self.console.print("[dim]/bio GeneAgent <genes>[/dim] - Analyze gene set (e.g., TP53,BRCA1,EGFR)")
            self.console.print("[dim]/bio GeneAgent <genes> --analysis_type <type>[/dim] - Specific analysis type")
            self.console.print("[dim]/bio GeneAgent <genes> --output_format <format>[/dim] - Output format")
            self.console.print("\n[dim]Analysis Types:[/dim]")
            self.console.print("[dim]  comprehensive - Full analysis (functional, pathways, interactions, clinical)[/dim]")
            self.console.print("[dim]  functional - Biological functions and processes[/dim]")
            self.console.print("[dim]  enrichment - GO/KEGG enrichment analysis[/dim]")
            self.console.print("[dim]  interactions - Protein-protein interactions[/dim]")
            self.console.print("[dim]  clinical - Disease associations and drug targets[/dim]")
            self.console.print("[dim]  custom - Answer custom questions[/dim]")
            self.console.print("\n[dim]Output Formats:[/dim]")
            self.console.print("[dim]  detailed - Full detailed analysis (default)[/dim]")
            self.console.print("[dim]  summary - Concise summary[/dim]")
            self.console.print("[dim]  structured - JSON-structured output[/dim]")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio GeneAgent TP53,BRCA1,EGFR                          # Basic analysis[/dim]")
            self.console.print("[dim]  /bio GeneAgent MYC,JUN,FOS --analysis_type interactions  # Interaction analysis[/dim]")
            self.console.print("[dim]  /bio GeneAgent CD4,CD8A,IL2 --output_format summary      # Summary output[/dim]")
            self.console.print("[dim]  /bio GeneAgent IFNG,TNF,IL6 --save_results true          # Save results[/dim]")
            self.console.print()
            return None
        
        # Extract gene list from first parameter
        genes = parts[2]
        
        # Parse additional parameters
        analysis_type = "comprehensive"
        output_format = "detailed"
        save_results = False
        custom_questions = []
        
        # Simple parameter parsing
        i = 3
        while i < len(parts):
            if parts[i] == "--analysis_type" and i + 1 < len(parts):
                analysis_type = parts[i + 1]
                i += 2
            elif parts[i] == "--output_format" and i + 1 < len(parts):
                output_format = parts[i + 1]
                i += 2
            elif parts[i] == "--save_results" and i + 1 < len(parts):
                save_results = parts[i + 1].lower() in ["true", "1", "yes"]
                i += 2
            elif parts[i] == "--custom_questions":
                # Collect all remaining parts as questions
                questions = parts[i + 1:]
                custom_questions = [q.strip('"\'') for q in questions if q.strip()]
                break
            else:
                i += 1
        
        # Build the gene agent command
        self.console.print(f"\n[bold cyan]🧬 Starting GeneAgent Analysis[/bold cyan]")
        self.console.print(f"[dim]Genes: {genes}[/dim]")
        self.console.print(f"[dim]Analysis type: {analysis_type}[/dim]")
        self.console.print(f"[dim]Output format: {output_format}[/dim]")
        if custom_questions:
            self.console.print(f"[dim]Custom questions: {len(custom_questions)} questions[/dim]")
        self.console.print("[dim]Powered by Pantheon-CLI's built-in Agent capabilities...[/dim]\n")
        
        # Format command for the gene agent toolset
        cmd_parts = [f"gene_agent GeneAgent {genes}"]
        cmd_parts.append(f"--analysis_type {analysis_type}")
        cmd_parts.append(f"--output_format {output_format}")
        if save_results:
            cmd_parts.append(f"--save_results {save_results}")
        if custom_questions:
            questions_str = " ".join([f'"{q}"' for q in custom_questions])
            cmd_parts.append(f"--custom_questions {questions_str}")
        
        return " ".join(cmd_parts)
    
    def _handle_rna_command(self, parts) -> str:
        """Handle RNA-specific commands with special logic"""
        
        if len(parts) == 2:
            # Just /bio rna - show RNA help
            self.console.print("\n[bold]🧬 RNA-seq Analysis Helper[/bold]")
            self.console.print("[dim]/bio rna init[/dim] - Enter RNA-seq analysis mode")
            self.console.print("[dim]/bio rna upstream <folder>[/dim] - Run upstream RNA-seq analysis on folder")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio rna init                     # Enter RNA mode[/dim]")
            self.console.print("[dim]  /bio rna upstream ./fastq_data   # Analyze FASTQ data[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter RNA mode - simple mode activation without automation
            self.console.print("\n[bold cyan]🧬 Entering RNA-seq Analysis Mode[/bold cyan]")
            
            # Clear all existing todos when entering RNA mode
            clear_message = """
RNA INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any rna.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  RNA init ready • todos={len(todos)}
"""
            
            self.console.print("[dim]Clearing existing todos and preparing RNA environment...[/dim]")
            self.console.print("[dim]Ready for RNA-seq analysis assistance...[/dim]")
            self.console.print("[dim]RNA-seq mode activated. You can now use RNA tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio rna init - Enter RNA mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio rna upstream <folder> - Run upstream RNA analysis on specific folder[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "upstream":
            # Handle upstream analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a folder path[/red]")
                self.console.print("[dim]Usage: /bio rna upstream <folder_path>[/dim]")
                self.console.print("[dim]Example: /bio rna upstream ./fastq_data[/dim]")
                return None
                
            try:
                from ..cli.prompt.rna_bulk_upstream import generate_rna_analysis_message
                
                folder_path = parts[3]
                self.console.print(f"\n[bold cyan]🧬 Starting RNA-seq Analysis[/bold cyan]")
                self.console.print(f"[dim]Target folder: {folder_path}[/dim]")
                self.console.print("[dim]Preparing RNA-seq analysis pipeline...[/dim]\n")
                
                # Generate the analysis message with folder
                rna_message = generate_rna_analysis_message(folder_path=folder_path)
                
                self.console.print("[dim]Sending RNA-seq analysis request...[/dim]\n")
                
                return rna_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: RNA module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing analysis: {str(e)}[/red]")
                return None
        
        else:
            # Handle other RNA commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_rna_{command} {params}"
            else:
                return f"bio_rna_{command}"
    
    def _handle_dock_command(self, parts) -> str:
        """Handle molecular docking commands"""
        
        if len(parts) == 2:
            # Just /bio dock - show dock help
            self.console.print("\n[bold]🧬 Molecular Docking Analysis Helper[/bold]")
            self.console.print("[dim]/bio dock init[/dim] - Initialize molecular docking project")
            self.console.print("[dim]/bio dock check[/dim] - Check dependencies (meeko, vina, pymol)")
            self.console.print("[dim]/bio dock run_dock[/dim] - Interactive docking workflow (no path needed)")
            self.console.print("[dim]/bio dock run <folder>[/dim] - Run batch docking analysis on specific folder")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio dock init                   # Initialize docking project[/dim]")
            self.console.print("[dim]  /bio dock check                  # Check and install dependencies[/dim]")
            self.console.print("[dim]  /bio dock run_dock              # Interactive docking workflow[/dim]")
            self.console.print("[dim]  /bio dock run ./docking_data    # Run docking analysis on folder[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter molecular docking mode - simple mode activation
            self.console.print("\n[bold cyan]🧬 Initializing Molecular Docking Mode[/bold cyan]")
            
            # Clear all existing todos when entering dock mode
            clear_message = """
DOCK INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any dock.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  Dock init ready • todos={len(todos)}
"""
            
            self.console.print("[dim]Clearing existing todos and preparing docking environment...[/dim]")
            self.console.print("[dim]Ready for molecular docking analysis assistance...[/dim]")
            self.console.print("[dim]Dock mode activated. You can now use docking tools directly.[/dim]")
            self.console.print()
            self.console.print("[dim]The command structure is now clean:[/dim]")
            self.console.print("[dim]  - /bio dock init - Enter dock mode (simple prompt loading)[/dim]")
            self.console.print("[dim]  - /bio dock run <folder> - Run docking analysis on specific folder[/dim]")
            self.console.print()
            
            return clear_message
        
        elif command == "check":
            # Check dependencies
            self.console.print("\n[bold cyan]🔍 Checking Molecular Docking Dependencies[/bold cyan]")
            check_message = "Execute dock.Dock_Workflow('check_dependencies') to verify installations"
            self.console.print("[dim]Checking meeko, vina, pymol installations...[/dim]")
            return check_message
        
        elif command == "run_dock":
            # Handle interactive docking workflow (no path required)
            self.console.print("\n[bold cyan]🧬 Starting Interactive Molecular Docking Workflow[/bold cyan]")
            self.console.print("[dim]Preparing interactive docking pipeline...[/dim]\n")
            
            try:
                # Try importing the module first
                try:
                    from ..cli.prompt.dock_interactive import generate_interactive_dock_message
                except ImportError:
                    # Fallback: execute the file directly
                    import os
                    dock_interactive_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
                                                        'cli', 'prompt', 'dock_interactive.py')
                    namespace = {}
                    exec(open(dock_interactive_path).read(), namespace)
                    generate_interactive_dock_message = namespace['generate_interactive_dock_message']
                
                # Generate the interactive workflow message
                dock_message = generate_interactive_dock_message()
                
                self.console.print("[dim]Sending interactive docking workflow request...[/dim]\n")
                
                return dock_message
                
            except Exception as e:
                self.console.print(f"[red]Error generating interactive dock message: {e}[/red]")
                return None
        
        elif command == "run":
            # Handle docking analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify a folder path[/red]")
                self.console.print("[dim]Usage: /bio dock run <folder_path>[/dim]")
                self.console.print("[dim]Example: /bio dock run ./docking_data[/dim]")
                return None
                
            try:
                from ..cli.prompt.dock_molecular import generate_dock_analysis_message
                
                folder_path = parts[3]
                message = generate_dock_analysis_message(folder_path)
                
                self.console.print(f"\n[bold cyan]🧬 Molecular Docking Analysis Initiated[/bold cyan]")
                self.console.print(f"[green]Target folder: {folder_path}[/green]")
                self.console.print("[dim]Processing docking workflow...[/dim]")
                
                return message
                
            except ImportError as e:
                self.console.print(f"[red]Error importing dock prompt: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error generating dock message: {e}[/red]")
                return None
        
        else:
            # Handle other dock commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_dock_{command} {params}"
            else:
                return f"bio_dock_{command}"
    
    def _handle_hic_command(self, parts) -> str:
        """Handle Hi-C analysis commands"""
        
        if len(parts) == 2:
            # Just /bio hic - show hic help
            self.console.print("\n[bold]🧬 Hi-C Analysis Helper[/bold]")
            self.console.print("[dim]/bio hic init[/dim] - Initialize Hi-C project")
            self.console.print("[dim]/bio hic upstream <folder>[/dim] - Run Hi-C upstream analysis")
            self.console.print("[dim]/bio hic analysis <folder>[/dim] - Run Hi-C downstream analysis")
            self.console.print("\n[dim]Examples:[/dim]")
            self.console.print("[dim]  /bio hic init                   # Initialize Hi-C project[/dim]")
            self.console.print("[dim]  /bio hic upstream ./hic_data   # Run upstream processing[/dim]")
            self.console.print("[dim]  /bio hic analysis ./hic_data   # Run TAD/compartment analysis[/dim]")
            self.console.print()
            return None
        
        command = parts[2]
        
        if command == "init":
            # Enter Hi-C mode - strict init mode like dock
            self.console.print("\n[bold cyan]🧬 Entering Hi-C Analysis Mode[/bold cyan]")
            
            # Clear all existing todos when entering Hi-C mode
            clear_message = """
HiC INIT MODE — STRICT

Goal: ONLY clear TodoList and report the new status. Do NOT create or execute anything.

Allowed tools (whitelist):
  - clear_all_todos()
  - show_todos()

Hard bans (do NOT call under any circumstance in init):
  - add_todo(), mark_task_done(), execute_current_task()
  - any hic.* analysis tools

Steps:
  1) clear_all_todos()
  2) todos = show_todos()

Response format (single line):
  Hi-C init ready • todos={len(todos)}
"""
            return clear_message
        
        elif command == "upstream":
            # Handle Hi-C upstream analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify folder path[/red]")
                self.console.print("[dim]Usage: /bio hic upstream <folder>[/dim]")  
                self.console.print("[dim]Example: /bio hic upstream ./hic_data[/dim]")
                return None
            
            folder_path = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Hi-C Upstream Analysis Initiated[/bold cyan]")
            self.console.print(f"[dim]Target folder: {folder_path}[/dim]")
            self.console.print("[dim]Processing Hi-C data: QC, alignment, matrix building...[/dim]")
            self.console.print()
            
            try:
                # Import and use Hi-C upstream analysis message
                from ..cli.prompt.hic_analysis import generate_hic_analysis_message
                
                # Generate the analysis message with folder
                hic_message = generate_hic_analysis_message(folder_path=folder_path)
                
                self.console.print("[dim]Sending Hi-C upstream analysis request...[/dim]")
                
                return hic_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: Hi-C upstream module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing Hi-C upstream analysis: {str(e)}[/red]")
                return None
        
        elif command == "analysis":
            # Handle Hi-C downstream analysis
            if len(parts) < 4:
                self.console.print("[red]Error: Please specify folder path[/red]")
                self.console.print("[dim]Usage: /bio hic analysis <folder>[/dim]")
                self.console.print("[dim]Example: /bio hic analysis ./hic_data[/dim]")
                return None
            
            folder_path = parts[3]
            self.console.print(f"\n[bold cyan]🧬 Hi-C Analysis Initiated[/bold cyan]")
            self.console.print(f"[dim]Target folder: {folder_path}[/dim]")
            self.console.print("[dim]Processing Hi-C analysis: TADs, compartments, loops...[/dim]")
            self.console.print()
            
            try:
                # Import and use Hi-C analysis message
                from ..cli.prompt.hic_analysis import generate_hic_analysis_message
                
                # Generate the analysis message with folder
                hic_message = generate_hic_analysis_message(folder_path=folder_path)
                
                self.console.print("[dim]Sending Hi-C analysis request...[/dim]")
                
                return hic_message
                
            except ImportError as e:
                self.console.print(f"[red]Error: Hi-C analysis module not available: {e}[/red]")
                return None
            except Exception as e:
                self.console.print(f"[red]Error preparing Hi-C analysis: {str(e)}[/red]")
                return None
        
        else:
            # Handle other hic commands generically
            params = " ".join(parts[3:]) if len(parts) > 3 else ""
            if params:
                return f"bio_hic_{command} {params}"
            else:
                return f"bio_hic_{command}"
    
    async def handle_deprecated_atac_command(self, command: str) -> str:
        """
        Handle deprecated /atac commands with migration and auto-conversion
        
        Returns:
            str: Converted bio command message, or None if no conversion
        """
        parts = command.split(maxsplit=2)
        
        # Show deprecation warning
        self.console.print("\n[bold yellow]⚠️  Command Migration Notice[/bold yellow]")
        self.console.print("[yellow]ATAC commands have moved to the unified bio interface![/yellow]")
        
        if len(parts) == 1:
            # Just /atac - show migration help
            self._show_atac_migration_help()
            return None
        
        # Auto-convert old commands to new bio commands
        if parts[1] == "init":
            self.console.print("\n[bold cyan]→ Auto-converting to: /bio atac init[/bold cyan]")
            return "bio_atac_init"
        
        elif parts[1] == "upstream":
            # Auto-convert upstream command
            if len(parts) < 3:
                self.console.print("[red]Error: Please specify a folder path[/red]")
                self.console.print("[dim]New usage: /bio atac upstream <folder_path>[/dim]")
                self.console.print("[dim]Example: /bio atac upstream ./fastq_data[/dim]")
                return None
            
            folder_path = parts[2]
            self.console.print(f"\n[bold cyan]→ Auto-converting to: /bio atac upstream {folder_path}[/bold cyan]")
            return f"bio_atac_upstream {folder_path}"
        
        else:
            self.console.print(f"[red]Unknown ATAC command: {parts[1]}[/red]")
            self.console.print("[dim]Please use the new bio interface instead:[/dim]")
            self.console.print("[dim]  /bio atac init - Initialize ATAC project[/dim]")
            self.console.print("[dim]  /bio atac upstream <folder> - Run upstream ATAC analysis[/dim]")
            return None
    
    def _show_atac_migration_help(self):
        """Show ATAC migration help"""
        self.console.print("\n[dim]Old commands → New commands:[/dim]")
        self.console.print("[dim]/atac init → /bio atac init[/dim]")
        self.console.print("[dim]/atac upstream <folder> → /bio atac upstream <folder>[/dim]")
        self.console.print("\n[bold cyan]🧬 Available Bio Commands[/bold cyan]")
        self.console.print("[dim]/bio list[/dim] - List all available bio tools")
        self.console.print("[dim]/bio atac init[/dim] - Initialize ATAC-seq project")
        self.console.print("[dim]/bio atac upstream <folder>[/dim] - Run upstream ATAC analysis")
        self.console.print("")


# Command mapping for easy extension
BIO_COMMAND_MAP = {
    # Direct bio manager commands
    'list': 'bio list',
    'help': 'bio help',
    'info': 'bio info',
    
    # ATAC-seq commands
    'atac_init': 'bio_atac_init',
    'atac_upstream': 'bio_atac_upstream',
    'atac_check_dependencies': 'bio_atac_check_dependencies',
    'atac_setup_genome_resources': 'bio_atac_setup_genome_resources',
    'atac_auto_align_fastq': 'bio_atac_auto_align_fastq',
    'atac_call_peaks_macs2': 'bio_atac_call_peaks_macs2',
    'atac_generate_atac_qc_report': 'bio_atac_generate_atac_qc_report',
    
    # scRNA-seq commands  
    'scrna_init': 'bio_scrna_init',
    'scrna_load_data': 'bio_scrna_load_data',
    'scrna_run_quality_control': 'bio_scrna_run_quality_control',
    'scrna_run_preprocessing': 'bio_scrna_run_preprocessing', 
    'scrna_run_pca': 'bio_scrna_run_pca',
    'scrna_run_clustering': 'bio_scrna_run_clustering',
    'scrna_run_cell_type_annotation': 'bio_scrna_run_cell_type_annotation',
    
    
    # RNA-seq commands
    'rna_init': 'bio_rna_init',
    'rna_upstream': 'bio_rna_upstream',
    'rna_check_dependencies': 'bio_rna_check_dependencies',
    'rna_setup_genome_resources': 'bio_rna_setup_genome_resources',
    'rna_align_star': 'bio_rna_align_star',
    'rna_differential_expression': 'bio_rna_differential_expression',
    
    # ChIP-seq commands (for future use)
    'chipseq_init': 'bio_chipseq_init',
    'chipseq_call_peaks': 'bio_chipseq_call_peaks',
    'chipseq_find_motifs': 'bio_chipseq_find_motifs',
}

# Deprecated command conversions
DEPRECATED_ATAC_MAP = {
    '/atac init': '/bio atac init',
    '/atac upstream': '/bio atac upstream',
}

def get_bio_command_suggestions() -> list:
    """Get list of available bio command suggestions for autocomplete"""
    suggestions = [
        '/bio list',
        '/bio help',
        '/bio info atac',
        '/bio atac init',
        '/bio atac upstream',
        '/bio atac check_dependencies',
        '/bio atac setup_genome_resources',
        '/bio scatac init',
        '/bio scatac upstream',
        '/bio scrna init',
        '/bio scrna analysis',
        '/bio scrna qc',
        '/bio scrna preprocess',
        '/bio scrna annotate',
        '/bio rna init',
        '/bio rna upstream',
        '/bio GeneAgent',
        '/bio GeneAgent TP53,BRCA1,EGFR',
        '/bio chipseq init',  # Future
    ]
    return suggestions