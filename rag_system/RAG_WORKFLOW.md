# Pantheon-CLI RAG System Workflow

## Overview Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        USER INTERACTION FLOW                             │
└─────────────────────────────────────────────────────────────────────────┘

1. User Input in CLI
   ┌─────────────────────┐
   │ $ pantheon /h5ad    │
   │ > Analyze my.h5ad   │
   │ > "Find marker      │
   │   genes"            │
   └──────────┬──────────┘
              │
              ▼
2. Pantheon CLI (TypeScript)
   ┌─────────────────────┐
   │ h5adCommand.ts      │
   │ ├─ Parse input      │
   │ └─ Call H5AD Tool   │
   └──────────┬──────────┘
              │
              ▼
3. H5AD Analyzer Tool
   ┌─────────────────────┐
   │ h5ad-analyzer.ts    │
   │ ├─ Validate h5ad    │
   │ └─ Spawn Python     │
   └──────────┬──────────┘
              │
              ▼

┌─────────────────────────────────────────────────────────────────────────┐
│                          RAG SYSTEM (Python)                             │
└─────────────────────────────────────────────────────────────────────────┘

4. Query Enhancement
   ┌─────────────────────┐
   │ h5ad_rag_query.py   │
   │ ├─ Load task        │
   │ └─ Enhance query    │
   └──────────┬──────────┘
              │
              ▼
5. Vector Store Search
   ┌─────────────────────────────────────┐
   │     FAISS Vector Store              │
   │ ┌─────────────────────────────────┐ │
   │ │ Embedded Knowledge Base:        │ │
   │ │ • t_deg_annotated.py           │ │
   │ │ • t_cluster_annotated.py       │ │
   │ │ • t_preprocess_annotated.py    │ │
   │ │ • ... (50+ tutorials)          │ │
   │ └─────────────────────────────────┘ │
   │                                     │
   │ Similarity Search (k=46 docs) ────► │
   └──────────┬──────────────────────────┘
              │ Retrieved Context
              ▼
6. LLM Generation
   ┌─────────────────────┐
   │ Gemini API          │
   │ ├─ Context + Query  │
   │ └─ Generate Code    │
   └──────────┬──────────┘
              │
              ▼
7. Response Processing
   ┌─────────────────────┐
   │ Format Response     │
   │ ├─ Python code      │
   │ ├─ Explanations     │
   │ └─ Best practices   │
   └──────────┬──────────┘
              │
              ▼

┌─────────────────────────────────────────────────────────────────────────┐
│                         OUTPUT TO USER                                   │
└─────────────────────────────────────────────────────────────────────────┘

8. Display in CLI
   ┌─────────────────────────────────────┐
   │ Generated Analysis Code:            │
   │ ```python                           │
   │ import scanpy as sc                 │
   │ import omicverse as ov              │
   │                                     │
   │ # Load h5ad file                    │
   │ adata = sc.read_h5ad('my.h5ad')    │
   │                                     │
   │ # Find marker genes                 │
   │ sc.tl.rank_genes_groups(...)        │
   │ ...                                 │
   │ ```                                 │
   │                                     │
   │ Sources: 3 relevant tutorials       │
   └─────────────────────────────────────┘
```

## Detailed Component Flow

### 1. Knowledge Base Preparation (One-time setup)
```
OmicVerse Tutorials     Build Vector Store
     (*.py)       ────►  (FAISS + HuggingFace)
                              │
                              ▼
                    ┌─────────────────────┐
                    │ index.faiss         │
                    │ index.pkl           │
                    └─────────────────────┘
```

### 2. Query Processing Pipeline
```
User Query ──► Enhanced Query ──► Vector Search ──► LLM ──► Code Output
   │                │                   │             │          │
   │                │                   │             │          │
"Find         "For h5ad file,    Find top 46    Gemini    Executable
marker        find marker genes   similar code   generates  Python code
genes"        with scanpy..."     examples       solution   + explanation
```

### 3. Component Interactions
```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│  TypeScript  │────►│    Python    │────►│   Gemini     │
│   Frontend   │     │  RAG System  │     │     API      │
└──────────────┘     └──────────────┘     └──────────────┘
      │                     │                      │
      │                     ├─ Embeddings         │
      │                     ├─ Vector Store       │
      │                     └─ Query Enhancement  │
      │                                           │
      └───────────── JSON Response ───────────────┘
```

## Key Features

1. **Progress Tracking**: Shows embedding progress, search progress, and generation status
2. **Source Attribution**: Returns which tutorials were used as context
3. **Error Handling**: Validates h5ad files and provides clear error messages
4. **Flexible Output**: Can save generated code to files or display in CLI
5. **API Integration**: Uses same GEMINI_API_KEY as main Pantheon CLI

## Example Commands

```bash
# Basic usage
$ pantheon
> /h5ad
> Enter h5ad file path: data/sample.h5ad
> Enter analysis task: perform quality control and filtering

# Direct command
$ pantheon /h5ad analyze data/sample.h5ad "find differentially expressed genes"

# With parameters
$ pantheon /h5ad analyze data/sample.h5ad "cluster cells" --params "resolution=0.5"
```
