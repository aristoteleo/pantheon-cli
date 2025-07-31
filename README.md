# Pantheon CLI

<img width="1256" height="300" alt="Image" src="https://private-user-images.githubusercontent.com/46667721/472695511-ab7a5b43-0531-4870-a8a2-93f6709fb67b.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTM5MDc2MzksIm5iZiI6MTc1MzkwNzMzOSwicGF0aCI6Ii80NjY2NzcyMS80NzI2OTU1MTEtYWI3YTViNDMtMDUzMS00ODcwLWE4YTItOTNmNjcwOWZiNjdiLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA3MzAlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNzMwVDIwMjg1OVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTM5MzU1NTAzNTBkMDZmNmQ1ZmM4NDgxMWUxNzYyZWMzYTFlNDcwMDU2NzUyZWQ1ZjNmYzUyNTkyZDFiMjU3NGQmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.duSGtnaa7ONSzdNvRoqd-h6x5j4AFM1UQ9S7do5Vk-8" />

## Quickstart

You have two options to install Pantheon CLI.

### With Node

1. **Prerequisites:** Ensure you have [Node.js version 20](https://nodejs.org/en/download) or higher installed.
2. **Run the CLI:** Execute the following command in your terminal:

   ```bash
   npm run build
   ```

   Then, run the CLI from anywhere:

   ```bash
   npm run start
   ```

### Common Configuration steps

3. **Pick a color theme**
4. **Authenticate:** When prompted, sign in with your personal Google account. This will grant you up to 60 model requests per minute and 1,000 model requests per day using Gemini.

You are now ready to use the Pantheon CLI!

### Use a Gemini API key:

The Gemini API provides a free tier with [100 requests per day](https://ai.google.dev/gemini-api/docs/rate-limits#free-tier) using Gemini 2.5 Pro, control over which model you use, and access to higher rate limits (with a paid plan):

1. Generate a key from [Google AI Studio](https://aistudio.google.com/apikey).
2. Set it as an environment variable in your terminal. Replace `YOUR_API_KEY` with your generated key.

   ```bash
   export GEMINI_API_KEY="YOUR_API_KEY"
   ```

3. (Optionally) Upgrade your Gemini API project to a paid plan on the API key page (will automatically unlock [Tier 1 rate limits](https://ai.google.dev/gemini-api/docs/rate-limits#tier-1))

### Use a Vertex AI API key:

The Vertex AI API provides a [free tier](https://cloud.google.com/vertex-ai/generative-ai/docs/start/express-mode/overview) using express mode for Gemini 2.5 Pro, control over which model you use, and access to higher rate limits with a billing account:

1. Generate a key from [Google Cloud](https://cloud.google.com/vertex-ai/generative-ai/docs/start/api-keys).
2. Set it as an environment variable in your terminal. Replace `YOUR_API_KEY` with your generated key and set GOOGLE_GENAI_USE_VERTEXAI to true

   ```bash
   export GOOGLE_API_KEY="YOUR_API_KEY"
   export GOOGLE_GENAI_USE_VERTEXAI=true
   ```

3. (Optionally) Add a billing account on your project to get access to [higher usage limits](https://cloud.google.com/vertex-ai/generative-ai/docs/quotas)

For other authentication methods, including Google Workspace accounts, see the [authentication](./docs/cli/authentication.md) guide.

## Multi-Agent Architecture

Pantheon CLI now includes a powerful multi-agent architecture that automatically discovers and integrates agents from reference folders. The system supports both direct agent invocation and intelligent intent-to-function mapping.

### Available Commands

#### Agent Management

```bash
# List all available agents
pantheon-cli agents

# Show detailed agent information
pantheon-cli agents --verbose
```

#### Natural Language Processing

```bash
# Examples of natural language commands that map to functions
pantheon-cli "Normalize counts.h5ad"
pantheon-cli "Map reads.fq to hg38"
pantheon-cli "Log transform the expression data"
pantheon-cli "Cluster the cells"
pantheon-cli "Find marker genes"
```

#### RAG Knowledge Base Management

```bash
# Update RAG knowledge base with notebook content
pantheon-cli update-rag --notebook-dirs tutorials/experiments docs/ipynb

# Update from multiple directories
pantheon-cli update-rag --notebook-dirs dir1 dir2 dir3
```

### Agent Discovery

The system automatically scans the following locations for agent definitions:

- `pantheon-single-cell/` - Single-cell analysis agents
- `pantheon-agents/` - General-purpose multi-agent framework
- `packages/core/src/agents/` - Built-in agents

Agents are discovered by scanning Python files for `Agent(...)` class instantiations and extracting their metadata.

### Intent Mapping

When no specific agent matches a query, the system falls back to intent-based function mapping:

| Intent Pattern                          | Function Handler              | Description                              |
| --------------------------------------- | ----------------------------- | ---------------------------------------- |
| `normalize.*counts?\.h5ad`              | `scanpy.pp.normalize_total`   | Normalize single-cell RNA-seq counts     |
| `log.*transform`                        | `scanpy.pp.log1p`             | Log transform expression data            |
| `map.*reads?.*to.*(hg38\|mm10\|genome)` | `star_aligner`                | Map sequencing reads to reference genome |
| `cluster.*cells?`                       | `scanpy.tl.leiden`            | Cluster single cells                     |
| `find.*marker.*genes?`                  | `scanpy.tl.rank_genes_groups` | Find marker genes for clusters           |
| `umap.*visualization`                   | `scanpy.tl.umap`              | Generate UMAP visualization              |

### Runtime Behavior

The agent system is designed to work even after reference folders are removed:

1. Agent definitions are scanned and cached at build/initialization time
2. Only metadata is stored, not the full agent implementations
3. Reference folders can be safely deleted after the initial scan
4. The CLI will continue to recognize agent patterns and intent mappings

## Acknowledgments

This project is based on [Google Gemini CLI](https://github.com/google-gemini/gemini-cli). We acknowledge and appreciate the excellent work of the Gemini CLI team. Our main contribution focuses on parser-level adaptations to better support Qwen-Coder models and multi-agent architecture integration.
