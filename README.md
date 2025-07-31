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


## Development Guide

Follow these steps to add a new tool (`h5ad-analyzer`) to the project:

1. **Create the tool and its tests**  
   - Add a new tool implementation in:  
     ```
     packages/core/src/tools/h5ad-analyzer.ts
     ```  
   - Add a corresponding Jest test file in:  
     ```
     packages/core/src/tools/h5ad-analyzer.test.ts
     ```

2. **Register the tool in the core config**  
   - Open `packages/core/src/config/index.ts` (or the appropriate config file) and add:  
     ```ts
     import { H5adAnalyzerTool } from '../tools/h5ad-analyzer.js';
     registerCoreTool(H5adAnalyzerTool);
     ```

3. **Add a CLI command for the tool’s UI**  
   - Create a new command definition in the CLI’s UI layer:  
     ```
     packages/cli/src/ui/commands/h5adCommand.ts
     ```  
   - Register the command in the built-in command loader by editing:  
     ```ts
     // In packages/cli/src/services/BuiltinCommandLoader.ts
     import { h5adCommand } from '../ui/commands/h5adCommand.js';
     // …then ensure `h5adCommand` is included in the loader’s command list
     ```


### Runtime Behavior

The agent system is designed to work even after reference folders are removed:

1. Agent definitions are scanned and cached at build/initialization time
2. Only metadata is stored, not the full agent implementations
3. Reference folders can be safely deleted after the initial scan
4. The CLI will continue to recognize agent patterns and intent mappings

## Acknowledgments

This project is based on [Google Gemini CLI](https://github.com/google-gemini/gemini-cli). We acknowledge and appreciate the excellent work of the Gemini CLI team. Our main contribution focuses on parser-level adaptations to better support Pantheon-Coder models.
