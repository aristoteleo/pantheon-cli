import os
import platform
from typing import Dict, Optional
from rag_components import (
    load_documents, split_documents, initialize_embeddings, create_vector_store,
    load_vector_store, save_vector_store
)
from rag_chain_gemini import setup_rag_chain_gemini, query_rag_system_gemini


class GeminiSystemCheck:
    def __init__(self):
        self.system = platform.system()
        self.required_packages = [
            "langchain", "langchain-community", "langchain-google-genai",
            "sentence-transformers", "numpy",
        ]

    def check_gemini_api_key(self) -> bool:
        """Check if GEMINI_API_KEY is available"""
        api_key = os.getenv("GEMINI_API_KEY")
        return api_key is not None and len(api_key.strip()) > 0

    def check_python_packages(self) -> Dict[str, bool]:
        package_status = {}
        for package in self.required_packages:
            try:
                __import__(package.replace("-", "_"))
                package_status[package] = True
            except ImportError:
                package_status[package] = False
        return package_status

    def check_system_resources(self) -> Dict[str, bool]:
        """Basic system checks - less stringent than Ollama requirements"""
        import psutil
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')
        return {
            "memory": memory.available > 2 * 1024 * 1024 * 1024,  # 2GB (reduced from 4GB)
            "disk": disk.free > 5 * 1024 * 1024 * 1024,  # 5GB (reduced from 10GB)
        }


def initialize_rag_system_gemini(data_directory, vector_store_path="omicverse_vectorstore_gemini", model_name="gemini-2.5-pro"):
    """
    Initialize RAG system using Gemini API instead of Ollama.
    Much simpler setup - no local model management needed.
    """
    checker = GeminiSystemCheck()

    print("🔍 Checking system requirements for Gemini RAG...")
    
    # Check system resources (less stringent)
    system_resources = checker.check_system_resources()
    if not all(system_resources.values()):
        print("⚠️  Insufficient system resources:")
        for resource, status in system_resources.items():
            print(f"- {resource}: {'✓' if status else '✗'}")
        return None

    print("📦 Checking required Python packages...")
    package_status = checker.check_python_packages()
    missing_packages = [pkg for pkg, status in package_status.items() if not status]
    if missing_packages:
        print("❌ Missing required packages:")
        print(f"Please install: pip install {' '.join(missing_packages)}")
        return None

    print("🔑 Checking Gemini API key...")
    if not checker.check_gemini_api_key():
        print("❌ GEMINI_API_KEY environment variable not found.")
        print("Please set your Gemini API key:")
        print("export GEMINI_API_KEY='your-api-key-here'")
        return None

    try:
        print(f"🚀 Initializing RAG system with {model_name}...")
        embeddings = initialize_embeddings()

        if os.path.exists(vector_store_path):
            print("📚 Loading existing vector store...")
            vector_store = load_vector_store(vector_store_path, embeddings)
        else:
            print("🔨 Creating new vector store...")
            documents = load_documents(data_directory)
            texts = split_documents(documents)
            vector_store = create_vector_store(texts, embeddings)
            save_vector_store(vector_store, vector_store_path)

        qa_chain = setup_rag_chain_gemini(vector_store, model_name=model_name)
        print("✅ Gemini RAG system initialized successfully!")
        return qa_chain

    except Exception as e:
        print(f"❌ Error initializing Gemini RAG system: {str(e)}")
        return None


def process_query_gemini(query, data_directory):
    """
    Process query using Gemini-powered RAG system.
    """
    qa_chain = initialize_rag_system_gemini(data_directory)
    if qa_chain is None:
        return None

    try:
        print(f"🧠 Processing query with Gemini: {query}")
        result = query_rag_system_gemini(qa_chain, query)
        return result
    except Exception as e:
        print(f"❌ Error processing query: {str(e)}")
        return None


def demo_gemini_rag():
    """
    Demo function to test Gemini RAG integration
    """
    # Use current directory structure
    data_dir = "Converted_Scripts_Annotated"  # OV_Agent annotated scripts
    
    # Example queries for testing
    test_queries = [
        "How do I preprocess single-cell RNA-seq data?",
        "What's the best way to cluster cells in scanpy?",
        "How to find marker genes for each cluster?",
        "Show me how to create a UMAP visualization",
        "What quality control steps should I perform?"
    ]
    
    print("🧪 Gemini RAG System Demo")
    print("=" * 50)
    
    for i, query in enumerate(test_queries, 1):
        print(f"\n📝 Query {i}: {query}")
        result = process_query_gemini(query, data_dir)
        
        if result:
            print(f"🤖 Gemini Response:")
            print(result['answer'][:500] + "..." if len(result['answer']) > 500 else result['answer'])
            print(f"\n📚 Sources: {len(result['source_documents'])} documents")
        else:
            print("❌ No response generated")
        print("-" * 50)


if __name__ == "__main__":
    # Check if running as demo or with specific query
    import sys
    
    if len(sys.argv) > 1:
        # Process specific query
        query = " ".join(sys.argv[1:])
        data_dir = "Converted_Scripts_Annotated"
        result = process_query_gemini(query, data_dir)
        
        if result:
            print(f"\n🤖 **Gemini Response:**")
            print(result['answer'])
            print(f"\n📚 **Sources:**")
            for doc in result['source_documents']:
                print(f"- {doc.metadata.get('source', 'Unknown source')}")
        else:
            print("❌ Failed to generate response")
    else:
        # Run demo
        demo_gemini_rag()