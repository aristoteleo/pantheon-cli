import os
import platform
import subprocess
import time
import requests
import psutil
from typing import Dict, Optional
from rag_components import (
    load_documents, split_documents, initialize_embeddings, create_vector_store,
    load_vector_store, save_vector_store
)
from rag_chain import setup_rag_chain, query_rag_system
from rag_chain_gemini import setup_rag_chain_gemini, query_rag_system_gemini


class OllamaSystemCheck:
    def __init__(self):
        self.system = platform.system()
        self.ollama_port = 11434
        self.required_packages = [
            "langchain", "langchain-community",
            "sentence-transformers", "numpy"
        ]

    def check_ollama_installation(self) -> bool:
        try:
            if self.system == "Windows":
                result = subprocess.run(["where", "ollama"], capture_output=True, text=True)
            else:
                result = subprocess.run(["which", "ollama"], capture_output=True, text=True)
            return result.returncode == 0
        except Exception:
            return False

    def is_ollama_running(self) -> bool:
        try:
            response = requests.get(f"http://localhost:{self.ollama_port}/api/tags")
            return response.status_code == 200
        except requests.exceptions.RequestException:
            return False

    def start_ollama(self) -> bool:
        if self.is_ollama_running():
            return True

        try:
            if self.system == "Windows":
                subprocess.Popen(["ollama", "serve"], creationflags=subprocess.CREATE_NO_WINDOW)
            else:
                subprocess.Popen(["ollama", "serve"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            for _ in range(10):
                time.sleep(1)
                if self.is_ollama_running():
                    return True
            return False
        except Exception as e:
            print(f"Error starting Ollama: {str(e)}")
            return False

    def check_model_availability(self, model_name: str = "llama3.2:3b") -> bool:
        try:
            response = requests.get(f"http://localhost:{self.ollama_port}/api/tags")
            if response.status_code == 200:
                models = response.json().get("models", [])
                return any(model["name"] == model_name for model in models)
            return False
        except requests.exceptions.RequestException:
            return False

    def pull_model(self, model_name: str = "llama3.2:3b") -> bool:
        try:
            subprocess.run(["ollama", "pull", model_name], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False

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
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')
        return {
            "memory": memory.available > 4 * 1024 * 1024 * 1024,  # 4GB
            "disk": disk.free > 10 * 1024 * 1024 * 1024,  # 10GB
        }


def initialize_rag_system(data_directory, vector_store_path="omicverse_vectorstore", model_name="llama3.2:3b"):
    checker = OllamaSystemCheck()

    print("Checking system requirements...")
    system_resources = checker.check_system_resources()
    if not all(system_resources.values()):
        print("Insufficient system resources:")
        for resource, status in system_resources.items():
            print(f"- {resource}: {'✓' if status else '✗'}")
        return None

    print("Checking required Python packages...")
    package_status = checker.check_python_packages()
    missing_packages = [pkg for pkg, status in package_status.items() if not status]
    if missing_packages:
        print("Missing required packages:")
        print(f"Please install: pip install {' '.join(missing_packages)}")
        return None

    print("Checking Ollama installation...")
    if not checker.check_ollama_installation():
        print("Ollama is not installed. Please install it first:")
        print("Visit: https://ollama.com/download")
        return None

    print("Starting Ollama service...")
    if not checker.start_ollama():
        print("Failed to start Ollama service")
        return None

    print(f"Checking availability of {model_name}...")
    if not checker.check_model_availability(model_name):
        print(f"Pulling {model_name}...")
        if not checker.pull_model(model_name):
            print(f"Failed to pull {model_name}")
            return None

    try:
        print("Initializing RAG system...")
        embeddings = initialize_embeddings()

        if os.path.exists(vector_store_path):
            print("Loading existing vector store...")
            vector_store = load_vector_store(vector_store_path, embeddings)
        else:
            print("Creating new vector store...")
            documents = load_documents(data_directory)
            texts = split_documents(documents)
            vector_store = create_vector_store(texts, embeddings)
            save_vector_store(vector_store, vector_store_path)

        qa_chain = setup_rag_chain(vector_store, model_name=model_name)
        print("RAG system initialized successfully!")
        return qa_chain

    except Exception as e:
        print(f"Error initializing RAG system: {str(e)}")
        return None



def process_query(query, data_directory, use_gemini=True):
    """
    Process query using either Gemini API (default) or Ollama local model.
    Gemini API is faster, more reliable, and doesn't require local model management.
    """
    if use_gemini:
        # Use Gemini API - simpler and more reliable
        from main_gemini import process_query_gemini
        return process_query_gemini(query, data_directory)
    else:
        # Use Ollama local model - requires more setup
        qa_chain = initialize_rag_system(data_directory)
        if qa_chain is None:
            return None

        try:
            result = query_rag_system(qa_chain, query)
            return result
        except Exception as e:
            print(f"Error processing query: {str(e)}")
            return None



if __name__ == "__main__":
    import sys
    
    # Use current directory structure for Pantheon CLI
    data_dir = "Converted_Scripts_Annotated"  # OV_Agent annotated scripts
    
    # Check command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--ollama":
            # Use Ollama local model
            query = " ".join(sys.argv[2:]) if len(sys.argv) > 2 else "How do I preprocess single-cell RNA-seq data?"
            print("🦙 Using Ollama local model...")
            result = process_query(query, data_dir, use_gemini=False)
        else:
            # Use Gemini API (default)
            query = " ".join(sys.argv[1:])
            print("✨ Using Gemini API...")
            result = process_query(query, data_dir, use_gemini=True)
    else:
        # Default demo with Gemini API
        query = "How do I preprocess single-cell RNA-seq data using OmicVerse?"
        print("✨ Using Gemini API (default)...")
        result = process_query(query, data_dir, use_gemini=True)
    
    if result:
        print(f"\n🤖 **Answer:** {result['answer']}")
        if 'source_documents' in result:
            print(f"\n📚 **Sources:**")
            for doc in result['source_documents']:
                print(f"- {doc.metadata.get('source', 'Unknown source')}")
    else:
        print("❌ Failed to generate response")