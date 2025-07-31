#!/usr/bin/env python3
"""
H5AD RAG Query System - 专门为H5AD分析工具设计的RAG查询系统
支持进度跟踪和详细的结果返回
"""

import os
import sys
import json
import argparse
from typing import Dict, Any, List
from rag_components_enhanced import (
    ProgressTracker,
    load_documents_with_progress,
    split_documents_with_progress,
    initialize_embeddings_with_progress,
    create_vector_store_with_progress,
    load_vector_store_with_progress,
    save_vector_store_with_progress,
    query_vector_store_with_progress
)
from rag_chain_gemini import setup_rag_chain_gemini, query_rag_system_gemini

class H5ADRagQuerySystem:
    """H5AD分析专用的RAG查询系统"""
    
    def __init__(self, data_directory: str = "Converted_Scripts_Annotated", 
                 vector_store_path: str = "omicverse_vectorstore_gemini",
                 model_name: str = "gemini-2.5-pro"):
        # Use absolute paths relative to this script's directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.data_directory = os.path.join(script_dir, data_directory) if not os.path.isabs(data_directory) else data_directory
        self.vector_store_path = os.path.join(script_dir, vector_store_path) if not os.path.isabs(vector_store_path) else vector_store_path
        self.model_name = model_name
        self.progress_tracker = ProgressTracker(self._print_progress)
        self.qa_chain = None
        self.vector_store = None
        
    def _print_progress(self, message: str):
        """打印进度信息到stderr，避免与结果混淆"""
        print(message, file=sys.stderr, flush=True)
    
    def initialize_rag_system(self) -> bool:
        """初始化RAG系统，返回是否成功"""
        try:
            self.progress_tracker.update("🚀 Initializing H5AD RAG Query System...")
            
            # Check API key
            if not os.getenv("GEMINI_API_KEY"):
                self.progress_tracker.update("❌ GEMINI_API_KEY environment variable not found")
                return False
            
            # Initialize embeddings
            embeddings = initialize_embeddings_with_progress(self.progress_tracker)
            
            # Load or create vector store
            if os.path.exists(self.vector_store_path):
                self.progress_tracker.update("📚 Loading existing vector store...")
                self.vector_store = load_vector_store_with_progress(
                    self.vector_store_path, embeddings, self.progress_tracker
                )
            else:
                self.progress_tracker.update("🔨 Creating new vector store...")
                
                # Load and process documents
                documents = load_documents_with_progress(self.data_directory, self.progress_tracker)
                if not documents:
                    self.progress_tracker.update(f"❌ No documents found in {self.data_directory}")
                    return False
                
                texts = split_documents_with_progress(documents, self.progress_tracker)
                self.vector_store = create_vector_store_with_progress(texts, embeddings, self.progress_tracker)
                save_vector_store_with_progress(self.vector_store, self.vector_store_path, self.progress_tracker)
            
            # Setup QA chain
            self.progress_tracker.update(f"🤖 Setting up QA chain with {self.model_name}...")
            self.qa_chain = setup_rag_chain_gemini(self.vector_store, model_name=self.model_name)
            
            self.progress_tracker.update("✅ H5AD RAG Query System initialized successfully!")
            return True
            
        except Exception as e:
            self.progress_tracker.update(f"❌ Error initializing RAG system: {str(e)}")
            return False
    
    def query_for_h5ad_analysis(self, analysis_task: str, h5ad_path: str = "", 
                               parameters: str = "") -> Dict[str, Any]:
        """为H5AD分析查询RAG系统"""
        
        if not self.qa_chain:
            return {
                "success": False,
                "error": "RAG system not initialized"
            }
        
        try:
            # 构建专门的查询
            enhanced_query = f"""
            Single-cell RNA-seq Analysis Task: {analysis_task}
            
            File: {h5ad_path if h5ad_path else "H5AD format file"}
            Parameters: {parameters if parameters else "Standard parameters"}
            
            Please provide:
            1. Complete Python code using scanpy and/or omicverse
            2. Step-by-step explanation of the analysis
            3. Parameter validation and error handling
            4. Expected outputs and file saving
            5. Best practices for this specific analysis type
            
            Focus on:
            - Proper h5ad file loading and handling
            - Quality control steps if needed
            - Appropriate preprocessing for the analysis type
            - Clear variable names and comments
            - Saving results with descriptive filenames
            
            Return executable Python code that can run independently.
            """
            
            self.progress_tracker.update(f"🧠 Querying RAG system for: {analysis_task}")
            
            # Get relevant documents first
            relevant_docs = query_vector_store_with_progress(
                self.vector_store, enhanced_query, k=10, progress_tracker=self.progress_tracker
            )
            
            # Query the full RAG system
            result = query_rag_system_gemini(self.qa_chain, enhanced_query)
            
            self.progress_tracker.update("✅ RAG query completed successfully")
            
            return {
                "success": True,
                "answer": result["answer"],
                "source_documents": [
                    {
                        "content": doc.page_content[:500] + "..." if len(doc.page_content) > 500 else doc.page_content,
                        "metadata": doc.metadata,
                        "source": doc.metadata.get("source", "unknown")
                    }
                    for doc in result["source_documents"]
                ],
                "source_count": len(result["source_documents"]),
                "model_used": result.get("model_used", self.model_name),
                "query": enhanced_query
            }
            
        except Exception as e:
            self.progress_tracker.update(f"❌ Error querying RAG system: {str(e)}")
            return {
                "success": False,
                "error": str(e)
            }

def main():
    """主函数，处理命令行调用"""
    parser = argparse.ArgumentParser(description="H5AD RAG Query System")
    parser.add_argument("analysis_task", help="Analysis task description")
    parser.add_argument("--h5ad-path", default="", help="Path to H5AD file")
    parser.add_argument("--parameters", default="", help="Additional parameters")
    parser.add_argument("--data-dir", default="Converted_Scripts_Annotated", 
                       help="Data directory for RAG documents")
    parser.add_argument("--vector-store", default="omicverse_vectorstore_gemini",
                       help="Vector store path")
    parser.add_argument("--model", default="gemini-2.5-pro", help="Gemini model to use")
    
    args = parser.parse_args()
    
    # Create and initialize RAG system
    rag_system = H5ADRagQuerySystem(
        data_directory=args.data_dir,
        vector_store_path=args.vector_store,
        model_name=args.model
    )
    
    # Initialize system
    if not rag_system.initialize_rag_system():
        sys.exit(1)
    
    # Query for analysis
    result = rag_system.query_for_h5ad_analysis(
        analysis_task=args.analysis_task,
        h5ad_path=args.h5ad_path,
        parameters=args.parameters
    )
    
    # Output result as JSON to stdout
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()