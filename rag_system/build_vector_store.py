#!/usr/bin/env python3
"""
Build vector store for H5AD RAG system
"""

import os
import sys
from rag_components_enhanced import (
    ProgressTracker,
    load_documents_with_progress,
    split_documents_with_progress,
    initialize_embeddings_with_progress,
    create_vector_store_with_progress,
    save_vector_store_with_progress
)

def main():
    """Build the vector store from scratch"""
    print("🚀 Building H5AD RAG Vector Store...")
    print("=" * 50)
    
    # Initialize progress tracker
    progress_tracker = ProgressTracker()
    
    # Configuration
    data_directory = "Converted_Scripts_Annotated"
    vector_store_path = "omicverse_vectorstore_gemini"
    
    try:
        # Step 1: Load documents
        documents = load_documents_with_progress(data_directory, progress_tracker)
        
        if not documents:
            print("❌ No documents found to process!")
            sys.exit(1)
        
        # Step 2: Split documents
        texts = split_documents_with_progress(documents, progress_tracker)
        
        # Step 3: Initialize embeddings
        embeddings = initialize_embeddings_with_progress(progress_tracker)
        
        # Step 4: Create vector store
        vector_store = create_vector_store_with_progress(texts, embeddings, progress_tracker)
        
        # Step 5: Save vector store
        save_vector_store_with_progress(vector_store, vector_store_path, progress_tracker)
        
        print("✅ Vector store built successfully!")
        print(f"📊 Processing summary:")
        print(f"   - Documents loaded: {len(documents)}")
        print(f"   - Text chunks created: {len(texts)}")
        print(f"   - Vector store saved to: {vector_store_path}")
        
    except Exception as e:
        print(f"❌ Error building vector store: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()