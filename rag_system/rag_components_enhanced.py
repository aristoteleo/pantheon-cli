"""
Enhanced RAG components with progress tracking for H5AD analyzer integration
"""
import os
import sys
import json
from typing import List, Dict, Any, Optional, Callable
from tqdm import tqdm
from langchain_community.document_loaders import TextLoader, DirectoryLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document

# Disable tokenizers parallelism warning
os.environ["TOKENIZERS_PARALLELISM"] = "false"

class ProgressTracker:
    """进度跟踪器，用于显示嵌入和RAG处理进度"""
    
    def __init__(self, progress_callback: Optional[Callable[[str], None]] = None):
        self.progress_callback = progress_callback or print
    
    def update(self, message: str):
        """更新进度信息"""
        self.progress_callback(message)

def load_documents_with_progress(directory_path: str, progress_tracker: ProgressTracker) -> List[Document]:
    """加载文档并显示进度"""
    progress_tracker.update(f"📂 Loading documents from {directory_path}...")
    
    # Load Python files with annotated OmicVerse scripts
    loader = DirectoryLoader(
        directory_path, 
        glob="t_*_annotated.py", 
        loader_cls=TextLoader,
        show_progress=True
    )
    
    documents = loader.load()
    progress_tracker.update(f"✅ Loaded {len(documents)} documents")
    
    return documents

def split_documents_with_progress(documents: List[Document], progress_tracker: ProgressTracker) -> List[Document]:
    """分割文档并显示进度"""
    progress_tracker.update("📄 Splitting documents into chunks...")
    
    text_splitter = RecursiveCharacterTextSplitter(
        chunk_size=1000, 
        chunk_overlap=200, 
        length_function=len, 
        separators=["\n\n", "\n", " ", ""]
    )
    
    # Split documents with progress
    all_texts = []
    with tqdm(total=len(documents), desc="Splitting documents") as pbar:
        for doc in documents:
            texts = text_splitter.split_documents([doc])
            all_texts.extend(texts)
            pbar.update(1)
    
    progress_tracker.update(f"✅ Created {len(all_texts)} text chunks")
    return all_texts

def initialize_embeddings_with_progress(progress_tracker: ProgressTracker) -> HuggingFaceEmbeddings:
    """初始化嵌入模型并显示进度"""
    progress_tracker.update("🤖 Initializing embedding model (all-MiniLM-L6-v2)...")
    
    try:
        embeddings = HuggingFaceEmbeddings(
            model_name="all-MiniLM-L6-v2",
            model_kwargs={'device': 'cpu'}  # Force CPU to avoid GPU issues
            # Remove encode_kwargs to avoid parameter conflicts
        )
        progress_tracker.update("✅ Embedding model initialized successfully")
        return embeddings
    except Exception as e:
        progress_tracker.update(f"❌ Error initializing embeddings: {str(e)}")
        raise

def create_vector_store_with_progress(
    texts: List[Document], 
    embeddings: HuggingFaceEmbeddings, 
    progress_tracker: ProgressTracker
) -> FAISS:
    """创建向量存储并显示嵌入进度"""
    progress_tracker.update(f"🔄 Creating vector store from {len(texts)} text chunks...")
    progress_tracker.update("⚡ Starting embedding process (this may take a few minutes)...")
    
    try:
        # Create embeddings with progress tracking
        with tqdm(total=len(texts), desc="Creating embeddings") as pbar:
            # Process texts in batches to show progress
            batch_size = 50
            vector_store = None
            
            for i in range(0, len(texts), batch_size):
                batch = texts[i:i + batch_size]
                
                if vector_store is None:
                    # Create initial vector store
                    vector_store = FAISS.from_documents(batch, embeddings)
                else:
                    # Add to existing vector store
                    batch_store = FAISS.from_documents(batch, embeddings)
                    vector_store.merge_from(batch_store)
                
                pbar.update(len(batch))
                progress_tracker.update(f"📊 Processed {min(i + batch_size, len(texts))}/{len(texts)} chunks")
        
        progress_tracker.update("✅ Vector store created successfully")
        return vector_store
        
    except Exception as e:
        progress_tracker.update(f"❌ Error creating vector store: {str(e)}")
        raise

def load_vector_store_with_progress(
    vector_store_path: str, 
    embeddings: HuggingFaceEmbeddings, 
    progress_tracker: ProgressTracker
) -> FAISS:
    """加载向量存储并显示进度"""
    progress_tracker.update(f"📚 Loading existing vector store from {vector_store_path}...")
    
    try:
        vector_store = FAISS.load_local(
            vector_store_path, 
            embeddings, 
            allow_dangerous_deserialization=True
        )
        progress_tracker.update("✅ Vector store loaded successfully")
        return vector_store
    except Exception as e:
        progress_tracker.update(f"❌ Error loading vector store: {str(e)}")
        raise

def save_vector_store_with_progress(
    vectorstore: FAISS, 
    vector_store_path: str, 
    progress_tracker: ProgressTracker
):
    """保存向量存储并显示进度"""
    progress_tracker.update(f"💾 Saving vector store to {vector_store_path}...")
    
    try:
        # Ensure directory exists
        os.makedirs(vector_store_path, exist_ok=True)
        vectorstore.save_local(vector_store_path)
        progress_tracker.update("✅ Vector store saved successfully")
    except Exception as e:
        progress_tracker.update(f"❌ Error saving vector store: {str(e)}")
        raise

def query_vector_store_with_progress(
    vector_store: FAISS,
    query: str,
    k: int = 5,
    progress_tracker: ProgressTracker = None
) -> List[Document]:
    """查询向量存储并显示进度"""
    if progress_tracker:
        progress_tracker.update(f"🔍 Searching for relevant documents (k={k})...")
    
    try:
        # Perform similarity search
        relevant_docs = vector_store.similarity_search(query, k=k)
        
        if progress_tracker:
            progress_tracker.update(f"✅ Found {len(relevant_docs)} relevant documents")
        
        return relevant_docs
    except Exception as e:
        if progress_tracker:
            progress_tracker.update(f"❌ Error querying vector store: {str(e)}")
        raise

# Backward compatibility functions
def load_documents(directory_path: str) -> List[Document]:
    """向后兼容的文档加载函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    return load_documents_with_progress(directory_path, progress_tracker)

def split_documents(documents: List[Document]) -> List[Document]:
    """向后兼容的文档分割函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    return split_documents_with_progress(documents, progress_tracker)

def initialize_embeddings() -> HuggingFaceEmbeddings:
    """向后兼容的嵌入初始化函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    return initialize_embeddings_with_progress(progress_tracker)

def create_vector_store(texts: List[Document], embeddings: HuggingFaceEmbeddings) -> FAISS:
    """向后兼容的向量存储创建函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    return create_vector_store_with_progress(texts, embeddings, progress_tracker)

def load_vector_store(vector_store_path: str, embeddings: HuggingFaceEmbeddings) -> FAISS:
    """向后兼容的向量存储加载函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    return load_vector_store_with_progress(vector_store_path, embeddings, progress_tracker)

def save_vector_store(vectorstore: FAISS, vector_store_path: str):
    """向后兼容的向量存储保存函数"""
    progress_tracker = ProgressTracker(lambda x: print(x, file=sys.stderr))
    save_vector_store_with_progress(vectorstore, vector_store_path, progress_tracker)