import os
from langchain.document_loaders import TextLoader, DirectoryLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.embeddings import HuggingFaceEmbeddings
from langchain.vectorstores import FAISS

# Disable tokenizers parallelism warning
os.environ["TOKENIZERS_PARALLELISM"] = "false"

def load_documents(directory_path):
    # Load Python files with annotated OmicVerse scripts
    loader = DirectoryLoader(directory_path, glob="t_*_annotated.py", loader_cls=TextLoader)
    return loader.load()

def split_documents(documents):
    text_splitter = RecursiveCharacterTextSplitter(
        chunk_size=1000, chunk_overlap=200, length_function=len, separators=["\n\n", "\n", " ", ""]
    )
    return text_splitter.split_documents(documents)

def initialize_embeddings():
    return HuggingFaceEmbeddings(model_name="sentence-transformers/all-MiniLM-L6-v2")

def create_vector_store(texts, embeddings):
    return FAISS.from_documents(texts, embeddings)

def load_vector_store(vector_store_path, embeddings):
    return FAISS.load_local(vector_store_path, embeddings, allow_dangerous_deserialization=True)


def save_vector_store(vectorstore, vector_store_path):
    vectorstore.save_local(vector_store_path)