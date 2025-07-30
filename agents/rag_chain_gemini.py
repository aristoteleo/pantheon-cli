from langchain.chains import RetrievalQA
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain.callbacks.manager import CallbackManager
from langchain.callbacks.streaming_stdout import StreamingStdOutCallbackHandler
from langchain.callbacks.base import BaseCallbackHandler
import os
import time

class RetryCallbackHandler(BaseCallbackHandler):
    def on_llm_error(self, error: Exception, **kwargs) -> None:
        print(f"Error occurred: {error}")
        time.sleep(1)

def setup_rag_chain_gemini(vector_store, model_name="gemini-2.5-pro", k=46):
    """
    Set up RAG chain using Gemini API instead of Ollama.
    Uses the same GEMINI_API_KEY that Pantheon CLI uses.
    """
    callback_manager = CallbackManager([StreamingStdOutCallbackHandler(), RetryCallbackHandler()])
    
    # Get API key from environment (same as Pantheon CLI)
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        raise ValueError("GEMINI_API_KEY environment variable not found. Please set it to use Gemini API.")
    
    # Initialize Gemini LLM
    llm = ChatGoogleGenerativeAI(
        model=model_name,
        google_api_key=api_key,
        temperature=0.3,  # Lower temperature for more consistent bioinformatics code
        callback_manager=callback_manager
    )
    
    qa_chain = RetrievalQA.from_chain_type(
        llm=llm,
        chain_type="stuff",
        retriever=vector_store.as_retriever(search_kwargs={"k": k}),
        return_source_documents=True
    )
    return qa_chain

def query_rag_system_gemini(qa_chain, query):
    """
    Query the RAG system with enhanced prompting for bioinformatics tasks.
    """
    # Enhance the query with bioinformatics context
    enhanced_query = f"""
    You are an expert bioinformatics assistant specializing in single-cell RNA-seq analysis.
    
    User Query: {query}
    
    Please provide:
    1. Clear, executable Python code using scanpy/omicverse
    2. Brief explanation of each step
    3. Best practices and parameter recommendations
    4. Expected outputs and next steps
    
    Focus on practical, working code that follows bioinformatics conventions.
    """
    
    result = qa_chain({"query": enhanced_query})
    return {
        "answer": result["result"], 
        "source_documents": result["source_documents"],
        "model_used": "gemini-api"
    }

# Backward compatibility functions
def setup_rag_chain(vector_store, model_name="gemini-2.5-pro", k=46):
    """Backward compatible function that now uses Gemini API"""
    return setup_rag_chain_gemini(vector_store, model_name, k)

def query_rag_system(qa_chain, query):
    """Backward compatible function that now uses Gemini API"""
    return query_rag_system_gemini(qa_chain, query)