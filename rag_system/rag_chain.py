from langchain.chains import RetrievalQA
from langchain.llms import Ollama
from langchain.callbacks.manager import CallbackManager
from langchain.callbacks.streaming_stdout import StreamingStdOutCallbackHandler
from langchain.callbacks.base import BaseCallbackHandler
import time

class RetryCallbackHandler(BaseCallbackHandler):
    def on_llm_error(self, error: Exception, **kwargs) -> None:
        print(f"Error occurred: {error}")
        time.sleep(1)

def setup_rag_chain(vector_store, model_name="llama3.2:3b", k=46):  # Use a reasonable default k value
    callback_manager = CallbackManager([StreamingStdOutCallbackHandler(), RetryCallbackHandler()])

    llm = Ollama(model=model_name, callback_manager=callback_manager, temperature=0.5)
    qa_chain = RetrievalQA.from_chain_type(
        llm=llm,
        chain_type="stuff",
        retriever=vector_store.as_retriever(search_kwargs={"k": k}),
        return_source_documents=True
    )
    return qa_chain

def query_rag_system(qa_chain, query):
    result = qa_chain({"query": query})
    return {"answer": result["result"], "source_documents": result["source_documents"]}