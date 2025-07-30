#!/usr/bin/env python3
"""
Test script to verify Gemini API integration works for OV_Agent RAG system.
This tests the integration without requiring all langchain dependencies.
"""

import os
import sys

def test_gemini_api_key():
    """Test if GEMINI_API_KEY is available"""
    api_key = os.getenv("GEMINI_API_KEY")
    if api_key:
        print("✅ GEMINI_API_KEY found in environment")
        print(f"   Key length: {len(api_key)} characters")
        return True
    else:
        print("❌ GEMINI_API_KEY not found in environment")
        print("   Please set: export GEMINI_API_KEY='your-api-key-here'")
        return False

def test_data_directory():
    """Test if the annotated scripts directory exists"""
    data_dir = "Converted_Scripts_Annotated"
    if os.path.exists(data_dir):
        files = [f for f in os.listdir(data_dir) if f.startswith('t_') and f.endswith('_annotated.py')]
        print(f"✅ Data directory found: {data_dir}")
        print(f"   Contains {len(files)} annotated Python scripts")
        print(f"   Sample files: {files[:3]}...")
        return True
    else:
        print(f"❌ Data directory not found: {data_dir}")
        return False

def test_required_packages():
    """Test if key packages can be imported"""
    required_packages = [
        "langchain",
        "langchain_google_genai", 
        "sentence_transformers",
        "numpy",
    ]
    
    available = {}
    for package in required_packages:
        try:
            if package == "faiss":
                __import__("faiss")
            else:
                __import__(package.replace("-", "_"))
            available[package] = True
            print(f"✅ {package} available")
        except ImportError:
            available[package] = False
            print(f"❌ {package} not available")
    
    return all(available.values())

def test_pantheon_cli_integration():
    """Test integration with Pantheon CLI structure"""
    pantheon_paths = [
        "../packages/cli/src/config/config.ts",
        "../package.json",
        "../bundle/pantheon.js"
    ]
    
    for path in pantheon_paths:
        if os.path.exists(path):
            print(f"✅ Pantheon CLI file found: {path}")
        else:
            print(f"⚠️  Pantheon CLI file not found: {path}")

def main():
    print("🧪 Testing Gemini API Integration for OV_Agent RAG System")
    print("=" * 60)
    
    tests = [
        ("API Key", test_gemini_api_key),
        ("Data Directory", test_data_directory),
        ("Required Packages", test_required_packages),
        ("Pantheon CLI Integration", test_pantheon_cli_integration)
    ]
    
    results = {}
    for test_name, test_func in tests:
        print(f"\n📋 Testing {test_name}...")
        results[test_name] = test_func()
    
    print("\n" + "=" * 60)
    print("📊 Test Results Summary:")
    
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"   {test_name}: {status}")
    
    if results["API Key"] and results["Data Directory"]:
        print("\n🎉 Core requirements met! Gemini RAG integration should work.")
        print("\n💡 To install missing packages:")
        print("   pip install -r requirements.txt")
        print("\n🚀 To test the system:")
        print("   python main_gemini.py 'How do I preprocess single-cell data?'")
    else:
        print("\n⚠️  Core requirements not met. Please address the failed tests above.")

if __name__ == "__main__":
    main()