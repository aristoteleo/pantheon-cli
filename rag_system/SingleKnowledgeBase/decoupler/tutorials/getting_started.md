# Getting Started with Decoupler

## Overview

`decoupler` is a python package containing different enrichment statistical methods to extract biologically driven scores from omics data within a unified framework. It is part of the scverse® ecosystem and provides a faster and memory efficient implementation for various enrichment analysis methods.

## Installation

You need to have Python 3.10 or newer installed on your system.

### Installation Options

1. **Minimal installation** from PyPI:
```bash
pip install decoupler
```

2. **Full installation** with extra dependencies:
```bash
pip install decoupler[full]
```

3. **Conda installation** (note the `-py` suffix):
```bash
mamba create -n=dcp conda-forge::decoupler-py
```

4. **Development version**:
```bash
pip install git+https://github.com/scverse/decoupler.git@main
```

## Basic Usage

Import decoupler as:

```python
import decoupler as dc
```

## Key Features

Decoupler provides multiple modules for different aspects of enrichment analysis:

- **Methods (`dc.mt`)**: Core enrichment statistical methods
- **Benchmarking (`dc.bm`)**: Tools for comparing method performance
- **Data structures (`dc.ds`)**: Data handling utilities
- **Omnipath (`dc.op`)**: Access to biological databases
- **Preprocessing (`dc.pp`)**: Data preprocessing functions
- **Tools (`dc.tl`)**: Additional analysis tools
- **Plotting (`dc.pl`)**: Visualization functions

## Next Steps

- Explore the [API documentation](../api/index.md) for detailed function references
- Check the [methods overview](methods_overview.md) for available enrichment methods
- Review the [changelog](../changelog.md) for recent updates