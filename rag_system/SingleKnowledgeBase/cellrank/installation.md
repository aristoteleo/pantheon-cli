# Installation
`cellrank` requires Python version >= 3.8 to run. We recommend using Mamba_ to manage Python environments.
If you encounter any problems, feel free to open an issue_.

## Conda
CellRank can also be installed via `conda-forge <https://anaconda.org/conda-forge/cellrank>`_ as::

    conda install -c conda-forge cellrank

This installation method is preferred because it also contains PETSc_ and SLEPc_,
libraries for large-scale linear algebra problems CellRank relies on (e.g., when computing macrostates or
fate probabilities).

## PyPI
CellRank is also available on PyPI::

    pip install cellrank

## Development Version
To stay up-to-date with the newest version, run::

    pip install git+https://github.com/theislab/cellrank.git@main

