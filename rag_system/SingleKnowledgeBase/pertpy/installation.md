```{highlight} shell

```

# Installation

## Stable release

### PyPI

To install pertpy, run this command in your terminal:

```console
pip install pertpy
```

This is the preferred method to install pertpy, as it will always install the most recent stable release.
If you don't have [pip] installed, this [Python installation guide] can guide you through the process.

### conda-forge

Alternatively, you can install pertpy from [conda-forge]:

```console
conda install -c conda-forge pertpy
```

### Additional dependency groups

#### Differential gene expression interface

The DGE interface of pertpy requires additional dependencies that can be installed by running:

```console
pip install pertpy[de]
```

Note that edger in pertpy requires edger and rpy2 to be installed:

```R
BiocManager::install("edgeR")
```

```console
pip install rpy2
```

#### milo

milo requires either the "de" extra for the "pydeseq2" solver:

```console
pip install 'pertpy[de]'
```

or, edger, statmod, and rpy2 for the "edger" solver:

```R
BiocManager::install("edgeR")
BiocManager::install("statmod")
```

```console
pip install rpy2
```

#### tascCODA

TascCODA requires an additional set of dependencies (ete4, pyqt6, and toytree) that can be installed by running:

```console
pip install pertpy[tcoda]
```

## From sources

The sources for pertpy can be downloaded from the [Github repo].

You can either clone the public repository:

```console
$ git clone git://github.com/scverse/pertpy
```

Or download the [tarball]:

```console
$ curl -OJL https://github.com/scverse/pertpy/tarball/master
```

[github repo]: https://github.com/scverse/pertpy
[pip]: https://pip.pypa.io
[conda-forge]: https://anaconda.org/conda-forge/pertpy
[python installation guide]: http://docs.python-guide.org/en/latest/starting/installation/
[tarball]: https://github.com/scverse/pertpy/tarball/master
[Homebrew]: https://brew.sh/
