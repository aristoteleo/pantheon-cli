# About Pertpy

Pertpy is an end-to-end framework for the analysis of large-scale single-cell perturbation experiments.
It provides access to harmonized perturbation datasets and metadata databases along with numerous fast and user-friendly implementations of both established and novel methods such as automatic metadata annotation or perturbation distances to efficiently analyze perturbation data.
As part of the scverse ecosystem, pertpy interoperates with existing single-cell analysis libraries and is designed to be easily extended.
If you find pertpy useful for your research, please check out {doc}`cite`.

## Design principles

Our framework is based on three key principles: `Modularity`, `Flexibility`, and `Scalability`.

### Modularity

Pertpy includes modules for analysis of single and combinatorial perturbations covering diverse types of perturbation data including genetic knockouts, drug screens, and disease states.
The framework is designed for flexibility, offering more than 100 composable and interoperable analysis functions organized in modules which further ease downstream interpretation and visualization.
These modules host fundamental building blocks for implementation and methods that share functionality and can be chained into custom pipelines.

A typical Pertpy workflow consists of several steps:

* Initial **data transformation** such as guide RNA assignment for CRISPR screens
* **Quality control** to address confounding factors and technical variation
* **Metadata annotation** against ontologies and enrichment from databases
* **Perturbation space analysis** to learn biologically interpretable embeddings
* **Downstream analysis** including differential expression, compositional analysis, and distance calculation

This modular approach yields a powerful and flexible framework as many analysis steps can be independently applied or chained together.

### Flexibility

Pertpy is purpose-built to organize, analyze, and visualize complex perturbation datasets.
It is flexible and can be applied to datasets of different assays, data types, sizes, and perturbations, thereby unifying previous data-type- or assay-specific single-problem approaches.
Designed to integrate external metadata with measured data, it enables unprecedented contextualization of results through swiftly built, experiment-specific pipelines, leading to more robust outcomes.

The inputs to a typical analysis with pertpy are unimodal scRNA-seq or multimodal perturbation readouts stored in AnnData or MuData objects.
While pertpy is primarily designed to explore perturbations such as genetic modifications, drug treatments, exposure to pathogens, and other environmental conditions, its utility extends to various other perturbation settings, including diverse disease states where experimental perturbations have not been applied.

### Scalability

Pertpy addresses a wide array of use-cases and different types of growing datasets through its sparse and memory-efficient implementations, which leverage the parallelization and GPU acceleration library Jax, and numba, thereby making them substantially faster than original implementations.
The framework can be applied to datasets ranging from thousands to millions of cells.

For example, when analyzing CRISPR screens, Pertpy's implementation of Mixscape is optimized using PyNNDescent for nearest neighbor search during the calculation of perturbation signatures.
Other methods such as scCODA and tascCODA are accelerated by replacing the Hamiltonian Monte Carlo algorithm in TensorFlow with the no-U-turn sampler from numpyro.
CINEMA-OT is optimized with ott-jax to make the implementation portable across hardware, enabling GPU acceleration.

## Why is it called "Pertpy"?

Pertpy is named for its core purpose: The analysis of **pert**urbations in **Py**thon.
The framework unifies perturbation analysis approaches across different data types and experimental designs, providing a comprehensive solution for understanding cellular responses to various stimuli.
