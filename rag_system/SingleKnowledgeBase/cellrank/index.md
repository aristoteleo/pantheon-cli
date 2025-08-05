|PyPI| |Downloads| |CI| |Docs| |Codecov| |Discourse|

# CellRank 2: Unified fate mapping in multiview single-cell data
    :width: 600px
    :align: center
    :class: only-light

    :width: 600px
    :align: center
    :class: only-dark

**CellRank** [lange:22,weiler:24] is a modular framework to study cellular dynamics based on Markov state modeling of
multi-view single-cell data. See [about CellRank <about/index>](about CellRank <about/index>.md) to learn more and [our citation guide <about/cite>](our citation guide <about/cite>.md) for guidance on
citing our work correctly. Two peer-reviewed publications accompany our software:

- [CellRank for directed single-cell fate mapping](https://doi.org/10.1038/s41592-021-01346-6)
- [CellRank 2: unified fate mapping in multiview single-cell data](https://doi.org/10.1038/s41592-024-02303-9)

    Please refer to [our citation guide <about/cite>](our citation guide <about/cite>.md) to cite our software correctly.

CellRank scales to large cell numbers, is fully compatible with the `scverse`_ ecosystem, and is easy to use. In the
backend, it is powered by the `pyGPCCA package <https://github.com/msmdev/pyGPCCA>`_ [reuter:19,reuter:22]. Feel
free to open an `issue`_ if you encounter a bug, need our help or just want to make a comment/suggestion.

    If you're moving from CellRank 1 to CellRank 2, check out [../about/version2](../about/version2.md).

## CellRank's key applications
- Estimate differentiation direction based on a varied number of biological priors, including
  [pseudotime <notebooks/tutorials/kernels/300_pseudotime>](pseudotime <notebooks/tutorials/kernels/300_pseudotime>.md),
  [developmental potential <notebooks/tutorials/kernels/400_cytotrace>](developmental potential <notebooks/tutorials/kernels/400_cytotrace>.md),
  [RNA velocity <notebooks/tutorials/kernels/200_rna_velocity>](RNA velocity <notebooks/tutorials/kernels/200_rna_velocity>.md),
  [experimental time points <notebooks/tutorials/kernels/500_real_time>](experimental time points <notebooks/tutorials/kernels/500_real_time>.md), and `more <cellrank.kernels>`.
- Compute initial, terminal and intermediate [macrostates <notebooks/tutorials/estimators/600_initial_terminal>](macrostates <notebooks/tutorials/estimators/600_initial_terminal>.md)
  [reuter:19,reuter:22].
- Infer [fate probabilities and driver genes <notebooks/tutorials/estimators/700_fate_probabilities>](fate probabilities and driver genes <notebooks/tutorials/estimators/700_fate_probabilities>.md).
- Visualize and cluster [gene expression trends <notebooks/tutorials/estimators/800_gene_trends>](gene expression trends <notebooks/tutorials/estimators/800_gene_trends>.md).
- ... and much more, check out our [API <api/index>](API <api/index>.md).

## Getting started with CellRank
We have [notebooks/tutorials/index](notebooks/tutorials/index.md) to help you getting started. To see CellRank in action, explore our
manuscripts [lange:22,weiler:24] in Nature Methods.

## Contributing
We actively encourage any contribution! To get started, please check out the [contributing](contributing.md).

    :caption: General
    :maxdepth: 3
    :hidden:

    installation
    api/index
    notebooks/tutorials/index
    release_notes
    contributing
    references

    :caption: About
    :maxdepth: 3
    :hidden:

    about/index
    about/version2
    about/cite
    about/team
    GitHub <https://github.com/theislab/cellrank>
    Discourse <https://discourse.scverse.org/c/ecosytem/cellrank/40>

    :target: https://pypi.org/project/cellrank
    :alt: PyPI

    :target: https://pepy.tech/project/cellrank
    :alt: Downloads

    :target: https://discourse.scverse.org/c/ecosystem/cellrank/
    :alt: Discourse

    :target: https://github.com/theislab/cellrank/actions
    :alt: CI

    :target: https://cellrank.readthedocs.io/
    :alt: Documentation

    :target: https://codecov.io/gh/theislab/cellrank
    :alt: Coverage

