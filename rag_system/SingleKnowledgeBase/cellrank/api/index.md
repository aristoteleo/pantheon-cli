
# API
Import CellRank as::

    import cellrank as cr

CellRank has a modular API, organized around kernels and estimators:

- `cellrank.kernels` compute cell-cell transition matrices using various input data modalities,
  including RNA velocity, any pseudotime, a developmental potential, experimental time points, and more.
- `cellrank.estimators` use the cell-cell transition matrix to derive insights about cellular dynamics,
  for example, they compute initial and terminal states, fate probabilities, and driver genes. Our recommended
  estimator is the `~cellrank.estimators.GPCCA` estimator.
- In addition, there are `cellrank.models` for gene trend fitting, `cellrank.pl`
  for visualization, and `cellrank.datasets` that help you getting started with CellRank.

    :caption: API
    :maxdepth: 1

    kernels
    estimators
    models
    plotting
    datasets
    developer
