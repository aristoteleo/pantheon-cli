
# Kernels
Kernels compute cell-cell transition probabilities based on various input data modalities, including molecular
similarity, RNA velocity [manno:18], experimental time points, and many more. They come with methods
for high-level, qualitative visualization, including vector field or random walk plots. For quantitative analysis of
kernel-computed transition matrices, we recommend taking look at the [estimators <estimators>](estimators <estimators>.md).

We don't use the term "kernel" in the way it is used in
`mathematics <https://en.wikipedia.org/wiki/Kernel_(algebra)>`_, but rather colloquially, to refer to a class that
takes in multi-view single-cell data and outputs a cell-cell transition matrix.


    :toctree: _autosummary/kernels

    kernels.VelocityKernel
    kernels.ConnectivityKernel
    kernels.PseudotimeKernel
    kernels.CytoTRACEKernel
    kernels.RealTimeKernel
    kernels.PrecomputedKernel
