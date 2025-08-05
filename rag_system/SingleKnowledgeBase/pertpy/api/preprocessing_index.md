```{eval-rst}
.. currentmodule:: pertpy
```

# Preprocessing

## Guide Assignment

Guide assignment is essential for quality control in single-cell Perturb-seq data, ensuring accurate mapping of guide RNAs to cells for reliable interpretation of gene perturbation effects.
pertpy provides a simple function to assign guides based on thresholds and a Gaussian mixture model {cite}`Replogle2022`.

```{eval-rst}
.. autosummary::
    :toctree: preprocessing
    :nosignatures:

    preprocessing.GuideAssignment
```

Example implementation:

```python
import pertpy as pt
import scanpy as sc

mdata = pt.dt.papalexi_2021()
gdo = mdata.mod["gdo"]
gdo.layers["counts"] = gdo.X.copy()
sc.pp.log1p(gdo)

ga = pt.pp.GuideAssignment()
ga.assign_by_threshold(gdo, 5, layer="counts", output_layer="assigned_guides")

ga.plot_heatmap(gdo, layer="assigned_guides")
```

See [guide assignment tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/guide_rna_assignment.html).
