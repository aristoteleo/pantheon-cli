```{eval-rst}
.. currentmodule:: pertpy
```

# Tools

## Differential gene expression

Differential gene expression involves the quantitative comparison of gene expression levels between two or more groups,
such as different cell types, tissues, or conditions to discern genes that are significantly up- or downregulated in response to specific biological contexts or stimuli.
Pertpy enables differential gene expression tests through a common interface that supports complex designs.

```{eval-rst}
.. autosummary::
    :toctree: preprocessing
    :nosignatures:

    tools.PyDESeq2
    tools.EdgeR
    tools.WilcoxonTest
    tools.TTest
    tools.Statsmodels
```

## Pooled CRISPR screens

### Perturbation assignment - Mixscape

CRISPR based screens can suffer from off-target effects but also limited efficacy of the guide RNAs.
When analyzing CRISPR screen data, it is vital to know which perturbations were successful and which ones were not
to accurately determine the effect of perturbations.

[Mixscape](https://www.nature.com/articles/s41588-021-00778-2) is a pipeline that aims to determine and remove unsuccessfully perturbed cells {cite}`Papalexi2021`.
First tries to remove confounding sources of variation such as cell cycle or replicate effect by calculating a perturbation signature
Next, it determines which targeted cells were affected by the genetic perturbation (=KO) and which targeted cells were not (=NP) with the use of mixture models.
Finally, it visualizes similarities and differences across different perturbations.

See [Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens](https://www.nature.com/articles/s41588-021-00778-2) for more details on the pipeline.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Mixscape
```

Example implementation:

```python
import pertpy as pt

mdata = pt.dt.papalexi_2021()
ms = pt.tl.Mixscape()
ms.perturbation_signature(mdata["rna"], "perturbation", "NT", "replicate")
ms.mixscape(adata=mdata["rna"], control="NT", labels="gene_target", layer="X_pert")
ms.lda(adata=mdata["rna"], labels="gene_target", layer="X_pert", control="NT")
ms.plot_lda(adata=mdata["rna"], control="NT")
```

See [mixscape tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/mixscape.html).

## Compositional analysis

Compositional data analysis focuses on identifying and quantifying variations in cell type composition across
different conditions or samples to uncover biological differences driven by changes in cellular makeup.

Generally, there's two ways of approaching this question:

1. Without labeled groups using graph based approaches
2. With labeled groups using pure statistical approaches

For a more in-depth explanation we refer to the corresponding [sc-best-practices compositional chapter](https://www.sc-best-practices.org/conditions/compositional.html).

### Without labeled groups - Milo

[Milo](https://www.nature.com/articles/s41587-021-01033-z) enables the exploration of differential abundance of cell types across different biological conditions or spatial locations {cite}`Dann2022`.
It employs a neighborhood-testing approach to statistically assess variations in cell type compositions, providing insights into the microenvironmental and functional heterogeneity within and across samples.

See [Differential abundance testing on single-cell data using k-nearest neighbor graphs](https://www.nature.com/articles/s41587-021-01033-z) for details on the statistical framework.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Milo
```

Example implementation:

```python
import pertpy as pt
import scanpy as sc

adata = pt.dt.stephenson_2021_subsampled()
adata.obs["COVID_severity"] = adata.obs["Status_on_day_collection_summary"].copy()
adata.obs[["patient_id", "COVID_severity"]].drop_duplicates()
adata = adata[adata.obs["Status"] != "LPS"].copy()

milo = pt.tl.Milo()
mdata = milo.load(adata)
sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150, n_pcs=10)
milo.make_nhoods(mdata["rna"], prop=0.1)
mdata = milo.count_nhoods(mdata, sample_col="patient_id")
mdata["rna"].obs["Status"] = (
    mdata["rna"].obs["Status"].cat.reorder_categories(["Healthy", "Covid"])
)
milo.da_nhoods(mdata, design="~Status")
```

See [milo tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/milo.html).

### With labeled groups - scCODA and tascCODA

[scCODA](https://www.nature.com/articles/s41467-021-27150-6) is designed to identify differences in cell type compositions from single-cell sequencing data across conditions for labeled groups {cite}`Büttner2021`.
It employs a Bayesian hierarchical model and Dirichlet-multinomial distribution, using Markov chain Monte Carlo (MCMC) for inference, to detect significant shifts in cell type composition across conditions.

[tascCODA](https://www.frontiersin.org/articles/10.3389/fgene.2021.766405/full) extends scCODA to analyze compositional count data from single-cell sequencing studies, incorporating hierarchical tree information and experimental covariates {cite}`Ostner2021`.
By integrating spike-and-slab Lasso penalization with latent tree-based parameters, tascCODA identifies differential abundance across hierarchical levels, offering parsimonious and predictive insights into compositional changes in cell populations.

See [scCODA is a Bayesian model for compositional single-cell data analysis](https://www.nature.com/articles/s41467-021-27150-6) and [tascCODA: Bayesian Tree-Aggregated Analysis of Compositional Amplicon and Single-Cell Data](https://www.frontiersin.org/articles/10.3389/fgene.2021.766405/full) for more details.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Sccoda
    tools.Tasccoda
```

Example implementation:

```python
import pertpy as pt

haber_cells = pt.dt.haber_2017_regions()
sccoda = pt.tl.Sccoda()
sccoda_data = sccoda.load(
    haber_cells,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)
sccoda_data.mod["coda_salm"] = sccoda_data["coda"][
    sccoda_data["coda"].obs["condition"].isin(["Control", "Salmonella"])
].copy()

sccoda_data = sccoda.prepare(
    sccoda_data,
    modality_key="coda_salm",
    formula="condition",
    reference_cell_type="Goblet",
)
sccoda.run_nuts(sccoda_data, modality_key="coda_salm")
sccoda.summary(sccoda_data, modality_key="coda_salm")
sccoda.plot_effects_barplot(
    sccoda_data, modality_key="coda_salm", parameter="Final Parameter"
)
```

See [sccoda tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda.html), [extended sccoda tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda_extended.html) and [tasccoda tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/tasccoda.html).

## Multicellular and gene programs

Multicellular programs are organized interactions and coordinated activities among different cell types within a tissue,
forming complex functional units that drive tissue-specific functions, responses to environmental changes, and pathological states.
These programs enable a higher level of biological organization by integrating signaling pathways, gene expression,
and cellular behaviors across the cellular community to maintain homeostasis and execute collective responses.

### Multicellular programs - DIALOGUE

[DIALOGUE](https://www.nature.com/articles/s41587-022-01288-0) identifies latent multicellular programs by mapping the data into
a feature space where the cell type specific representations are correlated across different samples and environments {cite}`JerbyArnon2022`.
Next, DIALOGUE employs multi-level hierarchical modeling to identify genes that comprise the latent features.

This is a **work in progress (!)** Python implementation of DIALOGUE for the discovery of multicellular programs.

See [DIALOGUE maps multicellular programs in tissue from single-cell or spatial transcriptomics data](https://www.nature.com/articles/s41587-022-01288-0) for more details on the methodology.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Dialogue
```

Example implementation:

```python
import pertpy as pt
import scanpy as sc

adata = pt.dt.dialogue_example()
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)


dl = pt.tl.Dialogue(
    sample_id="clinical.status",
    celltype_key="cell.subtypes",
    n_counts_key="nCount_RNA",
    n_mpcs=3,
)
adata, mcps, ws, ct_subs = dl.calculate_multifactor_PMD(adata, normalize=True)
all_results, new_mcps = dl.multilevel_modeling(
    ct_subs=ct_subs,
    mcp_scores=mcps,
    ws_dict=ws,
    confounder="gender",
)
```

See [DIALOGUE tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/dialogue.html).

### Enrichment

Enrichment tests for single-cell data assess whether specific biological pathways or gene sets are overrepresented in the expression profiles of individual cells,
aiding in the identification of functional characteristics and cellular states.
While pathway enrichment is a well-studied and commonly applied approach in single-cell RNA-seq, other data sources such as genes targeted by drugs can also be enriched.
Drug2cell performs such enrichment tests and is available in pertpy {cite}`Kanemaru2023`.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Enrichment
```

Example implementation:

```python
import pertpy as pt
import scanpy as sc

adata = sc.datasets.pbmc3k_processed()

pt_enricher = pt.tl.Enrichment()
pt_enricher.score(adata)
```

See [enrichment tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/enrichment.html).

## Distances and permutation tests

In settings where many perturbations are applied, it is often times unclear which perturbations had a strong effect and should be investigated further.
Differential gene expression poses one option to get candidate genes and p-values.
Determining statistical distances between the perturbations and applying a permutation test is another option {cite}`Peidli2024`.

For more details, we refer to [scPerturb: harmonized single-cell perturbation data](https://www.nature.com/articles/s41592-023-02144-y).

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Distance
    tools.DistanceTest
```

Example implementation:

```python
import pertpy as pt

adata = pt.dt.distance_example()

# Pairwise distances
distance = pt.tl.Distance(metric="edistance", obsm_key="X_pca")
pairwise_edistance = distance.pairwise(adata, groupby="perturbation")

# E-test (Permutation test using E-distance)
etest = pt.tl.PermutationTest(
    metric="edistance", obsm_key="X_pca", correction="holm-sidak"
)
tab = etest(adata, groupby="perturbation", contrast="control")
```

See [distance tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/distances.html)
and [distance tests tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/distance_tests.html).

## Response prediction

Response prediction describes computational models that predict how individual cells or cell populations will respond to
specific treatments, conditions, or stimuli based on their gene expression profiles, enabling insights into cellular behaviors and potential therapeutic strategies.
Such approaches can also order perturbations by their effect on groups of cells.

### Rank perturbations - Augur

[Augur](https://doi.org/10.1038/s41587-020-0605-1) aims to rank or prioritize cell types according to their response to experimental perturbations {cite}`Skinnider2021`.
Cells that respond strongly to perturbations are more easily distinguishable as treated or control in molecular space.
Augur quantifies this by training a classifier to predict experimental labels within each cell type across cross-validation runs.
Cell types are ranked by model accuracy—using AUC for categorical labels and concordance correlation for continuous ones—as a proxy for perturbation response.

For more details we refer to [Cell type prioritization in single-cell data](https://doi.org/10.1038/s41587-020-0605-1).

Example implementation:

```python
import pertpy as pt

adata = pt.dt.sc_sim_augur()
ag = pt.tl.Augur(estimator="random_forest_classifier")
adata = ag.load(adata)
adata, results = ag.predict(adata)

# metrics for each cell type
results["summary_metrics"]
```

See [augur tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/augur.html).

```{eval-rst}
.. autosummary::
    :toctree: tools
    :nosignatures:

    tools.Augur
```

### Gene expression prediction with scGen

scGen is a deep generative model that leverages autoencoders and adversarial training to integrate single-cell RNA sequencing data from different conditions or tissues,
enabling the generation of synthetic single-cell data for cross-condition analysis and predicting cell-type-specific responses to perturbations {cite}`Lotfollahi2019`.
See [scGen predicts single-cell perturbation responses](https://www.nature.com/articles/s41592-019-0494-8) for more details.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Scgen
```

Example implementation:

```python
import pertpy as pt

train = pt.dt.kang_2018()

train_new = train[
    ~((train.obs["cell_type"] == "CD4T") & (train.obs["condition"] == "stimulated"))
]
train_new = train_new.copy()

pt.tl.Scgen.setup_anndata(train_new, batch_key="condition", labels_key="cell_type")
scgen = pt.tl.Scgen(train_new)
scgen.train(max_epochs=100, batch_size=32)

pred, delta = scgen.predict(
    ctrl_key="control", stim_key="stimulated", celltype_to_predict="CD4T"
)
pred.obs["condition"] = "pred"
```

See [scgen tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/scgen_perturbation_prediction.html).

### Causal perturbation analysis with CINEMA-OT

CINEMA-OT is a causal framework for perturbation effect analysis to identify individual treatment effects and synergy at the single cell level {cite}`Dong2023`.
CINEMA-OT separates confounding sources of variation from perturbation effects to obtain an optimal transport matching that reflects counterfactual cell pairs.
These cell pairs represent causal perturbation responses permitting a number of novel analyses, such as individual treatment effect analysis, response clustering, attribution analysis, and synergy analysis.

See [Causal identification of single-cell experimental perturbation effects with CINEMA-OT](https://www.biorxiv.org/content/10.1101/2022.07.31.502173v3.abstract) for more details.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.Cinemaot
```

Example implementation:

```python
import pertpy as pt

adata = pt.dt.cinemaot_example()

model = pt.tl.Cinemaot()
de = model.causaleffect(
    adata,
    pert_key="perturbation",
    control="No stimulation",
    return_matching=True,
    thres=0.5,
    smoothness=1e-5,
    eps=1e-3,
    solver="Sinkhorn",
    preweight_label="cell_type0528",
)
```

See [CINEMA-OT tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/cinemaot.html).

## Perturbation space

Perturbation spaces depart from the individualistic perspective of cells and instead organizes cells into cohesive ensembles.
This specialized space enables comprehending the collective impact of perturbations on cells.
Pertpy offers various modules for calculating and evaluating perturbation spaces that are either based on summary statistics or clusters.

```{eval-rst}
.. autosummary::
    :toctree: tools

    tools.MLPClassifierSpace
    tools.LRClassifierSpace
    tools.CentroidSpace
    tools.DBSCANSpace
    tools.KMeansSpace
    tools.PseudobulkSpace
```

Example implementation:

```python
import pertpy as pt

mdata = pt.dt.papalexi_2021()
ps = pt.tl.PseudobulkSpace()
ps_adata = ps.compute(
    mdata["rna"],
    target_col="gene_target",
    groups_col="gene_target",
    mode="mean",
)
```

See [perturbation space tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/perturbation_space.html).
