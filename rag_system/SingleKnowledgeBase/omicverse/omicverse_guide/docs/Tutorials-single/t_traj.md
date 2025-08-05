# t_traj
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_traj.ipynb*

# Basic Trajectory Inference

Diffusion maps were introduced by Ronald Coifman and Stephane Lafon, and the underlying idea is to assume that the data are samples from a diffusion process.

Palantir is an algorithm to align cells along differentiation trajectories. Palantir models differentiation as a stochastic process where stem cells differentiate to terminally differentiated cells by a series of steps through a low dimensional phenotypic manifold. Palantir effectively captures the continuity in cell states and the stochasticity in cell fate determination.

Note that both methods require the input of cells in their initial state, and we will introduce other methods that do not require the input of artificial information, such as pyVIA, in subsequent analyses.


## Preprocess data

As an example, we apply differential kinetic analysis to dentate gyrus neurogenesis, which comprises multiple heterogeneous subpopulations.

```pythonimport scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import omicverse as ov
ov.plot_set()```
*Output:*
```
   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

Version: 1.6.3, Tutorials: https://omicverse.readthedocs.io/
```

```pythonimport scvelo as scv
adata=scv.datasets.dentategyrus()
adata```
*Output:*
```AnnData object with n_obs × n_vars = 2930 × 13913
    obs: 'clusters', 'age(days)', 'clusters_enlarged'
    uns: 'clusters_colors'
    obsm: 'X_umap'
    layers: 'ambiguous', 'spliced', 'unspliced'```

```pythonadata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=3000,)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
ov.pp.scale(adata)
ov.pp.pca(adata,layer='scaled',n_pcs=50)```
*Output:*
```Begin robust gene identification
After filtration, 13264/13913 genes are kept. Among 13264 genes, 13189 genes are robust.
End of robust gene identification.
Begin size normalization: shiftlog and HVGs selection pearson
normalizing counts per cell. The following highly-expressed genes are not considered during normalization factor computation:
['Hba-a1', 'Malat1', 'Ptgds', 'Hbb-bt']
    finished (0:00:00)
extracting highly variable genes
--> added
    'highly_variable', boolean vector (adata.var)
    'highly_variable_rank', float vector (adata.var)
    'highly_variable_nbatches', int vector (adata.var)
    'highly_variable_intersection', boolean vector (adata.var)
    'means', float vector (adata.var)
    'variances', float vector (adata.var)
    'residual_variances', float vector (adata.var)
Time to analyze data in cpu: 1.5488979816436768 seconds.
End of size normalization: shiftlog and HVGs selection pearson
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
```

Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells. In our experience, often a rough estimate of the number of PCs does fine.

```pythonov.utils.plot_pca_variance_ratio(adata)```
*Output:*
```<Figure size 320x320 with 1 Axes>```

## Trajectory inference with diffusion map

Here, we used `ov.single.TrajInfer` to construct a Trajectory Inference object.

```pythonTraj=ov.single.TrajInfer(adata,basis='X_umap',groupby='clusters',
                         use_rep='scaled|original|X_pca',n_comps=50,)
Traj.set_origin_cells('nIPC')```

```pythonTraj.inference(method='diffusion_map')```
*Output:*
```computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:04)
computing Diffusion Maps using n_comps=15(=n_dcs)
computing transitions
    finished (0:00:00)
    eigenvalues of transition matrix
    [1.         0.99972177 0.9991691  0.9990498  0.99902    0.99856234
     0.9955029  0.9947408  0.99265915 0.9906794  0.9799195  0.9787254
     0.97795963 0.9760296  0.969529  ]
    finished: added
    'X_diffmap', diffmap coordinates (adata.obsm)
    'diffmap_evals', eigenvalues of transition matrix (adata.uns) (0:00:00)
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:00)
drawing single-cell graph using layout 'fa'
    finished: added
    'X_draw_graph_fa', graph_drawing coordinates (adata.obsm) (0:00:10)
computing Diffusion Pseudotime using n_dcs=10
    finished: added
    'dpt_pseudotime', the pseudotime (adata.obs) (0:00:00)
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:00)
```

```pythonov.utils.embedding(adata,basis='X_umap',
                   color=['clusters','dpt_pseudotime'],
                   frameon='small',cmap='Reds')```
*Output:*
```<Figure size 772.8x320 with 3 Axes>```

PAGA graph abstraction has benchmarked as top-performing method for trajectory inference. It provides a graph-like map of the data topology with weighted edges corresponding to the connectivity between two clusters. 

Here, PAGA is extended by neighbor directionality.

```pythonov.utils.cal_paga(adata,use_time_prior='dpt_pseudotime',vkey='paga',
                 groups='clusters')```
*Output:*
```running PAGA using priors: ['dpt_pseudotime']
    finished
added
    'paga/connectivities', connectivities adjacency (adata.uns)
    'paga/connectivities_tree', connectivities subtree (adata.uns)
    'paga/transitions_confidence', velocity transitions (adata.uns)
```

```pythonov.utils.plot_paga(adata,basis='umap', size=50, alpha=.1,title='PAGA LTNN-graph',
            min_edge_width=2, node_size_scale=1.5,show=False,legend_loc=False)```
*Output:*
```<AxesSubplot: title={'center': 'PAGA LTNN-graph'}>```
```<Figure size 320x320 with 1 Axes>```

## Trajectory inference with Slingshot

Provides functions for inferring continuous, branching lineage structures in low-dimensional data. Slingshot was designed to model developmental trajectories in single-cell RNA sequencing data and serve as a component in an analysis pipeline after dimensionality reduction and clustering. It is flexible enough to handle arbitrarily many branching events and allows for the incorporation of prior knowledge through supervised graph construction.

```pythonTraj=ov.single.TrajInfer(adata,basis='X_umap',groupby='clusters',
                         use_rep='scaled|original|X_pca',n_comps=50)
Traj.set_origin_cells('nIPC')
#Traj.set_terminal_cells(["Granule mature","OL","Astrocytes"])```

If you only need the proposed timing and not the lineage of the process, then you can leave the debug_axes parameter unset.

```pythonTraj.inference(method='slingshot',num_epochs=1)```

else, you can set `debug_axes` to visualize the lineage

```pythonfig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
Traj.inference(method='slingshot',num_epochs=1,debug_axes=axes)```
*Output:*
```Lineages: [Lineage[13, 9, 5, 2, 7, 3], Lineage[13, 9, 5, 2, 8], Lineage[13, 9, 5, 4, 1], Lineage[13, 9, 5, 6, 11, 10], Lineage[13, 12, 0]]
```
```  0%|          | 0/1 [00:00<?, ?it/s]```
```Reversing from leaf to root
Averaging branch @2 with lineages: [0, 1] [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6de380>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7feda8b98850>]
Averaging branch @5 with lineages: [0, 1, 2, 3] [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda115300>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6df2b0>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6deef0>]
Averaging branch @13 with lineages: [0, 1, 2, 3, 4] [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6dec20>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6de8f0>]
Shrinking branch @13 with curves: [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6dec20>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6de8f0>]
Shrinking branch @5 with curves: [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda115300>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6df2b0>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6deef0>]
Shrinking branch @2 with curves: [<pcurvepy2.pcurve.PrincipalCurve object at 0x7fedda6de380>, <pcurvepy2.pcurve.PrincipalCurve object at 0x7feda8b98850>]
```
```<Figure size 640x640 with 4 Axes>```

```pythonov.utils.embedding(adata,basis='X_umap',
                   color=['clusters','slingshot_pseudotime'],
                   frameon='small',cmap='Reds')```
*Output:*
```<Figure size 772.8x320 with 3 Axes>```

```pythonsc.pp.neighbors(adata,use_rep='scaled|original|X_pca')
ov.utils.cal_paga(adata,use_time_prior='slingshot_pseudotime',vkey='paga',
                 groups='clusters')```
*Output:*
```computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:00)
running PAGA using priors: ['slingshot_pseudotime']
    finished
added
    'paga/connectivities', connectivities adjacency (adata.uns)
    'paga/connectivities_tree', connectivities subtree (adata.uns)
    'paga/transitions_confidence', velocity transitions (adata.uns)
```

```pythonov.utils.plot_paga(adata,basis='umap', size=50, alpha=.1,title='PAGA Slingshot-graph',
            min_edge_width=2, node_size_scale=1.5,show=False,legend_loc=False)```
*Output:*
```<AxesSubplot: title={'center': 'PAGA Slingshot-graph'}>```
```<Figure size 320x320 with 1 Axes>```

## Trajectory inference with Palantir

Palantir can be run by specifying an approxiate early cell.

Palantir can automatically determine the terminal states as well. In this dataset, we know the terminal states and we will set them using the terminal_states parameter

Here, we used `ov.single.TrajInfer` to construct a Trajectory Inference object.

```pythonTraj=ov.single.TrajInfer(adata,basis='X_umap',groupby='clusters',
                         use_rep='scaled|original|X_pca',n_comps=50)
Traj.set_origin_cells('nIPC')
Traj.set_terminal_cells(["Granule mature","OL","Astrocytes"])```

```pythonTraj.inference(method='palantir',num_waypoints=500)```
*Output:*
```Time for shortest paths: 0.21468504269917807 minutes
Iteratively refining the pseudotime...
Correlation at iteration 1: 0.9998
Correlation at iteration 2: 0.9999
Entropy and branch probabilities...
Markov chain construction...
Computing fundamental matrix and absorption probabilities...
Project results to all cells...
```
```<omicverse.external.palantir.presults.PResults at 0x7fedd8267e50>```

Palantir results can be visualized on the tSNE or UMAP using the plot_palantir_results function

```pythonTraj.palantir_plot_pseudotime(embedding_basis='X_umap',cmap='RdBu_r',s=3)```
*Output:*
```<Figure size 960x480 with 10 Axes>```

Once the cells are selected, it's often helpful to visualize the selection on the pseudotime trajectory to ensure we've isolated the correct cells for our specific trend. We can do this using the plot_branch_selection function:

```pythonTraj.palantir_cal_branch(eps=0)```
*Output:*
```<Figure size 1200x1200 with 6 Axes>```

```pythonov.external.palantir.plot.plot_trajectory(adata, "Granule mature",
                                cell_color="palantir_entropy",
                                n_arrows=10,
                                color="red",
                                scanpy_kwargs=dict(cmap="RdBu_r"),
                                )```
*Output:*
```[2024-06-28 14:34:58,879] [INFO    ] Using sparse Gaussian Process since n_landmarks (50) < n_samples (197) and rank = 1.0.
[2024-06-28 14:34:58,880] [INFO    ] Using covariance function Matern52(ls=1.8777667989311575).
[2024-06-28 14:34:58,930] [INFO    ] Computing 50 landmarks with k-means clustering.
```
```<AxesSubplot: title={'center': 'Branch: Granule mature'}, xlabel='UMAP1', ylabel='UMAP2'>```
```<Figure size 400x400 with 2 Axes>```

Palantir uses Mellon Function Estimator to determine the gene expression trends along different lineages. The marker trends can be determined using the following snippet. This computes the trends for all lineages. A subset of lineages can be used using the lineages parameter.

```pythongene_trends = Traj.palantir_cal_gene_trends(
    layers="MAGIC_imputed_data",
)```
*Output:*
```Granule mature
[2024-06-28 14:35:04,131] [INFO    ] Using non-sparse Gaussian Process since n_landmarks (500) >= n_samples (197) and rank = 1.0.
[2024-06-28 14:35:04,132] [INFO    ] Using covariance function Matern52(ls=1.0).
Astrocytes
[2024-06-28 14:35:05,461] [INFO    ] Using non-sparse Gaussian Process since n_landmarks (500) >= n_samples (177) and rank = 1.0.
[2024-06-28 14:35:05,462] [INFO    ] Using covariance function Matern52(ls=1.0).
OL
[2024-06-28 14:35:06,012] [INFO    ] Using non-sparse Gaussian Process since n_landmarks (500) >= n_samples (50) and rank = 1.0.
[2024-06-28 14:35:06,013] [INFO    ] Using covariance function Matern52(ls=1.0).
```

```pythongenes = ['Cdca3','Rasl10a','Mog','Aqp4']
Traj.palantir_plot_gene_trends(genes)
plt.show()```
*Output:*
```<Figure size 560x960 with 4 Axes>```

We can also use paga to visualize the cell stages

```pythonov.utils.cal_paga(adata,use_time_prior='palantir_pseudotime',vkey='paga',
                 groups='clusters')```
*Output:*
```running PAGA using priors: ['palantir_pseudotime']
    finished
added
    'paga/connectivities', connectivities adjacency (adata.uns)
    'paga/connectivities_tree', connectivities subtree (adata.uns)
    'paga/transitions_confidence', velocity transitions (adata.uns)
```

```pythonov.utils.plot_paga(adata,basis='umap', size=50, alpha=.1,title='PAGA LTNN-graph',
            min_edge_width=2, node_size_scale=1.5,show=False,legend_loc=False)```
*Output:*
```<AxesSubplot: title={'center': 'PAGA LTNN-graph'}>```
```<Figure size 320x320 with 1 Axes>```

