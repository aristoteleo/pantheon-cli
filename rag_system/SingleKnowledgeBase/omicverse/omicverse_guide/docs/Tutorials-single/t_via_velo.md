# t_via_velo
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_via_velo.ipynb*

# Trajectory Inference with VIA and scVelo

When scRNA-velocity is available, it can be used to guide the trajectory inference and automate initial state prediction. However, because RNA velocitycan be misguided by(Bergen 2021) boosts in expression, variable transcription rates and data capture scope limited to steady-state populations only, users might find it useful to adjust the level of influence the RNA-velocity data should exercise on the inferred TI.

Paper: [Generalized and scalable trajectory inference in single-cell omics data with VIA](https://www.nature.com/articles/s41467-021-25773-3)

Code: https://github.com/ShobiStassen/VIA

Colab_Reproducibility：https://colab.research.google.com/drive/1MtGr3e9uUb_BWOzKlcbOTiCYsZpljEyF?usp=sharing

```pythonimport omicverse as ov
import scanpy as sc
import scvelo as scv
import cellrank as cr
ov.utils.ov_plot_set()
```

## Data loading and preprocessing

We use a familiar endocrine-genesis dataset (Bastidas-Ponce et al. (2019).) to demonstrate initial state prediction at the EP Ngn3 low cells and automatic captures of the 4 differentiated islets (alpha, beta, delta and epsilon). As mentioned, it us useful to control the level of influence of RNA-velocity relative to gene-gene distance and this is done using the velo_weight parameter.

```pythonadata = cr.datasets.pancreas()
adata```
*Output:*
```AnnData object with n_obs × n_vars = 2531 × 27998
    obs: 'day', 'proliferation', 'G2M_score', 'S_score', 'phase', 'clusters_coarse', 'clusters', 'clusters_fine', 'louvain_Alpha', 'louvain_Beta', 'palantir_pseudotime'
    var: 'highly_variable_genes'
    uns: 'clusters_colors', 'clusters_fine_colors', 'day_colors', 'louvain_Alpha_colors', 'louvain_Beta_colors', 'neighbors', 'pca'
    obsm: 'X_pca', 'X_umap'
    layers: 'spliced', 'unspliced'
    obsp: 'connectivities', 'distances'```

```pythonn_pcs = 30
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)
sc.tl.pca(adata, n_comps = n_pcs)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
scv.tl.velocity(adata, mode='stochastic') # good results acheived with mode = 'stochastic' too
```
*Output:*
```Filtered out 22024 genes that are detected 20 counts (shared).
Normalized count data: X, spliced, unspliced.
Extracted 5000 highly variable genes.
Logarithmized X.
computing PCA
    on highly variable genes
    with n_comps=30
    finished (0:00:01)
computing moments based on connectivities
    finished (0:00:00) --> added 
    'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)
computing velocities
    finished (0:00:00) --> added 
    'velocity', velocity vectors for each individual cell (adata.layers)
```

## Initialize and run VIA



```pythonv0 = ov.single.pyVIA(adata=adata,adata_key='X_pca',adata_ncomps=n_pcs, basis='X_umap',
                         clusters='clusters',knn=20, root_user=None,
                         dataset='', random_seed=42,is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=50,
                         piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,
                         velocity_matrix=adata.layers['velocity'],gene_matrix=adata.X.todense(),velo_weight=0.5,
                         edgebundle_pruning_twice=False, edgebundle_pruning=0.15, pca_loadings = adata.varm['PCs']
                        )

v0.run()```
*Output:*
```2023-04-08 16:32:57.976202	Running VIA over input data of 2531 (samples) x 30 (features)
2023-04-08 16:32:57.976308	Knngraph has 20 neighbors
2023-04-08 16:32:58.800324	Finished global pruning of 20-knn graph used for clustering at level of 0.15. Kept 49.3 % of edges. 
2023-04-08 16:32:58.807060	Number of connected components used for clustergraph  is 1
2023-04-08 16:32:59.274711	Commencing community detection
2023-04-08 16:32:59.296206	Finished running Leiden algorithm. Found 67 clusters.
2023-04-08 16:32:59.296818	Merging 50 very small clusters (<10)
2023-04-08 16:32:59.297396	Finished detecting communities. Found 17 communities
2023-04-08 16:32:59.297494	Making cluster graph. Global cluster graph pruning level: 0.15
2023-04-08 16:32:59.300366	Graph has 1 connected components before pruning
2023-04-08 16:32:59.301039	Graph has 1 connected components after pruning
2023-04-08 16:32:59.301105	Graph has 1 connected components after reconnecting
2023-04-08 16:32:59.301317	0.0% links trimmed from local pruning relative to start
2023-04-08 16:32:59.301323	57.4% links trimmed from global pruning relative to start
2023-04-08 16:32:59.302247	Starting make edgebundle viagraph...
2023-04-08 16:32:59.302252	Make via clustergraph edgebundle
2023-04-08 16:33:00.142548	Hammer dims: Nodes shape: (17, 2) Edges shape: (52, 3)
size velocity matrix 2531 (2531, 5000)
2023-04-08 16:33:00.173534	Looking for initial states
2023-04-08 16:33:00.175672	Stationary distribution normed [0.086 0.038 0.181 0.009 0.084 0.027 0.003 0.11  0.028 0.037 0.103 0.097
 0.049 0.031 0.051 0.017 0.048]
2023-04-08 16:33:00.175849	Top 5 candidates for root: [ 6  3 15  5  8 13  9  1 16 12] with stationary prob (%) [0.323 0.932 1.742 2.733 2.78  3.051 3.681 3.762 4.82  4.902]
2023-04-08 16:33:00.175966	Top 5 candidates for terminal: [ 2  7 10 11  0]
2023-04-08 16:33:00.176867	component number 0 out of  [0]
2023-04-08 16:33:00.186598	Using the RNA velocity graph, A top3 candidate for initial state is 6 comprising predominantly of Ngn3 low EP cells
2023-04-08 16:33:00.186997	Using the RNA velocity graph, A top3 candidate for initial state is 3 comprising predominantly of Ngn3 high EP cells
2023-04-08 16:33:00.187360	Using the RNA velocity graph, A top3 candidate for initial state is 15 comprising predominantly of Delta cells
2023-04-08 16:33:00.187716	Using the RNA velocity graph, A top3 candidate for initial state is 5 comprising predominantly of Fev+ cells
2023-04-08 16:33:00.188051	Using the RNA velocity graph, A top3 candidate for initial state is 8 comprising predominantly of Ngn3 high EP cells
2023-04-08 16:33:00.188384	Using the RNA velocity graph, A top3 candidate for initial state is 13 comprising predominantly of Ngn3 high EP cells
2023-04-08 16:33:00.188723	Using the RNA velocity graph, A top3 candidate for initial state is 9 comprising predominantly of Ngn3 high EP cells
2023-04-08 16:33:00.189079	Using the RNA velocity graph, A top3 candidate for initial state is 1 comprising predominantly of Alpha cells
2023-04-08 16:33:00.189405	Using the RNA velocity graph, A top3 candidate for initial state is 16 comprising predominantly of Epsilon cells
2023-04-08 16:33:00.189749	Using the RNA velocity graph, A top3 candidate for initial state is 12 comprising predominantly of Epsilon cells
2023-04-08 16:33:00.189752	Using the RNA velocity graph, the suggested initial root state is 6 comprising predominantly of Ngn3 low EP cells
2023-04-08 16:33:00.189759	Computing lazy-teleporting expected hitting times
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:33:12.895333	Identifying terminal clusters corresponding to unique lineages...
2023-04-08 16:33:12.895502	Closeness:[0, 1, 2, 3, 4, 5, 6, 7, 10]
2023-04-08 16:33:12.895513	Betweenness:[0, 1, 4, 5, 6, 10, 12, 13, 15]
2023-04-08 16:33:12.895516	Out Degree:[0, 1, 2, 3, 4, 5, 6, 7, 10, 12, 15]
2023-04-08 16:33:12.895627	Cluster 0 had 3 or more neighboring terminal states [2, 4, 5, 7] and so we removed cluster 5
2023-04-08 16:33:12.895647	Cluster 2 had 3 or more neighboring terminal states [0, 7, 10] and so we removed cluster 0
2023-04-08 16:33:12.895758	Terminal clusters corresponding to unique lineages in this component are [1, 2, 4, 7, 10, 12] 
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:33:24.271222	From root 6,  the Terminal state 1 is reached 266 times.
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:33:34.962528	From root 6,  the Terminal state 2 is reached 640 times.
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:33:45.675023	From root 6,  the Terminal state 4 is reached 467 times.
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:33:56.019835	From root 6,  the Terminal state 7 is reached 628 times.
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:34:06.545564	From root 6,  the Terminal state 10 is reached 645 times.
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
2023-04-08 16:34:17.597466	From root 6,  the Terminal state 12 is reached 237 times.
2023-04-08 16:34:17.614112	Terminal clusters corresponding to unique lineages are {1: 'Alpha', 2: 'Beta', 4: 'Alpha', 7: 'Beta', 10: 'Beta', 12: 'Epsilon'}
2023-04-08 16:34:17.614139	Begin projection of pseudotime and lineage likelihood
2023-04-08 16:34:17.913029	Transition matrix with weight of 0.5 on RNA velocity
2023-04-08 16:34:17.913694	Graph has 1 connected components before pruning
2023-04-08 16:34:17.914918	Graph has 1 connected components after pruning
2023-04-08 16:34:17.914987	Graph has 1 connected components after reconnecting
2023-04-08 16:34:17.915221	40.4% links trimmed from local pruning relative to start
2023-04-08 16:34:17.915231	65.4% links trimmed from global pruning relative to start
2023-04-08 16:34:17.917529	Start making edgebundle milestone...
2023-04-08 16:34:17.917562	Start finding milestones
2023-04-08 16:34:18.189857	End milestones
2023-04-08 16:34:18.189946	Will use via-pseudotime for edges, otherwise consider providing a list of numeric labels (single cell level) or via_object
2023-04-08 16:34:18.193933	Recompute weights
2023-04-08 16:34:18.203565	pruning milestone graph based on recomputed weights
2023-04-08 16:34:18.204081	Graph has 1 connected components before pruning
2023-04-08 16:34:18.204384	Graph has 1 connected components after pruning
2023-04-08 16:34:18.204457	Graph has 1 connected components after reconnecting
2023-04-08 16:34:18.204894	63.9% links trimmed from global pruning relative to start
2023-04-08 16:34:18.204906	regenerate igraph on pruned edges
2023-04-08 16:34:18.207778	Setting numeric label as single cell pseudotime for coloring edges
2023-04-08 16:34:18.213883	Making smooth edges
2023-04-08 16:34:18.577366	Time elapsed 80.3 seconds
```

```pythonfig, ax, ax1 = v0.plot_piechart_graph(clusters='clusters',cmap='Reds',dpi=80,
                                   show_legend=False,ax_text=False,fontsize=4)
fig.set_size_inches(8,4)```
*Output:*
```/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/Pyomic/single/_via.py:914: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/Pyomic/single/_via.py:1011: MatplotlibDeprecationWarning: Auto-removal of grids by pcolor() and pcolormesh() is deprecated since 3.5 and will be removed two minor releases later; please call grid(False) first.
```
```<Figure size 640x320 with 20 Axes>```

## Visualize trajectory and cell progression

Fine grained vector fields

```pythonv0.plot_trajectory_gams(basis='X_umap',clusters='clusters',draw_all_curves=False)```
*Output:*
```100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
100% (11 of 11) |########################| Elapsed Time: 0:00:00 Time:  0:00:00
```
```2023-04-08 16:35:44.495793	Super cluster 1 is a super terminal with sub_terminal cluster 1
2023-04-08 16:35:44.495818	Super cluster 2 is a super terminal with sub_terminal cluster 2
2023-04-08 16:35:44.495824	Super cluster 4 is a super terminal with sub_terminal cluster 4
2023-04-08 16:35:44.495830	Super cluster 7 is a super terminal with sub_terminal cluster 7
2023-04-08 16:35:44.495836	Super cluster 10 is a super terminal with sub_terminal cluster 10
2023-04-08 16:35:44.495841	Super cluster 12 is a super terminal with sub_terminal cluster 12
```
```(<Figure size 640x320 with 3 Axes>,
 <AxesSubplot: title={'center': 'True Labels: ncomps:30. knn:20'}>,
 <AxesSubplot: title={'center': 'Pseudotime'}>)```
```<Figure size 640x320 with 3 Axes>```

```pythonv0.plot_stream(basis='X_umap',clusters='clusters',
               density_grid=0.8, scatter_size=30, scatter_alpha=0.3, linewidth=0.5)```
*Output:*
```(<Figure size 320x320 with 1 Axes>,
 <AxesSubplot: title={'center': 'Streamplot'}>)```
```<Figure size 320x320 with 1 Axes>```

## Draw lineage likelihoods

These indicate potential pathways corresponding to the 4 islets (two types of Beta islets Lineage 5 and 12)

```pythonv0.plot_lineage_probability()```
*Output:*
```2023-04-08 16:36:18.551874	Marker_lineages: [1, 2, 4, 7, 10, 12]
2023-04-08 16:36:18.561868	The number of components in the original full graph is 1
2023-04-08 16:36:18.561916	For downstream visualization purposes we are also constructing a low knn-graph 
2023-04-08 16:36:19.459300	Check sc pb [0.08  0.261 0.209 0.194 0.18  0.076]
2023-04-08 16:36:19.494450	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 1: [6, 3, 8, 9, 16, 12, 1]
2023-04-08 16:36:19.494481	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 2: [6, 3, 8, 9, 16, 11, 7, 2]
2023-04-08 16:36:19.494486	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 4: [6, 3, 8, 13, 5, 0, 4]
2023-04-08 16:36:19.494491	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 7: [6, 3, 8, 9, 16, 11, 7]
2023-04-08 16:36:19.494495	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 10: [6, 3, 8, 9, 16, 11, 7, 2, 10]
2023-04-08 16:36:19.494499	Cluster path on clustergraph starting from Root Cluster 6 to Terminal Cluster 12: [6, 3, 8, 9, 16, 12]
2023-04-08 16:36:19.805794	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 1 along path: [6, 6, 6, 6, 3, 9, 16, 12, 1, 1, 1]
2023-04-08 16:36:19.816800	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 2 along path: [6, 6, 6, 6, 3, 9, 16, 2, 2, 2, 2, 2, 2]
2023-04-08 16:36:19.827095	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 4 along path: [6, 6, 6, 6, 3, 9, 16, 11, 4, 4, 4]
2023-04-08 16:36:19.837298	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 7 along path: [6, 6, 6, 6, 3, 9, 16, 12, 4, 7]
2023-04-08 16:36:19.847061	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 10 along path: [6, 6, 6, 6, 3, 9, 16, 2, 10, 10, 10, 10, 10]
2023-04-08 16:36:19.856823	Revised Cluster level path on sc-knnGraph from Root Cluster 6 to Terminal Cluster 12 along path: [6, 6, 6, 6, 3, 9, 16, 12, 12, 12, 12]
```
```(<Figure size 640x320 with 14 Axes>,
 array([[<AxesSubplot: title={'center': 'Lineage: 1-Alpha'}>,
         <AxesSubplot: title={'center': 'Lineage: 2-Beta'}>,
         <AxesSubplot: title={'center': 'Lineage: 4-Alpha'}>,
         <AxesSubplot: title={'center': 'Lineage: 7-Beta'}>],
        [<AxesSubplot: title={'center': 'Lineage: 10-Beta'}>,
         <AxesSubplot: title={'center': 'Lineage: 12-Epsilon'}>,
         <AxesSubplot: >, <AxesSubplot: >]], dtype=object))```
```<Figure size 640x320 with 14 Axes>```

