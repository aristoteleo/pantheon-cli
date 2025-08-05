# t_scdeg
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_scdeg.ipynb*

# Differential expression analysis [Meta Cell]

Sometimes we need to compare differentially expressed genes or differentially expressed features between two cell types on single cell data, but existing methods focus more on cell-specific gene analysis. Researchers need to transfer bulk RNA-seq analysis to single-cell analysis, which involves interaction between different programming languages or programming tools, adding significantly to the workload of the researcher.

Here, we use omicverse's bulk RNA-seq pyDEG method to complete differential expression analysis at the single cell level. We will present two different perspectives, one from the perspective of all cells and one from the perspective of the metacellular.

Colab_Reproducibility：https://colab.research.google.com/drive/12faBRh0xT7v6KSy8NCSRqbegF_AEoDXr?usp=sharing

```pythonimport omicverse as ov
import scanpy as sc
import scvelo as scv

ov.utils.ov_plot_set()```
*Output:*
```
   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

Version: 1.5.5, Tutorials: https://omicverse.readthedocs.io/
```

## Data preprocessed 

We need to normalized and scale the data at first.

```pythonadata = scv.datasets.pancreas()
adata```
*Output:*
```AnnData object with n_obs × n_vars = 3696 × 27998
    obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score'
    var: 'highly_variable_genes'
    uns: 'clusters_coarse_colors', 'clusters_colors', 'day_colors', 'neighbors', 'pca'
    obsm: 'X_pca', 'X_umap'
    layers: 'spliced', 'unspliced'
    obsp: 'distances', 'connectivities'```

```pythonadata.X.max()```
*Output:*
```2286.0```

We found that the max value of anndata object larger than 10 and type is int. We need to normalize and log1p it

```python#quantity control
adata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.05, 'nUMIs': 500, 'detected_genes': 250})
#normalize and high variable genes (HVGs) calculated
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,)

#save the whole genes and filter the non-HVGs
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]

#scale the adata.X
ov.pp.scale(adata)

#Dimensionality Reduction
ov.pp.pca(adata,layer='scaled',n_pcs=50)```
*Output:*
```Calculate QC metrics
End calculation of QC metrics.
Original cell number: 3696
Begin of post doublets removal and QC plot
Running Scrublet
filtered out 12261 genes that are detected in less than 3 cells
normalizing counts per cell
    finished (0:00:00)
extracting highly variable genes
    finished (0:00:00)
--> added
    'highly_variable', boolean vector (adata.var)
    'means', float vector (adata.var)
    'dispersions', float vector (adata.var)
    'dispersions_norm', float vector (adata.var)
normalizing counts per cell
    finished (0:00:00)
normalizing counts per cell
    finished (0:00:00)
Embedding transcriptomes using PCA...
Automatically set threshold at doublet score = 0.35
Detected doublet rate = 0.2%
Estimated detectable doublet fraction = 56.2%
Overall doublet rate:
	Expected   = 5.0%
	Estimated  = 0.4%
    Scrublet finished (0:00:02)
Cells retained after scrublet: 3688, 8 removed.
End of post doublets removal and QC plots.
Filters application (seurat or mads)
Lower treshold, nUMIs: 500; filtered-out-cells: 0
Lower treshold, n genes: 250; filtered-out-cells: 0
Lower treshold, mito %: 0.05; filtered-out-cells: 0
Filters applicated.
Total cell filtered out with this last --mode seurat QC (and its chosen options): 0
Cells retained after scrublet and seurat filtering: 3688, 8 removed.
filtered out 12263 genes that are detected in less than 3 cells
Begin robust gene identification
After filtration, 15735/15735 genes are kept. Among 15735 genes, 15735 genes are robust.
End of robust gene identification.
Begin size normalization: shiftlog and HVGs selection pearson
normalizing counts per cell The following highly-expressed genes are not considered during normalization factor computation:
['Ghrl']
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
End of size normalization: shiftlog and HVGs selection pearson
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
```

```pythonadata.X.max()```
*Output:*
```11.971903```

## Different expression in total level

We then select the target cells to be analysed, including `Alpha` and `Beta`, derive the expression matrix using `to_df()` and build the differential expression analysis module using `pyDEG`

```pythontest_adata=adata[adata.obs['clusters'].isin(['Alpha','Beta'])]
test_adata```
*Output:*
```View of AnnData object with n_obs × n_vars = 1065 × 2000
    obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'highly_variable_genes', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'clusters_coarse_colors', 'clusters_colors', 'day_colors', 'neighbors', 'pca', 'scrublet', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues'
    obsm: 'X_pca', 'X_umap', 'scaled|original|X_pca'
    varm: 'scaled|original|pca_loadings'
    layers: 'spliced', 'unspliced', 'counts', 'scaled', 'lognorm'
    obsp: 'distances', 'connectivities'```

```pythondds=ov.bulk.pyDEG(test_adata.to_df(layer='lognorm').T)```

```pythondds.drop_duplicates_index()
print('... drop_duplicates_index success')```
*Output:*
```... drop_duplicates_index success
```

We also need to set up an experimental group and a control group, i.e. the two types of cells to be compared and analysed

```pythontreatment_groups=test_adata.obs[test_adata.obs['clusters']=='Alpha'].index.tolist()
control_groups=test_adata.obs[test_adata.obs['clusters']=='Beta'].index.tolist()
result=dds.deg_analysis(treatment_groups,control_groups,method='ttest')
```

```pythonresult.sort_values('qvalue').head()```
*Output:*
```         pvalue  qvalue  FoldChange  -log(pvalue)  -log(qvalue)  BaseMean  \
index                                                                       
Sytl4       0.0     0.0    0.030428           inf           inf  1.304178   
Slc25a5     0.0     0.0    1.098909           inf           inf  6.919838   
Etv1        0.0     0.0    5.125477           inf           inf  1.925031   
Ins2        0.0     0.0    0.096964           inf           inf  2.488212   
Adra2a      0.0     0.0    0.012644           inf           inf  1.188721   

         log2(BaseMean)    log2FC  abs(log2FC)      size  sig  
index                                                          
Sytl4          0.383141 -5.038435     5.038435  0.003043  sig  
Slc25a5        2.790738  0.136071     0.136071  0.109891  sig  
Etv1           0.944881  2.357686     2.357686  0.512548  sig  
Ins2           1.315110 -3.366402     3.366402  0.009696  sig  
Adra2a         0.249410 -6.305355     6.305355  0.001264  sig  ```

```python# -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)```
*Output:*
```... Fold change threshold: 1.7822449207305908
```

```pythondds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)```
*Output:*
```(<Figure size 320x320 with 1 Axes>,
 <AxesSubplot: title={'center': 'DEG Analysis'}, xlabel='$log_{2}FC$', ylabel='$-log_{10}(qvalue)$'>)```
```<Figure size 320x320 with 1 Axes>```

```pythondds.plot_boxplot(genes=['Irx1','Adra2a'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))```
*Output:*
```(<Figure size 160x240 with 1 Axes>,
 <AxesSubplot: title={'center': 'Gene Expression'}>)```
```<Figure size 160x240 with 1 Axes>```

```pythonov.utils.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=['clusters','Irx1','Adra2a'])```
*Output:*
```<Figure size 1159.2x320 with 5 Axes>```

## Different expression in Metacells level

Here, we calculated the metacells from the whole scRNA-seq datasets using SEACells, and the same analyze with total level.

### Constructing a metacellular object

We can use `ov.single.MetaCell` to construct a metacellular object to train the SEACells model, the arguments can be found in below.

- :param ad: (AnnData) annotated data matrix
- :param build_kernel_on: (str) key corresponding to matrix in ad.obsm which is used to compute kernel for metacells
                        Typically 'X_pca' for scRNA or 'X_svd' for scATAC
- :param n_SEACells: (int) number of SEACells to compute
- :param use_gpu: (bool) whether to use GPU for computation
- :param verbose: (bool) whether to suppress verbose program logging
- :param n_waypoint_eigs: (int) number of eigenvectors to use for waypoint initialization
- :param n_neighbors: (int) number of nearest neighbors to use for graph construction
- :param convergence_epsilon: (float) convergence threshold for Franke-Wolfe algorithm
- :param l2_penalty: (float) L2 penalty for Franke-Wolfe algorithm
- :param max_franke_wolfe_iters: (int) maximum number of iterations for Franke-Wolfe algorithm
- :param use_sparse: (bool) whether to use sparse matrix operations. Currently only supported for CPU implementation.

```pythonmeta_obj=ov.single.MetaCell(adata,use_rep='scaled|original|X_pca',n_metacells=150,
                           use_gpu=True)```
*Output:*
```Welcome to SEACells GPU!
Computing kNN graph using scanpy NN ...
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:01)
Computing radius for adaptive bandwidth kernel...
```
```  0%|          | 0/3688 [00:00<?, ?it/s]```
```Making graph symmetric...
Parameter graph_construction = union being used to build KNN graph...
Computing RBF kernel...
```
```  0%|          | 0/3688 [00:00<?, ?it/s]```
```Building similarity LIL matrix...
```
```  0%|          | 0/3688 [00:00<?, ?it/s]```
```Constructing CSR matrix...
```

```pythonmeta_obj.initialize_archetypes()```
*Output:*
```Building kernel on scaled|original|X_pca
Computing diffusion components from scaled|original|X_pca for waypoint initialization ... 
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:00)
Done.
Sampling waypoints ...
Done.
Selecting 135 cells from waypoint initialization.
Initializing residual matrix using greedy column selection
Initializing f and g...
```
```100%|██████████| 25/25 [00:00<00:00, 409.92it/s]```
```Selecting 15 cells from greedy initialization.
```
```
```

## Train and save the model

```pythonmeta_obj.train(min_iter=10, max_iter=50)```
*Output:*
```Randomly initialized A matrix.
Setting convergence threshold at 0.10173
Starting iteration 1.
Completed iteration 1.
Starting iteration 10.
Completed iteration 10.
Converged after 10 iterations.
```

```pythonmeta_obj.save('seacells/model.pkl')```

```pythonmeta_obj.load('seacells/model.pkl')```

## Predicted the metacells

we can use `predicted` to predicted the metacells of raw scRNA-seq data. There are two method can be selected, one is `soft`, the other is `hard`. 

In the `soft` method, Aggregates cells within each SEACell, summing over all raw data x assignment weight for all cells belonging to a SEACell. Data is un-normalized and pseudo-raw aggregated counts are stored in .layers['raw']. Attributes associated with variables (.var) are copied over, but relevant per SEACell attributes must be manually copied, since certain attributes may need to be summed, or averaged etc, depending on the attribute.

In the `hard` method, Aggregates cells within each SEACell, summing over all raw data for all cells belonging to a SEACell. Data is unnormalized and raw aggregated counts are stored .layers['raw']. Attributes associated with variables (.var) are copied over, but relevant per SEACell attributes must be manually copied, since certain attributes may need to be summed, or averaged etc, depending on the attribute.

```pythonad=meta_obj.predicted(method='soft',celltype_label='clusters',
                     summarize_layer='lognorm')```
*Output:*
```100%|██████████| 150/150 [00:06<00:00, 22.47it/s]
```

```pythonad.X.min(),ad.X.max()```
*Output:*
```(0.0, 10.963231738519204)```

```pythonimport matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4))
ov.utils.embedding(
    meta_obj.adata,
    basis="X_umap",
    color=['clusters'],
    frameon='small',
    title="Meta cells",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    size=10,
    ax=ax,
    alpha=0.2,
    #legend_loc='', 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=False,
    #palette=ov.utils.blue_color[:],
    #legend_fontweight='normal'
)
ov.single._metacell.plot_metacells(ax,meta_obj.adata,color='#CB3E35',
                                  )```
*Output:*
```<AxesSubplot: title={'center': 'Meta cells'}, xlabel='X_umap1', ylabel='X_umap2'>```
```<Figure size 320x320 with 1 Axes>```

### Differentially expressed analysis

Similar to total cells for differential expression analysis, we used metacells for differential expression in the same way.

```pythontest_adata=ad[ad.obs['celltype'].isin(['Alpha','Beta'])]
test_adata```
*Output:*
```View of AnnData object with n_obs × n_vars = 30 × 2000
    obs: 'Pseudo-sizes', 'celltype', 'celltype_purity'```

```pythondds_meta=ov.bulk.pyDEG(test_adata.to_df().T)```

```pythondds_meta.drop_duplicates_index()
print('... drop_duplicates_index success')```
*Output:*
```... drop_duplicates_index success
```

We also need to set up an experimental group and a control group, i.e. the two types of cells to be compared and analysed

```pythontreatment_groups=test_adata.obs[test_adata.obs['celltype']=='Alpha'].index.tolist()
control_groups=test_adata.obs[test_adata.obs['celltype']=='Beta'].index.tolist()
result=dds_meta.deg_analysis(treatment_groups,control_groups,method='ttest')```

```pythonresult.sort_values('qvalue').head()```
*Output:*
```               pvalue        qvalue  FoldChange  -log(pvalue)  -log(qvalue)  \
index                                                                         
Pdx1     1.077365e-15  3.715051e-14    0.281001     14.967637     13.430035   
Ctxn2    7.831602e-15  2.654780e-13    8.963869     14.106149     12.575971   
Nkx6-1   4.272626e-14  1.424209e-12    0.404696     13.369305     11.846426   
Gng12    2.126609e-13  6.972490e-12    0.294999     12.672312     11.156612   
Smarca1  2.937360e-13  9.475355e-12    3.167775     12.532043     11.023405   

         BaseMean  log2(BaseMean)    log2FC  abs(log2FC)      size  sig  
index                                                                    
Pdx1     3.454246        1.788371 -1.831352     1.831352  0.028100  sig  
Ctxn2    1.357059        0.440484  3.164122     3.164122  0.896387  sig  
Nkx6-1   3.753654        1.908296 -1.305089     1.305089  0.040470  sig  
Gng12    3.523417        1.816975 -1.761220     1.761220  0.029500  sig  
Smarca1  2.741508        1.454970  1.663470     1.663470  0.316777  sig  ```

```python# -1 means automatically calculates
dds_meta.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)```
*Output:*
```... Fold change threshold: 2.36044418611735
```

```pythondds_meta.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)```
*Output:*
```(<Figure size 320x320 with 1 Axes>,
 <AxesSubplot: title={'center': 'DEG Analysis'}, xlabel='$log_{2}FC$', ylabel='$-log_{10}(qvalue)$'>)```
```<Figure size 320x320 with 1 Axes>```

```pythondds_meta.plot_boxplot(genes=['Ctxn2','Mnx1'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))```
*Output:*
```(<Figure size 160x240 with 1 Axes>,
 <AxesSubplot: title={'center': 'Gene Expression'}>)```
```<Figure size 160x240 with 1 Axes>```

```pythonov.utils.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=['clusters','Ctxn2','Mnx1'])```
*Output:*
```<Figure size 1159.2x320 with 5 Axes>```

