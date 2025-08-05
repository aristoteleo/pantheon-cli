# t_preprocess
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_preprocess.ipynb*

# Preprocessing the data of scRNA-seq with omicverse

The count table, a numeric matrix of genes × cells, is the basic input data structure in the analysis of single-cell RNA-sequencing data. A common preprocessing step is to adjust the counts for variable sampling efficiency and to transform them so that the variance is similar across the dynamic range. 

Suitable methods to preprocess the scRNA-seq is important. Here, we introduce some preprocessing step to help researchers can perform downstream analysis easyier.

User can compare our tutorial with [scanpy'tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) to learn how to use omicverse well

Colab_Reproducibility：https://colab.research.google.com/drive/1DXLSls_ppgJmAaZTUvqazNC_E7EDCxUe?usp=sharing

```pythonimport omicverse as ov
import scanpy as sc
ov.ov_plot_set()
```

The data consist of 3k PBMCs from a Healthy Donor and are freely available from 10x Genomics ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)). On a unix system, you can uncomment and run the following to download and unpack the data. The last line creates a directory for writing processed data.

```python# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write```

```pythonadata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata```
*Output:*
```... reading from cache file cache/data-filtered_gene_bc_matrices-hg19-matrix.h5ad
```
```AnnData object with n_obs × n_vars = 2700 × 32738
    var: 'gene_ids'```

```pythonadata.var_names_make_unique()
adata.obs_names_make_unique()```

## Preprocessing

### Quantity control

For single-cell data, we require quality control prior to analysis, including the removal of cells containing double cells, low-expressing cells, and low-expressing genes. In addition to this, we need to filter based on mitochondrial gene ratios, number of transcripts, number of genes expressed per cell, cellular Complexity, etc. For a detailed description of the different QCs please see the document: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

```pythonadata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.05, 'nUMIs': 500, 'detected_genes': 250})
adata```
*Output:*
```Calculate QC metrics
End calculation of QC metrics.
Original cell number: 2700
Begin of post doublets removal and QC plot
Running Scrublet
filtered out 19024 genes that are detected in less than 3 cells
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
Automatically set threshold at doublet score = 0.31
Detected doublet rate = 1.4%
Estimated detectable doublet fraction = 35.1%
Overall doublet rate:
	Expected   = 5.0%
	Estimated  = 4.0%
    Scrublet finished (0:00:02)
Cells retained after scrublet: 2662, 38 removed.
End of post doublets removal and QC plots.
Filters application (seurat or mads)
Lower treshold, nUMIs: 500; filtered-out-cells: 0
Lower treshold, n genes: 250; filtered-out-cells: 3
Lower treshold, mito %: 0.05; filtered-out-cells: 56
Filters applicated.
Total cell filtered out with this last --mode seurat QC (and its chosen options): 59
Cells retained after scrublet and seurat filtering: 2603, 97 removed.
filtered out 19107 genes that are detected in less than 3 cells
```
```AnnData object with n_obs × n_vars = 2603 × 13631
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells'
    uns: 'scrublet'```

### High variable Gene Detection

Here we try to use Pearson's method to calculate highly variable genes. This is the method that is proposed to be superior to ordinary normalisation. See [Article](https://www.nature.com/articles/s41592-023-01814-1#MOESM3) in *Nature Method* for details.


Sometimes we need to recover the original counts for some single-cell calculations, but storing them in the layer layer may result in missing data, so we provide two functions here, a store function and a release function, to save the original data.

We set `layers=counts`, the counts will be stored in `adata.uns['layers_counts']`

```pythonov.utils.store_layers(adata,layers='counts')
adata```
*Output:*
```......The X of adata have been stored in counts
```
```AnnData object with n_obs × n_vars = 2603 × 13631
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells'
    uns: 'scrublet', 'layers_counts'```

normalize|HVGs：We use | to control the preprocessing step, | before for the normalisation step, either `shiftlog` or `pearson`, and | after for the highly variable gene calculation step, either `pearson` or `seurat`. Our default is `shiftlog|pearson`.

- if you use `mode`=`shiftlog|pearson` you need to set `target_sum=50*1e4`, more people like to se `target_sum=1e4`, we test the result think 50*1e4 will be better
- if you use `mode`=`pearson|pearson`, you don't need to set `target_sum`

<div class="admonition warning">
  <p class="admonition-title">Note</p>
  <p>
    if the version of `omicverse` lower than `1.4.13`, the mode can only be set between `scanpy` and `pearson`.
  </p>
</div>


```pythonadata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,)
adata```
*Output:*
```Begin robust gene identification
After filtration, 13631/13631 genes are kept. Among 13631 genes, 13631 genes are robust.
End of robust gene identification.
Begin size normalization: shiftlog and HVGs selection pearson
normalizing counts per cell The following highly-expressed genes are not considered during normalization factor computation:
[]
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
```
```AnnData object with n_obs × n_vars = 2603 × 13631
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg'
    layers: 'counts'```

Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

```pythonadata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata```
*Output:*
```View of AnnData object with n_obs × n_vars = 2603 × 2000
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg'
    layers: 'counts'```

We find that the adata.X matrix is normalized at this point, including the data in raw, but we want to get the unnormalized data, so we can use the retrieve function `ov.utils.retrieve_layers`

```pythonadata_counts=adata.copy()
ov.utils.retrieve_layers(adata_counts,layers='counts')
print('normalize adata:',adata.X.max())
print('raw count adata:',adata_counts.X.max())```
*Output:*
```......The X of adata have been stored in raw
......The layers counts of adata have been retreved
normalize adata: 11.381063
raw count adata: 419.0
```

```pythonadata_counts```
*Output:*
```AnnData object with n_obs × n_vars = 2603 × 2000
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'layers_raw'
    layers: 'counts'```

If we wish to recover the original count matrix at the whole gene level, we can try the following code

```pythonadata_counts=adata.raw.to_adata().copy()
ov.utils.retrieve_layers(adata_counts,layers='counts')
print('normalize adata:',adata.X.max())
print('raw count adata:',adata_counts.X.max())
adata_counts```
*Output:*
```......The X of adata have been stored in raw
......The layers counts of adata have been retreved
normalize adata: 11.381063
raw count adata: 419.0
```
```AnnData object with n_obs × n_vars = 2603 × 13631
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'layers_raw'```

## Principal component analysis

In contrast to scanpy, we do not directly scale the variance of the original expression matrix, but store the results of the variance scaling in the layer, due to the fact that scale may cause changes in the data distribution, and we have not found scale to be meaningful in any scenario other than a principal component analysis

```pythonov.pp.scale(adata)
adata```
*Output:*
```... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
```
```AnnData object with n_obs × n_vars = 2603 × 2000
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg'
    layers: 'counts', 'scaled'```

If you want to perform pca in normlog layer, you can set `layer`=`normlog`, but we think scaled is necessary in PCA.

```pythonov.pp.pca(adata,layer='scaled',n_pcs=50)
adata```
*Output:*
```AnnData object with n_obs × n_vars = 2603 × 2000
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues'
    obsm: 'scaled|original|X_pca'
    varm: 'scaled|original|pca_loadings'
    layers: 'counts', 'scaled', 'lognorm'```

```pythonadata.obsm['X_pca']=adata.obsm['scaled|original|X_pca']
ov.utils.embedding(adata,
                  basis='X_pca',
                  color='CST3',
                  frameon='small')```
*Output:*
```<Figure size 320x320 with 2 Axes>```

## Embedding the neighborhood graph

We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below. It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:

```pythonsc.pp.neighbors(adata, n_neighbors=15, n_pcs=50,
               use_rep='scaled|original|X_pca')```
*Output:*
```computing neighbors
```
```OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
```
```    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:01)
```

To visualize the PCA’s embeddings, we use the `pymde` package wrapper in omicverse. This is an alternative to UMAP that is GPU-accelerated.

```pythonadata.obsm["X_mde"] = ov.utils.mde(adata.obsm["scaled|original|X_pca"])
adata```
*Output:*
```AnnData object with n_obs × n_vars = 2603 × 2000
    obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
    var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
    uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues', 'neighbors'
    obsm: 'scaled|original|X_pca', 'X_pca', 'X_mde'
    varm: 'scaled|original|pca_loadings'
    layers: 'counts', 'scaled', 'lognorm'
    obsp: 'distances', 'connectivities'```

```pythonov.utils.embedding(adata,
                basis='X_mde',
                color='CST3',
                frameon='small')```
*Output:*
```<Figure size 320x320 with 2 Axes>```

You also can use `umap` to visualize the neighborhood graph

```pythonsc.tl.umap(adata)```
*Output:*
```computing UMAP
    finished: added
    'X_umap', UMAP coordinates (adata.obsm) (0:00:02)
```

```pythonov.utils.embedding(adata,
                basis='X_umap',
                color='CST3',
                frameon='small')```
*Output:*
```<Figure size 320x320 with 2 Axes>```

## Clustering the neighborhood graph

As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018). Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

```pythonsc.tl.leiden(adata)```
*Output:*
```running Leiden clustering
    finished: found 11 clusters and added
    'leiden', the cluster labels (adata.obs, categorical) (0:00:00)
```

We redesigned the visualisation of embedding to distinguish it from scanpy's embedding by adding the parameter `fraemon='small'`, which causes the axes to be scaled with the colourbar

```pythonov.utils.embedding(adata,
                basis='X_mde',
                color=['leiden', 'CST3', 'NKG7'],
                frameon='small')```
*Output:*
```<Figure size 1159.2x320 with 5 Axes>```

We also provide a boundary visualisation function `ov.utils.plot_ConvexHull` to visualise specific clusters.

Arguments: 
- color: if None will use the color of clusters
- alpha: default is 0.2

```pythonimport matplotlib.pyplot as plt
fig,ax=plt.subplots( figsize = (4,4))

ov.utils.embedding(adata,
                basis='X_mde',
                color=['leiden'],
                frameon='small',
                show=False,
                ax=ax)

ov.utils.plot_ConvexHull(adata,
                basis='X_mde',
                cluster_key='leiden',
                hull_cluster='0',
                ax=ax)
```
*Output:*
```leiden_colors
```
```<AxesSubplot: title={'center': 'leiden'}, xlabel='X_mde1', ylabel='X_mde2'>```
```<Figure size 320x320 with 1 Axes>```

If you have too many labels, e.g. too many cell types, and you are concerned about cell overlap, then consider trying the `ov.utils.gen_mpl_labels` function, which improves text overlap.
In addition, we make use of the `patheffects` function, which makes our text have outlines

- adjust_kwargs: it could be found in package `adjusttext`
- text_kwargs: it could be found in class `plt.texts`

```pythonfrom matplotlib import patheffects
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4))

ov.utils.embedding(adata,
                  basis='X_mde',
                  color=['leiden'],
                   show=False, legend_loc=None, add_outline=False, 
                   frameon='small',legend_fontoutline=2,ax=ax
                 )

ov.utils.gen_mpl_labels(
    adata,
    'leiden',
    exclude=("None",),  
    basis='X_mde',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 12 ,weight='bold',
                     path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)```
*Output:*
```<Figure size 320x320 with 1 Axes>```

```pythonmarker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']```

```pythonsc.pl.dotplot(adata, marker_genes, groupby='leiden',
             standard_scale='var');```
*Output:*
```<Figure size 623.2x388 with 4 Axes>```

## Finding marker genes

Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.

```pythonsc.tl.dendrogram(adata,'leiden',use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(adata, 'leiden', use_rep='scaled|original|X_pca',
                        method='t-test',use_raw=False,key_added='leiden_ttest')
sc.pl.rank_genes_groups_dotplot(adata,groupby='leiden',
                                cmap='Spectral_r',key='leiden_ttest',
                                standard_scale='var',n_genes=3)```
*Output:*
```Storing dendrogram info using `.uns['dendrogram_leiden']`
ranking genes
    finished: added to `.uns['leiden_ttest']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
```
```<Figure size 1160.8x388 with 6 Axes>```

cosg is also considered to be a better algorithm for finding marker genes. Here, omicverse provides the calculation of cosg

Paper: [Accurate and fast cell marker gene identification with COSG](https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext)

Code: https://github.com/genecell/COSG


```pythonsc.tl.rank_genes_groups(adata, groupby='leiden', 
                        method='t-test',use_rep='scaled|original|X_pca',)
ov.single.cosg(adata, key_added='leiden_cosg', groupby='leiden')
sc.pl.rank_genes_groups_dotplot(adata,groupby='leiden',
                                cmap='Spectral_r',key='leiden_cosg',
                                standard_scale='var',n_genes=3)```
*Output:*
```ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
**finished identifying marker genes by COSG**
```
```<Figure size 1160.8x388 with 6 Axes>```

## Other plotting

Next, let's try another chart, which we call the Stacked Volcano Chart. We need to prepare two dictionaries, a `data_dict` and a `color_dict`, both of which have the same key requirements.

For `data_dict`. we require the contents within each key to be a DataFrame containing ['names','logfoldchanges','pvals_adj'], where names stands for gene names, logfoldchanges stands for differential expression multiplicity, pvals_adj stands for significance p-value


```pythondata_dict={}
for i in adata.obs['leiden'].cat.categories:
    data_dict[i]=sc.get.rank_genes_groups_df(adata, group=i, key='leiden_ttest',
                                            pval_cutoff=None,log2fc_min=None)```

```pythondata_dict.keys()```
*Output:*
```dict_keys(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])```

```pythondata_dict[i].head()```
*Output:*
```   names     scores  logfoldchanges         pvals     pvals_adj
0    PF4  92.918243       17.017431  4.248822e-15  7.122920e-15
1   PPBP  85.572166       17.316359  6.224002e-15  1.041674e-14
2   SDPR  72.907913       16.093300  5.214614e-14  8.662149e-14
3  GNG11  72.813034       16.640303  5.858609e-14  9.723832e-14
4   GPX1  62.821903        9.493515  2.581554e-24  4.875457e-24```

For `color_dict`, we require that the colour to be displayed for the current key is stored within each key.`

```pythontype_color_dict=dict(zip(adata.obs['leiden'].cat.categories,
                         adata.uns['leiden_colors']))
type_color_dict```
*Output:*
```{'0': '#1f77b4',
 '1': '#ff7f0e',
 '2': '#279e68',
 '3': '#d62728',
 '4': '#aa40fc',
 '5': '#8c564b',
 '6': '#e377c2',
 '7': '#b5bd61',
 '8': '#17becf',
 '9': '#aec7e8',
 '10': '#ffbb78'}```

There are a number of parameters available here for us to customise the settings. Note that when drawing stacking_vol with omicverse version less than 1.4.13, there is a bug that the vertical coordinate is constant at [-15,15], so we have added some code in this tutorial for visualisation.

- data_dict: dict, in each key, there is a dataframe with columns of ['logfoldchanges','pvals_adj','names']
- color_dict: dict, in each key, there is a color for each omic
- pval_threshold: float, pvalue threshold for significant genes
- log2fc_threshold: float, log2fc threshold for significant genes
- figsize: tuple, figure size
- sig_color: str, color for significant genes
- normal_color: str, color for non-significant genes
- plot_genes_num: int, number of genes to plot
- plot_genes_fontsize: int, fontsize for gene names
- plot_genes_weight: str, weight for gene names

```pythonfig,axes=ov.utils.stacking_vol(data_dict,type_color_dict,
            pval_threshold=0.01,
            log2fc_threshold=2,
            figsize=(8,4),
            sig_color='#a51616',
            normal_color='#c7c7c7',
            plot_genes_num=2,
            plot_genes_fontsize=6,
            plot_genes_weight='bold',
            )

#The following code will be removed in future
y_min,y_max=0,0
for i in data_dict.keys():
    y_min=min(y_min,data_dict[i]['logfoldchanges'].min())
    y_max=max(y_max,data_dict[i]['logfoldchanges'].max())
for i in adata.obs['leiden'].cat.categories:
    axes[i].set_ylim(y_min,y_max)
plt.suptitle('Stacking_vol',fontsize=12)   ```
*Output:*
```0 2
2 4
4 6
6 8
8 10
10 12
12 14
14 16
16 18
18 20
20 22
```
```Text(0.5, 0.98, 'Stacking_vol')```
```<Figure size 640x320 with 11 Axes>```

