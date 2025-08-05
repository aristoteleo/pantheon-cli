# t_cellfate_genesets
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_cellfate_genesets.ipynb*

# Timing-associated geneset analysis with cellfategenie

In our single-cell analysis, we analyse the underlying temporal state in the cell, which we call pseudotime. and identifying the genes associated with pseudotime becomes the key to unravelling models of gene dynamic regulation. In traditional analysis, we would use correlation coefficients, or gene dynamics model fitting. The correlation coefficient approach will have a preference for genes at the beginning and end of the time series, and the gene dynamics model requires RNA velocity information. Unbiased identification of chronosequence-related genes, as well as the need for no additional dependency information, has become a challenge in current chronosequence analyses.

Here, we developed CellFateGenie, which first removes potential noise from the data through metacells, and then constructs an adaptive ridge regression model to find the minimum set of genes needed to satisfy the timing fit.CellFateGenie has similar accuracy to gene dynamics models while eliminating preferences for the start and end of the time series.

We provided the AUCell to evaluate the geneset of adata

Colab_Reproducibility：https://colab.research.google.com/drive/1upcKKZHsZMS78eOliwRAddbaZ9ICXSrc?usp=sharing

```pythonimport omicverse as ov
import scvelo as scv
import matplotlib.pyplot as plt
ov.ov_plot_set()```

## Data preprocessed

We using dataset of dentategyrus in scvelo to demonstrate the timing-associated genes analysis. Firstly, We use `ov.pp.qc` and `ov.pp.preprocess` to preprocess the dataset.

Then we use `ov.pp.scale` and `ov.pp.pca` to analysis the principal component of the data

```pythonadata=ov.read('data/tutorial_meta_den.h5ad')
adata=adata.raw.to_adata()
adata```
*Output:*
```AnnData object with n_obs × n_vars = 200 × 13118
    obs: 'Pseudo-sizes', 'celltype', 'celltype_purity', 'pt_via'
    uns: 'celltype_colors', 'hvg', 'neighbors', 'scaled|original|cum_sum_eigenvalues', 'scaled|original|pca_var_ratios', 'umap'
    obsm: 'X_umap', 'scaled|original|X_pca'
    obsp: 'connectivities', 'distances'```

## Genesets evaluata

```pythonimport omicverse as ov
pathway_dict=ov.utils.geneset_prepare('../placenta/genesets/GO_Biological_Process_2021.txt',organism='Mouse')
len(pathway_dict.keys())```
*Output:*
```6036```

```python##Assest all pathways
adata_aucs=ov.single.pathway_aucell_enrichment(adata,
                                                pathways_dict=pathway_dict,
                                                num_workers=8)```

```pythonadata_aucs.obs=adata[adata_aucs.obs.index].obs
adata_aucs.obsm=adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata[adata_aucs.obs.index].obsp
adata_aucs.uns=adata[adata_aucs.obs.index].uns

adata_aucs```
*Output:*
```AnnData object with n_obs × n_vars = 200 × 6036
    obs: 'Pseudo-sizes', 'celltype', 'celltype_purity', 'pt_via'
    uns: 'celltype_colors', 'hvg', 'neighbors', 'scaled|original|cum_sum_eigenvalues', 'scaled|original|pca_var_ratios', 'umap'
    obsm: 'X_umap', 'scaled|original|X_pca'
    obsp: 'connectivities', 'distances'```

## Timing-associated genes analysis

We have encapsulated the cellfategenie algorithm into omicverse, and we can simply use omicverse to analysis.

```pythoncfg_obj=ov.single.cellfategenie(adata_aucs,pseudotime='pt_via')
cfg_obj.model_init()```
*Output:*
```$MSE|RMSE|MAE|R^2$:0.0057|0.075|0.058|0.94
```
```                                                        coef  abs(coef)  \
Regulon                                                                   
kidney morphogenesis (GO:0060993)                  -0.073373   0.073373   
protein de-ADP-ribosylation (GO:0051725)            0.060204   0.060204   
coronary artery morphogenesis (GO:0060982)         -0.055837   0.055837   
isoleucine metabolic process (GO:0006549)          -0.055651   0.055651   
fatty acid derivative catabolic process (GO:190... -0.052161   0.052161   
...                                                      ...        ...   
negative regulation of peptidase activity (GO:0...  0.000000   0.000000   
negative regulation of pathway-restricted SMAD ...  0.000000   0.000000   
negative regulation of oxidoreductase activity ...  0.000000   0.000000   
negative regulation of oxidative stress-induced...  0.000000   0.000000   
zymogen inhibition (GO:0097341)                     0.000000   0.000000   

                                                      values  
Regulon                                                       
kidney morphogenesis (GO:0060993)                   0.131669  
protein de-ADP-ribosylation (GO:0051725)            0.182627  
coronary artery morphogenesis (GO:0060982)          0.073891  
isoleucine metabolic process (GO:0006549)           0.269442  
fatty acid derivative catabolic process (GO:190...  0.423015  
...                                                      ...  
negative regulation of peptidase activity (GO:0...  0.000000  
negative regulation of pathway-restricted SMAD ...  0.000000  
negative regulation of oxidoreductase activity ...  0.000000  
negative regulation of oxidative stress-induced...  0.000000  
zymogen inhibition (GO:0097341)                     0.000000  

[6036 rows x 3 columns]```

We used Adaptive Threshold Regression to calculate the minimum number of gene sets that would have the same accuracy as the regression model constructed for all genes.

```pythoncfg_obj.ATR(stop=500)```
*Output:*
```coef_threshold:0.01839457079768181, r2:0.939279664515779
```
```     coef_threshold        r2
0          0.060204  0.455190
1          0.055837  0.467974
2          0.055651  0.613062
3          0.052161  0.689825
4          0.050659  0.709170
..              ...       ...
995        0.008238  0.942300
996        0.008232  0.942246
997        0.008226  0.942244
998        0.008220  0.942232
999        0.008216  0.942131

[1000 rows x 2 columns]```

```pythonfig,ax=cfg_obj.plot_filtering(color='#5ca8dc')
ax.set_title('Dentategyrus Metacells\nCellFateGenie')```
*Output:*
```Text(0.5, 1.0, 'Dentategyrus Metacells\nCellFateGenie')```
```<Figure size 240x240 with 1 Axes>```

```pythonres=cfg_obj.model_fit()```
*Output:*
```$MSE|RMSE|MAE|R^2$:0.0068|0.083|0.061|0.93
```

## Visualization

We prepared a series of function to visualize the result. we can use `plot_color_fitting` to observe the different cells how to transit with the pseudotime.

```pythoncfg_obj.plot_color_fitting(type='raw',cluster_key='celltype')```
*Output:*
```(<Figure size 240x240 with 1 Axes>,
 <AxesSubplot: title={'center': 'Regression Genes\nDimension: 6036'}, xlabel='True pseudotime', ylabel='Predicted pseudotime'>)```
```<Figure size 240x240 with 1 Axes>```

```pythoncfg_obj.plot_color_fitting(type='filter',cluster_key='celltype')```
*Output:*
```(<Figure size 240x240 with 1 Axes>,
 <AxesSubplot: title={'center': 'Regression Genes\nDimension: 333'}, xlabel='True pseudotime', ylabel='Predicted pseudotime'>)```
```<Figure size 240x240 with 1 Axes>```

## Kendalltau test

We can further narrow down the set of genes that satisfy the maximum regression coefficient. We used the kendalltau test to calculate the trend significance for each gene.

```pythonkt_filter=cfg_obj.kendalltau_filter()
kt_filter.head()```
*Output:*
```                                                    kendalltau_sta  \
kidney morphogenesis (GO:0060993)                        -0.506842   
chromatin silencing at telomere (GO:0006348)             -0.428416   
apical protein localization (GO:0045176)                 -0.102495   
fatty acid derivative catabolic process (GO:190...       -0.213859   
positive regulation of the force of heart contr...        0.566918   

                                                          pvalue  
kidney morphogenesis (GO:0060993)                   3.358608e-25  
chromatin silencing at telomere (GO:0006348)        3.625202e-19  
apical protein localization (GO:0045176)            3.229373e-02  
fatty acid derivative catabolic process (GO:190...  7.950218e-06  
positive regulation of the force of heart contr...  1.637204e-30  ```

```pythonvar_name=kt_filter.loc[kt_filter['pvalue']<kt_filter['pvalue'].mean()].index.tolist()
gt_obj=ov.single.gene_trends(adata_aucs,'pt_via',var_name)
gt_obj.calculate(n_convolve=10)```

```pythonprint(f"Dimension: {len(var_name)}")```
*Output:*
```Dimension: 269
```

```pythonfig,ax=gt_obj.plot_trend(color=ov.utils.blue_color[3])
ax.set_title(f'Dentategyrus meta\nCellfategenie',fontsize=13)```
*Output:*
```Text(0.5, 1.0, 'Dentategyrus meta\nCellfategenie')```
```<Figure size 240x240 with 1 Axes>```

```pythong=ov.utils.plot_heatmap(adata_aucs,var_names=var_name,
                  sortby='pt_via',col_color='celltype',
                 n_convolve=10,figsize=(1,6),show=False)

g.fig.set_size_inches(2, 6)
g.fig.suptitle('CellFateGenie',x=0.25,y=0.83,
               horizontalalignment='left',fontsize=12,fontweight='bold')
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),fontsize=12)
plt.show()```
*Output:*
```<Figure size 160x480 with 4 Axes>```

```pythongw_obj1=ov.utils.geneset_wordcloud(adata=adata_aucs[:,var_name],
                                  cluster_key='celltype',pseudotime='pt_via',figsize=(3,6))
gw_obj1.get()```
*Output:*
```Granule mature 300 26
Granule immature 300 4
Neuroblast 300 22
Cajal Retzius 300 37
GABA 300 37
Microglia 300 115
Mossy 300 44
nIPC 300 91
Astrocytes 300 53
OPC 300 17
Radial Glia-like 300 13
OL 300 69
Cck-Tox 300 64
```
```{'Granule mature': <wordcloud.wordcloud.WordCloud at 0x7fdbe79be880>,
 'Granule immature': <wordcloud.wordcloud.WordCloud at 0x7fdb7f2797f0>,
 'Neuroblast': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4a30>,
 'Cajal Retzius': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4d00>,
 'GABA': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4b50>,
 'Microglia': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4ee0>,
 'Mossy': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4880>,
 'nIPC': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4c40>,
 'Astrocytes': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4c70>,
 'OPC': <wordcloud.wordcloud.WordCloud at 0x7fdb7f06b040>,
 'Radial Glia-like': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4e20>,
 'OL': <wordcloud.wordcloud.WordCloud at 0x7fdb7efc4be0>,
 'Cck-Tox': <wordcloud.wordcloud.WordCloud at 0x7fdb7efd7430>}```

```pythong=gw_obj1.plot_heatmap(figwidth=6,cmap='RdBu_r')
plt.suptitle('CellFateGenie',x=0.18,y=0.95,
               horizontalalignment='left',fontsize=12,fontweight='bold')
```
*Output:*
```Text(0.18, 0.95, 'CellFateGenie')```
```<Figure size 480x480 with 14 Axes>```

