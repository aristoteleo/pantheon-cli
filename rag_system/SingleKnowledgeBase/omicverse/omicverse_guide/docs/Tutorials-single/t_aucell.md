# t_aucell
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_aucell.ipynb*

# Pathway analysis with AUCell

Single-cell RNA sequencing (scRNA-seq) is a powerful tool for exploring variations in cell types between conditions, tissue types, species, and individuals. When conducting scRNA-seq analysis, the differential gene expression (DEG) analysis of the single-cell data is almost always followed by gene set enrichment analysis. The aim of this analysis is to identify gene programs and biological processes, gene ontologies, or regulatory pathways that are overrepresented in a case group compared to a control group.

Paper: [AUCell: Identifying cells with active gene sets](https://bioconductor.org/packages/AUCell)

Code: https://github.com/aertslab/AUCell

Colab_Reproducibility：https://colab.research.google.com/drive/1Rk7Zopil-Ve1WFRQLV_AZOwCHRyhsOXk?usp=sharing

<div class="admonition warning">
  <p class="admonition-title">Warning</p>
  <p>
    There are many methods to determine the enrichment pathway between two groups, and the choice of method can significantly impact the conclusion.
  </p>
</div>


This tutorial focuses on using AUCell to complete the gene set enrichment in scRNA-seq data.


## Part.1 The Mathematical Principles of AUCell

AUCell uses the “Area Under the Curve” (AUC) to calculate whether a critical subset of the input gene set is enriched within the expressed genes for each cell. The distribution of AUC scores across all the cells allows exploring the relative expression of the signature. Since the scoring method is ranking-based, AUCell is independent of the gene expression units and the normalization procedure.

 In brief, the scoring method is based on a recovery analysis where the x-axis  is the ranking of all genes based on expression level (genes with the same expression value, e.g., '0', are randomly sorted); and the y-axis is the number of genes recovered from the input set. AUCell then uses the AUC to **calculate whether a critical subset of the input gene set is enriched at the top of the ranking for each cell**. In this way, the AUC represents the proportion of expressed genes in the signature and their relative expression values compared to the other genes within the cell. The output of this step is a matrix with the AUC score for each gene set in each cell.


## Part.2 Data preprocess

In this part, we load a test data and perform preliminary processing of the data, such as normalization and logarithmization, in order to make the data more interpretable



```pythonimport omicverse as ov
import scanpy as sc
import scvelo as scv

ov.plot_set()```
*Output:*
```Using numpy backend - consider installing ott-jax for GPU acceleration: pip install ott-jax
🔬 Starting plot initialization...
🧬 Detecting CUDA devices…
✅ [GPU 0] NVIDIA H100 80GB HBM3
    • Total memory: 79.1 GB
    • Compute capability: 9.0

   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

🔖 Version: 1.7.2rc1   📚 Tutorials: https://omicverse.readthedocs.io/
✅ plot_set complete.

```

```pythonov.utils.download_pathway_database()
ov.utils.download_geneid_annotation_pair()```
*Output:*
```......Pathway Geneset download start: GO_Biological_Process_2021
......Downloading dataset save to genesets/GO_Biological_Process_2021.txt
......[GO_Biological_Process_2021 Size of file]: 0.15 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.649472713470459.2f s
......Pathway Geneset download start: GO_Cellular_Component_2021
......Downloading dataset save to genesets/GO_Cellular_Component_2021.txt
......[GO_Cellular_Component_2021 Size of file]: 0.03 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.330129623413086.2f s
......Pathway Geneset download start: GO_Molecular_Function_2021
......Downloading dataset save to genesets/GO_Molecular_Function_2021.txt
......[GO_Molecular_Function_2021 Size of file]: 0.03 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.3111729621887207.2f s
......Pathway Geneset download start: WikiPathway_2021_Human
......Downloading dataset save to genesets/WikiPathway_2021_Human.txt
......[WikiPathway_2021_Human Size of file]: 0.02 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.3150238990783691.2f s
......Pathway Geneset download start: WikiPathways_2019_Mouse
......Downloading dataset save to genesets/WikiPathways_2019_Mouse.txt
......[WikiPathways_2019_Mouse Size of file]: 0.01 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.0507233142852783.2f s
......Pathway Geneset download start: Reactome_2022
......Downloading dataset save to genesets/Reactome_2022.txt
......[Reactome_2022 Size of file]: 0.07 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！1.4951279163360596.2f s
......Pathway Geneset download finished!
......Other Genesets can be dowload in `https://maayanlab.cloud/Enrichr/#libraries`
......Geneid Annotation Pair download start: pair_GRCm39
......Loading dataset from genesets/pair_GRCm39.tsv
......Geneid Annotation Pair download start: pair_T2TCHM13
......Loading dataset from genesets/pair_T2TCHM13.tsv
......Geneid Annotation Pair download start: pair_GRCh38
......Loading dataset from genesets/pair_GRCh38.tsv
......Geneid Annotation Pair download start: pair_GRCh37
......Loading dataset from genesets/pair_GRCh37.tsv
......Geneid Annotation Pair download start: pair_danRer11
......Loading dataset from genesets/pair_danRer11.tsv
......Geneid Annotation Pair download start: pair_danRer7
......Loading dataset from genesets/pair_danRer7.tsv
......Geneid Annotation Pair download finished!
```

The dataset used here is the mouse pancreas dataset that comes with scvelo

```pythonadata = scv.datasets.pancreas()
adata```
*Output:*
```try downloading from url
https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad
... this may take a while but only happens once
creating directory data/Pancreas/ for saving data
```
```  0%|          | 0.00/50.0M [00:00<?, ?B/s]```
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

```pythonsc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)```
*Output:*
```normalizing counts per cell
    finished (0:00:00)
```

```pythonadata.X.max()```
*Output:*
```7.785111```

Now the max value of anndata object is 7.78

## Part.3 Pathway anaylsis

In this part, we will demonstrate how to utilize the existing gene set to conduct enrichment analysis and how to create a gene set based on our own ideas for enrichment analysis in the test dataset.

First, we need to download data sets, such as the Gene Ontology(GO) or the Kyoto Encyclopedia of Genes and Genomes(KEGG). It should be noted that here we need to select the correct species, such as 'Human'or 'Mouse'.



```pythonpathway_dict=ov.utils.geneset_prepare('genesets/GO_Biological_Process_2021.txt',organism='Mouse')```

When working with existing datasets, it is possible to use the `ov.single.geneset_aucell` to calculate the activity of a gene set that corresponds to a particular signaling pathway within the dataset. 

Additionally, we can use the `sc.pl.embedding` function to visualize the distribution of gene set activity. By doing so, we can gain insights into the behavior of the gene set within the dataset and how it relates to the signaling pathway of interest.

```python##Assest one geneset
geneset_name='response to vitamin (GO:0033273)'
ov.single.geneset_aucell(adata,
                            geneset_name=geneset_name,
                            geneset=pathway_dict[geneset_name])
sc.pl.embedding(adata,
                basis='umap',
          color=["{}_aucell".format(geneset_name)])```
*Output:*
```ctxcore have been install version: 0.2.0
```
```100%|██████████| 3696/3696 [00:01<00:00, 2199.20it/s]
```
```<Figure size 320x320 with 2 Axes>```

We also can calculate the AUCell of more than one geneset

```python##Assest more than one geneset
geneset_names=['response to vitamin (GO:0033273)','response to vitamin D (GO:0033280)']
ov.single.pathway_aucell(adata,
                            pathway_names=geneset_names,
                            pathways_dict=pathway_dict)
sc.pl.embedding(adata,
                basis='umap',
                color=[i+'_aucell' for i in geneset_names])```
*Output:*
```ctxcore have been install version: 0.2.0
```
```100%|██████████| 3696/3696 [00:01<00:00, 2179.88it/s]
```
```<Figure size 772.8x320 with 4 Axes>```

In certain situations, the pathway we wish to investigate may not be available in the database. In such cases, we can manually define the gene set and its corresponding genes, calculate the AUCell activity of the gene set, and then visualize the results. By doing so, we can gain insights into the behavior of the gene set and its associated pathway, even in the absence of a pre-existing database entry.


```python##Assest test geneset
ov.single.geneset_aucell(adata,
                            geneset_name='Sox',
                            geneset=['Sox17', 'Sox4', 'Sox7', 'Sox18', 'Sox5'])
sc.pl.embedding(adata,
                basis='umap',
          color=["Sox_aucell"])```
*Output:*
```ctxcore have been install version: 0.2.0
```
```100%|██████████| 3696/3696 [00:01<00:00, 2167.21it/s]
```
```<Figure size 320x320 with 2 Axes>```

Occasionally, we may wish to examine clusters-specific signaling pathways from a comprehensive perspective. In such cases, we can compute the AUCell scores of all signaling pathways in the database and store the results as an anndata file. Once saved, we can utilize differential expression gene calculation functions and visualization functions from Scanpy to conduct downstream analyses. By taking this approach, we can investigate signaling pathways that are specific to certain clusters and gain insights into the underlying biological mechanisms.


```python##Assest all pathways
adata_aucs=ov.single.pathway_aucell_enrichment(adata,
                                                  pathways_dict=pathway_dict,
                                                  num_workers=8)```

We can use the anndata class to store the resulting gene set pathway and visualize it using the functions provided by the anndata class.

```pythonadata_aucs.obs=adata[adata_aucs.obs.index].obs
adata_aucs.obsm=adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata[adata_aucs.obs.index].obsp
adata_aucs```
*Output:*
```AnnData object with n_obs × n_vars = 3696 × 6036
    obs: 'clusters_coarse', 'clusters', 'S_score', 'G2M_score', 'response to vitamin (GO:0033273)_aucell', 'response to vitamin D (GO:0033280)_aucell', 'Sox_aucell'
    obsm: 'X_pca', 'X_umap'
    obsp: 'distances', 'connectivities'```

```pythonadata_aucs.write_h5ad('data/pancreas_auce.h5ad',compression='gzip')```

We can load the AUCell result of scRNA-seq and visualize the geneset we test above

```pythonadata_aucs=sc.read('data/pancreas_auce.h5ad')```

```pythonsc.pl.embedding(adata_aucs,
                basis='umap',
          color=geneset_names)```
*Output:*
```<Figure size 772.8x320 with 4 Axes>```

## Part4. Visualize differential enrichment pathways between different cell clusters.

We first read the AUCell scores files of all previously calculated signal pathways。

Given that the AUCell score is stored in the anndata structure, and is roughly equivalent to the level of gene expression, we can utilize the algorithm that calculates differential expression genes across clusters from Scanpy  to determine clusters-specific signaling pathways, such as `sc.tl.rank_genes_groups`. By employing this approach, we can identify the pathways that are most distinctive to certain clusters and gain a better understanding of their underlying biology.


```python#adata_aucs.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata_aucs, 'clusters', method='t-test',n_genes=100)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='clusters',
                                cmap='Spectral_r',
                                standard_scale='var',n_genes=3)```
*Output:*
```ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:02)
WARNING: dendrogram data not found (using key=dendrogram_clusters). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
    using 'X_pca' with n_pcs = 50
Storing dendrogram info using `.uns['dendrogram_clusters']`
```
```<Figure size 894.4x304 with 6 Axes>```

Furthermore, we can extract specific signaling pathways that are unique to particular clusters and then visualize them. By doing so, we can explore the differences between these clusters and gain insights into the molecular mechanisms that underlie these distinctions. 

```pythondegs = sc.get.rank_genes_groups_df(adata_aucs, group='Beta', key='rank_genes_groups', log2fc_min=2, 
                                    pval_cutoff=0.05)['names'].squeeze()
degs```
*Output:*
```0               insulin metabolic process (GO:1901142)
1       amylin receptor signaling pathway (GO:0097647)
2    calcitonin family receptor signaling pathway (...
3    regulation of osteoclast differentiation (GO:0...
Name: names, dtype: object```

```pythonimport matplotlib.pyplot as plt
#fig, axes = plt.subplots(4,3,figsize=(12,9))
axes=sc.pl.embedding(adata_aucs,ncols=3,
                basis='umap',show=False,return_fig=True,wspace=0.55,hspace=0.65,
                color=['clusters']+degs.values.tolist(),
                title=[ov.utils.plot_text_set(i,3,20)for i in ['clusters']+degs.values.tolist()])

axes.tight_layout()```
*Output:*
```<Figure size 1488x640 with 9 Axes>```

## Part.5 Enrichment of geneset in scRNA-seq

In addition to using AUCell to assess gene sets, we can also use the cell's marker genes to find functions specific to each cell type. This is also commonly known as enrichment analysis.

```pythonadata.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata, 'clusters', method='t-test',n_genes=100)```
*Output:*
```ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
```

Additionally, we can calculate and visualize specific signaling pathways for all clusters with the function `ov.single.pathway_enrichment`.




```pythonres=ov.single.pathway_enrichment(adata,pathways_dict=pathway_dict,organism='Mouse',
                                     group_by='clusters',plot=True)```
*Output:*
```2025-06-13 16:29:57,818 [WARNING] Downloading mmusculus_gene_ensembl for the first time. It might take a couple of miniutes.
```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```
```<Figure size 960x320 with 1 Axes>```

To complete our analysis, we can use heatmaps to visualize the specific signaling pathways of each clusters with the function `ov.single.pathway_enrichment_plot`, with the depth of color reflecting the AUCell score for each respective clusters. This approach enables us to easily visualize and compare the activity of various signaling pathways across different clusters. By examining the heatmaps, we can identify patterns and differences between the clusters, providing us with a more in-depth understanding of the underlying biological mechanisms. Overall, this approach can be a powerful tool for identifying novel therapeutic targets and guiding the development of personalized treatments.





```pythonax=ov.single.pathway_enrichment_plot(res,plot_title='Enrichment',cmap='Reds',
                                         xticklabels=True,cbar=False,square=True,vmax=10,
                                         yticklabels=True,cbar_kws={'label': '-log10(qvalue)','shrink': 0.5,})```
*Output:*
```<Figure size 240x800 with 1 Axes>```

