# t_anno_trans
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_anno_trans.ipynb*

# Celltype annotation transfer in multi-omics

In the field of multi-omics research, transferring cell type annotations from one data modality to another is a crucial step. For instance, when annotating cell types in single-cell ATAC sequencing (scATAC-seq) data, it's often desirable to leverage the cell type labels already annotated in single-cell RNA sequencing (scRNA-seq) data. This process involves integrating information from both scRNA-seq and scATAC-seq data modalities.

GLUE is a prominent algorithm used for cross-modality integration, allowing researchers to combine data from different omics modalities effectively. However, GLUE does not inherently provide a method for transferring cell type labels from scRNA-seq to scATAC-seq data. To address this limitation, an approach was implemented in the omicverse platform using K-nearest neighbor (KNN) graphs.

The KNN graph-based approach likely involves constructing KNN graphs separately for scRNA-seq and scATAC-seq data. In these graphs, each cell is connected to its K nearest neighbors based on certain similarity metrics, which could be calculated using gene expression profiles in scRNA-seq and accessibility profiles in scATAC-seq. Once these graphs are constructed, the idea is to transfer the cell type labels from the scRNA-seq side to the scATAC-seq side by assigning labels to scATAC-seq cells based on the labels of their KNN neighbors in the scRNA-seq graph.

Colab_Reproducibility：https://colab.research.google.com/drive/1aIMmSgyIw-PGjJ65WvMgz4Ob3EtoK_UV?usp=sharing

```pythonimport omicverse as ov
import matplotlib.pyplot as plt
import scanpy as sc
ov.ov_plot_set()```

## Loading the data preprocessed with GLUE

Here, we use two output files from the GLUE cross-modal integration, and their common feature is that they both have the `obsm['X_glue']` layer. And the rna have been annotated.

```pythonrna=sc.read("data/analysis_lymph/rna-emb.h5ad")
atac=sc.read("data/analysis_lymph/atac-emb.h5ad")```

We can visualize the intergrated effect of GLUE with UMAP

```pythonimport scanpy as sc
combined=sc.concat([rna,atac],merge='same')
combined```
*Output:*
```AnnData object with n_obs × n_vars = 68415 × 0
    obs: 'balancing_weight', 'domain'
    var: 'chromStart', 'chromEnd', 'highly_variable'
    obsm: 'X_glue'
    varm: 'X_glue'```

```pythoncombined.obsm['X_mde']=ov.utils.mde(combined.obsm['X_glue'])```

We can see that the two layers are correctly aligned

```pythonov.utils.embedding(combined,
               basis='X_mde',
               color='domain',
                title='Layers',
                show=False,
                palette=ov.utils.red_color,
                frameon='small'
               )```
*Output:*
```<AxesSubplot: title={'center': 'Layers'}, xlabel='X_mde1', ylabel='X_mde2'>```
```<Figure size 320x320 with 1 Axes>```

And the RNA modality has an already annotated cell type label on it

```pythonov.utils.embedding(rna,
               basis='X_mde',
               color='major_celltype',
                title='Cell type',
                show=False,
                #palette=ov.utils.red_color,
                frameon='small'
               )```
*Output:*
```<AxesSubplot: title={'center': 'Cell type'}, xlabel='X_mde1', ylabel='X_mde2'>```
```<Figure size 320x320 with 1 Axes>```

## Celltype transfer

We train a knn nearest neighbour classifier using `X_glue` features

```pythonknn_transformer=ov.utils.weighted_knn_trainer(
    train_adata=rna,
    train_adata_emb='X_glue',
    n_neighbors=15,
)```
*Output:*
```Weighted KNN with n_neighbors = 15 ... ```

```pythonlabels,uncert=ov.utils.weighted_knn_transfer(
    query_adata=atac,
    query_adata_emb='X_glue',
    label_keys='major_celltype',
    knn_model=knn_transformer,
    ref_adata_obs=rna.obs,
)```
*Output:*
```finished!
```

We migrate the training results of the KNN classifier to atac. `unc` stands for uncertainty, with higher uncertainty demonstrating lower migration accuracy, suggesting that the cell in question may be a double-fate signature or some other type of cell.

```pythonatac.obs["transf_celltype"]=labels.loc[atac.obs.index,"major_celltype"]
atac.obs["transf_celltype_unc"]=uncert.loc[atac.obs.index,"major_celltype"]```

```pythonatac.obs["major_celltype"]=atac.obs["transf_celltype"].copy()```

```pythonov.utils.embedding(atac,
               basis='X_umap',
               color=['transf_celltype_unc','transf_celltype'],
                #title='Cell type Un',
                show=False,
                palette=ov.palette()[11:],
                frameon='small'
               )```
*Output:*
```[<AxesSubplot: title={'center': 'transf_celltype_unc'}, xlabel='X_umap1', ylabel='X_umap2'>,
 <AxesSubplot: title={'center': 'transf_celltype'}, xlabel='X_umap1', ylabel='X_umap2'>]```
```<Figure size 772.8x320 with 3 Axes>```

## Visualization

We can merge atac and rna after migration annotation and observe on the umap plot whether the cell types are consistent after merging the modalities.

```pythonimport scanpy as sc
combined1=sc.concat([rna,atac],merge='same')
combined1```
*Output:*
```AnnData object with n_obs × n_vars = 68415 × 0
    obs: 'major_celltype', 'balancing_weight', 'domain'
    var: 'chromStart', 'chromEnd', 'highly_variable'
    obsm: 'X_glue'
    varm: 'X_glue'```

```pythoncombined1.obsm['X_mde']=ov.utils.mde(combined1.obsm['X_glue'])```

We found that the annotation was better, suggesting that the KNN nearest-neighbour classifier we constructed can effectively migrate cell type labels from RNA to ATAC.

```pythonov.utils.embedding(combined1,
               basis='X_mde',
               color=['domain','major_celltype'],
                title=['Layers','Cell type'],
                show=False,
                palette=ov.palette()[11:],
                frameon='small'
               )```
*Output:*
```[<AxesSubplot: title={'center': 'Layers'}, xlabel='X_mde1', ylabel='X_mde2'>,
 <AxesSubplot: title={'center': 'Cell type'}, xlabel='X_mde1', ylabel='X_mde2'>]```
```<Figure size 772.8x320 with 2 Axes>```

