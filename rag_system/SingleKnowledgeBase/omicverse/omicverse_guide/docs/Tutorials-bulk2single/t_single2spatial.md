# t_single2spatial
*Converted from: omicverse/omicverse_guide/docs/Tutorials-bulk2single/t_single2spatial.ipynb*

# Single RNA-seq to Spatial RNA-seq

Single2Spatial is used for scRNA-seq transformat to stRNA-seq. We extracted the deep-forest part of the Bulk2Space algorithm and constructed an algorithm that can calculate the spatial RNA-seq from scRNA-seq. In addition, we have redesigned the input and output of the data so that it can be more compatible with the analysis conventions in the Python environment.

Paper: [De novo analysis of bulk RNA-seq data at spatially resolved single-cell resolution](https://www.nature.com/articles/s41467-022-34271-z)

Code: https://github.com/ZJUFanLab/bulk2space

Colab_Reproducibility：https://colab.research.google.com/drive/1E3TlGUW2bxEew8JYCLE9qA6saYbIjBnY?usp=sharing

This tutorial walks through how to read, set-up and train the model from scRNA-seq and reference stRNA-seq data. We use the pdac datasets as example

```pythonimport scanpy as sc
import pandas as pd
import numpy as np
import omicverse as ov
import matplotlib.pyplot as plt

ov.utils.ov_plot_set()```

## loading data

Single2Spatial calcualted the the spot of celltype in spatial RNA-seq from scRNA-seq by deep-forest tree, we need an average normalize, log and scale the expression matrix of scRNA-seq, 

We also need an anndata object of stRNA-seq as reference to calculate. The data can be downloaded from [here](https://drive.google.com/file/d/1xB-Gk_KLxQA320-tycJp4CFHA66zF3LE/view?usp=sharing). The processed the PDAC scRNA-seq data and ST data (GSE111672)

```pythonimport anndata
raw_data=pd.read_csv('data/pdac/sc_data.csv', index_col=0)
single_data=anndata.AnnData(raw_data.T)
single_data.obs = pd.read_csv('data/pdac/sc_meta.csv', index_col=0)[['Cell_type']]
single_data```
*Output:*
```AnnData object with n_obs × n_vars = 1926 × 19104
    obs: 'Cell_type'```

```pythonraw_data=pd.read_csv('data/pdac/st_data.csv', index_col=0)
spatial_data=anndata.AnnData(raw_data.T)
spatial_data.obs = pd.read_csv('data/pdac/st_meta.csv', index_col=0)
spatial_data```
*Output:*
```AnnData object with n_obs × n_vars = 428 × 19104
    obs: 'Spot', 'xcoord', 'ycoord'```

## set up, training, saving, and loading

We can now set up the Single2Spatial object, which will ensure everything the model needs is in place for training. We need to specify the cell type of the scRNA-seq and the spot_key of the stRNA-seq. And specify the number of marker genes for each cell type for training.

```pythonst_model=ov.bulk2single.Single2Spatial(single_data=single_data,
                                           spatial_data=spatial_data,
                                           celltype_key='Cell_type',
                                           spot_key=['xcoord','ycoord'],
                                           #gpu='mps',
                                          )```
*Output:*
```...loading data
```

Now we can start to train our Single2Spatial model.


```python%%time
sp_adata=st_model.train(spot_num=500,
                cell_num=10,
                df_save_dir='data/pdac/predata_net/save_model',
                df_save_name='pdac_df',
                k=10,num_epochs=1000,batch_size=1000,predicted_size=32)```
*Output:*
```ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:03)
sucessfully create positive data
sucessfully create negative data
select top 500 marker genes of each cell type...
ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:03)
...The model size is 14478
...Net training
```
```Epoch [1000/1000], Loss: 0.5809: 100%|██████| 1000/1000 [05:31<00:00,  3.02it/s]
```
```...Net training done!
...save trained net in data/pdac/predata_net/save_model/pdac_df.pth.
Calculating scores...
```
```Now calculating scores for spot 428/428: 100%|██| 14/14 [01:07<00:00,  4.84s/it]
```
```Calculating scores done.
...save trained net in data/pdac/predata_net/save_model/pdac_df.pth.
CPU times: user 7min 42s, sys: 3min 1s, total: 10min 43s
Wall time: 7min 7s
```

We can also load our previously trained model directly

```python%%time
sp_adata=st_model.load(modelsize=14478,df_load_dir='data/pdac/predata_net/save_model/pdac_df.pth',
                k=10,predicted_size=32)```
*Output:*
```select top 500 marker genes of each cell type...
ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:03)
Calculating scores...
```
```Now calculating scores for spot 428/428: 100%|██| 14/14 [01:02<00:00,  4.43s/it]
```
```Calculating scores done.
CPU times: user 36.7 s, sys: 51.4 s, total: 1min 28s
Wall time: 1min 13s
```

## Investigate deconvolution results at spot-level resolution

We also evaluated the spatial deconvolution results at spot-level resolution. We firstly calculated cell type proportion per spot and then aggregated gene expression of cells per spot.

```pythonsp_adata_spot=st_model.spot_assess()```

## Plot spatial expression patterns of marker genes with single cell resolution

The spatial expression patterns of marker genes were also consistent with the spatial distribution of cell types.

```pythonsc.pl.embedding(
    sp_adata,
    basis="X_spatial",
    color=['REG1A', 'CLDN1', 'KRT16', 'MUC5B'],
    frameon=False,
    ncols=4,
    #save='_figure1_gene.png',
    show=False,
)```
*Output:*
```[<Axes: title={'center': 'REG1A'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'CLDN1'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'KRT16'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'MUC5B'}, xlabel='X_spatial1', ylabel='X_spatial2'>]```
```<Figure size 1545.6x320 with 8 Axes>```

We also can plot gene in spot level

```pythonsc.pl.embedding(
    sp_adata_spot,
    basis="X_spatial",
    color=['REG1A', 'CLDN1', 'KRT16', 'MUC5B'],
    frameon=False,
    ncols=4,
    #save='_figure1_gene.png',
    show=False,
)```
*Output:*
```[<Axes: title={'center': 'REG1A'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'CLDN1'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'KRT16'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'MUC5B'}, xlabel='X_spatial1', ylabel='X_spatial2'>]```
```<Figure size 1545.6x320 with 8 Axes>```

## Plot cell type proportion per spot

We can visualize the cell type proportion per spot by `sc.pl.embedding`

```pythonsc.pl.embedding(
    sp_adata_spot,
    basis="X_spatial",
    color=['Acinar cells','Cancer clone A','Cancer clone B','Ductal'],
    frameon=False,
    ncols=4,
    show=False,
    #save='_figure1_celltype_spot.png',
)```
*Output:*
```[<Axes: title={'center': 'Acinar cells'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'Cancer clone A'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'Cancer clone B'}, xlabel='X_spatial1', ylabel='X_spatial2'>,
 <Axes: title={'center': 'Ductal'}, xlabel='X_spatial1', ylabel='X_spatial2'>]```
```<Figure size 1545.6x320 with 8 Axes>```

```pythonsc.pl.embedding(
    sp_adata,
    basis="X_spatial",
    color=['Cell_type'],
    frameon=False,
    ncols=4,
    show=False,
    palette=ov.utils.ov_palette()[11:]
)```
*Output:*
```<Axes: title={'center': 'Cell_type'}, xlabel='X_spatial1', ylabel='X_spatial2'>```
```<Figure size 320x320 with 1 Axes>```

