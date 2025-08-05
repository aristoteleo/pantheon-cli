# t_slat
*Converted from: omicverse/omicverse_guide/docs/Tutorials-space/t_slat.ipynb*

# Single cell spatial alignment tools

SLAT (Spatially-Linked Alignment Tool), a graph-based algorithm for efficient and effective alignment of spatial slices. Adopting a graph adversarial matching strategy, SLAT is the first algorithm capable of aligning heterogenous spatial data across distinct technologies and modalities. 

We made two improvements in integrating the STT algorithm in OmicVerse:

- **Fix the running error in alignment**: We fixed some issues with the scSLAT package on pypi.
- **Added more downstream analysis**: We have expanded on the original tutorial by combining the tutorial and reproduce code given by the authors for downstream analysis.

If you found this tutorial helpful, please cite SLAT and OmicVerse: 

- Xia, CR., Cao, ZJ., Tu, XM. et al. Spatial-linked alignment tool (SLAT) for aligning heterogenous slices. Nat Commun 14, 7236 (2023). https://doi.org/10.1038/s41467-023-43105-5

```pythonimport omicverse as ov
import os

import scanpy as sc
import numpy as np
import pandas as pd
import torch
ov.plot_set()```
*Output:*
```
   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

Version: 1.6.0, Tutorials: https://omicverse.readthedocs.io/
```

```python#import scSLAT
from omicverse.external.scSLAT.model import load_anndatas, Cal_Spatial_Net, run_SLAT, scanpy_workflow, spatial_match
from omicverse.external.scSLAT.viz import match_3D_multi, hist, Sankey, match_3D_celltype, Sankey,Sankey_multi,build_3D
from omicverse.external.scSLAT.metrics import region_statistics```

## Preprocess Data

adata1.h5ad: E11.5 mouse embryo dataset, download from [here](https://drive.google.com/uc?export=download&id=1KkuJt6aSlKS1AJzFZjE_odypY-GINRuD)

adata2.h5ad: E12.5 mouse embryo dataset, download from [here](https://drive.google.com/uc?export=download&id=1YIiEmjGfHxcDbGn4nv2kzmTHUo3_q5hJ)

```pythonadata1 = sc.read_h5ad('data/E115_Stereo.h5ad')
adata2 = sc.read_h5ad('data/E125_Stereo.h5ad')```

```pythonadata1.obs['week']='E11.5'
adata2.obs['week']='E12.5'```

```pythonsc.pl.spatial(adata1, color='annotation', spot_size=3)
sc.pl.spatial(adata2, color='annotation', spot_size=3)```
*Output:*
```<Figure size 320x320 with 1 Axes>```
```<Figure size 320x320 with 1 Axes>```

## Run SLAT

Then we run SLAT as usual

```pythonCal_Spatial_Net(adata1, k_cutoff=20, model='KNN')
Cal_Spatial_Net(adata2, k_cutoff=20, model='KNN')
edges, features = load_anndatas([adata1, adata2], feature='DPCA', check_order=False)```
*Output:*
```Calculating spatial neighbor graph ...
The graph contains 218282 edges, 10000 cells.
21.8282 neighbors per cell on average.
Calculating spatial neighbor graph ...
The graph contains 219259 edges, 10001 cells.
21.923707629237075 neighbors per cell on average.
Use DPCA feature to format graph
If you pass `n_top_genes`, all cutoffs are ignored.
extracting highly variable genes
--> added
    'highly_variable', boolean vector (adata.var)
    'highly_variable_rank', float vector (adata.var)
    'means', float vector (adata.var)
    'variances', float vector (adata.var)
    'variances_norm', float vector (adata.var)
normalizing counts per cell
    finished (0:00:00)
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
Warning! Dual PCA is using GPU, which may lead to OUT OF GPU MEMORY in big dataset!
```

```pythonembd0, embd1, time = run_SLAT(features, edges, LGCN_layer=5)```
*Output:*
```Choose GPU:0 as device
Running
---------- epochs: 1 ----------
---------- epochs: 2 ----------
---------- epochs: 3 ----------
---------- epochs: 4 ----------
---------- epochs: 5 ----------
---------- epochs: 6 ----------
Training model time: 1.21
```

```pythonbest, index, distance = spatial_match([embd0, embd1], reorder=False, adatas=[adata1,adata2])```

```pythonmatching = np.array([range(index.shape[0]), best])
best_match = distance[:,0]
region_statistics(best_match, start=0.5, number_of_interval=10)```
*Output:*
```0.500~0.550         5          0.050%
0.550~0.600        48          0.480%
0.600~0.650       162          1.620%
0.650~0.700       425          4.250%
0.700~0.750      1004         10.039%
0.750~0.800      1578         15.778%
0.800~0.850      1900         18.998%
0.850~0.900      1701         17.008%
0.900~0.950      2141         21.408%
0.950~1.000      1036         10.359%
```

## Visualization of alignment

```pythonimport matplotlib.pyplot as plt
matching_list=[matching]
model = build_3D([adata1,adata2], matching_list,subsample_size=300, )
ax=model.draw_3D(hide_axis=True, line_color='#c2c2c2', height=1, size=[6,6], line_width=1)```
*Output:*
```Mapping 0th layer 
```
```<Figure size 480x480 with 1 Axes>```

Then we check the alignment quality of the whole slide

```pythonadata2.obs['low_quality_index']= best_match
adata2.obs['low_quality_index'] = adata2.obs['low_quality_index'].astype(float)```

```pythonadata2.obsm['spatial']```
*Output:*
```array([[-118.07557033,  391.48711306],
       [ -80.27430014,  356.96083362],
       [-130.37495374,  164.18395609],
       ...,
       [ -91.13328399,  258.15252187],
       [   7.42631645,  302.86275741],
       [ -98.93455418,  344.64032554]])```

```pythonsc.pl.spatial(adata2, color='low_quality_index', spot_size=3, title='Quality')```
*Output:*
```<Figure size 320x320 with 2 Axes>```

We use a Sankey diagram to show the correspondence between cell types at different stages of development

```pythonfig=Sankey_multi(adata_li=[adata1,adata2],
             prefix_li=['E11.5','E12.5'],
             matching_li=[matching],
                clusters='annotation',filter_num=10,
             node_opacity = 0.8,
             link_opacity = 0.2,
                layout=[800,500],
           font_size=12,
           font_color='Black',
           save_name=None,
           format='png',
           width=1200,
           height=1000,
           return_fig=True)
fig.show()```
*Output:*

```pythonfig.write_html("slat_sankey.html")```

## Focus on developing Kidney

We highlighted the “Kidney” cells in E12.5 and their aligned precursor cells in E11.5 in alignment results. Consistent with our biological priors, the precursors of the kidney are the mesonephros and the metanephros

Then we focus on another organ: ‘Ovary’, and found ovary only has single spatial origin. It is interesting that precursors of ovary are spatially close to the mesonephros (see Kidney part), because mammalian ovary originates from the regressed mesonephros.

```pythoncolor_dict1=dict(zip(adata1.obs['annotation'].cat.categories,
                    adata1.uns['annotation_colors'].tolist()))
adata1_df = pd.DataFrame({'index':range(embd0.shape[0]),
                          'x': adata1.obsm['spatial'][:,0],
                          'y': adata1.obsm['spatial'][:,1],
                          'celltype':adata1.obs['annotation'],
                         'color':adata1.obs['annotation'].map(color_dict1)
                         }
                        )
color_dict2=dict(zip(adata2.obs['annotation'].cat.categories,
                    adata2.uns['annotation_colors'].tolist()))
adata2_df = pd.DataFrame({'index':range(embd1.shape[0]),
                          'x': adata2.obsm['spatial'][:,0],
                          'y': adata2.obsm['spatial'][:,1],
                          'celltype':adata2.obs['annotation'],
                         'color':adata2.obs['annotation'].map(color_dict2)
                         }
                        )```

```pythonkidney_align = match_3D_celltype(adata1_df, adata2_df, matching, meta='celltype', 
                                 highlight_celltype = [['Urogenital ridge'],['Kidney','Ovary']],
                                 subsample_size=10000, highlight_line = ['blue'], scale_coordinate = True )
kidney_align.draw_3D(size= [6, 6], line_width =0.8, point_size=[0.6,0.6], hide_axis=True)```
*Output:*
```dataset1: 18 cell types; dataset2: 22 cell types; 
                    Total :29 celltypes; Overlap: 11 cell types 
                    Not overlap :[['Dermomyotome', 'Surface ectoderm', 'Sclerotome', 'Inner ear', 'Spinal cord', 'Head mesenchyme', 'Branchial arch', 'Sympathetic nerve', 'Ovary', 'Pancreas', 'Mucosal epithelium', 'Muscle', 'Mesentery', 'Choroid plexus', 'Kidney', 'Epidermis', 'Cartilage primordium', 'Mesothelium']]
Subsampled 10000 pairs from 10001
```
```<Figure size 480x480 with 1 Axes>```

We can get the lineage of the query's cells and mappings using the following function

```pythondef cal_matching_cell(target_adata,query_adata,matching,query_cell,clusters='annotation',):
    adata1_df = pd.DataFrame({'index':range(target_adata.shape[0]),
                          'x': target_adata.obsm['spatial'][:,0],
                          'y': target_adata.obsm['spatial'][:,1],
                          'celltype':target_adata.obs[clusters]})
    adata2_df = pd.DataFrame({'index':range(query_adata.shape[0]),
                              'x': query_adata.obsm['spatial'][:,0],
                              'y': query_adata.obsm['spatial'][:,1],
                              'celltype':query_adata.obs[clusters]})
    query_adata = target_adata[matching[1,adata2_df.loc[adata2_df.celltype==query_cell,'index'].values],:]
    #adata2_df['target_celltype'] = adata1_df.iloc[matching[1,:],:]['celltype'].to_list()
    #adata2_df['target_obs_names'] = adata1_df.iloc[matching[1,:],:].index.to_list()
    
    #query_obs=adata2_df.loc[adata2_df['celltype']==query_cell,'target_obs_names'].tolist()
    return query_adata
    
```

We find that maps mapped on 3D also show up well on 2D

```pythonquery_adata=cal_matching_cell(target_adata=adata1,
                              query_adata=adata2,
                              matching=matching,
                              query_cell='Kidney',clusters='annotation')
query_adata```
*Output:*
```View of AnnData object with n_obs × n_vars = 59 × 26854
    obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'annotation', 'Regulon - A1cf', 'kidney_c0', 'kidney_c1', 'kidney_c2', 'kidney_c3', 'kidney_anno', 'week'
    var: 'n_cells', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts'
    uns: 'Spatial_Net', 'annotation_colors', 'kidney_c0_colors', 'kidney_c1_colors', 'kidney_c2_colors', 'kidney_c3_colors', 'kidney_anno_colors'
    obsm: 'spatial'
    varm: 'PCs'```

```pythonadata1.obs['kidney_anno']=''
adata1.obs.loc[query_adata.obs.index,'kidney_anno']=query_adata.obs['annotation']```

```pythonsc.pl.spatial(adata1, color='kidney_anno', spot_size=3,
             palette=['#F5F5F5','#ff7f0e', 'green',])```
*Output:*
```<Figure size 320x320 with 1 Axes>```

We are concerned with Kidney lineage development, so we integrated the cells corresponding to the Kidney lineage on the two sections of E11 and E12, and then we could use the method of difference analysis to study the dynamic process of Kidney lineage development.

```pythonkidney_lineage_ad=sc.concat([query_adata,adata2[adata2.obs['annotation']=='Kidney']],merge='same')
kidney_lineage_ad=ov.pp.preprocess(kidney_lineage_ad,mode='shiftlog|pearson',n_HVGs=3000,target_sum=1e4)
kidney_lineage_ad.raw = kidney_lineage_ad
kidney_lineage_ad = kidney_lineage_ad[:, kidney_lineage_ad.var.highly_variable_features]
ov.pp.scale(kidney_lineage_ad)
ov.pp.pca(kidney_lineage_ad)
ov.pp.neighbors(kidney_lineage_ad,use_rep='scaled|original|X_pca',metric="cosine")
ov.utils.cluster(kidney_lineage_ad,method='leiden',resolution=1)
ov.pp.umap(kidney_lineage_ad)```
*Output:*
```Begin robust gene identification
After filtration, 14823/26436 genes are kept. Among 14823 genes, 14823 genes are robust.
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
Time to analyze data in cpu: 0.07181000709533691 seconds.
End of size normalization: shiftlog and HVGs selection pearson
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:06)
running Leiden clustering
    finished: found 6 clusters and added
    'leiden', the cluster labels (adata.obs, categorical) (0:00:00)
computing UMAP
    finished: added
    'X_umap', UMAP coordinates (adata.obsm) (0:00:00)
```

```pythonov.pl.embedding(kidney_lineage_ad,basis='X_umap',
               color=['annotation','week','leiden'],
               frameon='small')```
*Output:*
```<Figure size 1159.2x320 with 3 Axes>```

```python# Nphs1 https://www.nature.com/articles/s41467-021-22266-1
sc.pl.dotplot(kidney_lineage_ad,{'nephron progenitors':['Wnt9b','Osr1','Nphs1','Lhx1','Pax2','Pax8'],
                         'metanephric':['Eya1','Shisa3','Foxc1'], 
                         'kidney':['Wt1','Wnt4','Nr2f2','Dach1','Cd44']} ,
              'leiden',dendrogram=False,colorbar_title='Expression')```
*Output:*
```<Figure size 534.4x248 with 5 Axes>```

```pythonkidney_lineage_ad.obs['re_anno'] = 'Unknown'
kidney_lineage_ad.obs.loc[kidney_lineage_ad.obs.leiden.isin(['4']),'re_anno'] = 'Nephron progenitors (E11.5)'
kidney_lineage_ad.obs.loc[kidney_lineage_ad.obs.leiden.isin(['2','3','1','5']),'re_anno'] = 'Metanephron progenitors (E11.5)'
kidney_lineage_ad.obs.loc[kidney_lineage_ad.obs.leiden=='0','re_anno'] = 'Kidney (E12.5)'```

```python# kidney_all = kidney_all[kidney_all.obs.leiden!='3',:]
kidney_lineage_ad.obs.leiden = list(kidney_lineage_ad.obs.leiden)
ov.pl.embedding(kidney_lineage_ad,basis='X_umap',
               color=['annotation','re_anno'],
               frameon='small')```
*Output:*
```<Figure size 772.8x320 with 2 Axes>```

```pythonadata1.obs['kidney_anno']=''
adata1.obs.loc[kidney_lineage_ad[kidney_lineage_ad.obs['week']=='E11.5'].obs.index,'kidney_anno']=kidney_lineage_ad[kidney_lineage_ad.obs['week']=='E11.5'].obs['re_anno']```

```pythonimport matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
sc.pl.spatial(adata1, color='kidney_anno', spot_size=1.5,
             palette=['#F5F5F5','#ff7f0e', 'green',],show=False,ax=ax)```
*Output:*
```[<AxesSubplot: title={'center': 'kidney_anno'}, xlabel='spatial1', ylabel='spatial2'>]```
```<Figure size 640x640 with 1 Axes>```

We can also differentially analyse Kidney's developmental pedigree to find different marker genes, and we can analyse transcription factors and thus find the regulatory units involved.

```pythontest_adata=kidney_lineage_ad
dds=ov.bulk.pyDEG(test_adata.to_df(layer='lognorm').T)
dds.drop_duplicates_index()
print('... drop_duplicates_index success')
treatment_groups=test_adata.obs[test_adata.obs['week']=='E12.5'].index.tolist()
control_groups=test_adata.obs[test_adata.obs['week']=='E11.5'].index.tolist()
result=dds.deg_analysis(treatment_groups,control_groups,method='ttest')
# -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)```
*Output:*
```... drop_duplicates_index success
... Fold change threshold: 2.1686935424804688
```

```pythondds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)```
*Output:*
```<AxesSubplot: title={'center': 'DEG Analysis'}, xlabel='$log_{2}FC$', ylabel='$-log_{10}(qvalue)$'>```
```<Figure size 320x320 with 1 Axes>```

```pythonup_gene=dds.result.loc[dds.result['sig']=='up'].sort_values('qvalue')[:3].index.tolist()
down_gene=dds.result.loc[dds.result['sig']=='down'].sort_values('qvalue')[:3].index.tolist()
deg_gene=up_gene+down_gene```

```pythonsc.pl.dotplot(kidney_lineage_ad,deg_gene,
             groupby='re_anno')```
*Output:*
```<Figure size 297.6x200 with 4 Axes>```

In addition to analysing directly using differential expression, we can also look for weekly marker genes by considering weeks as categories.

```pythonsc.tl.dendrogram(kidney_lineage_ad,'re_anno',use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(kidney_lineage_ad, 're_anno', use_rep='scaled|original|X_pca',
                        method='t-test',use_raw=False,key_added='re_anno_ttest')
sc.pl.rank_genes_groups_dotplot(kidney_lineage_ad,groupby='re_anno',
                                cmap='RdBu_r',key='re_anno_ttest',
                                standard_scale='var',n_genes=3)```
*Output:*
```Storing dendrogram info using `.uns['dendrogram_re_anno']`
ranking genes
    finished: added to `.uns['re_anno_ttest']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
```
```<Figure size 450.4x200 with 6 Axes>```

