# t_mofa_glue
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_mofa_glue.ipynb*

# Multi omics analysis by MOFA and GLUE
MOFA is a factor analysis model that provides a general framework for the integration of multi-omic data sets in an unsupervised fashion.

Most of the time, however, we did not get paired cells in the multi-omics analysis. Here, we can pair cells using GLUE, a dimensionality reduction algorithm that can integrate different histological layers, and it can efficiently merge data from different histological layers.

This tutorial focuses on how to perform mofa in multi-omics using GLUE.

Paper: [MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1) and [Multi-omics single-cell data integration and regulatory inference with graph-linked embedding](https://www.nature.com/articles/s41587-022-01284-4)

Code: https://github.com/bioFAM/mofapy2 and https://github.com/gao-lab/GLUE

Colab_Reproducibility：https://colab.research.google.com/drive/1zlakFf20IoBdyuOQDocwFQHu8XOVizRL?usp=sharing

We used the result anndata object `rna-emb.h5ad` and `atac.emb.h5ad` from [GLUE'tutorial](https://scglue.readthedocs.io/en/latest/training.html)

```pythonimport omicverse as ov
ov.utils.ov_plot_set()
```
*Output:*
```2023-05-19 03:10:24.139649: I tensorflow/core/platform/cpu_feature_guard.cc:193] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  AVX2 FMA
To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.
2023-05-19 03:10:24.226148: W tensorflow/compiler/xla/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libcudart.so.11.0'; dlerror: libcudart.so.11.0: cannot open shared object file: No such file or directory
2023-05-19 03:10:24.226183: I tensorflow/compiler/xla/stream_executor/cuda/cudart_stub.cc:29] Ignore above cudart dlerror if you do not have a GPU set up on your machine.
2023-05-19 03:10:24.729457: W tensorflow/compiler/xla/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libnvinfer.so.7'; dlerror: libnvinfer.so.7: cannot open shared object file: No such file or directory
2023-05-19 03:10:24.729553: W tensorflow/compiler/xla/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libnvinfer_plugin.so.7'; dlerror: libnvinfer_plugin.so.7: cannot open shared object file: No such file or directory
2023-05-19 03:10:24.729559: W tensorflow/compiler/tf2tensorrt/utils/py_utils.cc:38] TF-TRT Warning: Cannot dlopen some TensorRT libraries. If you would like to use Nvidia GPU with TensorRT, please make sure the missing libraries mentioned above are installed properly.
```
```/mnt/data/env/pyomic/lib/python3.8/site-packages/phate/__init__.py
```

## Load the data

We use `ov.utils.read` to read the `h5ad` files

```pythonrna=ov.utils.read("chen_rna-emb.h5ad")
atac=ov.utils.read("chen_atac-emb.h5ad")```

## Pair the omics 

Each cell in our rna and atac data has a feature vector, X_glue, based on which we can calculate the Pearson correlation coefficient to perform cell matching.

```pythonpair_obj=ov.single.GLUE_pair(rna,atac)
pair_obj.correlation()```
*Output:*
```......Extract GLUE layer from obs
......Prepare for pair
......Start to calculate the Pearson coef
```
```Now Pearson block is 1/2: 100%|██████████████████████████████████████████████████████████████████████████████████| 2/2 [00:01<00:00,  1.45it/s]
Now rna_index is 4999/9999, all is 9190: 100%|█████████████████████████████████████████████████████████████| 5000/5000 [01:00<00:00, 82.59it/s]
```
```Now epoch is 0, 5000/9190
```
```Now Pearson block is 1/2: 100%|██████████████████████████████████████████████████████████████████████████████████| 2/2 [00:01<00:00,  1.70it/s]
Now rna_index is 9189/13379, all is 9190: 100%|████████████████████████████████████████████████████████████| 4190/4190 [00:58<00:00, 71.04it/s]```
```Now epoch is 1, 9190/9190
```
```
```

We counted the top 50 highly correlated cells in another histology layer for each cell in one of the histology layers to avoid missing data due to one cell being highly correlated with multiple cells. The default minimum threshold for high correlation is 0.9. We can obtain more paired cells by increasing the depth, but note that increasing the depth may lead to higher errors in cell matching

```pythonres_pair=pair_obj.find_neighbor_cell(depth=20)
res_pair.to_csv('models/chen_pair_res.csv')```
*Output:*
```Now depth is 19/20: 100%|██████████████████████████████████████████████████████████████████████████████████████| 20/20 [00:01<00:00, 19.49it/s]
```

We filter to get paired cells

```pythonrna1=rna[res_pair['omic_1']]
atac1=atac[res_pair['omic_2']]
rna1.obs.index=res_pair.index
atac1.obs.index=res_pair.index
rna1,atac1```
*Output:*
```(View of AnnData object with n_obs × n_vars = 5802 × 28930
     obs: 'domain', 'cell_type', 'balancing_weight'
     var: 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'mean', 'std', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts', 'gene_id', 'gene_type', 'mgi_id', 'havana_gene', 'tag'
     uns: '__scglue__', 'cell_type_colors', 'hvg', 'log1p', 'neighbors', 'pca', 'umap'
     obsm: 'X_glue', 'X_pca', 'X_umap'
     varm: 'PCs', 'X_glue'
     layers: 'counts'
     obsp: 'connectivities', 'distances',
 View of AnnData object with n_obs × n_vars = 5802 × 241757
     obs: 'domain', 'cell_type', 'balancing_weight'
     var: 'chrom', 'chromStart', 'chromEnd', 'highly_variable'
     uns: '__scglue__', 'cell_type_colors', 'neighbors', 'umap'
     obsm: 'X_glue', 'X_lsi', 'X_umap'
     varm: 'X_glue'
     obsp: 'connectivities', 'distances')```

We can use mudata to store the multi-omics 

```pythonfrom mudata import MuData

mdata = MuData({'rna': rna1, 'atac': atac1})
mdata```
*Output:*
```MuData object with n_obs × n_vars = 5802 × 270687
  var:	'highly_variable', 'chrom', 'chromStart', 'chromEnd'
  2 modalities
    rna:	5802 x 28930
      obs:	'domain', 'cell_type', 'balancing_weight'
      var:	'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'mean', 'std', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts', 'gene_id', 'gene_type', 'mgi_id', 'havana_gene', 'tag'
      uns:	'__scglue__', 'cell_type_colors', 'hvg', 'log1p', 'neighbors', 'pca', 'umap'
      obsm:	'X_glue', 'X_pca', 'X_umap'
      varm:	'PCs', 'X_glue'
      layers:	'counts'
      obsp:	'connectivities', 'distances'
    atac:	5802 x 241757
      obs:	'domain', 'cell_type', 'balancing_weight'
      var:	'chrom', 'chromStart', 'chromEnd', 'highly_variable'
      uns:	'__scglue__', 'cell_type_colors', 'neighbors', 'umap'
      obsm:	'X_glue', 'X_lsi', 'X_umap'
      varm:	'X_glue'
      obsp:	'connectivities', 'distances'```

```pythonmdata.write("chen_mu.h5mu",compression='gzip')```

## MOFA prepare

In the MOFA analysis, we only need to use highly variable genes, for which we perform one filter

```pythonrna1=mdata['rna']
rna1=rna1[:,rna1.var['highly_variable']==True]
atac1=mdata['atac']
atac1=atac1[:,atac1.var['highly_variable']==True]
rna1.obs.index=res_pair.index
atac1.obs.index=res_pair.index```

```pythonimport random
random_obs_index=random.sample(list(rna1.obs.index),5000)```

```pythonfrom sklearn.metrics import adjusted_rand_score as ari
ari_random=ari(rna1[random_obs_index].obs['cell_type'], atac1[random_obs_index].obs['cell_type'])
ari_raw=ari(rna1.obs['cell_type'], atac1.obs['cell_type'])
print('raw ari:{}, random ari:{}'.format(ari_raw,ari_random))```
*Output:*
```raw ari:0.7497091922584931, random ari:0.7491244560264247
```

```python#rna1=rna1[random_obs_index]
#atac1=atac1[random_obs_index]```

## MOFA analysis

In this part, we construct a model of mofa by scRNA-seq and scATAC-seq

```pythontest_mofa=ov.single.pyMOFA(omics=[rna1,atac1],
                             omics_name=['RNA','ATAC'])```

```pythontest_mofa.mofa_preprocess()
test_mofa.mofa_run(outfile='models/chen_rna_atac.hdf5')```
*Output:*
```
        #########################################################
        ###           __  __  ____  ______                    ### 
        ###          |  \/  |/ __ \|  ____/\    _             ### 
        ###          | \  / | |  | | |__ /  \ _| |_           ### 
        ###          | |\/| | |  | |  __/ /\ \_   _|          ###
        ###          | |  | | |__| | | / ____ \|_|            ###
        ###          |_|  |_|\____/|_|/_/    \_\              ###
        ###                                                   ### 
        ######################################################### 
       
 
        
Groups names not provided, using default naming convention:
- group1, group2, ..., groupG

Successfully loaded view='RNA' group='group0' with N=5802 samples and D=2000 features...
Successfully loaded view='ATAC' group='group0' with N=5802 samples and D=25488 features...


Warning: 15 features(s) in view 0 have zero variance, consider removing them before training the model...

Warning: 204 features(s) in view 1 have zero variance, consider removing them before training the model...

Model options:
- Automatic Relevance Determination prior on the factors: True
- Automatic Relevance Determination prior on the weights: True
- Spike-and-slab prior on the factors: False
- Spike-and-slab prior on the weights: True
Likelihoods:
- View 0 (RNA): gaussian
- View 1 (ATAC): gaussian



GPU mode is activated



######################################
## Training the model with seed 112 ##
######################################


ELBO before training: -2115272814.25 

Iteration 1: time=19.35, ELBO=191213094.17, deltaELBO=2306485908.422 (109.03964221%), Factors=19
Iteration 2: time=16.58, ELBO=196818457.51, deltaELBO=5605363.336 (0.26499482%), Factors=18
Iteration 3: time=15.64, ELBO=196891451.69, deltaELBO=72994.187 (0.00345082%), Factors=17
Iteration 4: time=14.82, ELBO=196995605.43, deltaELBO=104153.733 (0.00492389%), Factors=16
Iteration 5: time=14.06, ELBO=197054469.84, deltaELBO=58864.412 (0.00278283%), Factors=15
Iteration 6: time=13.33, ELBO=197108040.24, deltaELBO=53570.400 (0.00253255%), Factors=14
Iteration 7: time=12.57, ELBO=197166806.53, deltaELBO=58766.287 (0.00277819%), Factors=13
Iteration 8: time=11.82, ELBO=197221227.29, deltaELBO=54420.761 (0.00257275%), Factors=12
Iteration 9: time=11.02, ELBO=197274435.32, deltaELBO=53208.035 (0.00251542%), Factors=11
Iteration 10: time=10.27, ELBO=197321461.12, deltaELBO=47025.800 (0.00222316%), Factors=10
Iteration 11: time=9.44, ELBO=197307721.09, deltaELBO=-13740.036 (0.00064956%), Factors=9
Warning, lower bound is decreasing...
Iteration 12: time=8.66, ELBO=197345997.01, deltaELBO=38275.920 (0.00180950%), Factors=8
Iteration 13: time=8.12, ELBO=197354473.50, deltaELBO=8476.491 (0.00040073%), Factors=8
Iteration 14: time=8.14, ELBO=197361554.07, deltaELBO=7080.568 (0.00033474%), Factors=8

Converged!



#######################
## Training finished ##
#######################


Saving model in models/chen_rna_atac.hdf5...
```

## MOFA Visualization

In this part, we provide a series of function to visualize the result of mofa.

```pythonpymofa_obj=ov.single.pyMOFAART(model_path='models/chen_rna_atac.hdf5')```

```pythonpymofa_obj.get_factors(rna1)
rna1```
*Output:*
```......Add factors to adata and store to adata.obsm["X_mofa"]
```
```AnnData object with n_obs × n_vars = 5802 × 2000
    obs: 'domain', 'cell_type', 'balancing_weight', 'factor1', 'factor2', 'factor3', 'factor4', 'factor5', 'factor6', 'factor7', 'factor8'
    var: 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'mean', 'std', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts', 'gene_id', 'gene_type', 'mgi_id', 'havana_gene', 'tag'
    uns: '__scglue__', 'cell_type_colors', 'hvg', 'log1p', 'neighbors', 'pca', 'umap'
    obsm: 'X_glue', 'X_pca', 'X_umap', 'X_mofa'
    varm: 'PCs', 'X_glue'
    layers: 'counts'
    obsp: 'connectivities', 'distances'```

### Visualize the varience of each view

```pythonpymofa_obj.plot_r2()```
*Output:*
```(<Figure size 160x240 with 2 Axes>,
 <AxesSubplot: title={'center': 'Varience'}, xlabel='View', ylabel='Factor'>)```
```<Figure size 160x240 with 2 Axes>```

```pythonpymofa_obj.get_r2()```
*Output:*
```        RNA      ATAC
0  1.531975  0.020675
1 -0.036567  1.056309
2  0.751953  0.039987
3  0.060023  0.674215
4  0.364469  0.251295
5  0.141428  0.051458
6  0.114254  0.053613
7  0.119133  0.034315```

### Visualize the correlation between factor and celltype

```pythonpymofa_obj.plot_cor(rna1,'cell_type',figsize=(4,6))```
*Output:*
```(<Figure size 320x480 with 2 Axes>,
 <AxesSubplot: title={'center': 'Correlation'}, xlabel='Factor', ylabel='cell_type'>)```
```<Figure size 320x480 with 2 Axes>```

```pythonpymofa_obj.get_cor(rna1,'cell_type')```
*Output:*
```               factor1     factor2     factor3     factor4     factor5  \
InV           0.665011   20.657980    9.248108   17.089494   13.422602   
Clau          2.826121    5.808516    0.811940    0.268856    0.263320   
E5Parm1       7.532620   17.957307    1.763991   27.139456   84.494624   
OliM          6.513966   17.629268  500.000000   24.141743   10.914886   
E4Thsd7a     16.938637  237.595015   10.175490   13.044958   45.402276   
Mic          18.721245    0.863540   97.697177    1.925362    1.950863   
InP           3.293589   38.579983   26.957391   25.219145   15.585743   
OliI          2.847773    3.391799  141.002571    4.200228    2.875140   
E5Tshz2       2.709568    7.122963    1.657839   18.795893    2.700763   
E5Galnt14    14.255253    4.202316    9.685683   15.230670   46.138543   
OPC          23.475918   18.455308   32.199756   14.478467   10.162953   
Peri          0.658162    1.418508    0.884916    1.978547    3.519710   
InN           1.333910   14.820248   14.356910   10.069938    6.115415   
E3Rmst       10.885583   20.024174   11.668498    4.168142   17.267523   
E3Rorb       36.750410   38.807227   48.466332   81.247901  285.221588   
E4Il1rapl2   18.191082   80.250141   11.056312    1.086077    6.790889   
E5Sulf1       5.240678    1.483262    3.790759    6.827697  135.848121   
E6Tle4       29.553202   39.878014   32.630424  104.177836  500.000000   
InS           3.389144   31.927599    6.973416   17.450485    6.059428   
Endo          0.638567    1.657187    0.245722    1.829868    0.430014   
E2Rasgrf2    68.612084   10.432200   43.132372  500.000000   99.122213   
Ast         500.000000   72.170444    9.236814   67.929809   42.193385   

               factor6     factor7     factor8  
InV          32.132796   24.124465   44.255535  
Clau          0.261479    0.105971    0.908162  
E5Parm1       0.502745    0.651795    0.132044  
OliM         23.621465   25.395877    4.926897  
E4Thsd7a      4.448403    1.102767    5.383185  
Mic           1.756524    0.071363    0.533441  
InP         161.657573  250.170196  141.766482  
OliI        228.418422    0.743355    6.744291  
E5Tshz2       0.392799    0.462671    2.154666  
E5Galnt14     0.958545    1.906899    3.441217  
OPC          82.261963    5.433134    3.726258  
Peri          0.687741    0.444232    0.211925  
InN          45.999768   31.705707   51.326725  
E3Rmst        0.959733    4.561186   10.930971  
E3Rorb       14.008255   36.506786   69.163363  
E4Il1rapl2    2.651834    1.754884    5.952761  
E5Sulf1       1.140704    0.166260    0.283728  
E6Tle4       16.817520    5.015146    0.432723  
InS          81.276278   74.919020  125.230026  
Endo          0.024523    0.465449    0.896513  
E2Rasgrf2    23.616769    7.072289    4.607651  
Ast           5.108831    0.439351    4.559738  ```

```pythonpymofa_obj.plot_factor(rna1,'cell_type','Ast',figsize=(3,3),
                    factor1=1,factor2=3,)```
*Output:*
```<Figure size 240x240 with 1 Axes>```
```(<Figure size 240x240 with 1 Axes>,
 <AxesSubplot: title={'center': 'Ast'}, xlabel='X_mofa1', ylabel='X_mofa3'>)```

### Visualize the factor in UMAP

To visualize the GLUE’s learned embeddings, we use the pymde package wrapperin scvi-tools. This is an alternative to UMAP that is GPU-accelerated.

You can use `sc.tl.umap` insteaded.

```pythonfrom scvi.model.utils import mde
import scanpy as sc
sc.pp.neighbors(rna1, use_rep="X_glue", metric="cosine")
rna1.obsm["X_mde"] = mde(rna1.obsm["X_glue"])```
*Output:*
```computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:12)
```

```pythonsc.pl.embedding(
    rna1,
    basis="X_mde",
    color=["factor1","factor3","cell_type"],
    frameon=False,
    ncols=3,
    #palette=ov.utils.pyomic_palette(),
    show=False,
    cmap='Greens',
    vmin=0,
)
```
*Output:*
```[<AxesSubplot: title={'center': 'factor1'}, xlabel='X_mde1', ylabel='X_mde2'>,
 <AxesSubplot: title={'center': 'factor3'}, xlabel='X_mde1', ylabel='X_mde2'>,
 <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde1', ylabel='X_mde2'>]```
```<Figure size 1159.2x320 with 5 Axes>```

### Weights ranked
A visualization of factor weights familiar to MOFA and MOFA+ users is implemented with some modifications in `plot_weight_gene_d1`, `plot_weight_gene_d2`, and `plot_weights`.

```pythonpymofa_obj.plot_weight_gene_d1(view='RNA',factor1=1,factor2=3,)```
*Output:*
```(<Figure size 240x240 with 1 Axes>,
 <AxesSubplot: xlabel='factor_1', ylabel='factor_3'>)```
```<Figure size 240x240 with 1 Axes>```

```pythonpymofa_obj.plot_weights(view='RNA',factor=1,
                        ascending=False)```
*Output:*
```(<Figure size 240x320 with 1 Axes>,
 <AxesSubplot: title={'center': 'factor_1'}, xlabel='Feature rank', ylabel='Weight'>)```
```<Figure size 240x320 with 1 Axes>```

### Weights heatmap

While trying to annotate factors, a global overview of top features defining them could be helpful.

```pythonpymofa_obj.plot_top_feature_heatmap(view='RNA')```
*Output:*
```ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
WARNING: dendrogram data not found (using key=dendrogram_Factor). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
WARNING: You’re trying to run this on 2000 dimensions of `.X`, if you really want this, set `use_rep='X'`.
         Falling back to preprocessing with `sc.pp.pca` and default params.
computing PCA
    with n_comps=15
    finished (0:00:00)
Storing dendrogram info using `.uns['dendrogram_Factor']`
```
```{'mainplot_ax': <AxesSubplot: >,
 'group_extra_ax': <AxesSubplot: >,
 'gene_group_ax': <AxesSubplot: >,
 'color_legend_ax': <AxesSubplot: title={'center': 'Mean expression\nin group'}>}```
```<Figure size 894.4x304 with 5 Axes>```

