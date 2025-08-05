# t_deg_single
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_deg_single.ipynb*

# Differential expression and celltype analysis [All Cell]

1. Differential gene expression (DGE) analysis identifies genes that show statistically significant differences in expression levels across distinct cell populations or conditions. This analysis helps in identifying which cell types are most affected by a condition of interest such as a disease, and characterizing their functional signatures.

2. Differential Compositional analysis identifies Quantifies changes in the relative abundances of each cell type across conditions (e.g., case vs. control, time points, treatment groups). Reveals expansions or contractions of specific populations that may not be captured by gene-level analyses alone.

Here, we introduced `omicverse.single.DEG` and `omicverse.single.DCT` to performed these two analysis.

```pythonimport scanpy as sc
import pertpy as pt
import omicverse as ov
ov.plot_set()```
*Output:*
```🔬 Starting plot initialization...
🧬 Detecting CUDA devices…
✅ [GPU 0] NVIDIA TITAN Xp
    • Total memory: 11.9 GB
    • Compute capability: 6.1

   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

🔖 Version: 1.6.11   📚 Tutorials: https://omicverse.readthedocs.io/
✅ plot_set complete.

```

## Data Preprocess

The data we use in the following example comes from [Haber et al., 2017]. It contains samples from the small intestinal epithelium of mice with different conditions. We first load the raw cell-level data. The dataset contains gene expressions of 9842 cells. They are annotated with their sample identifier (`batch`), the condition of the subjects and the type of each cell (`cell_label`).

```pythonadata = pt.dt.haber_2017_regions()
adata```
*Output:*
```AnnData object with n_obs × n_vars = 9842 × 15215
    obs: 'batch', 'barcode', 'condition', 'cell_label'```

For our first example, we want to look at how the `Salmonella` infection influences the cell composition. Therefore, we create a subset of our compositional data that only contains the `healthy` and `Salmonella-infected` samples as a new data modality.

```python# Select control and salmonella data
adata = adata[
    adata.obs["condition"].isin(["Control", "Salmonella"])
].copy()
print(adata)```
*Output:*
```AnnData object with n_obs × n_vars = 5010 × 15215
    obs: 'batch', 'barcode', 'condition', 'cell_label'
```

```pythonadata.obs["condition"].unique()```
*Output:*
```['Control', 'Salmonella']
Categories (2, object): ['Control', 'Salmonella']```

## DEG with wilcoxon/t-test

In omicverse, we only need one function `ov.single.DEG` to perfrom all DEG analysis tasks. First, I will introduce the nonparametric Wilcoxon test and the t-test—two widely used methods for differential expression analysis. Due to their high computational efficiency, we applied them in our DEG analysis.

We need to set the `ctrl_group` and `test_group` to perform the analysis. Besides, we also need to define the celltype to be explode. If you set `celltype_group` is None, all the celltype will be calculated.

```pythondeg_obj=ov.single.DEG(
    adata,
    condition='condition',
    ctrl_group='Control',
    test_group='Salmonella',
    method='wilcoxon',
)
deg_obj.run(
    celltype_key='cell_label',
    celltype_group=['TA'],
)
```
*Output:*
```✅ Differential expression analysis initialized
📊 DEG analysis using wilcoxon method
📊 Condition: condition, Control group: Control, Test group: Salmonella
📊 Celltype key: cell_label, Celltype group: ['TA']
Total cells: 533 will be used for DEG analysis
normalizing counts per cell
    finished (0:00:00)
ranking genes
    finished: added to `.uns['rank_genes_groups']`
    'names', sorted np.recarray to be indexed by group ids
    'scores', sorted np.recarray to be indexed by group ids
    'logfoldchanges', sorted np.recarray to be indexed by group ids
    'pvals', sorted np.recarray to be indexed by group ids
    'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
✅ wilcoxon DEG analysis completed
```

```pythonres_wilcoxon=deg_obj.get_results()
res_wilcoxon.head()```
*Output:*
```          log2FC        pvalue          padj        qvalue      size  sig  \
Reg3b   5.958403  2.899264e-58  4.411230e-54  4.411230e-54  0.595840  sig   
Reg3g   4.425004  1.576382e-51  1.199233e-47  1.199233e-47  0.442500  sig   
Apoa1   2.589844  2.025723e-40  1.027379e-36  1.027379e-36  0.258984  sig   
Guca2b  2.380329  1.732794e-21  5.272893e-18  5.272893e-18  0.238033  sig   
Zg16    2.077089  2.553716e-21  6.475798e-18  6.475798e-18  0.207709  sig   

        -log(pvalue)  -log(qvalue)  
Reg3b      57.537712     53.355440  
Reg3g      50.802339     46.921097  
Apoa1      39.693420     35.988269  
Guca2b     20.761253     17.277951  
Zg16       20.592827     17.188707  ```

We can use `sc.pl.violin` to compare the gene expression between different condition.

```pythoncelltypes_li=['TA']

sc.pl.violin(
    adata[adata.obs['cell_label'].isin(celltypes_li)],
    keys=['Reg3b','Reg3g','Apoa1'],
    groupby='condition'
)```
*Output:*
```<Figure size 1116.72x320 with 3 Axes>```

## DEG with memento

memento is a python package for performing differential mean, variability, and correlation in single-cell RNA sequencing data.

Memento, an end-to-end method that implements a hierarchical model for estimating mean, residual variance, and gene correlation from scRNA-seq data and provides a statistical framework for hypothesis testing of these parameters.


```pythondeg_obj=ov.single.DEG(
    adata,
    condition='condition',
    ctrl_group='Control',
    test_group='Salmonella',
    method='memento-de',
)
deg_obj.run(
    celltype_key='cell_label',
    celltype_group=['TA'],
    capture_rate=0.07, 
    num_cpus=12,
    num_boot=5000
)
```
*Output:*
```✅ Differential expression analysis initialized
📊 DEG analysis using memento-de method
📊 Condition: condition, Control group: Control, Test group: Salmonella
📊 Celltype key: cell_label, Celltype group: ['TA']
Total cells: 533 will be used for DEG analysis
```
```[Parallel(n_jobs=12)]: Using backend LokyBackend with 12 concurrent workers.
[Parallel(n_jobs=12)]: Done  26 tasks      | elapsed:    8.7s
[Parallel(n_jobs=12)]: Done 176 tasks      | elapsed:   11.8s
[Parallel(n_jobs=12)]: Done 528 tasks      | elapsed:   15.3s
[Parallel(n_jobs=12)]: Done 1228 tasks      | elapsed:   20.2s
[Parallel(n_jobs=12)]: Done 2128 tasks      | elapsed:   26.2s
[Parallel(n_jobs=12)]: Done 3228 tasks      | elapsed:   34.0s
[Parallel(n_jobs=12)]: Done 4521 out of 4544 | elapsed:   42.8s remaining:    0.2s
[Parallel(n_jobs=12)]: Done 4544 out of 4544 | elapsed:   42.9s finished
```
```✅ memento-de DEG analysis completed
```

```pythonres_memento=deg_obj.get_results()
res_memento.query('dv_coef > 1 & de_coef > 0').sort_values('dv_pval').head(5)```
*Output:*
```         gene    tx   de_coef     de_se   de_pval   dv_coef     dv_se  \
443      Btf3  stim  0.838806  0.230108  0.000213  2.940342  0.580788   
397     Birc5  stim  0.631178  0.185273  0.000661  1.198823  0.313596   
3421  Serinc2  stim  0.030935  0.238055  0.826321  2.288608  0.511458   
917     Dcaf8  stim  0.037631  0.316794  0.846642  1.833301  0.500972   
536     Ccnb2  stim  0.132247  0.234056  0.565606  1.126923  0.323989   

           dv_pval  
443   4.519054e-07  
397   8.280704e-05  
3421  1.560657e-04  
917   2.687823e-04  
536   3.477457e-04  ```

```pythoncelltypes_li=['TA']

sc.pl.violin(
    adata[adata.obs['cell_label'].isin(celltypes_li)],
    keys=['Btf3','Serinc2','Birc5'],
    groupby='condition'
)```
*Output:*
```<Figure size 1116.72x320 with 3 Axes>```

## DCT with scCODA

In omicverse, we only need one function `ov.single.DCT` to perfrom all DCT analysis tasks. We included `scCODA` and `milo` to perform the celltype abundance analysis.

Besides, you can also perform the analysis with `pertpy`'s api. `dct_obj.model` will be helpful.

More tutorial could be found in https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/sccoda.html

```pythondct_obj=ov.single.DCT(
    adata,
    condition='condition',
    ctrl_group='Control',
    test_group='Salmonella',
    cell_type_key='cell_label',
    method='sccoda',
    sample_key='batch',
)```
*Output:*
```[94m•[0m Automatic reference selection! Reference cell type set to Endocrine
```

No-U-turn HMC sampling is then initiated by calling sccoda_model.`run_nuts`().

We can use `help(dct_obj.model.run_nuts)` to obtain the argument as input.

```pythondct_obj.run(
    num_samples=5000, #number of sampled values after burn-in.
    num_warmup=500, #Number of burn-in (warmup) samples.
)```
*Output:*
```sample: 100%|██████████| 5500/5500 [03:17<00:00, 27.79it/s, 127 steps of size 2.91e-02. acc. prob=0.76]
```

```pythonres=dct_obj.get_results()
res.head()```
*Output:*
```                                             Final Parameter  HDI 3%  HDI 97%  \
Covariate             Cell Type                                                 
conditionT.Salmonella Endocrine                     0.000000   0.000    0.000   
                      Enterocyte                    1.367619   0.864    1.886   
                      Enterocyte.Progenitor         0.000000  -0.393    0.647   
                      Goblet                        0.000000  -0.228    1.037   
                      Stem                          0.000000  -0.750    0.251   

                                                SD  Inclusion probability  \
Covariate             Cell Type                                             
conditionT.Salmonella Endocrine              0.000                 0.0000   
                      Enterocyte             0.262                 1.0000   
                      Enterocyte.Progenitor  0.163                 0.3072   
                      Goblet                 0.298                 0.4688   
                      Stem                   0.199                 0.3560   

                                             Expected Sample  log2-fold change  
Covariate             Cell Type                                                 
conditionT.Salmonella Endocrine                    25.782816         -0.495162  
                      Enterocyte                  325.488661          1.477896  
                      Enterocyte.Progenitor       100.454836         -0.495162  
                      Goblet                       43.541227         -0.495162  
                      Stem                        120.266274         -0.495162  ```

### Adjusting the False discovery rate

scCODA selects credible effects based on their inclusion probability. The cutoff between credible and non-credible effects depends on the desired false discovery rate (FDR). A smaller FDR value will produce more conservative results, but might miss some effects, while a larger FDR value selects more effects at the cost of a larger number of false discoveries.

The desired FDR level can be easily set after inference via sim_results.set_fdr(). Per default, the value is 0.05, but we recommend to increase it up to 0.2 if no effects are found at a more conservative level.

In our example, setting a desired FDR of 0.4 reveals small effects on Endocrine and Tuft cells. Keep in mind that we chose this value only for instructive purposes, since there are no credible effects beside Enterocytes at lower FDR levels. In practice, expecting 40% of all credible effects to be false-positives is usually not recommended.

```pythondct_obj.model.set_fdr(dct_obj.sccoda_data, 
                      modality_key="coda", 
                      est_fdr=0.6)
res=dct_obj.get_results()
res.sort_values('Final Parameter',ascending=False).head()```
*Output:*
```                                             Final Parameter  HDI 3%  HDI 97%  \
Covariate             Cell Type                                                 
conditionT.Salmonella Enterocyte                    1.367619   0.864    1.886   
                      Goblet                        0.353104  -0.228    1.037   
                      Enterocyte.Progenitor         0.115878  -0.393    0.647   
                      TA.Early                      0.008501  -0.408    0.486   
                      Endocrine                     0.000000   0.000    0.000   

                                                SD  Inclusion probability  \
Covariate             Cell Type                                             
conditionT.Salmonella Enterocyte             0.262                 1.0000   
                      Goblet                 0.298                 0.4688   
                      Enterocyte.Progenitor  0.163                 0.3072   
                      TA.Early               0.121                 0.2748   
                      Endocrine              0.000                 0.0000   

                                             Expected Sample  log2-fold change  
Covariate             Cell Type                                                 
conditionT.Salmonella Enterocyte                  327.583710          1.487152  
                      Goblet                       62.378974          0.023516  
                      Enterocyte.Progenitor       113.522620         -0.318729  
                      TA.Early                    142.968607         -0.473642  
                      Endocrine                    25.948771         -0.485906  ```

```python#ov.plot_set()
dct_obj.model.plot_boxplots(
    dct_obj.sccoda_data, 
    modality_key="coda", 
    feature_name="condition", 
    add_dots=True,
    figsize=(4,4),
    dpi=80,
)
ov.plt.show()```
*Output:*
```<Figure size 320x320 with 1 Axes>```

```pythondct_obj.model.summary(
    dct_obj.sccoda_data, 
    modality_key="coda"
)```
*Output:*
```[3m                                          Compositional Analysis summary                                           [0m
┌──────────────────────────────────────────────┬──────────────────────────────────────────────────────────────────┐
│[1m [0m[1mName                                        [0m[1m [0m│[1m [0m[1mValue                                                           [0m[1m [0m│
├──────────────────────────────────────────────┼──────────────────────────────────────────────────────────────────┤
│[36m [0m[36mData                                        [0m[36m [0m│ Data: [1;36m6[0m samples, [1;36m8[0m cell types                                    │
│[36m [0m[36mReference cell type                         [0m[36m [0m│ Endocrine                                                        │
│[36m [0m[36mFormula                                     [0m[36m [0m│ condition                                                        │
└──────────────────────────────────────────────┴──────────────────────────────────────────────────────────────────┘
```
```┌─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│[1m [0m[1mIntercepts                                                                                                     [0m[1m [0m│
├─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│                        Final Parameter  Expected Sample                                                         │
│ Cell Type                                                                                                       │
│ Endocrine                  [1;36m1.183[0m            [1;36m36.340[0m                                                              │
│ Enterocyte                 [1;36m2.351[0m           [1;36m116.854[0m                                                              │
│ Enterocyte.Progenitor      [1;36m2.543[0m           [1;36m141.589[0m                                                              │
│ Goblet                     [1;36m1.707[0m            [1;36m61.370[0m                                                              │
│ Stem                       [1;36m2.723[0m           [1;36m169.513[0m                                                              │
│ TA                         [1;36m2.126[0m            [1;36m93.310[0m                                                              │
│ TA.Early                   [1;36m2.881[0m           [1;36m198.528[0m                                                              │
│ Tuft                       [1;36m0.452[0m            [1;36m17.495[0m                                                              │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```
```┌─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│[1m [0m[1mEffects                                                                                                        [0m[1m [0m│
├─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│                                              Final Parameter  Expected Sample  log2-fold change                 │
│ Covariate             Cell Type                                                                                 │
│ conditionT.Salmonella Endocrine                   [1;36m0.000[0m           [1;36m25.949[0m            [1;36m-0.486[0m                      │
│                       Enterocyte                  [1;36m1.368[0m          [1;36m327.584[0m             [1;36m1.487[0m                      │
│                       Enterocyte.Progenitor       [1;36m0.116[0m          [1;36m113.523[0m            [1;36m-0.319[0m                      │
│                       Goblet                      [1;36m0.353[0m           [1;36m62.379[0m             [1;36m0.024[0m                      │
│                       Stem                       [1;36m-0.229[0m           [1;36m96.260[0m            [1;36m-0.816[0m                      │
│                       TA                         [1;36m-0.212[0m           [1;36m53.916[0m            [1;36m-0.791[0m                      │
│                       TA.Early                    [1;36m0.009[0m          [1;36m142.969[0m            [1;36m-0.474[0m                      │
│                       Tuft                       [1;36m-0.006[0m           [1;36m12.421[0m            [1;36m-0.494[0m                      │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

## DCT with milo

Many biological conditions (disease, development, genetic KOs) can induce shifts in cell composition, where cells of a given state become enriched or depleted in response to a perturbation. With differential abundance analysis, we quantify consistent changes in cell composition across replicate samples. While differential abundance analysis can be performed on cell type clusters, it’s not always possible or practical to use precisely labeled clusters, especially when we are interested in studying transitional states, such as during developmental processes, or when we expect only a subpopulation of a cell type to be affected by the condition of interest. Milo is a method to detect compositional changes occurring in smaller subpopulations of cells, defined as neighbourhoods over the k-nearest neighbor (KNN) graph of cell-cell similarities.

Tutorials could be found in https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/milo.html

### Build KNN graph

We can use omicverse functions to build a KNN graph. We set the dimensionality and value for k to use in subsequent steps.

Here the value of k indicates the smallest possible size of neighbourhood in which we will quantify differential abundance (i.e. with k=50 the smallest neighbourhood will have 50 cells). Depending on the number of samples, you might want to use a high value of k for neighbourhood analysis, to have sufficient power to estimate abundance fold-changes. Since here we have data from > 100 patients, we set k=150 to have on average more than one cell per donor in each neighbourhood.

```pythonov.settings.cpu_gpu_mixed_init()```
*Output:*
```CPU-GPU mixed mode activated
```

```pythonadata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,
                       target_sum=50*1e4)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
ov.single.batch_correction(adata,batch_key='batch',
                                        methods='harmony',n_pcs=50)
ov.pp.neighbors(adata, n_neighbors=15, n_pcs=50,
               use_rep='X_harmony')
ov.pp.umap(adata)```
*Output:*
```Begin robust gene identification
After filtration, 15215/15215 genes are kept.     Among 15215 genes, 15215 genes are robust.
End of robust gene identification.
Begin size normalization: shiftlog and HVGs selection pearson
normalizing counts per cell. The following highly-expressed genes are not considered during normalization factor computation:
['Defa24', 'Fabp6', 'Gcg', 'Gip', 'Nts', 'Reg3b', 'Reg4', 'Sct', 'Spink4', 'Sst', 'Tff3', 'Zg16']
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
Time to analyze data in cpu: 7.434457540512085 seconds.
End of size normalization: shiftlog and HVGs selection pearson
...Begin using harmony to correct batch effect
🚀 Using GPU to calculate PCA...
📊 [GPU 0] [92m[90m------------------------------[0m 5/12288 MiB (0.0%)
computing PCA🔍
    with n_comps=50
    finished✅ (0:00:01)
```
```2025-05-24 03:28:39,033 - harmonypy - INFO - Computing initial centroids with sklearn.KMeans...
2025-05-24 03:28:40,914 - harmonypy - INFO - sklearn.KMeans initialization complete.
2025-05-24 03:28:40,966 - harmonypy - INFO - Iteration 1 of 10
2025-05-24 03:28:42,555 - harmonypy - INFO - Iteration 2 of 10
2025-05-24 03:28:44,190 - harmonypy - INFO - Iteration 3 of 10
2025-05-24 03:28:45,633 - harmonypy - INFO - Iteration 4 of 10
2025-05-24 03:28:46,905 - harmonypy - INFO - Iteration 5 of 10
2025-05-24 03:28:48,659 - harmonypy - INFO - Iteration 6 of 10
2025-05-24 03:28:49,316 - harmonypy - INFO - Converged after 6 iterations
```
```🖥️ Using Scanpy CPU to calculate neighbors...
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:02)
🔍 [2025-05-24 03:28:52] Running UMAP in 'cpu-gpu-mixed' mode...
🚀 Using torch GPU to calculate UMAP...
📊 [GPU 0] [92m[90m------------------------------[0m 5/12288 MiB (0.0%)
computing UMAP🚀
    finished ✅: added
    'X_umap', UMAP coordinates (adata.obsm)
    'umap', UMAP parameters (adata.uns) (0:00:05)
✅ UMAP completed successfully.
```

```pythonov.pl.embedding(
    adata,
    basis='X_umap',
    color=['batch','cell_label'],
)```
*Output:*
```<Figure size 772.8x320 with 2 Axes>```

### Differential abundance testing with GLM

Similar to the scCODA approach, we use `omicverse.single.DCT` to perform differential cell–abundance analysis, except that here we set the `method` to `milo`. In future releases, we will gradually remove the rpy2‐based edgeR dependency so that milo analyses can run entirely in a native Python environment.


```pythondct_obj=ov.single.DCT(
    adata,
    condition='condition',
    ctrl_group='Control',
    test_group='Salmonella',
    cell_type_key='cell_label',
    method='milo',
    sample_key='batch',
    use_rep='X_harmony'
)```
*Output:*
```computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:00)
```

```pythondct_obj.run()```

We can start inspecting the results of our DA analysis from a couple of standard diagnostic plots.

```pythonimport matplotlib.pyplot as plt
old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [6, 3]
plt.subplot(1, 2, 1)
plt.hist(dct_obj.mdata["milo"].var.PValue, bins=50)
plt.xlabel("P-Vals")
plt.subplot(1, 2, 2)
plt.plot(dct_obj.mdata["milo"].var.logFC, -ov.np.log10(dct_obj.mdata["milo"].var.SpatialFDR), ".")
plt.xlabel("log-Fold Change")
plt.ylabel("- log10(Spatial FDR)")
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize```
*Output:*
```<Figure size 480x240 with 2 Axes>```

We can see that for the majority of neighbourhoods, almost all cells have the same cell type label. We can rename neighbourhoods where less than 60% of the cells have the top label as “Mixed”

```pythonplt.hist(dct_obj.mdata["milo"].var["nhood_annotation_frac"], bins=30)
plt.xlabel("celltype fraction")```
*Output:*
```Text(0.5, 0, 'celltype fraction')```
```<Figure size 480x240 with 1 Axes>```

```pythonres=dct_obj.get_results(mix_threshold=0.6)
res.head()```
*Output:*
```                                        index_cell  kth_distance     logFC  \
0  B1_AATAAGCTAGAGAT_Control_Enterocyte.Progenitor     11.069503 -0.721935   
1                 B1_ACATGGTGCCGTTC_Control_Goblet     15.162609  0.157505   
2                   B1_AGCGGCACCAGAAA_Control_Stem      9.486914 -0.668452   
3                     B1_ATTCAAGATTCACT_Control_TA      9.663015 -0.571189   
4                     B1_CCTAGAGAATTCGG_Control_TA      9.019321 -0.288963   

      logCPM         F    PValue       FDR  SpatialFDR  Nhood_size  \
0  12.275702  9.473558  0.017160  0.091587    0.086902       333.0   
1  11.849096  0.154712  0.705412  0.791629    0.796892       242.0   
2  13.153671  7.922610  0.025113  0.092657    0.088924       599.0   
3  12.738083  9.086466  0.018799  0.091587    0.087399       460.0   
4  13.097860  2.213955  0.179021  0.309080    0.306139       599.0   

  nhood_annotation  nhood_annotation_frac  
0            Mixed               0.462462  
1           Goblet               0.991736  
2             Stem               0.796327  
3            Mixed               0.423913  
4            Mixed               0.400668  ```

### Visualization

This is my favorite Milo plot. First, we create a `color_dict` to specify the cell colors we want to visualize.


```pythoncolor_dict=dict(zip(
    adata.obs['cell_label'].cat.categories,
    ov.pl.green_color[:4]+ov.pl.purple_color
))
color_dict['Mixed']='#c2c2c2'
fig, ax = plt.subplots(figsize=(3, 3))
ov.pl.embedding(
    adata,
    basis='X_umap',
    color=['cell_label'],
    palette=color_dict,
    ax=ax,
    #fig_size=(3,3)
)```
*Output:*
```<Figure size 240x240 with 1 Axes>```

```python#fig, ax = plt.subplots(figsize=(3, 4))
dct_obj.model.plot_da_beeswarm(
    dct_obj.mdata, 
    alpha=0.1,
    return_fig=True,
    palette=color_dict,
)
ov.plt.xticks(fontsize=12)
ov.plt.yticks(fontsize=12)
ov.plt.ylabel('T/NK Cells',fontsize=13)
ov.plt.text(-1,-0.75,'Enriched in\nControl',fontsize=13)
ov.plt.text(1,-0.75,'Enriched in\nSalmonella',fontsize=13)
#fig```
*Output:*
```Text(1, -0.75, 'Enriched in\nSalmonella')```
```<Figure size 480x240 with 1 Axes>```

