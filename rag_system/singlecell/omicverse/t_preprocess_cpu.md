# Preprocessing the data of scRNA-seq with omicverse[CPU-GPU-mixed]

The count table, a numeric matrix of genes × cells, is the basic input data structure in the analysis of single-cell RNA-sequencing data. A common preprocessing step is to adjust the counts for variable sampling efficiency and to transform them so that the variance is similar across the dynamic range.

Suitable methods to preprocess the scRNA-seq is important. Here, we introduce some preprocessing step to help researchers can perform downstream analysis easyier.

User can compare our tutorial with [scanpy'tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) to learn how to use omicverse well

Colab_Reproducibility：https://colab.research.google.com/drive/1DXLSls_ppgJmAaZTUvqazNC_E7EDCxUe?usp=sharing

```python
import scanpy as sc
import omicverse as ov
ov.plot_set(font_path='Arial')

# Enable auto-reload for development
%load_ext autoreload
%autoreload 2
```

    🔬 Starting plot initialization...
    Using already downloaded Arial font from: /tmp/omicverse_arial.ttf
    Registered as: Arial
    🧬 Detecting CUDA devices…
    ✅ [GPU 0] Tesla V100-SXM2-16GB
        • Total memory: 15.8 GB
        • Compute capability: 7.0

       ____            _     _    __
      / __ \____ ___  (_)___| |  / /__  _____________
     / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \
    / /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/
    \____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/

    🔖 Version: 1.7.2rc1   📚 Tutorials: https://omicverse.readthedocs.io/
    ✅ plot_set complete.

<div class="admonition warning">
  <p class="admonition-title">Note</p>
  <p>
    “When OmicVerse is upgraded to version > 1.7.0, it supports CPU–GPU mixed acceleration without requiring `rapids_singlecell` as a dependency—enjoy faster single-cell analysis!”

  </p>
</div>

```python
#ov.settings.cpu_gpu_mixed_init()
```

    CPU-GPU mixed mode activated

The data consist of 3k PBMCs from a Healthy Donor and are freely available from 10x Genomics ([here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) from this [webpage](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)). On a unix system, you can uncomment and run the following to download and unpack the data. The last line creates a directory for writing processed data.

```python
# !mkdir data
#!wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
#!cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write
```

```python
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata
```

    ... reading from cache file cache/data-filtered_gene_bc_matrices-hg19-matrix.h5ad





    AnnData object with n_obs × n_vars = 2700 × 32738
        var: 'gene_ids'

```python
adata.var_names_make_unique()
adata.obs_names_make_unique()
```

## Preprocessing

### Quantity control

For single-cell data, we require quality control prior to analysis, including the removal of cells containing double cells, low-expressing cells, and low-expressing genes. In addition to this, we need to filter based on mitochondrial gene ratios, number of transcripts, number of genes expressed per cell, cellular Complexity, etc. For a detailed description of the different QCs please see the document: https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

<div class="admonition warning">
  <p class="admonition-title">Note</p>
  <p>
    if the version of `omicverse` larger than `1.6.4`, the `doublets_method` can be set between `scrublet` and `sccomposite`.
  </p>
</div>

COMPOSITE (COMpound POiSson multIplet deTEction model) is a computational tool for multiplet detection in both single-cell single-omics and multiomics settings. It has been implemented as an automated pipeline and is available as both a cloud-based application with a user-friendly interface and a Python package.

Hu, H., Wang, X., Feng, S. et al. A unified model-based framework for doublet or multiplet detection in single-cell multiomics data. Nat Commun 15, 5562 (2024). https://doi.org/10.1038/s41467-024-49448-x

```python
%%time
adata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.2, 'nUMIs': 500, 'detected_genes': 250},
               doublets_method='scrublet',
              batch_key=None)
adata
```

    [95m[1m⚙️ Using CPU/GPU mixed mode for QC...[0m
    [94mNVIDIA CUDA GPUs detected:[0m
    📊 [CUDA 0] Tesla V100-SXM2-16GB
        [92m||[90m----------------------------[0m 1410/16384 MiB (8.6%)

    [95m[1m🔍 Quality Control Analysis (CPU-GPU Mixed):[0m
       [96mDataset shape: [1m2,700 cells × 32,738 genes[0m
       [94mQC mode: [1mseurat[0m
       [94mDoublet detection: [1mscrublet[0m
       [94mMitochondrial genes: [1mMT-[0m

    [95m[1m📊 Step 1: Calculating QC Metrics[0m
       [96mMitochondrial genes (prefix 'MT-'): [1m13[0m[96m found[0m
       [92m✓ QC metrics calculated:[0m
         [94m• Mean nUMIs: [1m2367[0m[94m (range: 548-15844)[0m
         [94m• Mean genes: [1m847[0m[94m (range: 212-3422)[0m
         [94m• Mean mitochondrial %: [1m2.2%[0m[94m (max: 22.6%)[0m

    [95m[1m🔍 Step 2: Doublet Detection[0m
       [93m⚠️  Note: 'scrublet' detection is legacy and may not work optimally[0m
       [96m💡 Consider using 'doublets_method=sccomposite' for better results[0m
       [92m🔍 Running scrublet doublet detection...[0m
    Running Scrublet🔍
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
        using data matrix X directly
    Automatically set threshold at doublet score = 0.23
    Detected doublet rate = 1.5%
    Estimated detectable doublet fraction = 37.4%
    Overall doublet rate:
    	Expected   = 5.0%
    	Estimated  = 4.1%
        Scrublet finished✅ (0:00:22)
       [92m✓ Scrublet completed: [1m41[0m[92m doublets removed (1.5%)[0m

    [95m[1m🔧 Step 3: Quality Filtering (SEURAT)[0m
       [96mThresholds: mito≤0.2, nUMIs≥500, genes≥250[0m
       [94m📊 Seurat Filter Results:[0m
         [96m• nUMIs filter (≥500): [1m0[0m[96m cells failed (0.0%)[0m
         [96m• Genes filter (≥250): [1m3[0m[96m cells failed (0.1%)[0m
         [96m• Mitochondrial filter (≤0.2): [1m2[0m[96m cells failed (0.1%)[0m
       [92m✓ Combined QC filters: [1m5[0m[92m cells removed (0.2%)[0m

    [95m[1m🎯 Step 4: Final Filtering[0m
       [96mParameters: min_genes=200, min_cells=3[0m
       [96mRatios: max_genes_ratio=1, max_cells_ratio=1[0m
    filtered out 19094 genes that are detected in less than 3 cells
       [92m✓ Final filtering: [1m0[0m[92m cells, [1m19,094[0m[92m genes removed[0m

    [92m✅ Quality Control Analysis Completed![0m

    [95m[1m📈 Final Summary:[0m
       [96m📊 Original: [1m2,700[0m[96m cells × [1m32,738[0m[96m genes[0m
       [92m✓ Final: [1m2,654[0m[92m cells × [1m13,644[0m[92m genes[0m
       [94m📉 Total removed: [1m46[0m[94m cells (1.7%)[0m
       [94m📉 Total removed: [1m19,094[0m[94m genes (58.3%)[0m
       [92m💯 Quality: [1m98.3%[0m[92m retention (Excellent retention rate)[0m

    [96m────────────────────────────────────────────────────────────[0m
    CPU times: user 2min 44s, sys: 2min 12s, total: 4min 57s
    Wall time: 22.6 s





    AnnData object with n_obs × n_vars = 2654 × 13644
        obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
        var: 'gene_ids', 'mt', 'n_cells'
        uns: 'scrublet', 'status', 'status_args'

### High variable Gene Detection

Here we try to use Pearson's method to calculate highly variable genes. This is the method that is proposed to be superior to ordinary normalisation. See [Article](https://www.nature.com/articles/s41592-023-01814-1#MOESM3) in _Nature Method_ for details.

normalize|HVGs：We use | to control the preprocessing step, | before for the normalisation step, either `shiftlog` or `pearson`, and | after for the highly variable gene calculation step, either `pearson` or `seurat`. Our default is `shiftlog|pearson`.

- if you use `mode`=`shiftlog|pearson` you need to set `target_sum=50*1e4`, more people like to se `target_sum=1e4`, we test the result think 50\*1e4 will be better
- if you use `mode`=`pearson|pearson`, you don't need to set `target_sum`

<div class="admonition warning">
  <p class="admonition-title">Note</p>
  <p>
    if the version of `omicverse` lower than `1.4.13`, the mode can only be set between `scanpy` and `pearson`.
  </p>
</div>

```python
%%time
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=2000,
                       target_sum=50*1e4)
adata
```

    Begin robust gene identification
    After filtration, 13644/13644 genes are kept.     Among 13644 genes, 13644 genes are robust.
    End of robust gene identification.
    Begin size normalization: shiftlog and HVGs selection pearson
    normalizing counts per cell
    The following highly-expressed genes are not considered during normalization factor computation:
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
    Time to analyze data in cpu: 0.4788978099822998 seconds.
    End of size normalization: shiftlog and HVGs selection pearson
    CPU times: user 1.44 s, sys: 96 ms, total: 1.53 s
    Wall time: 536 ms





    AnnData object with n_obs × n_vars = 2654 × 13644
        obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
        var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'means', 'variances', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'status', 'status_args', 'log1p', 'hvg', 'REFERENCE_MANU'
        layers: 'counts'

You can use `recover_counts` to recover the raw counts after normalize and log1p

```python
adata[:,'CD3D'].to_df().T
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AAACATACAACCAC-1</th>
      <th>AAACATTGAGCTAC-1</th>
      <th>AAACATTGATCAGC-1</th>
      <th>AAACCGTGCTTCCG-1</th>
      <th>AAACCGTGTATGCG-1</th>
      <th>AAACGCACTGGTAC-1</th>
      <th>AAACGCTGACCAGT-1</th>
      <th>AAACGCTGGTTCTT-1</th>
      <th>AAACGCTGTAGCCA-1</th>
      <th>AAACGCTGTTTCTG-1</th>
      <th>...</th>
      <th>TTTCAGTGTCACGA-1</th>
      <th>TTTCAGTGTCTATC-1</th>
      <th>TTTCAGTGTGCAGT-1</th>
      <th>TTTCCAGAGGTGAG-1</th>
      <th>TTTCGAACACCTGA-1</th>
      <th>TTTCGAACTCTCAT-1</th>
      <th>TTTCTACTGAGGCA-1</th>
      <th>TTTCTACTTCCTCG-1</th>
      <th>TTTGCATGAGAGGC-1</th>
      <th>TTTGCATGCCTCAC-1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CD3D</th>
      <td>6.718757</td>
      <td>0.0</td>
      <td>7.371373</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5.447429</td>
      <td>6.132899</td>
      <td>6.499361</td>
      <td>5.974209</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>6.532146</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>6.224622</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 2654 columns</p>
</div>

```python
adata[:,'CD3D'].to_df(layer='counts').T
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AAACATACAACCAC-1</th>
      <th>AAACATTGAGCTAC-1</th>
      <th>AAACATTGATCAGC-1</th>
      <th>AAACCGTGCTTCCG-1</th>
      <th>AAACCGTGTATGCG-1</th>
      <th>AAACGCACTGGTAC-1</th>
      <th>AAACGCTGACCAGT-1</th>
      <th>AAACGCTGGTTCTT-1</th>
      <th>AAACGCTGTAGCCA-1</th>
      <th>AAACGCTGTTTCTG-1</th>
      <th>...</th>
      <th>TTTCAGTGTCACGA-1</th>
      <th>TTTCAGTGTCTATC-1</th>
      <th>TTTCAGTGTGCAGT-1</th>
      <th>TTTCCAGAGGTGAG-1</th>
      <th>TTTCGAACACCTGA-1</th>
      <th>TTTCGAACTCTCAT-1</th>
      <th>TTTCTACTGAGGCA-1</th>
      <th>TTTCTACTTCCTCG-1</th>
      <th>TTTGCATGAGAGGC-1</th>
      <th>TTTGCATGCCTCAC-1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CD3D</th>
      <td>4.0</td>
      <td>0.0</td>
      <td>10.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>2.0</td>
      <td>3.0</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>3.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.0</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 2654 columns</p>
</div>

```python
X_counts_recovered, size_factors_sub=ov.pp.recover_counts(adata.X, 50*1e4, 50*1e5, log_base=None,
                                                          chunk_size=10000)
adata.layers['recover_counts']=X_counts_recovered
adata[:,'CD3D'].to_df(layer='recover_counts').T
```

    100%|██████████| 2654/2654 [00:03<00:00, 831.68it/s]

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AAACATACAACCAC-1</th>
      <th>AAACATTGAGCTAC-1</th>
      <th>AAACATTGATCAGC-1</th>
      <th>AAACCGTGCTTCCG-1</th>
      <th>AAACCGTGTATGCG-1</th>
      <th>AAACGCACTGGTAC-1</th>
      <th>AAACGCTGACCAGT-1</th>
      <th>AAACGCTGGTTCTT-1</th>
      <th>AAACGCTGTAGCCA-1</th>
      <th>AAACGCTGTTTCTG-1</th>
      <th>...</th>
      <th>TTTCAGTGTCACGA-1</th>
      <th>TTTCAGTGTCTATC-1</th>
      <th>TTTCAGTGTGCAGT-1</th>
      <th>TTTCCAGAGGTGAG-1</th>
      <th>TTTCGAACACCTGA-1</th>
      <th>TTTCGAACTCTCAT-1</th>
      <th>TTTCTACTGAGGCA-1</th>
      <th>TTTCTACTTCCTCG-1</th>
      <th>TTTGCATGAGAGGC-1</th>
      <th>TTTGCATGCCTCAC-1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CD3D</th>
      <td>4</td>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 2654 columns</p>
</div>

Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

```python
%%time
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata
```

    CPU times: user 13.3 ms, sys: 2.15 ms, total: 15.4 ms
    Wall time: 13.4 ms





    View of AnnData object with n_obs × n_vars = 2654 × 2000
        obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
        var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'means', 'variances', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'status', 'status_args', 'log1p', 'hvg', 'REFERENCE_MANU'
        layers: 'counts', 'recover_counts'

## Principal component analysis

In contrast to scanpy, we do not directly scale the variance of the original expression matrix, but store the results of the variance scaling in the layer, due to the fact that scale may cause changes in the data distribution, and we have not found scale to be meaningful in any scenario other than a principal component analysis

```python
%%time
ov.pp.scale(adata)
adata
```

    CPU times: user 862 ms, sys: 54.7 ms, total: 917 ms
    Wall time: 871 ms





    AnnData object with n_obs × n_vars = 2654 × 2000
        obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
        var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'means', 'variances', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'status', 'status_args', 'log1p', 'hvg', 'REFERENCE_MANU'
        layers: 'counts', 'recover_counts', 'scaled'

If you want to perform pca in normlog layer, you can set `layer`=`normlog`, but we think scaled is necessary in PCA.

```python
%%time
ov.pp.pca(adata,layer='scaled',n_pcs=50)
adata
```

    🚀 Using GPU to calculate PCA...
    [94mNVIDIA CUDA GPUs detected:[0m
    📊 [CUDA 0] Tesla P100-PCIE-16GB
        [92m[90m------------------------------[0m 3/16384 MiB (0.0%)
    computing PCA🔍
        with n_comps=50
    [KeOps] Compiling cuda jit compiler engine ... OK
    [pyKeOps] Compiling nvrtc binder for python ... OK
        finished✅ (0:00:36)
    CPU times: user 3min 11s, sys: 2min 1s, total: 5min 13s
    Wall time: 36.6 s





    AnnData object with n_obs × n_vars = 2654 × 2000
        obs: 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'n_genes'
        var: 'gene_ids', 'mt', 'n_cells', 'percent_cells', 'robust', 'means', 'variances', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'status', 'status_args', 'log1p', 'hvg', 'REFERENCE_MANU', 'pca', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues'
        obsm: 'X_pca', 'scaled|original|X_pca'
        varm: 'PCs', 'scaled|original|pca_loadings'
        layers: 'counts', 'recover_counts', 'scaled'

```python
adata.obsm['X_pca']=adata.obsm['scaled|original|X_pca']
ov.pl.embedding(adata,
                  basis='X_pca',
                  color='CST3',
                  frameon='small')
```

![png](output_23_0.png)

## Embedding the neighborhood graph

We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below. It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:

```python
%%time
ov.pp.neighbors(adata, n_neighbors=15, n_pcs=50,
               use_rep='scaled|original|X_pca')
```

    🖥️ Using Scanpy CPU to calculate neighbors...
    computing neighbors
        finished: added to `.uns['neighbors']`
        `.obsp['distances']`, distances for each pair of neighbors
        `.obsp['connectivities']`, weighted adjacency matrix (0:00:11)
    CPU times: user 11 s, sys: 279 ms, total: 11.3 s
    Wall time: 11.2 s

You also can use `umap` to visualize the neighborhood graph

```python
%%time
ov.pp.umap(adata)
```

    🔍 [2025-06-24 19:57:24] Running UMAP in 'cpu-gpu-mixed' mode...
    🚀 Using torch GPU to calculate UMAP...
    [94mNVIDIA CUDA GPUs detected:[0m
    📊 [CUDA 0] Tesla P100-PCIE-16GB
        [92m[90m------------------------------[0m 3/16384 MiB (0.0%)
    computing UMAP🚀
        finished ✅: added
        'X_umap', UMAP coordinates (adata.obsm)
        'umap', UMAP parameters (adata.uns) (0:00:11)
    ✅ UMAP completed successfully.
    CPU times: user 4.15 s, sys: 1.37 s, total: 5.52 s
    Wall time: 11.7 s

```python
ov.pl.embedding(adata,
                basis='X_umap',
                color='CST3',
                frameon='small')
```

![png](output_28_0.png)

To visualize the PCA’s embeddings, we use the `pymde` package wrapper in omicverse. This is an alternative to UMAP that is GPU-accelerated.

```python
ov.pp.mde(adata,embedding_dim=2,n_neighbors=15, basis='X_mde',
          n_pcs=50, use_rep='scaled|original|X_pca',)
```

    computing neighbors
        finished: added to `.uns['neighbors']`
        `.obsm['X_mde']`, MDE coordinates
        `.obsp['neighbors_distances']`, distances for each pair of neighbors
        `.obsp['neighbors_connectivities']`, weighted adjacency matrix (0:00:03)

```python
ov.pl.embedding(adata,
                basis='X_mde',
                color='CST3',
                frameon='small')
```

![png](output_31_0.png)

## Score cell cyle

In OmicVerse, we store both G1M/S and G2M genes into the function (both human and mouse), so you can run cell cycle analysis without having to manually enter cycle genes!

```python
adata_raw=adata.raw.to_adata()
ov.pp.score_genes_cell_cycle(adata_raw,species='human')
```

    calculating cell cycle phase
    computing score 'S_score'
    WARNING: genes are not in var_names and ignored: Index(['DTL', 'UHRF1', 'MLF1IP', 'CDC6', 'EXO1', 'CASP8AP2', 'BRIP1', 'E2F8'], dtype='object')
        finished: added
        'S_score', score of gene set (adata.obs).
        644 total control genes are used. (0:00:00)
    computing score 'G2M_score'
    WARNING: genes are not in var_names and ignored: Index(['FAM64A', 'BUB1', 'HJURP', 'CDCA3', 'TTK', 'CDC25C', 'DLGAP5', 'CDCA2',
           'CDCA8', 'ANLN', 'NEK2', 'GAS2L3'],
          dtype='object')
        finished: added
        'G2M_score', score of gene set (adata.obs).
        815 total control genes are used. (0:00:00)
    -->     'phase', cell cycle phase (adata.obs)

```python
ov.pl.embedding(adata_raw,
                basis='X_mde',
                color='phase',
                frameon='small')
```

![png](output_34_0.png)

## Clustering the neighborhood graph

As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag _et al._ (2018). Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

```python
ov.pp.leiden(adata,resolution=1)
```

    🖥️ Using Scanpy CPU Leiden...
    running Leiden clustering
        finished: found 10 clusters and added
        'leiden', the cluster labels (adata.obs, categorical) (0:00:00)

We redesigned the visualisation of embedding to distinguish it from scanpy's embedding by adding the parameter `fraemon='small'`, which causes the axes to be scaled with the colourbar

```python
ov.pl.embedding(adata,
                basis='X_mde',
                color=['leiden', 'CST3', 'NKG7'],
                frameon='small')
```

![png](output_38_0.png)

We also provide a boundary visualisation function `ov.utils.plot_ConvexHull` to visualise specific clusters.

Arguments:

- color: if None will use the color of clusters
- alpha: default is 0.2

```python
import matplotlib.pyplot as plt
fig,ax=plt.subplots( figsize = (4,4))

ov.pl.embedding(adata,
                basis='X_mde',
                color=['leiden'],
                frameon='small',
                show=False,
                ax=ax)

ov.pl.ConvexHull(adata,
                basis='X_mde',
                cluster_key='leiden',
                hull_cluster='0',
                ax=ax)

```

    leiden_colors





    <Axes: title={'center': 'leiden'}, xlabel='X_mde1', ylabel='X_mde2'>

![png](output_40_2.png)

If you have too many labels, e.g. too many cell types, and you are concerned about cell overlap, then consider trying the `ov.utils.gen_mpl_labels` function, which improves text overlap.
In addition, we make use of the `patheffects` function, which makes our text have outlines

- adjust_kwargs: it could be found in package `adjusttext`
- text_kwargs: it could be found in class `plt.texts`

```python
from matplotlib import patheffects
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4))

ov.pl.embedding(adata,
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
)
```

![png](output_42_0.png)

```python
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
```

```python
ov.pl.dotplot(adata, marker_genes, groupby='leiden',
             standard_scale='var');
```

![png](output_44_0.png)

## Finding marker genes

Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.

```python
sc.tl.dendrogram(adata,'leiden',use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(adata, 'leiden', use_rep='scaled|original|X_pca',
                        method='t-test',use_raw=False,key_added='leiden_ttest')
ov.pl.rank_genes_groups_dotplot(adata,groupby='leiden',
                                cmap='Spectral_r',key='leiden_ttest',
                                standard_scale='var',n_genes=3,dendrogram=False)
```

    Storing dendrogram info using `.uns['dendrogram_leiden']`
    ranking genes
        finished: added to `.uns['leiden_ttest']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)

![png](output_46_1.png)

cosg is also considered to be a better algorithm for finding marker genes. Here, omicverse provides the calculation of cosg

Paper: [Accurate and fast cell marker gene identification with COSG](https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext)

Code: https://github.com/genecell/COSG

```python
sc.tl.rank_genes_groups(adata, groupby='leiden',
                        method='t-test',use_rep='scaled|original|X_pca',)
ov.single.cosg(adata, key_added='leiden_cosg', groupby='leiden')
ov.pl.rank_genes_groups_dotplot(adata,groupby='leiden',
                                cmap='Spectral_r',key='leiden_cosg',
                                standard_scale='var',n_genes=3,dendrogram=False)
```

    ranking genes
        finished: added to `.uns['rank_genes_groups']`
        'names', sorted np.recarray to be indexed by group ids
        'scores', sorted np.recarray to be indexed by group ids
        'logfoldchanges', sorted np.recarray to be indexed by group ids
        'pvals', sorted np.recarray to be indexed by group ids
        'pvals_adj', sorted np.recarray to be indexed by group ids (0:00:00)
    **finished identifying marker genes by COSG**

![png](output_48_1.png)
