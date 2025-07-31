# Data integration and batch correction

An important task of single-cell analysis is the integration of several samples, which we can perform with omicverse.

Here we demonstrate how to merge data using omicverse and perform a corrective analysis for batch effects. We provide a total of 4 methods for batch effect correction in omicverse, including harmony, scanorama and combat which do not require GPU, and SIMBA which requires GPU. if available, we recommend using GPU-based scVI and scANVI to get the best batch effect correction results.

```python
import omicverse as ov
#print(f"omicverse version: {ov.__version__}")
import scanpy as sc
#print(f"scanpy version: {sc.__version__}")
ov.utils.ov_plot_set()
```

## Data integration

First, we need to concat the data of scRNA-seq from different batch. We can use `sc.concat` to perform it。

The dataset we will use to demonstrate data integration contains several samples of bone marrow mononuclear cells. These samples were originally created for the Open Problems in Single-Cell Analysis NeurIPS Competition 2021.

We selected sample of `s1d3`, `s2d1` and `s3d7` to perform integrate. The individual data can be downloaded from figshare.

- s1d3:
- s2d1:
- s3d7:

```python
adata1=ov.read('neurips2021_s1d3.h5ad')
adata1.obs['batch']='s1d3'
adata2=ov.read('neurips2021_s2d1.h5ad')
adata2.obs['batch']='s2d1'
adata3=ov.read('neurips2021_s3d7.h5ad')
adata3.obs['batch']='s3d7'
```

```python
adata=sc.concat([adata1,adata2,adata3],merge='same')
adata
```

    AnnData object with n_obs × n_vars = 27423 × 13953
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train'
        var: 'feature_types', 'gene_id'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap'
        layers: 'counts'

We can see that there are now three elements in the batch

```python
adata.obs['batch'].unique()
```

    array(['s1d3', 's2d1', 's3d7'], dtype=object)

```python
import numpy as np
adata.X=adata.X.astype(np.int64)
```

## Data preprocess and Batch visualize

We first performed quality control of the data and normalisation with screening for highly variable genes. Then visualise potential batch effects in the data.

Here, we can set `batch_key=batch` to correct the doublet detectation and Highly variable genes identifcation.

```python
adata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.2, 'nUMIs': 500, 'detected_genes': 250},
              batch_key='batch')
adata
```

We can store the raw counts if we need the raw counts after filtered the HVGs.

```python
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',
                       n_HVGs=3000,batch_key=None)
adata
```

    Begin robust gene identification
    After filtration, 13953/13953 genes are kept. Among 13953 genes, 13953 genes are robust.
    End of robust gene identification.
    Begin size normalization: shiftlog and HVGs selection pearson
    normalizing counts per cell The following highly-expressed genes are not considered during normalization factor computation:
    ['IGKC', 'HBB', 'MALAT1', 'IGHA1', 'IGHM', 'HBA2', 'IGLC1', 'IGLC2', 'IGLC3']
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





    AnnData object with n_obs × n_vars = 26707 × 13953
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'layers_counts', 'log1p', 'hvg'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap'
        layers: 'counts'

```python
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata
```

    View of AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'layers_counts', 'log1p', 'hvg'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap'
        layers: 'counts'

We can save the pre-processed data.

```python
adata.write_h5ad('neurips2021_batch_normlog.h5ad',compression='gzip')
```

Similarly, we calculated PCA for HVGs and visualised potential batch effects in the data using pymde. pymde is GPU-accelerated UMAP.

```python
ov.pp.scale(adata)
ov.pp.pca(adata,layer='scaled',n_pcs=50)

adata.obsm["X_mde_pca"] = ov.utils.mde(adata.obsm["scaled|original|X_pca"])
```

    ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
    computing PCA
        with n_comps=50
        finished (0:00:13)

There is a very clear batch effect in the data

```python
ov.utils.embedding(adata,
                basis='X_mde_pca',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_pca1', ylabel='X_mde_pca2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_pca1', ylabel='X_mde_pca2'>]

![png](output_18_1.png)

## Harmony

Harmony is an algorithm for performing integration of single cell genomics datasets. Please check out manuscript on [Nature Methods](https://www.nature.com/articles/s41592-019-0619-0).

![harmony](https://portals.broadinstitute.org/harmony/articles/main.jpg)

The function `ov.single.batch_correction` can be set in three methods: `harmony`,`combat` and `scanorama`

```python
ov.single.batch_correction(adata,batch_key='batch',
                                        methods='harmony',n_pcs=50)
adata
```

    ...Begin using harmony to correct batch effect
    ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption


    2023-11-19 20:25:03,351 - harmonypy - INFO - Computing initial centroids with sklearn.KMeans...
    INFO:harmonypy:Computing initial centroids with sklearn.KMeans...
    2023-11-19 20:25:12,444 - harmonypy - INFO - sklearn.KMeans initialization complete.
    INFO:harmonypy:sklearn.KMeans initialization complete.
    2023-11-19 20:25:12,725 - harmonypy - INFO - Iteration 1 of 10
    INFO:harmonypy:Iteration 1 of 10
    2023-11-19 20:25:19,161 - harmonypy - INFO - Iteration 2 of 10
    INFO:harmonypy:Iteration 2 of 10
    2023-11-19 20:25:25,779 - harmonypy - INFO - Iteration 3 of 10
    INFO:harmonypy:Iteration 3 of 10
    2023-11-19 20:25:32,350 - harmonypy - INFO - Iteration 4 of 10
    INFO:harmonypy:Iteration 4 of 10
    2023-11-19 20:25:38,889 - harmonypy - INFO - Iteration 5 of 10
    INFO:harmonypy:Iteration 5 of 10
    2023-11-19 20:25:43,768 - harmonypy - INFO - Converged after 5 iterations
    INFO:harmonypy:Converged after 5 iterations





    AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'topic_0', 'topic_1', 'topic_2', 'topic_3', 'topic_4', 'topic_5', 'topic_6', 'topic_7', 'topic_8', 'topic_9', 'topic_10', 'topic_11', 'topic_12', 'topic_13', 'topic_14', 'LDA_cluster'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues', 'batch_colors', 'cell_type_colors', 'topic_dendogram'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap', 'scaled|original|X_pca', 'X_mde_pca', 'X_topic_compositions', 'X_umap_features', 'X_mde_mira', 'X_mde_mira_topic', 'X_mde_mira_feature', 'X_harmony'
        varm: 'scaled|original|pca_loadings', 'topic_feature_compositions', 'topic_feature_activations'
        layers: 'counts', 'scaled', 'lognorm'

```python
adata.obsm["X_mde_harmony"] = ov.utils.mde(adata.obsm["X_harmony"])
```

```python
ov.utils.embedding(adata,
                basis='X_mde_harmony',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_harmony1', ylabel='X_mde_harmony2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_harmony1', ylabel='X_mde_harmony2'>]

![png](output_23_1.png)

## Combat

combat is a batch effect correction method that is very widely used in bulk RNA-seq, and it works just as well on single-cell sequencing data.

```python
ov.single.batch_correction(adata,batch_key='batch',
                                        methods='combat',n_pcs=50)
adata
```

    ...Begin using combat to correct batch effect
    Standardizing Data across genes.

    Found 3 batches

    Found 0 numerical variables:


    Fitting L/S model and finding priors

    Finding parametric adjustments

    Adjusting data






    AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'topic_0', 'topic_1', 'topic_2', 'topic_3', 'topic_4', 'topic_5', 'topic_6', 'topic_7', 'topic_8', 'topic_9', 'topic_10', 'topic_11', 'topic_12', 'topic_13', 'topic_14', 'LDA_cluster'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues', 'batch_colors', 'cell_type_colors', 'topic_dendogram'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap', 'scaled|original|X_pca', 'X_mde_pca', 'X_topic_compositions', 'X_umap_features', 'X_mde_mira', 'X_mde_mira_topic', 'X_mde_mira_feature', 'X_harmony', 'X_mde_harmony', 'X_combat'
        varm: 'scaled|original|pca_loadings', 'topic_feature_compositions', 'topic_feature_activations'
        layers: 'counts', 'scaled', 'lognorm'

```python
adata.obsm["X_mde_combat"] = ov.utils.mde(adata.obsm["X_combat"])
```

```python
ov.utils.embedding(adata,
                basis='X_mde_combat',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_combat1', ylabel='X_mde_combat2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_combat1', ylabel='X_mde_combat2'>]

![png](output_27_1.png)

## scanorama

Integration of single-cell RNA sequencing (scRNA-seq) data from multiple experiments, laboratories and technologies can uncover biological insights, but current methods for scRNA-seq data integration are limited by a requirement for datasets to derive from functionally similar cells. We present Scanorama, an algorithm that identifies and merges the shared cell types among all pairs of datasets and accurately integrates heterogeneous collections of scRNA-seq data.

![scanorama](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41587-019-0113-3/MediaObjects/41587_2019_113_Fig1_HTML.png?as=webp)

```python
ov.single.batch_correction(adata,batch_key='batch',
                                        methods='scanorama',n_pcs=50)
adata
```

    ...Begin using scanorama to correct batch effect
    s1d3
    s2d1
    s3d7
    Found 3000 genes among all datasets
    [[0.         0.50093205 0.5758346 ]
     [0.         0.         0.60733037]
     [0.         0.         0.        ]]
    Processing datasets (1, 2)
    Processing datasets (0, 2)
    Processing datasets (0, 1)
    (26707, 50)





    AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'topic_0', 'topic_1', 'topic_2', 'topic_3', 'topic_4', 'topic_5', 'topic_6', 'topic_7', 'topic_8', 'topic_9', 'topic_10', 'topic_11', 'topic_12', 'topic_13', 'topic_14', 'LDA_cluster'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'scrublet', 'layers_counts', 'log1p', 'hvg', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues', 'batch_colors', 'cell_type_colors', 'topic_dendogram'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap', 'scaled|original|X_pca', 'X_mde_pca', 'X_topic_compositions', 'X_umap_features', 'X_mde_mira', 'X_mde_mira_topic', 'X_mde_mira_feature', 'X_harmony', 'X_mde_harmony', 'X_combat', 'X_mde_combat', 'X_scanorama'
        varm: 'scaled|original|pca_loadings', 'topic_feature_compositions', 'topic_feature_activations'
        layers: 'counts', 'scaled', 'lognorm'

```python
adata.obsm["X_mde_scanorama"] = ov.utils.mde(adata.obsm["X_scanorama"])
```

```python
ov.utils.embedding(adata,
                basis='X_mde_scanorama',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_scanorama1', ylabel='X_mde_scanorama2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_scanorama1', ylabel='X_mde_scanorama2'>]

![png](output_31_1.png)

## scVI

An important task of single-cell analysis is the integration of several samples, which we can perform with scVI. For integration, scVI treats the data as unlabelled. When our dataset is fully labelled (perhaps in independent studies, or independent analysis pipelines), we can obtain an integration that better preserves biology using scANVI, which incorporates cell type annotation information. Here we demonstrate this functionality with an integrated analysis of cells from the lung atlas integration task from the scIB manuscript. The same pipeline would generally be used to analyze any collection of scRNA-seq datasets.

```python
model=ov.single.batch_correction(adata,batch_key='batch',
                           methods='scVI',n_layers=2, n_latent=30, gene_likelihood="nb")
adata
```

    ...Begin using scVI to correct batch effect


    Global seed set to 0
    No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
    GPU available: True (cuda), used: True
    TPU available: False, using: 0 TPU cores
    IPU available: False, using: 0 IPUs
    HPU available: False, using: 0 HPUs
    LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]


    Epoch 300/300: 100%|██████████| 300/300 [05:51<00:00,  1.18s/it, loss=1.51e+03, v_num=1]

    `Trainer.fit` stopped: `max_epochs=300` reached.


    Epoch 300/300: 100%|██████████| 300/300 [05:51<00:00,  1.17s/it, loss=1.51e+03, v_num=1]





    AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'topic_0', 'topic_1', 'topic_2', 'topic_3', 'topic_4', 'topic_5', 'topic_6', 'topic_7', 'topic_8', 'topic_9', 'topic_10', 'topic_11', 'topic_12', 'topic_13', 'topic_14', 'LDA_cluster', '_scvi_batch', '_scvi_labels'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features'
        uns: 'batch_colors', 'cell_type_colors', 'hvg', 'layers_counts', 'log1p', 'scaled|original|cum_sum_eigenvalues', 'scaled|original|pca_var_ratios', 'scrublet', 'topic_dendogram', '_scvi_uuid', '_scvi_manager_uuid'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap', 'X_combat', 'X_harmony', 'X_mde_combat', 'X_mde_harmony', 'X_mde_mira', 'X_mde_mira_feature', 'X_mde_mira_topic', 'X_mde_pca', 'X_mde_scanorama', 'X_scanorama', 'X_topic_compositions', 'X_umap_features', 'scaled|original|X_pca', 'X_scVI'
        varm: 'scaled|original|pca_loadings', 'topic_feature_activations', 'topic_feature_compositions'
        layers: 'counts', 'lognorm', 'scaled'

```python
adata.obsm["X_mde_scVI"] = ov.utils.mde(adata.obsm["X_scVI"])
```

```python
ov.utils.embedding(adata,
                basis='X_mde_scVI',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_scVI1', ylabel='X_mde_scVI2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_scVI1', ylabel='X_mde_scVI2'>]

![png](output_35_1.png)

## CellANOVA

The integration of cells across samples to remove unwanted batch variation plays a critical role in single cell analyses. When the samples are expected to be biologically distinct, it is often unclear how aggressively the cells should be aligned across samples to achieve uniformity. CellANOVA is a Python package for batch integration with signal recovery in single cell data. It builds on existing single cell data integration methods, and uses a pool of control samples to quantify the batch effect and separate meaningful biological variation from unwanted batch variation. When used with an existing integration method, CellAnova allows the recovery of biological signals that are lost during integration.

In omicverse, you only need to prepare the `control_dict`(At least two samples are required!) when you want to try `CellANOVA`. When you're done running it, there are two outputs you need to be aware of:

1. the first one being: `adata.layers['denoised']`, which stores the matrix after the batch effect is removed.
2. The second is `adata.obsm['X_mde_cellANOVA']`, which stores the low-dimensional representation of the cell after removing the batch effect

- Zhang, Z., Mathew, D., Lim, T.L. et al. Recovery of biological signals lost in single-cell batch integration with CellANOVA. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02463-1

```python
## construct control pool
control_dict = {
    'pool1': ['s1d3','s2d1'],
}

ov.single.batch_correction(adata,batch_key='batch',n_pcs=50,
                           methods='CellANOVA',control_dict=control_dict)
adata
```

    ...Begin using CellANOVA to correct batch effect
    computing PCA
        with n_comps=70
        finished (0:01:05)


    2024-12-16 16:00:00,139 - harmonypy - INFO - Computing initial centroids with sklearn.KMeans...
    2024-12-16 16:00:26,367 - harmonypy - INFO - sklearn.KMeans initialization complete.
    2024-12-16 16:00:26,739 - harmonypy - INFO - Iteration 1 of 30
    2024-12-16 16:01:12,771 - harmonypy - INFO - Iteration 2 of 30
    2024-12-16 16:01:53,924 - harmonypy - INFO - Iteration 3 of 30
    2024-12-16 16:02:37,141 - harmonypy - INFO - Iteration 4 of 30
    2024-12-16 16:03:19,467 - harmonypy - INFO - Iteration 5 of 30
    2024-12-16 16:03:44,733 - harmonypy - INFO - Converged after 5 iterations


    computing PCA
        with n_comps=50
        finished (0:01:14)





    AnnData object with n_obs × n_vars = 26707 × 3000
        obs: 'GEX_n_genes_by_counts', 'GEX_pct_counts_mt', 'GEX_size_factors', 'GEX_phase', 'ADT_n_antibodies_by_counts', 'ADT_total_counts', 'ADT_iso_count', 'cell_type', 'batch', 'ADT_pseudotime_order', 'GEX_pseudotime_order', 'Samplename', 'Site', 'DonorNumber', 'Modality', 'VendorLot', 'DonorID', 'DonorAge', 'DonorBMI', 'DonorBloodType', 'DonorRace', 'Ethnicity', 'DonorGender', 'QCMeds', 'DonorSmoker', 'is_train', 'nUMIs', 'mito_perc', 'detected_genes', 'cell_complexity', 'n_genes', 'doublet_score', 'predicted_doublet', 'passing_mt', 'passing_nUMIs', 'passing_ngenes', 'topic_0', 'topic_1', 'topic_2', 'topic_3', 'topic_4', 'topic_5', 'topic_6', 'topic_7', 'topic_8', 'topic_9', 'topic_10', 'topic_11', 'topic_12', 'topic_13', 'topic_14', 'LDA_cluster', '_scvi_batch', '_scvi_labels'
        var: 'feature_types', 'gene_id', 'mt', 'n_cells', 'percent_cells', 'robust', 'mean', 'var', 'residual_variances', 'highly_variable_rank', 'highly_variable_features', 'highly_variable'
        uns: '_scvi_manager_uuid', '_scvi_uuid', 'batch_colors', 'cell_type_colors', 'hvg', 'layers_counts', 'log1p', 'scaled|original|cum_sum_eigenvalues', 'scaled|original|pca_var_ratios', 'scrublet', 'topic_dendogram', 'S_BE_basis', 'S_TE_basis', 'pca', 'denoised|original|pca_var_ratios', 'denoised|original|cum_sum_eigenvalues'
        obsm: 'ADT_X_pca', 'ADT_X_umap', 'ADT_isotype_controls', 'GEX_X_pca', 'GEX_X_umap', 'X_combat', 'X_harmony', 'X_mde_combat', 'X_mde_harmony', 'X_mde_mira', 'X_mde_mira_feature', 'X_mde_mira_topic', 'X_mde_pca', 'X_mde_scVI', 'X_mde_scanorama', 'X_scVI', 'X_scanorama', 'X_topic_compositions', 'X_umap_features', 'scaled|original|X_pca', 'Cmat', 'X_pca', 'denoised|original|X_pca', 'X_cellanova'
        varm: 'scaled|original|pca_loadings', 'topic_feature_activations', 'topic_feature_compositions', 'Mmat', 'V_BE_basis', 'W_TE_basis', 'PCs', 'denoised|original|pca_loadings'
        layers: 'counts', 'lognorm', 'scaled', 'main_effect', 'corrected', 'trt_effect', 'denoised'

```python
adata.obsm["X_mde_cellANOVA"] = ov.utils.mde(adata.obsm["X_cellanova"])
```

```python
ov.pl.embedding(adata,
                basis='X_mde_cellANOVA',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_cellANOVA1', ylabel='X_mde_cellANOVA2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_cellANOVA1', ylabel='X_mde_cellANOVA2'>]

![png](output_39_1.png)

## MIRA+CODAL

Topic modeling of batched single-cell data is challenging because these models cannot typically distinguish between biological and technical effects of the assay. CODAL (COvariate Disentangling Augmented Loss) uses a novel mutual information regularization technique to explicitly disentangle these two sources of variation.

```python
LDA_obj=ov.utils.LDA_topic(adata,feature_type='expression',
                  highly_variable_key='highly_variable_features',
                 layers='counts',batch_key='batch',learning_rate=1e-3)
```

    INFO:mira.adata_interface.topic_model:Predicting expression from genes from col: highly_variable_features


    mira have been install version: 2.1.0



    Gathering dataset statistics:   0%|          | 0/26707 [00:00<?, ?it/s]



    Learning rate range test:   0%|          | 0/98 [00:00<?, ?it/s]


    INFO:mira.topic_model.base:Set learning rates to: (0.0061957597093704065, 0.22248668375233174)

![png](output_41_5.png)

```python
LDA_obj.plot_topic_contributions(6)
```

    Gathering dataset statistics:   0%|          | 0/26707 [00:00<?, ?it/s]



    Epoch 0:   0%|          | 0/24 [00:00<?, ?it/s]



    Predicting latent vars:   0%|          | 0/105 [00:00<?, ?it/s]

![png](output_42_3.png)

```python
LDA_obj.predicted(15)
```

    INFO:mira.adata_interface.topic_model:Predicting expression from genes from col: highly_variable_features


    running LDA topic predicted



    Gathering dataset statistics:   0%|          | 0/26707 [00:00<?, ?it/s]



    Epoch 0:   0%|          | 0/24 [00:00<?, ?it/s]


    INFO:mira.topic_model.base:Moving model to device: cpu



    Predicting latent vars:   0%|          | 0/105 [00:00<?, ?it/s]


    INFO:mira.adata_interface.core:Added key to obsm: X_topic_compositions
    INFO:mira.adata_interface.core:Added key to obsm: X_umap_features
    INFO:mira.adata_interface.topic_model:Added cols: topic_0, topic_1, topic_2, topic_3, topic_4, topic_5, topic_6, topic_7, topic_8, topic_9, topic_10, topic_11, topic_12, topic_13, topic_14
    INFO:mira.adata_interface.core:Added key to varm: topic_feature_compositions
    INFO:mira.adata_interface.core:Added key to varm: topic_feature_activations
    INFO:mira.adata_interface.topic_model:Added key to uns: topic_dendogram


    finished: found 15 clusters and added
        'LDA_cluster', the cluster labels (adata.obs, categorical)

```python
adata.obsm["X_mde_mira_topic"] = ov.utils.mde(adata.obsm["X_topic_compositions"])
adata.obsm["X_mde_mira_feature"] = ov.utils.mde(adata.obsm["X_umap_features"])
```

```python
ov.utils.embedding(adata,
                basis='X_mde_mira_topic',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_mira_topic1', ylabel='X_mde_mira_topic2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_mira_topic1', ylabel='X_mde_mira_topic2'>]

![png](output_45_1.png)

```python
ov.utils.embedding(adata,
                basis='X_mde_mira_feature',frameon='small',
                color=['batch','cell_type'],show=False)
```

    [<AxesSubplot: title={'center': 'batch'}, xlabel='X_mde_mira_feature1', ylabel='X_mde_mira_feature2'>,
     <AxesSubplot: title={'center': 'cell_type'}, xlabel='X_mde_mira_feature1', ylabel='X_mde_mira_feature2'>]

![png](output_46_1.png)

## Benchmarking test

The methods demonstrated here are selected based on results from benchmarking experiments including the single-cell integration benchmarking project [Luecken et al., 2021]. This project also produced a software package called [scib](https://www.github.com/theislab/scib) that can be used to run a range of integration methods as well as the metrics that were used for evaluation. In this section, we show how to use this package to evaluate the quality of an integration.

```python
adata.write_h5ad('neurips2021_batch_all.h5ad',compression='gzip')
```

```python
adata=sc.read('neurips2021_batch_all.h5ad')
```

```python
adata.obsm['X_pca']=adata.obsm['scaled|original|X_pca'].copy()
adata.obsm['X_mira_topic']=adata.obsm['X_topic_compositions'].copy()
adata.obsm['X_mira_feature']=adata.obsm['X_umap_features'].copy()
```

```python
from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="cell_type",
    embedding_obsm_keys=["X_pca", "X_combat", "X_harmony",'X_cellanova',
                         'X_scanorama','X_mira_topic','X_mira_feature','X_scVI',],
    n_jobs=8,
)
bm.benchmark()
```

    computing PCA
        with n_comps=50
        finished (0:00:44)


    Computing neighbors: 100%|██████████| 8/8 [03:08<00:00, 23.52s/it]
    Embeddings:   0%|[32m          [0m| 0/8 [00:00<?, ?it/s]
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:07<01:08,  7.64s/it, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:07<01:08,  7.64s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:14<00:55,  6.93s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:14<00:55,  6.93s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:14<00:28,  4.07s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:14<00:28,  4.07s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:19<00:25,  4.23s/it, Bio conservation: clisi_knn][A
    Metrics:  40%|[34m████      [0m| 4/10 [00:19<00:25,  4.23s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [01:56<03:08, 37.75s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [01:56<03:08, 37.75s/it, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [01:56<01:40, 25.10s/it, Batch correction: ilisi_knn][A
    Metrics:  60%|[34m██████    [0m| 6/10 [01:56<01:40, 25.10s/it, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [03:57<02:48, 56.23s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [03:57<02:48, 56.23s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [03:57<01:16, 38.40s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [03:57<01:16, 38.40s/it, Batch correction: pcr_comparison]    [A
    Embeddings:  12%|[32m█▎        [0m| 1/8 [03:59<27:55, 239.40s/it]tch correction: pcr_comparison][A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.36it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.36it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:02<00:12,  1.58s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:02<00:12,  1.58s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:03<00:08,  1.16s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:03<00:08,  1.16s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:03<00:06,  1.16s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:04<00:03,  1.48it/s, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:04<00:03,  1.48it/s, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:04<00:02,  1.48it/s, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:58<01:17, 25.91s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:58<01:17, 25.91s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:58<00:39, 19.68s/it, Batch correction: graph_connectivity][A
    Embeddings:  25%|[32m██▌       [0m| 2/8 [05:58<16:50, 168.35s/it]tch correction: pcr_comparison]    [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.35it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.35it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:20,  2.51s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:20,  2.51s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:05<00:11,  1.67s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:05<00:11,  1.67s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:05<00:10,  1.67s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:05<00:04,  1.04it/s, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:05<00:04,  1.04it/s, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:06<00:03,  1.04it/s, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:57<01:16, 25.38s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:57<01:16, 25.38s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:57<00:38, 19.28s/it, Batch correction: graph_connectivity][A
    Embeddings:  38%|[32m███▊      [0m| 3/8 [07:55<12:05, 145.06s/it]tch correction: pcr_comparison]    [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.36it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.36it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:18,  2.29s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:18,  2.29s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.55s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.55s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:04<00:09,  1.55s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:05<00:04,  1.11it/s, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:05<00:04,  1.11it/s, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:05<00:03,  1.11it/s, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:13<00:47, 15.77s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:13<00:47, 15.77s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:14<00:23, 11.99s/it, Batch correction: graph_connectivity][A
    Embeddings:  50%|[32m█████     [0m| 4/8 [09:09<07:48, 117.05s/it]tch correction: pcr_comparison]    [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.31it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:06,  1.31it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:03<00:14,  1.76s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:03<00:14,  1.76s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:03<00:08,  1.26s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:03<00:08,  1.26s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:03<00:07,  1.26s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:04<00:03,  1.37it/s, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:04<00:03,  1.37it/s, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:04<00:02,  1.37it/s, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:02<00:40, 13.38s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:02<00:40, 13.38s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:02<00:20, 10.18s/it, Batch correction: graph_connectivity][A
    Embeddings:  62%|[32m██████▎   [0m| 5/8 [10:12<04:52, 97.48s/it] tch correction: pcr_comparison]    [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:07,  1.21it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:07,  1.21it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:18,  2.29s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:18,  2.29s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.50s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.50s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:04<00:08,  1.50s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:12<00:14,  2.96s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:12<00:14,  2.96s/it, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:12<00:11,  2.96s/it, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:28<00:55, 18.56s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:28<00:55, 18.56s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:28<00:28, 14.14s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:28<00:28, 14.14s/it, Batch correction: pcr_comparison]    [A
    Embeddings:  75%|[32m███████▌  [0m| 6/8 [11:42<03:09, 94.98s/it]atch correction: pcr_comparison][A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:07,  1.13it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:00<00:07,  1.13it/s, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:03<00:17,  2.15s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:03<00:17,  2.15s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:09,  1.41s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:09,  1.41s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:04<00:08,  1.41s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:14<00:16,  3.35s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:14<00:16,  3.35s/it, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:14<00:13,  3.35s/it, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:18<00:48, 16.33s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:18<00:48, 16.33s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:19<00:24, 12.42s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:19<00:24, 12.42s/it, Batch correction: pcr_comparison]    [A
    Embeddings:  88%|[32m████████▊ [0m| 7/8 [13:02<01:30, 90.13s/it]atch correction: pcr_comparison][A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s][A
                                                                                             [A
    Metrics:   0%|[34m          [0m| 0/10 [00:00<?, ?it/s, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:01<00:09,  1.07s/it, Bio conservation: isolated_labels][A
    Metrics:  10%|[34m█         [0m| 1/10 [00:01<00:09,  1.07s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:17,  2.20s/it, Bio conservation: nmi_ari_cluster_labels_kmeans][A
    Metrics:  20%|[34m██        [0m| 2/10 [00:04<00:17,  2.20s/it, Bio conservation: silhouette_label]             [A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.46s/it, Bio conservation: silhouette_label][A
    Metrics:  30%|[34m███       [0m| 3/10 [00:04<00:10,  1.46s/it, Bio conservation: clisi_knn]       [A
    Metrics:  40%|[34m████      [0m| 4/10 [00:04<00:08,  1.46s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:15<00:18,  3.62s/it, Batch correction: silhouette_batch][A
    Metrics:  50%|[34m█████     [0m| 5/10 [00:15<00:18,  3.62s/it, Batch correction: ilisi_knn]       [A
    Metrics:  60%|[34m██████    [0m| 6/10 [00:15<00:14,  3.62s/it, Batch correction: kbet_per_label][A

    [34mINFO    [0m NK CD158e1+ consists of a single batch or is too small. Skip.



    Metrics:  70%|[34m███████   [0m| 7/10 [01:19<00:48, 16.31s/it, Batch correction: kbet_per_label][A
    Metrics:  70%|[34m███████   [0m| 7/10 [01:19<00:48, 16.31s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:19<00:24, 12.41s/it, Batch correction: graph_connectivity][A
    Metrics:  80%|[34m████████  [0m| 8/10 [01:19<00:24, 12.41s/it, Batch correction: pcr_comparison]    [A
    Embeddings: 100%|[32m██████████[0m| 8/8 [14:23<00:00, 107.89s/it]tch correction: pcr_comparison][A

                                                                                             [A

```python
bm.plot_results_table(min_max_scale=False)
```

![png](output_52_0.png)

    <plottable.table.Table at 0x7f33707f6620>

We can find that harmony removes the batch effect the best of the three methods that do not use the GPU, scVI is method to remove batch effect using GPU.
