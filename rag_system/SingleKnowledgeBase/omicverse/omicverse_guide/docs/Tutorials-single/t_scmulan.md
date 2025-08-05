# t_scmulan
*Converted from: omicverse/omicverse_guide/docs/Tutorials-single/t_scmulan.ipynb*

## Using scMulan to annotate cell types in Heart, Lung, Liver, Bone marrow, Blood, Brain, and Thymus

In this study, the authors enrich the pre-training paradigm by integrating an abundance of metadata and a multiplicity of pre-training tasks, and obtain scMulan, a multitask generative pre-trained language model tailored for single-cell analysis. They represent a cell as a structured cell sentence (c-sentence) by encoding its gene expression, metadata terms, and target tasks as words of tuples, each consisting of entities and their corresponding values. They construct a unified generative framework to model the cell language on c-sentence and design three pretraining tasks to bridge the microscopic and macroscopic information within the c-sentences. They pre-train scMulan on 10 million single-cell transcriptomic data and their corresponding metadata, with 368 million parameters. As a single model, scMulan can accomplish tasks zero-shot for cell type annotation, batch integration, and conditional cell generation, guided by different task prompts.

#### we provide a liver dataset sampled (percentage of 20%) from Suo C, 2022 (doi/10.1126/science.abo0510)
**Paper:** [scMulan: a multitask generative pre-trained language model for single-cell analysis](https://www.biorxiv.org/content/10.1101/2024.01.25.577152v1)    
**Data download:** https://cloud.tsinghua.edu.cn/f/45a7fd2a27e543539f59/?dl=1   
**Pre-train model download:** https://cloud.tsinghua.edu.cn/f/2250c5df51034b2e9a85/?dl=1

If you found this tutorial helpful, please cite scMulan and OmicVerse:         
Bian H, Chen Y, Dong X, et al. scMulan: a multitask generative pre-trained language model for single-cell analysis[C]//International Conference on Research in Computational Molecular Biology. Cham: Springer Nature Switzerland, 2024: 479-482.

```pythonimport os
#os.environ["CUDA_VISIBLE_DEVICES"] = "-1" # if use CPU only
import scanpy as sc
import omicverse as ov
ov.plot_set()
#import scMulan
#from scMulan import GeneSymbolUniform```
*Output:*
```Dependency error: (pydeseq2 0.4.1 (/opt/miniforge/envs/omicverse/lib/python3.10/site-packages), Requirement.parse('pydeseq2<=0.4.0,>=0.3'))
```

## 1. load h5ad
You can download the liver dataset from the following link: https://cloud.tsinghua.edu.cn/f/45a7fd2a27e543539f59/?dl=1 

It's recommended that you use h5ad here with raw count (and after your QC)


```pythonadata = sc.read('./data/liver_test.h5ad')```

```pythonadata```
*Output:*
```AnnData object with n_obs × n_vars = 27436 × 43878
    obs: 'cid', 'seq_tech', 'donor_ID', 'donor_gender', 'donor_age', 'donor_status', 'original_name', 'organ', 'region', 'subregion', 'sample_status', 'treatment', 'ethnicity', 'cell_type', 'cell_id', 'study_id'
    var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable'
    obsm: 'umap'```

```python
from scipy.sparse import csc_matrix
adata.X = csc_matrix(adata.X)```

## 2. transform original h5ad with uniformed genes (42117 genes)

This step transform the genes in input adata to 42117 gene symbols and reserves the corresponding gene expression values. The gene symbols are the same as the pre-trained scMulan model.

```pythonadata_GS_uniformed = ov.external.scMulan.GeneSymbolUniform(input_adata=adata,
                                 output_dir="./data",
                                 output_prefix='liver')```
*Output:*
```{message}
The shape of query data is: (27436, 43878)
The length of reference gene_list is: 42117
Performing gene symbol uniform, this step may take several minutes
```
```Processing: 100%|██████████████████████████████████| 43878/43878 [00:28<00:00, 1565.30it/s]
```
```Building output data, this step may take several minutes
```
```Processing: 100%|██████████████████████████████████| 42117/42117 [00:12<00:00, 3304.21it/s]
```
```Shape of output data is (27436, 42117). It should have 42117 genes with cell number unchanged.
h5ad file saved in:/data/hulei/Projects/omicverse_scripts/tutorials/data/liver_uniformed.h5ad
report file saved in: /data/hulei/Projects/omicverse_scripts/tutorials/data/liver_report.csv
```

## 3. process uniformed data (simply norm and log1p)

```python## you can read the saved uniformed adata

adata_GS_uniformed=sc.read_h5ad('./data/liver_uniformed.h5ad')```

```pythonadata_GS_uniformed```
*Output:*
```AnnData object with n_obs × n_vars = 27436 × 42117
    obs: 'cid', 'seq_tech', 'donor_ID', 'donor_gender', 'donor_age', 'donor_status', 'original_name', 'organ', 'region', 'subregion', 'sample_status', 'treatment', 'ethnicity', 'cell_type', 'cell_id', 'study_id'```

```python# norm and log1p count matrix
# in some case, the count matrix is not normalized, and log1p is not applied.
# So we need to normalize the count matrix
if adata_GS_uniformed.X.max() > 10:
    sc.pp.normalize_total(adata_GS_uniformed, target_sum=1e4) 
    sc.pp.log1p(adata_GS_uniformed)```

## 4. load scMulan

```python# you should first download ckpt from https://cloud.tsinghua.edu.cn/f/2250c5df51034b2e9a85/?dl=1
# put it under .ckpt/ckpt_scMulan.pt
# by: wget https://cloud.tsinghua.edu.cn/f/2250c5df51034b2e9a85/?dl=1  -O ckpt/ckpt_scMulan.pt

ckp_path = './ckpt/ckpt_scMulan.pt'```

```pythonscml = ov.external.scMulan.model_inference(ckp_path, adata_GS_uniformed)
base_process = scml.cuda_count()```
*Output:*
```number of parameters: 368.80M
✅ adata passed check
👸 scMulan is ready
scMulan is currently available to 8 GPUs.
```

```pythonscml.get_cell_types_and_embds_for_adata(parallel=True, n_process = 1)
# scml.get_cell_types_and_embds_for_adata(parallel=False) # for only using CPU, but it is really slow.```
*Output:*
```⚡ Speed up by multiprocessing with 1 processes and 8 GPUs...
```
```/opt/miniforge/envs/omicverse/lib/python3.10/site-packages/torch_geometric/typing.py:68: UserWarning: An issue occurred while importing 'pyg-lib'. Disabling its usage. Stacktrace: /opt/miniforge/envs/omicverse/lib/python3.10/site-packages/libpyg.so: undefined symbol: _ZN3c1010Dispatcher17runRecordFunctionERN2at14RecordFunctionESt17reference_wrapperIKNS_14FunctionSchemaEENS_11DispatchKeyE
  warnings.warn(f"An issue occurred while importing 'pyg-lib'. "
2024-10-28 05:34:49.934679: E external/local_xla/xla/stream_executor/cuda/cuda_fft.cc:485] Unable to register cuFFT factory: Attempting to register factory for plugin cuFFT when one has already been registered
2024-10-28 05:34:49.948925: E external/local_xla/xla/stream_executor/cuda/cuda_dnn.cc:8454] Unable to register cuDNN factory: Attempting to register factory for plugin cuDNN when one has already been registered
2024-10-28 05:34:49.953191: E external/local_xla/xla/stream_executor/cuda/cuda_blas.cc:1452] Unable to register cuBLAS factory: Attempting to register factory for plugin cuBLAS when one has already been registered
2024-10-28 05:34:51.088012: W tensorflow/compiler/tf2tensorrt/utils/py_utils.cc:38] TF-TRT Warning: Could not find TensorRT
```
```
   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

Version: 1.6.8, Tutorials: https://omicverse.readthedocs.io/
```
```⏳ Generating cell type labels and embds for each cell on device 0: 100%|██████████| 27436/27436 [13:26<00:00, 34.02it/s]
```

The predicted cell types are stored in scml.adata.obs['cell_type_from_scMulan'], besides the cell embeddings (for multibatch integration) in scml.adata.obsm['X_scMulan'] (not used in this tutorial).

## 5. visualization

Here, we visualize the cell types predicted by scMulan. And we also visualize the original cell types in the dataset.

```pythonadata_mulan = scml.adata.copy()```

```python# calculated the 2-D embedding of the adata using pyMDE
ov.pp.scale(adata_mulan)
ov.pp.pca(adata_mulan)

#sc.pl.pca_variance_ratio(adata_mulan)
ov.pp.mde(adata_mulan,embedding_dim=2,n_neighbors=15, basis='X_mde',
          n_pcs=10, use_rep='scaled|original|X_pca',)```
*Output:*
```... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
computing PCA
    with n_comps=50
    finished (0:00:06)
computing neighbors
    finished: added to `.uns['neighbors']`
    `.obsm['X_mde']`, MDE coordinates
    `.obsp['neighbors_distances']`, distances for each pair of neighbors
    `.obsp['neighbors_connectivities']`, weighted adjacency matrix (0:00:31)
```

```python# Here, we can see the cell type annotation from scMulan
ov.pl.embedding(adata_mulan,basis='X_mde',
                color=["cell_type_from_scMulan",],
                ncols=1,frameon='small')```
*Output:*
```<Figure size 320x320 with 1 Axes>```

```pythonadata_mulan.obsm['X_umap']=adata_mulan.obsm['X_mde']```

```python# you can run smoothing function to filter the false positives
ov.external.scMulan.cell_type_smoothing(adata_mulan, threshold=0.1)```
*Output:*
```computing neighbors
    finished: added to `.uns['Smoothing']`
    `.obsp['Smoothing_distances']`, distances for each pair of neighbors
    `.obsp['Smoothing_connectivities']`, weighted adjacency matrix (0:00:03)
```
```100%|██████████████████████████████████████████████| 27436/27436 [00:12<00:00, 2257.40it/s]
```

```python# cell_type_from_mulan_smoothing: pred+smoothing
# cell_type: original annotations by the authors
ov.pl.embedding(adata_mulan,basis='X_mde',
                color=["cell_type_from_mulan_smoothing","cell_type"],
                ncols=1,frameon='small')```
*Output:*
```<Figure size 386.4x640 with 2 Axes>```

```pythonadata_mulan```
*Output:*
```AnnData object with n_obs × n_vars = 27436 × 2000
    obs: 'cid', 'seq_tech', 'donor_ID', 'donor_gender', 'donor_age', 'donor_status', 'original_name', 'organ', 'region', 'subregion', 'sample_status', 'treatment', 'ethnicity', 'cell_type', 'cell_id', 'study_id', 'cell_type_from_scMulan', 'cell_type_from_mulan_smoothing', 'smoothing_score', 'selected_celltype'
    uns: 'pca', 'scaled|original|pca_var_ratios', 'scaled|original|cum_sum_eigenvalues', 'neighbors', 'Smoothing', 'cell_type_from_scMulan_colors', 'cell_type_from_mulan_smoothing_colors', 'cell_type_colors', 'selected_celltype_colors'
    obsm: 'X_scMulan', 'X_pca', 'scaled|original|X_pca', 'X_mde', 'X_umap'
    varm: 'PCs', 'scaled|original|pca_loadings'
    layers: 'scaled', 'lognorm'
    obsp: 'distances', 'connectivities', 'Smoothing_distances', 'Smoothing_connectivities'```

```pythontop_celltypes = adata_mulan.obs.cell_type_from_scMulan.value_counts().index[:20]```

```python# you can select some cell types of interest (from scMulan's prediction) for visulization
# selected_cell_types = ["NK cell", "Kupffer cell", "Conventional dendritic cell 2"] # as example
selected_cell_types = top_celltypes
ov.external.scMulan.visualize_selected_cell_types(adata_mulan,selected_cell_types,smoothing=True)```
*Output:*
```selected_celltype
Kupffer cell                           9963
NK cell                                3029
Sinusoidal endothelial cell            2529
Hepatocyte                             2189
Conventional dendritic cell 2          1334
Erythroid progenitor cell              1234
Monocyte                               1111
Naive B cell                           1107
Pre-B cell                             1031
Immune                                  780
others                                  709
CD4 Treg                                645
Megakaryocyte                           500
Monocyte-dendritic progenitor (MDP)     418
Plasmacytoid dendritic cell             191
Dendritic cell                          173
CD4 T cell                              142
Classical monocyte                      135
Microglia                               124
Pro-B cell                               76
Memory CD8 T cell                        16
Name: count, dtype: int64
```
```<Figure size 320x320 with 1 Axes>```

