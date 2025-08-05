# t_tcga
*Converted from: omicverse/omicverse_guide/docs/Tutorials-bulk/t_tcga.ipynb*

# TCGA database preprocess

We often download patient survival data from the TCGA database for analysis in order to verify the importance of genes in cancer. However, the pre-processing of the TCGA database is often a headache. Here, we have introduced the TCGA module in ov, a way to quickly process the file formats we download from the TCGA database. We need to prepare 3 files as input:

- gdc_sample_sheet (`.tsv`): The `Sample Sheet` button of TCGA, and we can get `tsv` file from it
- gdc_download_files (`folder`): The `Download/Cart` button of TCGA, and we get `tar.gz` included all file you selected/
- clinical_cart (`folder`): The `Clinical` button of TCGA, and we can get `tar.gz` included all clinical of your files

```pythonimport omicverse as ov
import scanpy as sc
ov.plot_set()```
*Output:*
```
   ____            _     _    __                  
  / __ \____ ___  (_)___| |  / /__  _____________ 
 / / / / __ `__ \/ / ___/ | / / _ \/ ___/ ___/ _ \ 
/ /_/ / / / / / / / /__ | |/ /  __/ /  (__  )  __/ 
\____/_/ /_/ /_/_/\___/ |___/\___/_/  /____/\___/                                              

Version: 1.6.4, Tutorials: https://omicverse.readthedocs.io/
All dependencies are satisfied.
```

## TCGA counts read

Here, we use `ov.bulk.TCGA` to perform the `gdc_sample_sheet`, `gdc_download_files` and `clinical_cart` you download before. The raw count, fpkm and tpm matrix will be stored in anndata object

```python%%time
gdc_sample_sheep='data/TCGA_OV/gdc_sample_sheet.2024-07-05.tsv'
gdc_download_files='data/TCGA_OV/gdc_download_20240705_180129.081531'
clinical_cart='data/TCGA_OV/clinical.cart.2024-07-05'
aml_tcga=ov.bulk.pyTCGA(gdc_sample_sheep,gdc_download_files,clinical_cart)
aml_tcga.adata_init()```
*Output:*
```tcga module init success
...index init
... expression matrix init
...anndata construct
CPU times: user 1min 30s, sys: 4.16 s, total: 1min 34s
Wall time: 1min 34s
```

We can save the anndata object for the next use

```pythonaml_tcga.adata.write_h5ad('data/TCGA_OV/ov_tcga_raw.h5ad',compression='gzip')```

Note: Each time we read the anndata file, we need to initialize the TCGA object using three paths so that the subsequent TCGA functions such as survival analysis can be used properly

If you wish to create your own TCGA data, we have provided [sample data](https://figshare.com/ndownloader/files/47461946) here for download:

TCGA OV: https://figshare.com/ndownloader/files/47461946

```pythongdc_sample_sheep='data/TCGA_OV/gdc_sample_sheet.2024-07-05.tsv'
gdc_download_files='data/TCGA_OV/gdc_download_20240705_180129.081531'
clinical_cart='data/TCGA_OV/clinical.cart.2024-07-05'
aml_tcga=ov.bulk.pyTCGA(gdc_sample_sheep,gdc_download_files,clinical_cart)
aml_tcga.adata_read('data/TCGA_OV/ov_tcga_raw.h5ad')```
*Output:*
```tcga module init success
... anndata reading
```

## Meta init
As the TCGA reads the gene_id, we need to convert it to gene_name as well as adding basic information about the patient. Therefore we need to initialise the patient's meta information.

```pythonaml_tcga.adata_meta_init()```
*Output:*
```...anndata meta init ['gene_name', 'gene_type'] ['Case ID', 'Sample Type']
```
```AnnData object with n_obs × n_vars = 429 × 60664
    obs: 'Case ID', 'Sample Type'
    var: 'gene_name', 'gene_type', 'gene_id'
    layers: 'deseq_normalize', 'fpkm', 'tpm'```

## Survial init
We set up the path for Clinical earlier, but in fact we did not import the patient information in the previous process, we only initially determined the id of the patient's TCGA, so we attracted to initialize the clinical information

```pythonaml_tcga.survial_init()
aml_tcga.adata```
*Output:*
```AnnData object with n_obs × n_vars = 429 × 60664
    obs: 'Case ID', 'Sample Type', 'vital_status', 'days'
    var: 'gene_name', 'gene_type', 'gene_id'
    layers: 'deseq_normalize', 'fpkm', 'tpm'```

To visualize the gene you interested, we can use `survival_analysis` to finish it. 

```pythonaml_tcga.survival_analysis('MYC',layer='deseq_normalize',plot=True)```
*Output:*
```(0.9528907736380975, 0.328984565709021)```
```<Figure size 240x240 with 1 Axes>```

If you want to calculate the survival of all genes, you can also use the `survial_analysis_all` to finish it. It may calculate a lot of times.

```pythonaml_tcga.survial_analysis_all()
aml_tcga.adata```

Don't forget to save your result.

```pythonaml_tcga.adata.write_h5ad('data/TCGA_OV/ov_tcga_survial_all.h5ad',compression='gzip')```

