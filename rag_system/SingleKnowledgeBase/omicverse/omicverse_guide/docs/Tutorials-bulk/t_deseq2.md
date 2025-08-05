# t_deseq2
*Converted from: omicverse/omicverse_guide/docs/Tutorials-bulk/t_deseq2.ipynb*

# Different Expression Analysis with DEseq2

An important task of bulk rna-seq analysis is the different expression , which we can perform with omicverse. For different expression analysis, ov change the `gene_id` to `gene_name` of matrix first. 

Now we can use `PyDEseq2` to perform DESeq2 analysis like R

Paper: [PyDESeq2: a python package for bulk RNA-seq differential expression analysis](https://www.biorxiv.org/content/10.1101/2022.12.14.520412v1)

Code: https://github.com/owkin/PyDESeq2

Colab_Reproducibility：https://colab.research.google.com/drive/1fZS-v0zdIYkXrEoIAM1X5kPoZVfVvY5h?usp=sharing

```pythonimport omicverse as ov
ov.utils.ov_plot_set()```
*Output:*
```/Users/fernandozeng/miniforge3/envs/scbasset/lib/python3.8/site-packages/phate/__init__.py
```

Note that this dataset has not been processed in any way and is only exported by `featureCounts`, and Sequence alignment was performed from the genome file of CRCm39

```pythondata=ov.utils.read('https://raw.githubusercontent.com/Starlitnightly/Pyomic/master/sample/counts.txt',index_col=0,header=1)
#replace the columns `.bam` to `` 
data.columns=[i.split('/')[-1].replace('.bam','') for i in data.columns]
data.head()```
*Output:*
```                    1--1  1--2  2--1  2--2  3--1  3--2  4--1  4--2  4-3  4-4  \
Geneid                                                                         
ENSMUSG00000102628     0     0     0     0     5     0     0     0    0    0   
ENSMUSG00000100595     0     0     0     0     0     0     0     0    0    0   
ENSMUSG00000097426     5     0     0     0     0     0     0     1    0    0   
ENSMUSG00000104478     0     0     0     0     0     0     0     0    0    0   
ENSMUSG00000104385     0     0     0     0     0     0     0     0    0    0   

                    Blank-1  Blank-2  
Geneid                                
ENSMUSG00000102628        0        9  
ENSMUSG00000100595        0        0  
ENSMUSG00000097426        0        0  
ENSMUSG00000104478        0        0  
ENSMUSG00000104385        0        0  ```

## ID mapping

We performed the gene_id mapping by the mapping pair file `GRCm39` downloaded before.

```pythonov.utils.download_geneid_annotation_pair()```

```pythondata=ov.bulk.Matrix_ID_mapping(data,'genesets/pair_GRCm39.tsv')
data.head()```
*Output:*
```         1--1  1--2  2--1  2--2  3--1  3--2  4--1  4--2   4-3   4-4  Blank-1  \
Gm14845   115   116    84    86   133   170   130   105    91   127      124   
Hdc        97   579   123   172   571   119   106    28   217   156        2   
H2bu2      60    59    58    22    71    73    75   138    55    38       18   
Gm6693      0     0     0     0     0     0     0     0     0     0        0   
Rnd3     2423  2289  1996  1750  2304  2669  2952  2109  2030  2026      875   

         Blank-2  
Gm14845       99  
Hdc           51  
H2bu2         53  
Gm6693         0  
Rnd3        2555  ```

## Different expression analysis with ov

We can do differential expression analysis very simply by ov, simply by providing an expression matrix. To run DEG, we simply need to:

- Read the raw count by featureCount or any other qualify methods.
- Create an ov DEseq object.

```pythondds=ov.bulk.pyDEG(data)```

We notes that the gene_name mapping before exist some duplicates, we will process the duplicate indexes to retain only the highest expressed genes

```pythondds.drop_duplicates_index()
print('... drop_duplicates_index success')```
*Output:*
```... drop_duplicates_index success
```

Now we can calculate the different expression gene from matrix, we need to input the treatment and control groups

```pythontreatment_groups=['4-3','4-4']
control_groups=['1--1','1--2']
result=dds.deg_analysis(treatment_groups,control_groups,method='DEseq2')
```
*Output:*
```Fitting size factors...
... done in 0.00 seconds.

Fitting dispersions...
... done in 1.59 seconds.

Fitting dispersion trend curve...
... done in 2.82 seconds.

logres_prior=1.1538905878789707, sigma_prior=0.25
Fitting MAP dispersions...
... done in 1.57 seconds.

Fitting LFCs...
... done in 1.27 seconds.

Refitting 0 outliers.

Running Wald tests...
... done in 1.33 seconds.

Log2 fold change & Wald test p-value: condition Treatment vs Control
```
```            baseMean  log2FoldChange     lfcSE      stat    pvalue      padj
Gm14845   111.727600       -0.049168  0.470660 -0.104467  0.916799  0.975241
Hdc       258.120455       -0.809097  1.116541 -0.724646  0.468669  0.789482
H2bu2      52.656807       -0.323968  0.652995 -0.496127  0.619805  0.877166
Gm6693      0.000000             NaN       NaN       NaN       NaN       NaN
Rnd3     2180.318184       -0.183828  0.190533 -0.964809  0.334641  0.690369
...              ...             ...       ...       ...       ...       ...
Gm18244     0.000000             NaN       NaN       NaN       NaN       NaN
Gm50317     0.000000             NaN       NaN       NaN       NaN       NaN
Olfr516     0.000000             NaN       NaN       NaN       NaN       NaN
Gm37042     0.000000             NaN       NaN       NaN       NaN       NaN
Prelid1  3335.335908       -0.032464  0.190413 -0.170493  0.864622  0.958729

[54504 rows x 6 columns]```

One important thing is that we do not filter out low expression genes when processing DEGs, and in future versions I will consider building in the corresponding processing.

```pythonprint(result.shape)
result=result.loc[result['log2(BaseMean)']>1]
print(result.shape)```
*Output:*
```(54504, 14)
(23377, 14)
```

We also need to set the threshold of Foldchange, we prepare a method named `foldchange_set` to finish. This function automatically calculates the appropriate threshold based on the log2FC distribution, but you can also enter it manually.

```python# -1 means automatically calculates
dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=10)```
*Output:*
```... Fold change threshold: 1.6248531033651643
```

## Visualize the DEG result and specific genes

To visualize the DEG result, we use `plot_volcano` to do it. This fuction can visualize the gene interested or high different expression genes. There are some parameters you need to input:

- title: The title of volcano
- figsize: The size of figure
- plot_genes: The genes you interested
- plot_genes_num: If you don't have interested genes, you can auto plot it.

```pythondds.plot_volcano(title='DEG Analysis',figsize=(4,4),
                 plot_genes_num=8,plot_genes_fontsize=12,)```
*Output:*
```<Axes: title={'center': 'DEG Analysis'}, xlabel='$log_{2}FC$', ylabel='$-log_{10}(qvalue)$'>```
```<Figure size 320x320 with 1 Axes>```

To visualize the specific genes, we only need to use the `dds.plot_boxplot` function to finish it.

```pythondds.plot_boxplot(genes=['Ckap2','Lef1'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))```
*Output:*
```(<Figure size 160x240 with 1 Axes>,
 <Axes: title={'center': 'Gene Expression'}>)```
```<Figure size 160x240 with 1 Axes>```

```pythondds.plot_boxplot(genes=['Ckap2'],treatment_groups=treatment_groups,
                control_groups=control_groups,figsize=(2,3),fontsize=12,
                 legend_bbox=(2,0.55))```
*Output:*
```(<Figure size 160x240 with 1 Axes>,
 <Axes: title={'center': 'Gene Expression'}>)```
```<Figure size 160x240 with 1 Axes>```

## Pathway enrichment analysis by Pyomic

Here we use the `gseapy` package, which included the GSEA analysis and Enrichment. We have optimised the output of the package and given some better looking graph drawing functions

Similarly, we need to download the pathway/genesets first. Five genesets we prepare previously, you can use `Pyomic.utils.download_pathway_database()` to download automatically. Besides, you can download the pathway you interested from enrichr: https://maayanlab.cloud/Enrichr/#libraries

```pythonov.utils.download_pathway_database()```
*Output:*
```......Pathway Geneset download start: GO_Biological_Process_2021
......Loading dataset from genesets/GO_Biological_Process_2021.txt
......Pathway Geneset download start: GO_Cellular_Component_2021
......Loading dataset from genesets/GO_Cellular_Component_2021.txt
......Pathway Geneset download start: GO_Molecular_Function_2021
......Loading dataset from genesets/GO_Molecular_Function_2021.txt
......Pathway Geneset download start: WikiPathway_2021_Human
......Loading dataset from genesets/WikiPathway_2021_Human.txt
......Pathway Geneset download start: WikiPathways_2019_Mouse
......Loading dataset from genesets/WikiPathways_2019_Mouse.txt
......Pathway Geneset download start: Reactome_2022
......Loading dataset from genesets/Reactome_2022.txt
......Pathway Geneset download finished!
......Other Genesets can be dowload in `https://maayanlab.cloud/Enrichr/#libraries`
```

```pythonpathway_dict=ov.utils.geneset_prepare('genesets/WikiPathways_2019_Mouse.txt',organism='Mouse')```

To perform the GSEA analysis, we need to ranking the genes at first. Using `dds.ranking2gsea` can obtain a ranking gene's matrix sorted by -log10(padj).

$Metric=\frac{-log_{10}(padj)}{sign(log2FC)}$

```pythonrnk=dds.ranking2gsea()```

We used `ov.bulk.pyGSEA` to construst a GSEA object to perform enrichment.

```pythongsea_obj=ov.bulk.pyGSEA(rnk,pathway_dict)```

```pythonenrich_res=gsea_obj.enrichment()```
*Output:*
```2023-05-18 03:12:10,455 Input gene rankings contains NA values(gene name and ranking value), drop them all!
```

The results are stored in the `enrich_res` attribute.

```pythongsea_obj.enrich_res.head()```
*Output:*
```                                                 es       nes  pval       fdr  \
Term                                                                            
Complement and Coagulation Cascades WP449  0.732116  2.140070   0.0  0.000000   
Matrix Metalloproteinases WP441            0.879498  2.397240   0.0  0.000000   
TYROBP Causal Network WP3625               0.786372  2.358131   0.0  0.000000   
PPAR signaling pathway WP2316              0.681572  2.074737   0.0  0.003011   
Metapathway biotransformation WP1251       0.643937  1.991519   0.0  0.012042   

                                           geneset_size  matched_size  \
Term                                                                    
Complement and Coagulation Cascades WP449            62            56   
Matrix Metalloproteinases WP441                      29            27   
TYROBP Causal Network WP3625                         58            57   
PPAR signaling pathway WP2316                        81            69   
Metapathway biotransformation WP1251                141           120   

                                                                                       genes  \
Term                                                                                           
Complement and Coagulation Cascades WP449  Cfd;Masp1;F2r;C4b;Hc;Cfh;F7;F12;Pros1;Serping1...   
Matrix Metalloproteinases WP441            Mmp11;Mmp14;Mmp3;Mmp12;Timp4;Timp1;Mmp28;Mmp9;...   
TYROBP Causal Network WP3625               Itgax;Itgb2;Rgs1;Gpx1;Lhfpl2;Tcirg1;Cxcl16;Cd3...   
PPAR signaling pathway WP2316              Hmgcs2;Pck1;Slc27a1;Scd3;Acox3;Acsbg1;Scd1;Ang...   
Metapathway biotransformation WP1251       Cyp26b1;Cyp2e1;Fmo2;Gpx1;Cyp4b1;Cyp11a1;Mgst2;...   

                                                                                 ledge_genes  \
Term                                                                                           
Complement and Coagulation Cascades WP449  Cfd;Masp1;F2r;C4b;Hc;Cfh;F7;F12;Pros1;Serping1...   
Matrix Metalloproteinases WP441            Mmp11;Mmp14;Mmp3;Mmp12;Timp4;Timp1;Mmp28;Mmp9;...   
TYROBP Causal Network WP3625               Itgax;Itgb2;Rgs1;Gpx1;Lhfpl2;Tcirg1;Cxcl16;Cd3...   
PPAR signaling pathway WP2316              Hmgcs2;Pck1;Slc27a1;Scd3;Acox3;Acsbg1;Scd1;Ang...   
Metapathway biotransformation WP1251       Cyp26b1;Cyp2e1;Fmo2;Gpx1;Cyp4b1;Cyp11a1;Mgst2;...   

                                               logp      logc  num  fraction  \
Term                                                                           
Complement and Coagulation Cascades WP449  9.210340  2.140070   56  0.903226   
Matrix Metalloproteinases WP441            9.210340  2.397240   27  0.931034   
TYROBP Causal Network WP3625               9.210340  2.358131   57  0.982759   
PPAR signaling pathway WP2316              5.772960  2.074737   69  0.851852   
Metapathway biotransformation WP1251       4.411073  1.991519  120  0.851064   

                                                                                Term  \
Term                                                                                   
Complement and Coagulation Cascades WP449  Complement and Coagulation Cascades WP449   
Matrix Metalloproteinases WP441                      Matrix Metalloproteinases WP441   
TYROBP Causal Network WP3625                            TYROBP Causal Network WP3625   
PPAR signaling pathway WP2316                          PPAR signaling pathway WP2316   
Metapathway biotransformation WP1251            Metapathway biotransformation WP1251   

                                            P-value  
Term                                                 
Complement and Coagulation Cascades WP449  0.000000  
Matrix Metalloproteinases WP441            0.000000  
TYROBP Causal Network WP3625               0.000000  
PPAR signaling pathway WP2316              0.003011  
Metapathway biotransformation WP1251       0.012042  ```

To visualize the enrichment, we use `plot_enrichment` to do.
- num: The number of enriched terms to plot. Default is 10.
- node_size: A list of integers defining the size of nodes in the plot. Default is [5,10,15].
- cax_loc: The location of the colorbar on the plot. Default is 2.
- cax_fontsize: The fontsize of the colorbar label. Default is 12.
- fig_title: The title of the plot. Default is an empty string.
- fig_xlabel: The label of the x-axis. Default is 'Fractions of genes'.
- figsize: The size of the plot. Default is (2,4).
- cmap: The colormap to use for the plot. Default is 'YlGnBu'.

```pythongsea_obj.plot_enrichment(num=10,node_size=[10,20,30],
                        cax_fontsize=12,
                        fig_title='Wiki Pathway Enrichment',fig_xlabel='Fractions of genes',
                        figsize=(2,4),cmap='YlGnBu',
                        text_knock=2,text_maxsize=30,
                        cax_loc=[2.5, 0.45, 0.5, 0.02],
                          bbox_to_anchor_used=(-0.25, -13),node_diameter=10,)```
*Output:*
```<Axes: title={'center': 'Wiki Pathway Enrichment'}, xlabel='Fractions of genes'>```
```<Figure size 160x320 with 2 Axes>```

Not only the basic analysis, pyGSEA also help us to visualize the term with Ranked and Enrichment Score. 

We can select the number of term to plot, which stored in `gsea_obj.enrich_res.index`, the `0` is `Complement and Coagulation Cascades WP449` and the `1` is `Matrix Metalloproteinases WP441`

```pythongsea_obj.enrich_res.index[:5]```
*Output:*
```Index(['Complement and Coagulation Cascades WP449',
       'Matrix Metalloproteinases WP441', 'TYROBP Causal Network WP3625',
       'PPAR signaling pathway WP2316',
       'Metapathway biotransformation WP1251'],
      dtype='object', name='Term')```

We can set the `gene_set_title` to change the title of GSEA plot

```pythonfig=gsea_obj.plot_gsea(term_num=1,
                  gene_set_title='Matrix Metalloproteinases',
                  figsize=(3,4),
                  cmap='RdBu_r',
                  title_fontsize=14,
                  title_y=0.95)```
*Output:*
```<Figure size 240x320 with 4 Axes>```

