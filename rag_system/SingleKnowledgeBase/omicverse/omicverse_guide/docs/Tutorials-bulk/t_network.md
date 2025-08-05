# t_network
*Converted from: omicverse/omicverse_guide/docs/Tutorials-bulk/t_network.ipynb*

# Protein-Protein interaction (PPI) analysis by String-db
STRING is a database of known and predicted protein-protein interactions. The interactions include direct (physical) and indirect (functional) associations; they stem from computational prediction, from knowledge transfer between organisms, and from interactions aggregated from other (primary) databases.

Here we produce a tutorial that use python to construct protein-protein interaction network

Colab_Reproducibility：https://colab.research.google.com/drive/1ReLCFA5cNNcem_WaMXYN9da7W0GN4gzl?usp=sharing

```pythonimport omicverse as ov
ov.utils.ov_plot_set()```

## Prepare data

Here we use the example data of string-db to perform the analysis 

FAA4 and its ten most confident interactors.
FAA4 in yeast is a long chain fatty acyl-CoA synthetase; see it connected to other synthetases as well as regulators.

Saccharomyces cerevisiae
NCBI taxonomy Id: 4932
Other names: ATCC 18824, Candida robusta, NRRL Y-12632, S. cerevisiae, Saccharomyces capensis, Saccharomyces italicus, Saccharomyces oviformis, Saccharomyces uvarum var. melibiosus, lager beer yeast, yeast

```pythongene_list=['FAA4','POX1','FAT1','FAS2','FAS1','FAA1','OLE1','YJU3','TGL3','INA1','TGL5']```

Besides, we also need to set the gene's type and color. Here, we randomly set the top 5 genes named Type1, other named Type2

```pythongene_type_dict=dict(zip(gene_list,['Type1']*5+['Type2']*6))
gene_color_dict=dict(zip(gene_list,['#F7828A']*5+['#9CCCA4']*6))```

## STRING interaction analysis

The network API method also allows you to retrieve your STRING interaction network for one or multiple proteins in various text formats. It will tell you the combined score and all the channel specific scores for the set of proteins. You can also extend the network neighborhood by setting "add_nodes", which will add, to your network, new interaction partners in order of their confidence.

```pythonG_res=ov.bulk.string_interaction(gene_list,4932)
G_res.head()```
*Output:*
```     stringId_A    stringId_B preferredName_A preferredName_B ncbiTaxonId  \
0  4932.YBR041W  4932.YKL182W            FAT1            FAS1        4932   
1  4932.YBR041W  4932.YKL182W            FAT1            FAS1        4932   
2  4932.YBR041W  4932.YPL231W            FAT1            FAS2        4932   
3  4932.YBR041W  4932.YPL231W            FAT1            FAS2        4932   
4  4932.YBR041W  4932.YOR081C            FAT1            TGL5        4932   

   score nscore fscore pscore ascore escore dscore tscore  
0   0.69      0      0      0      0      0      0   0.69  
1   0.69      0      0      0      0      0      0   0.69  
2  0.692      0      0      0      0      0      0  0.692  
3  0.692      0      0      0      0      0      0  0.692  
4    0.7      0      0      0      0      0      0    0.7  ```

## STRING PPI network 

We also can use `ov.bulk.pyPPI` to get the PPI network of `gene_list`, we init it at first

```pythonppi=ov.bulk.pyPPI(gene=gene_list,
                      gene_type_dict=gene_type_dict,
                      gene_color_dict=gene_color_dict,
                      species=4932)
```
*Output:*
```<networkx.classes.graph.Graph at 0x1056abb80>```

Then we connect to string-db to calculate the protein-protein interaction

```pythonppi.interaction_analysis()```
*Output:*
```<networkx.classes.graph.Graph at 0x17aa86fd0>```

We provided a very simple function to plot the network, you can refer the `ov.utils.plot_network` to find out the parameter

```pythonppi.plot_network()```
*Output:*
```(<Figure size 320x320 with 1 Axes>, <AxesSubplot: >)```
```<Figure size 320x320 with 1 Axes>```

