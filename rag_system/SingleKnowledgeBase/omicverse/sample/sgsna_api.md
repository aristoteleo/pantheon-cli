# sgsna_api
*Converted from: omicverse/sample/sgsna_api.ipynb*

```pythonimport pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from scipy.stats import norm
from scipy import stats
import networkx as nx
import datetime
import seaborn as sns
import pandas as pd
import seaborn as sns  #用于绘制热图的工具包
from scipy.cluster import hierarchy  #用于进行层次聚类，话层次聚类图的工具包
from scipy import cluster   
import matplotlib.pyplot as plt
from sklearn import decomposition as skldec #用于主成分分析降维的包```

```python#from dynamicTree import cutreeHybrid
from scipy.spatial.distance import pdist
import numpy as np
from scipy.cluster.hierarchy import linkage,dendrogram```

```pythonimport ERgene
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
%matplotlib inline
from scipy import stats```

```pythondata=pd.read_csv('LiverFemale3600.csv')
data.dropna(inplace=True)
data.set_index(data.columns[0],inplace=True)
data.head()```
*Output:*
```                F2_2    F2_3     F2_14    F2_15    F2_19     F2_20    F2_23  \
substanceBXH                                                                  
MMT00000044  -0.0181  0.0642  0.000064 -0.05800  0.04830 -0.151974 -0.00129   
MMT00000046  -0.0773 -0.0297  0.112000 -0.05890  0.04430 -0.093800  0.09340   
MMT00000051  -0.0226  0.0617 -0.129000  0.08710 -0.11500 -0.065026  0.00249   
MMT00000080  -0.0487  0.0582 -0.048300 -0.03710  0.02510  0.085043  0.04450   
MMT00000102   0.1760 -0.1890 -0.065000 -0.00846 -0.00574 -0.018072 -0.12500   

                F2_24   F2_26    F2_37    ...       F2_324  F2_325  F2_326  \
substanceBXH                              ...                                
MMT00000044  -0.23600 -0.0307 -0.02610    ...     0.047700 -0.0488  0.0168   
MMT00000046   0.02690 -0.1330  0.07570    ...    -0.049200 -0.0350 -0.0738   
MMT00000051  -0.10200  0.1420 -0.10200    ...     0.000612  0.1210  0.0996   
MMT00000080   0.00167 -0.0680  0.00567    ...     0.113000 -0.0859 -0.1340   
MMT00000102  -0.06820  0.1250  0.00998    ...    -0.080000 -0.1200  0.1230   

              F2_327   F2_328  F2_329  F2_330  F2_332  F2_355    F2_357  
substanceBXH                                                             
MMT00000044  -0.0309  0.02740 -0.0310  0.0660 -0.0199 -0.0146  0.065000  
MMT00000046  -0.1730 -0.07380 -0.2010 -0.0820 -0.0939  0.0192 -0.049900  
MMT00000051   0.1090  0.02730  0.1200 -0.0629 -0.0395  0.1090  0.000253  
MMT00000080   0.0639  0.00731  0.1240 -0.0212  0.0870  0.0512  0.024300  
MMT00000102   0.1870  0.05410  0.0699  0.0708  0.1450 -0.0399  0.037500  

[5 rows x 135 columns]```

```python#Identification of similar co-expression gene/protein module
table=soft_corr_matrix(data,
                      method='pearson')
#Select module from actual corr matrix
module=select_module(table,
                  linkage_method='ward',
                  minClusterSize=30,
                  deepSplit=2,
                 )```
*Output:*
```...correlation coefficient matrix is being calculated
...direction correlation have been saved
...calculate time 0:00:05.085976
...indirect correlation matrix is being calculated
...indirection correlation have been saved
...calculate time 0:00:00.452784
...soft_threshold is being calculated
...appropriate soft_thresholds: 6.0
...distance have being calculated
...geneTree have being calculated
...dynamicMods have being calculated
..cutHeight not given, setting it to 448.9691031625521  ===>  99% of the (truncated) height range in dendro.
```
```D:\Anaconda\lib\site-packages\pandas\core\series.py:842: FutureWarning: 
Passing list-likes to .loc or [] with any missing label will raise
KeyError in the future, you can use .reindex() as an alternative.

See the documentation here:
https://pandas.pydata.org/pandas-docs/stable/indexing.html#deprecate-loc-reindex-listlike
  return self.loc[key]
```
```..done.
...total: 12
```
```D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
```
```<Figure size 2000x1000 with 1 Axes>```
```<Figure size 2000x1000 with 2 Axes>```
```<Figure size 1800x720 with 2 Axes>```

```pythonmol=select_module(table1)```
*Output:*

```python#Correlation analysis between gene modules and traits
cha=pd.read_csv('newnewcri.csv')
cha.set_index(cha.columns[0],inplace=True)
co_character(data,
             character=cha,
             module=module)```
*Output:*
```...PCA analysis have being done
...co-analysis have being done
```
```    weight_g  length_cm    ab_fat  other_fat
1   0.000579   0.130891  0.025539   0.208714
2   0.370098   0.084585  0.295754   0.351596
3   0.285441   0.023763  0.281423   0.213187
4   0.150483   0.122483  0.205454   0.047039
5   0.244029   0.152480  0.288425   0.038686
6   0.005157   0.067767  0.027740   0.002952
7   0.007820   0.058246  0.027043   0.138332
8   0.299864   0.143134  0.273537   0.134064
9   0.179886   0.137660  0.230412   0.029496
10  0.203550   0.031582  0.146253   0.207015
11  0.676308   0.173676  0.549321   0.534850
12  0.577942   0.111693  0.495211   0.447296```
```<Figure size 720x720 with 2 Axes>```

```pythondata1=pd.read_csv('GSE_121.csv')
data1.set_index(data1.columns[0],inplace=True)
data1.columns
eg=['DD18_old121_rep1',
       'DD19_old121_rep2', 'DD20_old121_rep3', 'DD21_old121_rep4',
       'DD22_old121_rep5', 'DD23_old121_rep6']
cg=['DD6_old_rep1', 'DD7_old_rep2', 'DD8_old_rep3', 'DD9_old_rep4',
       'DD10_old_rep5', 'DD11_old_rep6']```

```pythondata1=pd.read_csv('GSE_121.csv')
data1.set_index(data1.columns[0],inplace=True)
#Screening of differentially expressed genes/proteins
deg=find_DEG(data1,
             eg=eg,
             cg=cg,
             fold_threshold=0,
             pvalue_threshold=0.05
            )```
*Output:*
```up: 234
down: 236
```
```<Figure size 2000x1000 with 2 Axes>```
```<Figure size 720x720 with 4 Axes>```

```python#method='pearson'/'kendall'/'spearman'
def Trans_corr_matrix(data,method='pearson'):
    data_len=len(data)
    data_index=data.index
    
    #correlation coefficient
    start = datetime.datetime.now()
    print('...correlation coefficient matrix is being calculated')
    result=data.T.corr(method)
    end = datetime.datetime.now()
    result.to_csv('direction correlation matrix.csv')
    print("...direction correlation have been saved")
    print("...calculate time",end-start)
    
    #indirect correlation add
    start = datetime.datetime.now()
    print('...indirect correlation matrix is being calculated')
    np.fill_diagonal(result.values, 0)
    arr=abs(result)
    arr_len=len(arr)
    temp=np.zeros(arr_len)
    for i in range(1,3):    
        temp=temp+(arr**i)/i
    end = datetime.datetime.now()
    temp.to_csv('indirection correlation matrix.csv')
    print("...indirection correlation have been saved")
    print("...calculate time",end-start)
    
        
    #cal soft_threshold
    plt.rc('font', family='Arial')
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.weight"] = "12"
    plt.rcParams["axes.labelweight"] = "bold"
    my_dpi=300
    fig=plt.figure(figsize=(2000/my_dpi, 1000/my_dpi), dpi=my_dpi)
    start = datetime.datetime.now()
    print('...soft_threshold is being calculated')
    soft=6
    re1=pd.DataFrame(columns=['beta','r2','meank'])
    for j in range(1,12):
        result_i=np.float_power(temp,j)
        tt_0=np.sum(abs(result_i),axis=0)-1
        n=plt.hist(x = tt_0), #
        x=n[0][0]
        y=[]
        for i in range(len(n[0][1])-1):
            y.append((n[0][1][i]+n[0][1][i+1])/2)
        x=np.log10(x)
        y=np.log10(y)
        res=stats.linregress(x, y)
        r2=np.float_power(res.rvalue,2)
        k=tt_0.mean()
        re1=re1.append({'beta':j,'r2':r2,'meank':k},ignore_index=True)
    for i in re1['r2']:
        if i>0.85:
            soft=re1[re1['r2']==i]['beta'].iloc[0]
            break
    print('...appropriate soft_thresholds:',soft)
    plt.savefig('soft_threshold_hist.png',dpi=300)
    
    #select soft_threhold
    my_dpi=300
    fig=plt.figure(figsize=(2000/my_dpi, 1000/my_dpi), dpi=my_dpi)
    grid = plt.GridSpec(1, 4, wspace=1, hspace=0.1)
    #fig, (ax0, ax1) = plt.subplots(2, 1)
    plt.subplot(grid[0,0:2])
    cmap=sns.color_palette("seismic")
    plt.subplot(1,2, 1)
    p1=sns.regplot(x=re1["beta"], y=re1['r2'], fit_reg=False, marker="o", color=cmap[0])
    p1.axhline(y=0.9,ls=":",c=cmap[5])

    plt.subplot(grid[0,2:4])
    p1=sns.regplot(x=re1["beta"], y=re1['meank'], fit_reg=False, marker="o", color=cmap[4])
    plt.savefig('soft_threshold_select.png',dpi=300)
    
    test=np.float_power(temp,soft)
    np.fill_diagonal(test.values, 1.0)
    return test
    
    
    ```

```pythondef select_module(data,
                  linkage_method='ward',
                  minClusterSize=30,
                  deepSplit=2,
                 ):
    #distance
    print("...distance have being calculated")
    distances = pdist(1-data, "euclidean")
    
    #geneTree
    print("...geneTree have being calculated")
    geneTree=linkage(distances, linkage_method)
      
    #dynamicMods
    print("...dynamicMods have being calculated")
    dynamicMods=cutreeHybrid(geneTree,distM=distances,minClusterSize = minClusterSize,deepSplit = deepSplit, pamRespectsDendro = False)
    print("...total:",len(set(dynamicMods['labels'])))
    
    plt.rc('font', family='Arial')
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.weight"] = "12"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.figure(figsize=(25, 10))
    grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.1)
    #fig, (ax0, ax1) = plt.subplots(2, 1)
    plt.subplot(grid[0:2,0])
    hierarchy.set_link_color_palette(['#000000'])
    dn=hierarchy.dendrogram(geneTree,color_threshold=0, above_threshold_color='black')
    plt.tick_params( \
        axis='x',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off')
    
    #mol
    x=dn['ivl']
    y=[dynamicMods['labels'][int(x)] for x in dn['ivl']]
    yy=np.array([y])
    z=[data.index[int(x)] for x in dn['ivl']]
    mol=pd.DataFrame(columns=['ivl','module','name'])
    mol['ivl']=x
    mol['module']=y
    mol['name']=z
    
    plt.subplot(grid[2,0])
    ax1=plt.pcolor(yy,cmap='seismic')

    plt.savefig('module_tree.png',dpi=300)
    return mol
    ```

```pythondef Analysis_cocharacter(data,character,module):
    print("...PCA analysis have being done")
    pcamol=pd.DataFrame(columns=data.columns)
    set_index=set(module['module'])
    for j in set_index:
        newdata=pd.DataFrame(columns=data.columns)
        for i in list(module[module['module']==j].dropna()['name']):
            newdata=newdata.append(data[data.index==i])
        from sklearn.decomposition import PCA
        pca = PCA(n_components=1) 
        reduced_X = pca.fit_transform(newdata.T)
        tepcamol=pd.DataFrame(reduced_X.T,columns=data.columns)
        pcamol=pcamol.append(tepcamol,ignore_index=True)
    pcamol.index=set_index
    
    print("...co-analysis have being done")
    from scipy.stats import spearmanr,pearsonr,kendalltau
    # seed random number generator
    # calculate spearman's correlation
    result_1=pd.DataFrame(columns=character.columns)
    result_p=pd.DataFrame(columns=character.columns)
    for j in character.columns:
        co=[]
        pvv=[]
        for i in range(len(pcamol)):   
            tempcor=pd.DataFrame(columns=['x','y'])
            tempcor['x']=list(character[j])
            tempcor['y']=list(pcamol.iloc[i])
            tempcor=tempcor.dropna()
            coef,pv=pearsonr(tempcor['x'],tempcor['y'])
            co.append(coef)
            pvv.append(pv)
        result_1[j]=co
        result_p[j]=pvv
            #print(coef)
    result_1=abs(result_1)
    result_1.index=set_index
    
    plt.rc('font', family='Arial')
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.weight"] = "12"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.figure(figsize=(10,10))
    sns.heatmap(result_1,vmin=0, vmax=1,center=1,annot=True,square=True)
    plt.savefig('co_character.png',dpi=300)
    return result_1```

```pythonfont1 = {'family': 'Arial','weight' : 'bold','size'   : 15,}
data=data1
fold_change=-1
fold_threshold=0
pvalue_threshold=0.05
cmap="seismic"```

```python#pvalue
pvalue = []
for i in range(0, len(data)):
    ttest = stats.ttest_ind(list(data.iloc[i][eg].values), list(data.iloc[i][cg].values))
    pvalue.append(ttest[1])```

```pythonfrom statsmodels.stats.multitest import fdrcorrection```

```pythoncpvalue=fdrcorrection(np.array(pvalue), alpha=0.05, method='indep', is_sorted=False)```

```pythonpvalue[:20]```
*Output:*
```[0.6912767703980381,
 0.004595669513467789,
 0.15380129298258083,
 0.6339270973230824,
 0.005814769475742696,
 0.6554359366255393,
 0.44068138628043774,
 0.8323419402230038,
 0.21094281363692927,
 0.006812557462208787,
 0.5178967314347178,
 0.06601828313091539,
 0.15062003196331625,
 0.07423755381159991,
 0.17730428367014733,
 0.8754650581038927,
 0.0076893086202108225,
 0.2436311511624047,
 0.03640678474151587,
 0.18828740113876674]```

```pythoncpvalue[1][0:20]```
*Output:*
```array([0.84176081, 0.08908264, 0.40124096, 0.80914046, 0.09811343,
       0.82057703, 0.67806369, 0.91991951, 0.46687315, 0.10406133,
       0.73547297, 0.2642261 , 0.39666932, 0.27971328, 0.42927463,
       0.94277245, 0.10749819, 0.50278639, 0.20278749, 0.44242986])```

```python

#cal_mean
eg_mean=data[eg].mean(axis=1)
eg_mean.head()
cg_mean=data[cg].mean(axis=1)
cg_mean.head()

#cal_fold
my_dpi=300
fig=plt.figure(figsize=(2000/my_dpi, 1000/my_dpi), dpi=my_dpi)
grid = plt.GridSpec(1, 4, wspace=1, hspace=0.1)
#fig, (ax0, ax1) = plt.subplots(2, 1)
plt.subplot(grid[0,0:2])
fold=eg_mean/cg_mean
log2fold=-np.log2(fold)
foldp=plt.hist(log2fold,color="#384793")
if fold_change==-1:
    foldchange=(foldp[1][np.where(foldp[1]>0)[0][fold_threshold]]+foldp[1][np.where(foldp[1]>0)[0][fold_threshold+1]])/2
else:
    foldchange=fold_change
cmap1=sns.color_palette(cmap)
plt.title("log2fc",font1)
plt.ylabel('Density',font1)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig("fold_change.png",dpi=300,bbox_inches = 'tight')

#pvalue
pvalue = []
for i in range(0, len(data)):
    ttest = stats.ttest_ind(list(data.iloc[i][eg].values), list(data.iloc[i,2:4].values))
    pvalue.append(ttest[1])```

```pythondef find_DEG(data,eg,cg,
             fold_change=-1,
             fold_threshold=0,
             pvalue_threshold=0.05,
             cmap="seismic"
            ):
    #plt_set
    font1 = {
    'family': 'Arial',
    'weight' : 'bold',
    'size'   : 15,
    }
    
    
    #cal_mean
    eg_mean=data[eg].mean(axis=1)
    eg_mean.head()
    cg_mean=data[cg].mean(axis=1)
    cg_mean.head()
    
    #cal_fold
    my_dpi=300
    fig=plt.figure(figsize=(2000/my_dpi, 1000/my_dpi), dpi=my_dpi)
    grid = plt.GridSpec(1, 4, wspace=1, hspace=0.1)
    #fig, (ax0, ax1) = plt.subplots(2, 1)
    plt.subplot(grid[0,0:2])
    fold=eg_mean/cg_mean
    log2fold=-np.log2(fold)
    foldp=plt.hist(log2fold,color="#384793")
    if fold_change==-1:
        foldchange=(foldp[1][np.where(foldp[1]>0)[0][fold_threshold]]+foldp[1][np.where(foldp[1]>0)[0][fold_threshold+1]])/2
    else:
        foldchange=fold_change
    cmap1=sns.color_palette(cmap)
    plt.title("log2fc",font1)
    plt.ylabel('Density',font1)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.savefig("fold_change.png",dpi=300,bbox_inches = 'tight')
    
    #pvalue
    pvalue = []
    for i in range(0, len(data)):
        ttest = stats.ttest_ind(list(data.iloc[i][eg].values), list(data.iloc[i,2:4].values))
        pvalue.append(ttest[1])
    
    #result
    genearray = np.asarray(pvalue)
    result = pd.DataFrame({'pvalue':genearray,'FoldChange':fold})
    result['log(pvalue)'] = -np.log10(result['pvalue'])
    result['log2FC'] = -np.log2(result['FoldChange'])
    result['sig'] = 'normal'
    result['size']  =np.abs(result['FoldChange'])/10
    result.loc[(result.log2FC> foldchange )&(result.pvalue < pvalue_threshold),'sig'] = 'up'
    result.loc[(result.log2FC< 0-foldchange )&(result.pvalue < pvalue_threshold),'sig'] = 'down'
    result.to_csv("DEG_result.csv")
    print('up:',len(result[result['sig']=='up']))
    print('down:',len(result[result['sig']=='down']))
    
    #plt
    plt.subplot(grid[0,2:4])
    ax = sns.scatterplot(x="log2FC", y="log(pvalue)",
                      hue='sig',
                      hue_order = ('up','down','normal'),
                      palette=(cmap1[5],cmap1[0],"grey"),
                      size='sig',sizes=(50, 100),
                      data=result)
    ax.set_ylabel('-log(pvalue)',font1)                                    
    ax.set_xlabel('log2FC',font1)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.legend(bbox_to_anchor=(1.05,0), loc=3, borderaxespad=0)
    ax.tick_params(labelsize=15)
    plt.savefig("DER_fire.png",dpi=300,bbox_inches = 'tight')
    
    #snsclu
    fold_cutoff = foldchange
    pvalue_cutoff = 0.05
    filtered_ids = []
    for i in range(0, len(result)):
        if (abs(-np.log2(fold[i])) >= fold_cutoff) and (pvalue[i] <= pvalue_cutoff):
            filtered_ids.append(i)        
    filtered = data.iloc[filtered_ids,:]
    filtered.to_csv('fi.csv')
    a=sns.clustermap(filtered, cmap=cmap, standard_scale = 0)
    plt.savefig("sns2.png",dpi=300,bbox_inches = 'tight')
    return result```

```pythonimport gseapy as gp
from gseapy.plot import barplot, dotplot
def enrichment_KEGG(gene_list,
                    gene_sets=['KEGG_2019_Human'],
                    organism='Human',
                    description='test_name',
                    outdir='enrichment_kegg',
                    cutoff=0.5):
    enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=gene_sets,
                 organism=organism, # don't forget to set organism to the one you desired! e.g. Yeast
                 description=description,
                 outdir=outdir,
                 # no_plot=True,
                 cutoff=cutoff # test dataset, use lower value from range(0,1)
                )
    subp=dotplot(enr.res2d, title=description,cmap='viridis_r',ofname='kegg.png')
    print(subp)
    return enr.res2d```

```pythonhelp(gp.prerank)```
*Output:*
```Help on function prerank in module gseapy.gsea:

prerank(rnk, gene_sets, outdir='GSEA_Prerank', pheno_pos='Pos', pheno_neg='Neg', min_size=15, max_size=500, permutation_num=1000, weighted_score_type=1, ascending=False, processes=1, figsize=(6.5, 6), format='pdf', graph_num=20, no_plot=False, seed=None, verbose=False)
    Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.
    
    :param rnk: pre-ranked correlation table or pandas DataFrame. Same input with ``GSEA`` .rnk file.
    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
    :param outdir: results output directory.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param str weighted_score_type: Refer to :func:`algorithm.enrichment_score`. Default:1.
    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype.
    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.
    :param seed: Random seed. expect an integer. Default:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.
    
    :return: Return a Prerank obj. All results store to  a dictionary, obj.results,
             where contains::
    
                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set
                 |  ledge_genes: leading edge genes}

```

```pythondef enrichment_GO(gene_list,
                    go_mode='Bio',
                    organism='Human',
                    description='test_name',
                    outdir='enrichment_go',
                    cutoff=0.5):
    if(go_mode=='Bio'):
        geneset='GO_Biological_Process_2018'
    if(go_mode=='Cell'):
        geneset='GO_Cellular_Component_2018'
    if(go_mode=='Mole'):
        geneset='GO_Molecular_Function_2018'
    enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=geneset,
                 organism=organism, # don't forget to set organism to the one you desired! e.g. Yeast
                 description=description,
                 outdir=outdir,
                 # no_plot=True,
                 cutoff=cutoff # test dataset, use lower value from range(0,1)
                )
    subp=dotplot(enr.res2d, title=description,cmap='viridis_r',ofname='go_'+go_mode+'.png')
    print(subp)
    return enr.res2d```

```pythondef enrichment_GSEA(data,
                   gene_sets='KEGG_2016',
                   processes=4,
                   permutation_num=100,
                   outdir='prerank_report_kegg',
                   seed=6):
    rnk=pd.DataFrame(columns=['genename','FoldChange'])
    rnk['genename']=data.index
    rnk['FoldChange']=data['FoldChange'].tolist()
    rnk1=rnk.drop_duplicates(['genename'])
    rnk1=rnk1.sort_values(by='FoldChange', ascending=False)
    
    pre_res = gp.prerank(rnk=rnk1, gene_sets=gene_sets,
                     processes=processes,
                     permutation_num=permutation_num, # reduce number to speed up testing
                     outdir=outdir, format='png', seed=seed)
    pre_res.res2d.sort_index().to_csv('GSEA_result.csv')
    return pre_res```

```pythonresult[result['sig']=='up']```
*Output:*
```            pvalue  FoldChange  log(pvalue)    log2FC sig      size
SKAP2     0.026041    0.927901     1.584340  0.107957  up  0.092790
GABBR2    0.037239    0.921045     1.428999  0.118657  up  0.092104
SYT1      0.021706    0.925714     1.663422  0.111361  up  0.092571
CAMKK1    0.030513    0.920596     1.515512  0.119359  up  0.092060
MIF       0.040802    0.918903     1.389314  0.122016  up  0.091890
SNRPA     0.017790    0.911486     1.749830  0.133707  up  0.091149
NMT2      0.034770    0.911214     1.458801  0.134139  up  0.091121
RPS9      0.048605    0.912683     1.313321  0.131815  up  0.091268
PSME1     0.044375    0.904953     1.352863  0.144086  up  0.090495
ACAN      0.021108    0.890520     1.675544  0.167280  up  0.089052
ADAP1     0.012932    0.886673     1.888339  0.173526  up  0.088667
STMN1     0.039626    0.883142     1.402016  0.179283  up  0.088314
SCRN1     0.018609    0.885723     1.730278  0.175073  up  0.088572
NEDD8     0.038544    0.876894     1.414043  0.189526  up  0.087689
HOMER1    0.020396    0.871788     1.690450  0.197951  up  0.087179
BCAS3     0.043655    0.851541     1.359969  0.231852  up  0.085154
DNAJC11   0.040340    0.854242     1.394264  0.227284  up  0.085424
STX12     0.005186    0.841412     2.285185  0.249115  up  0.084141
CCDC149   0.013941    0.840104     1.855702  0.251360  up  0.084010
BTBD8     0.049152    0.816823     1.308462  0.291905  up  0.081682
HINT3     0.027380    0.827259     1.562573  0.273589  up  0.082726
CDIP1     0.048586    0.821218     1.313490  0.284162  up  0.082122
SLC38A10  0.005449    0.795782     2.263672  0.329555  up  0.079578
YJEFN3    0.024005    0.743570     1.619697  0.427460  up  0.074357```

```pythonresult=find_DEG(normdata,eg=['lab2','lab3'],cg=['non2','non3'],fold_change=0.1,cmap="seismic")
res=enrichment_GSEA(result)
GSEA_plot(res,0)```
*Output:*
```up: 24
down: 93
```
```<Figure size 2000x1000 with 2 Axes>```
```<Figure size 720x720 with 4 Axes>```
```<Figure size 432x396 with 4 Axes>```

```pythondef Plot_GSEA(data,num=0):
    terms = data.res2d.index
    from gseapy.plot import gseaplot
    # to save your figure, make sure that ofname is not None
    gseaplot(rank_metric=data.ranking, term=terms[num], **data.results[terms[num]])```

```pythonhelp(data.plot)```
*Output:*
```Help on FramePlotMethods in module pandas.plotting._core object:

class FramePlotMethods(BasePlotMethods)
 |  FramePlotMethods(data)
 |  
 |  DataFrame plotting accessor and method
 |  
 |  Examples
 |  --------
 |  >>> df.plot.line()
 |  >>> df.plot.scatter('x', 'y')
 |  >>> df.plot.hexbin()
 |  
 |  These plotting methods can also be accessed by calling the accessor as a
 |  method with the ``kind`` argument:
 |  ``df.plot(kind='line')`` is equivalent to ``df.plot.line()``
 |  
 |  Method resolution order:
 |      FramePlotMethods
 |      BasePlotMethods
 |      pandas.core.base.PandasObject
 |      pandas.core.base.StringMixin
 |      pandas.core.accessor.DirNamesMixin
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __call__(self, x=None, y=None, kind='line', ax=None, subplots=False, sharex=None, sharey=False, layout=None, figsize=None, use_index=True, title=None, grid=None, legend=True, style=None, logx=False, logy=False, loglog=False, xticks=None, yticks=None, xlim=None, ylim=None, rot=None, fontsize=None, colormap=None, table=False, yerr=None, xerr=None, secondary_y=False, sort_columns=False, **kwds)
 |      Make plots of DataFrame using matplotlib / pylab.
 |      
 |      *New in version 0.17.0:* Each plot kind has a corresponding method on the
 |      ``DataFrame.plot`` accessor:
 |      ``df.plot(kind='line')`` is equivalent to
 |      ``df.plot.line()``.
 |      
 |      Parameters
 |      ----------
 |      data : DataFrame
 |      x : label or position, default None
 |      y : label, position or list of label, positions, default None
 |          Allows plotting of one column versus another
 |      kind : str
 |          - 'line' : line plot (default)
 |          - 'bar' : vertical bar plot
 |          - 'barh' : horizontal bar plot
 |          - 'hist' : histogram
 |          - 'box' : boxplot
 |          - 'kde' : Kernel Density Estimation plot
 |          - 'density' : same as 'kde'
 |          - 'area' : area plot
 |          - 'pie' : pie plot
 |          - 'scatter' : scatter plot
 |          - 'hexbin' : hexbin plot
 |      ax : matplotlib axes object, default None
 |      subplots : boolean, default False
 |          Make separate subplots for each column
 |      sharex : boolean, default True if ax is None else False
 |          In case subplots=True, share x axis and set some x axis labels to
 |          invisible; defaults to True if ax is None otherwise False if an ax
 |          is passed in; Be aware, that passing in both an ax and sharex=True
 |          will alter all x axis labels for all axis in a figure!
 |      sharey : boolean, default False
 |          In case subplots=True, share y axis and set some y axis labels to
 |          invisible
 |      layout : tuple (optional)
 |          (rows, columns) for the layout of subplots
 |      figsize : a tuple (width, height) in inches
 |      use_index : boolean, default True
 |          Use index as ticks for x axis
 |      title : string or list
 |          Title to use for the plot. If a string is passed, print the string at
 |          the top of the figure. If a list is passed and `subplots` is True,
 |          print each item in the list above the corresponding subplot.
 |      grid : boolean, default None (matlab style default)
 |          Axis grid lines
 |      legend : False/True/'reverse'
 |          Place legend on axis subplots
 |      style : list or dict
 |          matplotlib line style per column
 |      logx : boolean, default False
 |          Use log scaling on x axis
 |      logy : boolean, default False
 |          Use log scaling on y axis
 |      loglog : boolean, default False
 |          Use log scaling on both x and y axes
 |      xticks : sequence
 |          Values to use for the xticks
 |      yticks : sequence
 |          Values to use for the yticks
 |      xlim : 2-tuple/list
 |      ylim : 2-tuple/list
 |      rot : int, default None
 |          Rotation for ticks (xticks for vertical, yticks for horizontal plots)
 |      fontsize : int, default None
 |          Font size for xticks and yticks
 |      colormap : str or matplotlib colormap object, default None
 |          Colormap to select colors from. If string, load colormap with that name
 |          from matplotlib.
 |      colorbar : boolean, optional
 |          If True, plot colorbar (only relevant for 'scatter' and 'hexbin' plots)
 |      position : float
 |          Specify relative alignments for bar plot layout.
 |          From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
 |      table : boolean, Series or DataFrame, default False
 |          If True, draw a table using the data in the DataFrame and the data will
 |          be transposed to meet matplotlib's default layout.
 |          If a Series or DataFrame is passed, use passed data to draw a table.
 |      yerr : DataFrame, Series, array-like, dict and str
 |          See :ref:`Plotting with Error Bars <visualization.errorbars>` for
 |          detail.
 |      xerr : same types as yerr.
 |      stacked : boolean, default False in line and
 |          bar plots, and True in area plot. If True, create stacked plot.
 |      sort_columns : boolean, default False
 |          Sort column names to determine plot ordering
 |      secondary_y : boolean or sequence, default False
 |          Whether to plot on the secondary y-axis
 |          If a list/tuple, which columns to plot on secondary y-axis
 |      mark_right : boolean, default True
 |          When using a secondary_y axis, automatically mark the column
 |          labels with "(right)" in the legend
 |      `**kwds` : keywords
 |          Options to pass to matplotlib plotting method
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or numpy.ndarray of them
 |      
 |      Notes
 |      -----
 |      
 |      - See matplotlib documentation online for more on this subject
 |      - If `kind` = 'bar' or 'barh', you can specify relative alignments
 |        for bar plot layout by `position` keyword.
 |        From 0 (left/bottom-end) to 1 (right/top-end). Default is 0.5 (center)
 |      - If `kind` = 'scatter' and the argument `c` is the name of a dataframe
 |        column, the values of that column are used to color each point.
 |      - If `kind` = 'hexbin', you can control the size of the bins with the
 |        `gridsize` argument. By default, a histogram of the counts around each
 |        `(x, y)` point is computed. You can specify alternative aggregations
 |        by passing values to the `C` and `reduce_C_function` arguments.
 |        `C` specifies the value at each `(x, y)` point and `reduce_C_function`
 |        is a function of one argument that reduces all the values in a bin to
 |        a single number (e.g. `mean`, `max`, `sum`, `std`).
 |  
 |  area(self, x=None, y=None, **kwds)
 |      Area plot
 |      
 |      Parameters
 |      ----------
 |      x, y : label or position, optional
 |          Coordinates for each point.
 |      `**kwds` : optional
 |          Additional keyword arguments are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or numpy.ndarray of them
 |  
 |  bar(self, x=None, y=None, **kwds)
 |      Vertical bar plot.
 |      
 |      A bar plot is a plot that presents categorical data with
 |      rectangular bars with lengths proportional to the values that they
 |      represent. A bar plot shows comparisons among discrete categories. One
 |      axis of the plot shows the specific categories being compared, and the
 |      other axis represents a measured value.
 |      
 |      Parameters
 |      ----------
 |      x : label or position, optional
 |          Allows plotting of one column versus another. If not specified,
 |          the index of the DataFrame is used.
 |      y : label or position, optional
 |          Allows plotting of one column versus another. If not specified,
 |          all numerical columns are used.
 |      **kwds
 |          Additional keyword arguments are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : matplotlib.axes.Axes or np.ndarray of them
 |          An ndarray is returned with one :class:`matplotlib.axes.Axes`
 |          per column when ``subplots=True``.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.plot.barh : Horizontal bar plot.
 |      pandas.DataFrame.plot : Make plots of a DataFrame.
 |      matplotlib.pyplot.bar : Make a bar plot with matplotlib.
 |      
 |      Examples
 |      --------
 |      Basic plot.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame({'lab':['A', 'B', 'C'], 'val':[10, 30, 20]})
 |          >>> ax = df.plot.bar(x='lab', y='val', rot=0)
 |      
 |      Plot a whole dataframe to a bar plot. Each column is assigned a
 |      distinct color, and each row is nested in a group along the
 |      horizontal axis.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
 |          >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
 |          >>> index = ['snail', 'pig', 'elephant',
 |          ...          'rabbit', 'giraffe', 'coyote', 'horse']
 |          >>> df = pd.DataFrame({'speed': speed,
 |          ...                    'lifespan': lifespan}, index=index)
 |          >>> ax = df.plot.bar(rot=0)
 |      
 |      Instead of nesting, the figure can be split by column with
 |      ``subplots=True``. In this case, a :class:`numpy.ndarray` of
 |      :class:`matplotlib.axes.Axes` are returned.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> axes = df.plot.bar(rot=0, subplots=True)
 |          >>> axes[1].legend(loc=2)  # doctest: +SKIP
 |      
 |      Plot a single column.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax = df.plot.bar(y='speed', rot=0)
 |      
 |      Plot only selected categories for the DataFrame.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax = df.plot.bar(x='lifespan', rot=0)
 |  
 |  barh(self, x=None, y=None, **kwds)
 |      Make a horizontal bar plot.
 |      
 |      A horizontal bar plot is a plot that presents quantitative data with
 |      rectangular bars with lengths proportional to the values that they
 |      represent. A bar plot shows comparisons among discrete categories. One
 |      axis of the plot shows the specific categories being compared, and the
 |      other axis represents a measured value.
 |      
 |      Parameters
 |      ----------
 |      x : label or position, default DataFrame.index
 |          Column to be used for categories.
 |      y : label or position, default All numeric columns in dataframe
 |          Columns to be plotted from the DataFrame.
 |      **kwds
 |          Keyword arguments to pass on to :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or numpy.ndarray of them.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.plot.bar: Vertical bar plot.
 |      pandas.DataFrame.plot : Make plots of DataFrame using matplotlib.
 |      matplotlib.axes.Axes.bar : Plot a vertical bar plot using matplotlib.
 |      
 |      Examples
 |      --------
 |      Basic example
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame({'lab':['A', 'B', 'C'], 'val':[10, 30, 20]})
 |          >>> ax = df.plot.barh(x='lab', y='val')
 |      
 |      Plot a whole DataFrame to a horizontal bar plot
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
 |          >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
 |          >>> index = ['snail', 'pig', 'elephant',
 |          ...          'rabbit', 'giraffe', 'coyote', 'horse']
 |          >>> df = pd.DataFrame({'speed': speed,
 |          ...                    'lifespan': lifespan}, index=index)
 |          >>> ax = df.plot.barh()
 |      
 |      Plot a column of the DataFrame to a horizontal bar plot
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
 |          >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
 |          >>> index = ['snail', 'pig', 'elephant',
 |          ...          'rabbit', 'giraffe', 'coyote', 'horse']
 |          >>> df = pd.DataFrame({'speed': speed,
 |          ...                    'lifespan': lifespan}, index=index)
 |          >>> ax = df.plot.barh(y='speed')
 |      
 |      Plot DataFrame versus the desired column
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> speed = [0.1, 17.5, 40, 48, 52, 69, 88]
 |          >>> lifespan = [2, 8, 70, 1.5, 25, 12, 28]
 |          >>> index = ['snail', 'pig', 'elephant',
 |          ...          'rabbit', 'giraffe', 'coyote', 'horse']
 |          >>> df = pd.DataFrame({'speed': speed,
 |          ...                    'lifespan': lifespan}, index=index)
 |          >>> ax = df.plot.barh(x='lifespan')
 |  
 |  box(self, by=None, **kwds)
 |      Make a box plot of the DataFrame columns.
 |      
 |      A box plot is a method for graphically depicting groups of numerical
 |      data through their quartiles.
 |      The box extends from the Q1 to Q3 quartile values of the data,
 |      with a line at the median (Q2). The whiskers extend from the edges
 |      of box to show the range of the data. The position of the whiskers
 |      is set by default to 1.5*IQR (IQR = Q3 - Q1) from the edges of the
 |      box. Outlier points are those past the end of the whiskers.
 |      
 |      For further details see Wikipedia's
 |      entry for `boxplot <https://en.wikipedia.org/wiki/Box_plot>`__.
 |      
 |      A consideration when using this chart is that the box and the whiskers
 |      can overlap, which is very common when plotting small sets of data.
 |      
 |      Parameters
 |      ----------
 |      by : string or sequence
 |          Column in the DataFrame to group by.
 |      **kwds : optional
 |          Additional keywords are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or numpy.ndarray of them
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.boxplot: Another method to draw a box plot.
 |      pandas.Series.plot.box: Draw a box plot from a Series object.
 |      matplotlib.pyplot.boxplot: Draw a box plot in matplotlib.
 |      
 |      Examples
 |      --------
 |      Draw a box plot from a DataFrame with four columns of randomly
 |      generated data.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> data = np.random.randn(25, 4)
 |          >>> df = pd.DataFrame(data, columns=list('ABCD'))
 |          >>> ax = df.plot.box()
 |  
 |  density = kde(self, bw_method=None, ind=None, **kwds)
 |  
 |  hexbin(self, x, y, C=None, reduce_C_function=None, gridsize=None, **kwds)
 |      Generate a hexagonal binning plot.
 |      
 |      Generate a hexagonal binning plot of `x` versus `y`. If `C` is `None`
 |      (the default), this is a histogram of the number of occurrences
 |      of the observations at ``(x[i], y[i])``.
 |      
 |      If `C` is specified, specifies values at given coordinates
 |      ``(x[i], y[i])``. These values are accumulated for each hexagonal
 |      bin and then reduced according to `reduce_C_function`,
 |      having as default the NumPy's mean function (:meth:`numpy.mean`).
 |      (If `C` is specified, it must also be a 1-D sequence
 |      of the same length as `x` and `y`, or a column label.)
 |      
 |      Parameters
 |      ----------
 |      x : int or str
 |          The column label or position for x points.
 |      y : int or str
 |          The column label or position for y points.
 |      C : int or str, optional
 |          The column label or position for the value of `(x, y)` point.
 |      reduce_C_function : callable, default `np.mean`
 |          Function of one argument that reduces all the values in a bin to
 |          a single number (e.g. `np.mean`, `np.max`, `np.sum`, `np.std`).
 |      gridsize : int or tuple of (int, int), default 100
 |          The number of hexagons in the x-direction.
 |          The corresponding number of hexagons in the y-direction is
 |          chosen in a way that the hexagons are approximately regular.
 |          Alternatively, gridsize can be a tuple with two elements
 |          specifying the number of hexagons in the x-direction and the
 |          y-direction.
 |      **kwds
 |          Additional keyword arguments are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      matplotlib.AxesSubplot
 |          The matplotlib ``Axes`` on which the hexbin is plotted.
 |      
 |      See Also
 |      --------
 |      DataFrame.plot : Make plots of a DataFrame.
 |      matplotlib.pyplot.hexbin : hexagonal binning plot using matplotlib,
 |          the matplotlib function that is used under the hood.
 |      
 |      Examples
 |      --------
 |      The following examples are generated with random data from
 |      a normal distribution.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> n = 10000
 |          >>> df = pd.DataFrame({'x': np.random.randn(n),
 |          ...                    'y': np.random.randn(n)})
 |          >>> ax = df.plot.hexbin(x='x', y='y', gridsize=20)
 |      
 |      The next example uses `C` and `np.sum` as `reduce_C_function`.
 |      Note that `'observations'` values ranges from 1 to 5 but the result
 |      plot shows values up to more than 25. This is because of the
 |      `reduce_C_function`.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> n = 500
 |          >>> df = pd.DataFrame({
 |          ...     'coord_x': np.random.uniform(-3, 3, size=n),
 |          ...     'coord_y': np.random.uniform(30, 50, size=n),
 |          ...     'observations': np.random.randint(1,5, size=n)
 |          ...     })
 |          >>> ax = df.plot.hexbin(x='coord_x',
 |          ...                     y='coord_y',
 |          ...                     C='observations',
 |          ...                     reduce_C_function=np.sum,
 |          ...                     gridsize=10,
 |          ...                     cmap="viridis")
 |  
 |  hist(self, by=None, bins=10, **kwds)
 |      Draw one histogram of the DataFrame's columns.
 |      
 |      A histogram is a representation of the distribution of data.
 |      This function groups the values of all given Series in the DataFrame
 |      into bins and draws all bins in one :class:`matplotlib.axes.Axes`.
 |      This is useful when the DataFrame's Series are in a similar scale.
 |      
 |      Parameters
 |      ----------
 |      by : str or sequence, optional
 |          Column in the DataFrame to group by.
 |      bins : int, default 10
 |          Number of histogram bins to be used.
 |      **kwds
 |          Additional keyword arguments are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : matplotlib.AxesSubplot histogram.
 |      
 |      See Also
 |      --------
 |      DataFrame.hist : Draw histograms per DataFrame's Series.
 |      Series.hist : Draw a histogram with Series' data.
 |      
 |      Examples
 |      --------
 |      When we draw a dice 6000 times, we expect to get each value around 1000
 |      times. But when we draw two dices and sum the result, the distribution
 |      is going to be quite different. A histogram illustrates those
 |      distributions.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame(
 |          ...     np.random.randint(1, 7, 6000),
 |          ...     columns = ['one'])
 |          >>> df['two'] = df['one'] + np.random.randint(1, 7, 6000)
 |          >>> ax = df.plot.hist(bins=12, alpha=0.5)
 |  
 |  kde(self, bw_method=None, ind=None, **kwds)
 |      Generate Kernel Density Estimate plot using Gaussian kernels.
 |      
 |      In statistics, `kernel density estimation`_ (KDE) is a non-parametric
 |      way to estimate the probability density function (PDF) of a random
 |      variable. This function uses Gaussian kernels and includes automatic
 |      bandwith determination.
 |      
 |      .. _kernel density estimation:
 |          https://en.wikipedia.org/wiki/Kernel_density_estimation
 |      
 |      Parameters
 |      ----------
 |      bw_method : str, scalar or callable, optional
 |          The method used to calculate the estimator bandwidth. This can be
 |          'scott', 'silverman', a scalar constant or a callable.
 |          If None (default), 'scott' is used.
 |          See :class:`scipy.stats.gaussian_kde` for more information.
 |      ind : NumPy array or integer, optional
 |          Evaluation points for the estimated PDF. If None (default),
 |          1000 equally spaced points are used. If `ind` is a NumPy array, the
 |          KDE is evaluated at the points passed. If `ind` is an integer,
 |          `ind` number of equally spaced points are used.
 |      **kwds : optional
 |          Additional keyword arguments are documented in
 |          :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : matplotlib.axes.Axes or numpy.ndarray of them
 |      
 |      See Also
 |      --------
 |      scipy.stats.gaussian_kde : Representation of a kernel-density
 |          estimate using Gaussian kernels. This is the function used
 |          internally to estimate the PDF.
 |      Series.plot.kde : Generate a KDE plot for a
 |          Series.
 |      
 |      Examples
 |      --------
 |      Given several Series of points randomly sampled from unknown
 |      distributions, estimate their PDFs using KDE with automatic
 |      bandwidth determination and plot the results, evaluating them at
 |      1000 equally spaced points (default):
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame({
 |          ...     'x': [1, 2, 2.5, 3, 3.5, 4, 5],
 |          ...     'y': [4, 4, 4.5, 5, 5.5, 6, 6],
 |          ... })
 |          >>> ax = df.plot.kde()
 |      
 |      A scalar bandwidth can be specified. Using a small bandwidth value can
 |      lead to overfitting, while using a large bandwidth value may result
 |      in underfitting:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax = df.plot.kde(bw_method=0.3)
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax = df.plot.kde(bw_method=3)
 |      
 |      Finally, the `ind` parameter determines the evaluation points for the
 |      plot of the estimated PDF:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax = df.plot.kde(ind=[1, 2, 3, 4, 5, 6])
 |  
 |  line(self, x=None, y=None, **kwds)
 |      Plot DataFrame columns as lines.
 |      
 |      This function is useful to plot lines using DataFrame's values
 |      as coordinates.
 |      
 |      Parameters
 |      ----------
 |      x : int or str, optional
 |          Columns to use for the horizontal axis.
 |          Either the location or the label of the columns to be used.
 |          By default, it will use the DataFrame indices.
 |      y : int, str, or list of them, optional
 |          The values to be plotted.
 |          Either the location or the label of the columns to be used.
 |          By default, it will use the remaining DataFrame numeric columns.
 |      **kwds
 |          Keyword arguments to pass on to :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or :class:`numpy.ndarray`
 |          Returns an ndarray when ``subplots=True``.
 |      
 |      See Also
 |      --------
 |      matplotlib.pyplot.plot : Plot y versus x as lines and/or markers.
 |      
 |      Examples
 |      --------
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          The following example shows the populations for some animals
 |          over the years.
 |      
 |          >>> df = pd.DataFrame({
 |          ...    'pig': [20, 18, 489, 675, 1776],
 |          ...    'horse': [4, 25, 281, 600, 1900]
 |          ...    }, index=[1990, 1997, 2003, 2009, 2014])
 |          >>> lines = df.plot.line()
 |      
 |      .. plot::
 |         :context: close-figs
 |      
 |         An example with subplots, so an array of axes is returned.
 |      
 |         >>> axes = df.plot.line(subplots=True)
 |         >>> type(axes)
 |         <class 'numpy.ndarray'>
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          The following example shows the relationship between both
 |          populations.
 |      
 |          >>> lines = df.plot.line(x='pig', y='horse')
 |  
 |  pie(self, y=None, **kwds)
 |      Generate a pie plot.
 |      
 |      A pie plot is a proportional representation of the numerical data in a
 |      column. This function wraps :meth:`matplotlib.pyplot.pie` for the
 |      specified column. If no column reference is passed and
 |      ``subplots=True`` a pie plot is drawn for each numerical column
 |      independently.
 |      
 |      Parameters
 |      ----------
 |      y : int or label, optional
 |          Label or position of the column to plot.
 |          If not provided, ``subplots=True`` argument must be passed.
 |      **kwds
 |          Keyword arguments to pass on to :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : matplotlib.axes.Axes or np.ndarray of them.
 |          A NumPy array is returned when `subplots` is True.
 |      
 |      See Also
 |      --------
 |      Series.plot.pie : Generate a pie plot for a Series.
 |      DataFrame.plot : Make plots of a DataFrame.
 |      
 |      Examples
 |      --------
 |      In the example below we have a DataFrame with the information about
 |      planet's mass and radius. We pass the the 'mass' column to the
 |      pie function to get a pie plot.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame({'mass': [0.330, 4.87 , 5.97],
 |          ...                    'radius': [2439.7, 6051.8, 6378.1]},
 |          ...                   index=['Mercury', 'Venus', 'Earth'])
 |          >>> plot = df.plot.pie(y='mass', figsize=(5, 5))
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> plot = df.plot.pie(subplots=True, figsize=(6, 3))
 |  
 |  scatter(self, x, y, s=None, c=None, **kwds)
 |      Create a scatter plot with varying marker point size and color.
 |      
 |      The coordinates of each point are defined by two dataframe columns and
 |      filled circles are used to represent each point. This kind of plot is
 |      useful to see complex correlations between two variables. Points could
 |      be for instance natural 2D coordinates like longitude and latitude in
 |      a map or, in general, any pair of metrics that can be plotted against
 |      each other.
 |      
 |      Parameters
 |      ----------
 |      x : int or str
 |          The column name or column position to be used as horizontal
 |          coordinates for each point.
 |      y : int or str
 |          The column name or column position to be used as vertical
 |          coordinates for each point.
 |      s : scalar or array_like, optional
 |          The size of each point. Possible values are:
 |      
 |          - A single scalar so all points have the same size.
 |      
 |          - A sequence of scalars, which will be used for each point's size
 |            recursively. For instance, when passing [2,14] all points size
 |            will be either 2 or 14, alternatively.
 |      
 |      c : str, int or array_like, optional
 |          The color of each point. Possible values are:
 |      
 |          - A single color string referred to by name, RGB or RGBA code,
 |            for instance 'red' or '#a98d19'.
 |      
 |          - A sequence of color strings referred to by name, RGB or RGBA
 |            code, which will be used for each point's color recursively. For
 |            intance ['green','yellow'] all points will be filled in green or
 |            yellow, alternatively.
 |      
 |          - A column name or position whose values will be used to color the
 |            marker points according to a colormap.
 |      
 |      **kwds
 |          Keyword arguments to pass on to :meth:`pandas.DataFrame.plot`.
 |      
 |      Returns
 |      -------
 |      axes : :class:`matplotlib.axes.Axes` or numpy.ndarray of them
 |      
 |      See Also
 |      --------
 |      matplotlib.pyplot.scatter : scatter plot using multiple input data
 |          formats.
 |      
 |      Examples
 |      --------
 |      Let's see how to draw a scatter plot using coordinates from the values
 |      in a DataFrame's columns.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame([[5.1, 3.5, 0], [4.9, 3.0, 0], [7.0, 3.2, 1],
 |          ...                    [6.4, 3.2, 1], [5.9, 3.0, 2]],
 |          ...                   columns=['length', 'width', 'species'])
 |          >>> ax1 = df.plot.scatter(x='length',
 |          ...                       y='width',
 |          ...                       c='DarkBlue')
 |      
 |      And now with the color determined by a column as well.
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> ax2 = df.plot.scatter(x='length',
 |          ...                       y='width',
 |          ...                       c='species',
 |          ...                       colormap='viridis')
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from BasePlotMethods:
 |  
 |  __init__(self, data)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.base.PandasObject:
 |  
 |  __sizeof__(self)
 |      Generates the total memory usage for an object that returns
 |      either a value or Series of values
 |  
 |  __unicode__(self)
 |      Return a string representation for a particular object.
 |      
 |      Invoked by unicode(obj) in py2 only. Yields a Unicode String in both
 |      py2/py3.
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.base.StringMixin:
 |  
 |  __bytes__(self)
 |      Return a string representation for a particular object.
 |      
 |      Invoked by bytes(obj) in py3 only.
 |      Yields a bytestring in both py2/py3.
 |  
 |  __repr__(self)
 |      Return a string representation for a particular object.
 |      
 |      Yields Bytestring in Py2, Unicode String in py3.
 |  
 |  __str__(self)
 |      Return a string representation for a particular Object
 |      
 |      Invoked by str(df) in both py2/py3.
 |      Yields Bytestring in Py2, Unicode String in py3.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from pandas.core.base.StringMixin:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.accessor.DirNamesMixin:
 |  
 |  __dir__(self)
 |      Provide method name lookup and completion
 |      Only provide 'public' methods

```

```pythondef density_norm(data,depth=2,legend=False,norm_by=0,xlim=-1,ylim=-1):
    font1 = {
    'family': 'Arial',
    'weight' : 'bold',
    'size'   : 15,
    }
    grid = plt.GridSpec(1, 4, wspace=2, hspace=0.1)
    ax1=plt.subplot(grid[0,0:2])
    data.plot(kind = 'density',legend=legend,fontsize=15,colormap='seismic',ax=ax1)
    plt.ylabel('Density',font1)
    plt.title('Raw Count',font1)
    if (xlim!=-1):
        plt.xlim(0-xlim,xlim)
    if (ylim!=-1):
        plt.ylim(0,ylim)
    ERlist=ERgene.FindERG(data,depth)
    data2=ERgene.normalizationdata(data,ERlist[norm_by])
    ax2=plt.subplot(grid[0,2:4])
    data2.plot(kind='density',legend=legend,fontsize=15,colormap='seismic',ax=ax2)
    plt.ylabel('',font1)
    plt.title('After Count',font1)
    if (xlim!=-1):
        plt.xlim(0-xlim,xlim)
    if (ylim!=-1):
        plt.ylim(0,ylim)
    return data2```

```pythonhelp(data.plot)```

```pythondata2=ERgene.normalizationdata(data,'MMT00040568')
data.plot(kind = 'density',legend=False)
plt.xlim(-1,1)
data2.plot(kind='density',legend=False)
plt.xlim(-1,1)```
*Output:*

```pythondata.dropna().head()```
*Output:*
```                    lab2        lab3       non2       non3
genename                                                  
Rbbp7           3.199028    0.000000   0.000000   1.210764
Pafah1b2       61.878472   14.433681  11.759722  20.620139
Col1a2         39.829861  100.621528  14.347917  46.003472
Pimreg         12.568750   17.022917   5.023264   8.711806
1700022I11Rik   4.738889    0.000000   0.000000   2.273854```

```pythonmapping.head()```
*Output:*
```   #queryIndex queryItem              stringId preferredName  \
0            1    Col1a2  9606.ENSP00000297268        COL1A2   
1            2    Pimreg  9606.ENSP00000250056        FAM64A   
2            3  Pafah1b2  9606.ENSP00000435289      PAFAH1B2   
3            4      Ece2  9606.ENSP00000384223          ECE2   
4            5       Ttr  9606.ENSP00000237014           TTR   

                                          annotation  
0  Collagen alpha-2(I) chain; Type I collagen is ...  
1  Protein PIMREG; During mitosis, may play a rol...  
2  Platelet-activating factor acetylhydrolase IB ...  
3  Endothelin-converting enzyme 2; Converts big e...  
4  Transthyretin; Thyroid hormone-binding protein...  ```

```pythondata1=ID_mapping(data.dropna(),mapping,'queryItem','preferredName')
data1.head()```
*Output:*
```                lab2        lab3        non2        non3
PAFAH1B2   61.878472   14.433681   11.759722   20.620139
COL1A2     39.829861  100.621528   14.347917   46.003472
FAM64A     12.568750   17.022917    5.023264    8.711806
ECE2       15.166667   15.785764    5.435069   12.391319
NDUFS8    366.145833  204.065972  250.118056  357.500000```

```pythondata=pd.read_csv('test.csv')
data=data.set_index(data.columns[0])

mapping=pd.read_csv('C:/Users/FernandoZeng/Downloads/string_mapping.tsv',sep='\t')
#queryItem=raw_index
#preferredName=matching name
data1=ID_mapping(data.dropna(),mapping,'queryItem','preferredName')
#Eliminate the sample batch effect
normdata=density_norm(data1,
                      depth=2,
                      legend=False,
                      norm_by=0,
                      xlim=1000,
                     )```
*Output:*
```calculate time:0.41s
['ATP1A3', 'SPTAN1', 'SPTBN1', 'MAP1A', 'CKB', 'DYNC1H1', 'INA', 'ANK2', 'CLTC', 'HSP90AB1', 'CAMK2A', 'ATP1A1', 'STXBP1', 'NEFM', 'DNM1', 'SLC1A2', 'CAMK2B', 'MBP', 'ATP1A2', 'GNAO1']
```
```<Figure size 432x288 with 2 Axes>```

```pythondef ID_mapping(raw_data,mapping_data,raw_col,map_col):
    raw_index=raw_data.index
    #raw_col='queryItem'
    raw_length=raw_data.columns.size
    #map_col='preferredName'
    map_index=[]
    map_data=pd.DataFrame(columns=raw_data.columns)
    mapping=mapping_data
    for i in sorted(set(list(raw_index)),key=list(raw_index).index):
        if(mapping[mapping[raw_col]==i][map_col].values.size==0):
            continue
        else:

            if(raw_data.loc[i].size==raw_length):
                map_index.append(mapping[mapping[raw_col]==i][map_col].values[0])
                map_data=map_data.append(raw_data.loc[i],ignore_index=True)
            else:
                map_index.append(mapping[mapping[raw_col]==i][map_col].values[0])
                testdata=raw_data.loc[i]
                map_data=map_data.append(testdata.iloc[np.where(testdata.mean(axis=1)==testdata.mean(axis=1).max())])
            #print(mapping[mapping[raw_col]==i][map_col].values[0],raw_data.loc[i])
    map_data.index=map_index
    return map_data```

```pythondef gene_expression_plot(data,
                         gene_list,
                         eg,
                         cg,
                         save_name="gene_expression_plot.png",
                        cmap='seismic'):
    ui=pd.DataFrame(columns=['gene','value','class'])
    gene=gene_list
    for ge in gene:
        for i in data.loc[ge][eg].values:
            ui=ui.append({'gene':ge,'value':i,'class':'experiment'},ignore_index=True)
        for i in data.loc[ge][cg].values:
            ui=ui.append({'gene':ge,'value':i,'class':'ctrl'},ignore_index=True)
    df=ui
    dyt=[]
    for i in df['gene']:
        if (i not in dyt):
            dyt.append(i)
    dyt

    #plt.figure(figsize=(10,5))
    plt.rc('font', family='Arial')
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["font.weight"] = "15"
    plt.rcParams["axes.labelweight"] = "bold"
    cmap1=sns.color_palette(cmap)
    ax = sns.violinplot(x="gene", y="value",hue="class", data=df,palette=[cmap1[0],cmap1[5]],
                        scale="area",split=True,
                        )
    ax.set_ylabel('Gene Expression',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    ax.set_xlabel('')
    for i in range(len(dyt)):
        ttest = stats.ttest_ind(df[df['gene']==dyt[i]][df[df['gene']==dyt[i]]['class']=='ctrl']['value'], df[df['gene']==dyt[i]][df[df['gene']==dyt[i]]['class']=='experiment']['value'])
        max=df[df['gene']==dyt[i]]['value'].max()
        if(ttest[1]<0.001):
            xing="***"   
        elif(ttest[1]<0.01):
            xing="**"
        elif(ttest[1]<0.05):
            xing="*"
        else:
            xing=' '
        print(ttest[1],xing)
        ax.text(i,max+0.5, xing,ha='center', va='bottom', fontsize=15)
    plt.savefig(save_name,dpi=300,bbox_inches = 'tight')```

```pythongene_expression_plot(normdata,['SKAP2','GABBR2'],eg=['lab2','lab3'],cg=['non2','non3'])```
*Output:*
```0.026041143752939457 *
0.03723929692928013 *
```
```<Figure size 432x288 with 1 Axes>```

