# t_visualize_colorsystem
*Converted from: omicverse/omicverse_guide/docs/Tutorials-plotting/t_visualize_colorsystem.ipynb*

# Color system

In OmicVerse, we offer a color system based on Eastern aesthetics, featuring 384 representative colors derived from the Forbidden City. We will utilize these colors in combination for future visualizations.

All color come from the book: "中国传统色：故宫里的色彩美学" (ISBN: 9787521716054)

```pythonimport omicverse as ov
import scanpy as sc
#import scvelo as scv
ov.plot_set()```

We utilized single-cell RNA-seq data (GEO accession: GSE95753) obtained from the dentate gyrus of the hippocampus in mouse.

```pythonadata = ov.read('data/DentateGyrus/10X43_1.h5ad')
adata```
*Output:*
```AnnData object with n_obs × n_vars = 2930 × 13913
    obs: 'clusters', 'age(days)', 'clusters_enlarged'
    uns: 'clusters_colors'
    obsm: 'X_umap'
    layers: 'ambiguous', 'spliced', 'unspliced'```

## Understanding the Color System

In OmicVerse, we offer a color system based on Eastern aesthetics, featuring 384 representative colors derived from the Forbidden City. We will utilize these colors in combination for future visualizations.

```pythonfb=ov.pl.ForbiddenCity()```

```pythonfrom IPython.display import HTML
HTML(fb.visual_color(loc_range=(0,384),
                    num_per_row=24))```
*Output:*
```<IPython.core.display.HTML object>```

we can get a color using `get_color`

```pythonfb.get_color(name='凝夜紫')```
*Output:*
```     num name       name_en     color_rgb color_html
161  161  凝夜紫  Noctilucence  (68, 36, 84)    #442454```

## Default Colormap

We have provided a range of default colors including `green`, `red`, `pink`, `purple`, `yellow`, `brown`, `blue`, and `grey`. Each of these colors comes with its own set of sub-colormaps, providing a more granular level of color differentiation.

Here's a breakdown of the sub-colormaps available for each default color:

- **Green**: 
  - `green1`: `Forbidden_Cmap(range(1, 19))`
  - `green2`: `Forbidden_Cmap(range(19, 41))`
  - `green3`: `Forbidden_Cmap(range(41, 62))`

- **Red**: 
  - `red1`: `Forbidden_Cmap(range(62, 77))`
  - `red2`: `Forbidden_Cmap(range(77, 104))`

- **Pink**: 
  - `pink1`: `Forbidden_Cmap(range(104, 134))`
  - `pink2`: `Forbidden_Cmap(range(134, 148))`

- **Purple**: 
  - `purple1`: `Forbidden_Cmap(range(148, 162))`
  - `purple2`: `Forbidden_Cmap(range(162, 176))`

- **Yellow**: 
  - `yellow1`: `Forbidden_Cmap(range(176, 196))`
  - `yellow2`: `Forbidden_Cmap(range(196, 207))`
  - `yellow3`: `Forbidden_Cmap(range(255, 276))`

- **Brown**: 
  - `brown1`: `Forbidden_Cmap(range(207, 228))`
  - `brown2`: `Forbidden_Cmap(range(228, 246))`
  - `brown3`: `Forbidden_Cmap(range(246, 255))`
  - `brown4`: `Forbidden_Cmap(range(276, 293))`

- **Blue**: 
  - `blue1`: `Forbidden_Cmap(range(293, 312))`
  - `blue2`: `Forbidden_Cmap(range(312, 321))`
  - `blue3`: `Forbidden_Cmap(range(321, 333))`
  - `blue4`: `Forbidden_Cmap(range(333, 339))`

- **Grey**: 
  - `grey1`: `Forbidden_Cmap(range(339, 356))`
  - `grey2`: `Forbidden_Cmap(range(356, 385))`

Each main color can be represented as a combination of its sub-colormaps:

- `green = green1 + green2 + green3`
- `red = red1 + red2`
- `pink = pink1 + pink2`
- `purple = purple1 + purple2`
- `yellow = yellow1 + yellow2 + yellow3`
- `brown = brown1 + brown2 + brown3 + brown4`
- `blue = blue1 + blue2 + blue3 + blue4`
- `grey = grey1 + grey2`

These colormaps can be utilized in various applications where color differentiation is necessary, providing flexibility in visual representation.

`palette` is the argument we need to revise

```pythonimport matplotlib.pyplot as plt
fig, axes = plt.subplots(1,3,figsize=(9,3)) 
ov.pl.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=["clusters"],
                   palette=fb.red[:],
                   ncols=3,
                show=False,
                legend_loc=None,
                    ax=axes[0])

ov.pl.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=["clusters"],
                   palette=fb.pink1[:],
                   ncols=3,show=False,
                legend_loc=None,
                    ax=axes[1])

ov.pl.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=["clusters"],
                   palette=fb.red1[:4]+fb.blue1,
                   ncols=3,show=False,
                    ax=axes[2])


```
*Output:*
```<AxesSubplot: title={'center': 'clusters'}, xlabel='X_umap1', ylabel='X_umap2'>```
```<Figure size 720x240 with 3 Axes>```

```pythoncolor_dict={'Astrocytes': '#e40414',
 'Cajal Retzius': '#ec5414',
 'Cck-Tox': '#ec4c2c',
 'Endothelial': '#d42c24',
 'GABA': '#2c5ca4',
 'Granule immature': '#acd4ec',
 'Granule mature': '#a4bcdc',
 'Microglia': '#8caccc',
 'Mossy': '#8cacdc',
 'Neuroblast': '#6c9cc4',
 'OL': '#6c94cc',
 'OPC': '#5c74bc',
 'Radial Glia-like': '#4c94c4',
 'nIPC': '#3474ac'}

ov.pl.embedding(adata,
                   basis='X_umap',
                    frameon='small',
                   color=["clusters"],
                   palette=color_dict,
                   ncols=3,show=False,
                    )

```
*Output:*
```<AxesSubplot: title={'center': 'clusters'}, xlabel='X_umap1', ylabel='X_umap2'>```
```<Figure size 320x320 with 1 Axes>```

## Segmented Colormap

When we need to create a continuous color gradient, we will use another function: `get_cmap_seg`, and we can combine the colors we need for visualization.

```pythoncolors=[
    fb.get_color_rgb('群青'),
    fb.get_color_rgb('半见'),
    fb.get_color_rgb('丹罽'),
]
fb.get_cmap_seg(colors)```
*Output:*
```<matplotlib.colors.LinearSegmentedColormap at 0x7f6012bc20e0>```

```pythoncolors=[
    fb.get_color_rgb('群青'),
    fb.get_color_rgb('山矾'),
    fb.get_color_rgb('丹罽'),
]
fb.get_cmap_seg(colors)```
*Output:*
```<matplotlib.colors.LinearSegmentedColormap at 0x7f6012c4c5e0>```

```pythoncolors=[
    fb.get_color_rgb('山矾'),
    fb.get_color_rgb('丹罽'),
]
fb.get_cmap_seg(colors)```
*Output:*
```<matplotlib.colors.LinearSegmentedColormap at 0x7f6012accb80>```

```pythonov.pl.embedding(adata,
                basis='X_umap',
                frameon='small',
                color=["Sox7"],
                cmap=fb.get_cmap_seg(colors),
                ncols=3,show=False,
                #vmin=-1,vmax=1
                )
```
*Output:*
```<AxesSubplot: title={'center': 'Sox7'}, xlabel='X_umap1', ylabel='X_umap2'>```
```<Figure size 320x320 with 2 Axes>```

