Supplemental Figure 2
================

## Setup

Load R libraries

``` r
library(tidyverse)
library(glue)
library(RColorBrewer)
library(reticulate)
use_python("/projects/home/nealpsmith/.conda/envs/updated_pegasus/bin/python")

setwd('/projects/home/ikernin/github_code/myocarditis/functions')
source('stacked_bar.R')
```

Load Python packages

``` python
import pegasus as pg
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions
```

Read in single-cell data

``` python
tissue_t = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_t.zarr')
```

``` python
tissue_b = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_b.zarr')
```

## Supplemental Figure 2A

``` python
tissue_t.obs['Condition'] = [x.capitalize() for x in tissue_t.obs['condition']]
tissue_t_anndata = tissue_t.to_anndata()
sc.tl.embedding_density(tissue_t_anndata, groupby='Condition')
sc.pl.embedding_density(tissue_t_anndata, basis='umap', key=f'umap_density_Condition', colorbar_loc=None)
```

    ## ... storing 'Condition' as categorical

<img src="supp_2_files/figure-gfm/supp_2a-1.png" width="1514" />

## Supplemental Figure 2B

``` python
stacked_bar_df = python_functions.get_stacked_bar_df(tissue_t, 't')
stacked_bar_order = tissue_t.obs['umap_name'].cat.categories.values
```

    ## Getting stacked bar info for: t

``` r
stacked_bar_order = py$stacked_bar_order[!str_detect(py$stacked_bar_order, 'Doublets')]
plot_clust_perc_by_donor(py$stacked_bar_df, 't', cluster_order = stacked_bar_order)
```

    ## Warning in py_to_r.pandas.core.frame.DataFrame(x): index contains duplicated
    ## values: row names not set

![](supp_2_files/figure-gfm/supp_2b_plot-3.png)<!-- -->

## Supplemental Figure 2E

``` python
python_functions.hex_featureplot(tissue_t, 'STMN1', cmap=python_functions.blues_cmap)
```

<img src="supp_2_files/figure-gfm/unnamed-chunk-1-1.png" width="576" />

## Supplemental Figure 2F

``` python
python_functions.plot_umap(tissue_b, 'Tissue: B and Plasma', python_functions.tissue_b_pal, marker_multiplier=5)
```

<img src="supp_2_files/figure-gfm/unnamed-chunk-2-3.png" width="960" />
