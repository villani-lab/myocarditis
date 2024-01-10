Supplemental Figure 8
================

## Set up

Load R libraries

``` r
# load packages
library(tidyverse)
library(rmarkdown)
library(rlang)
library(parameters)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(glue)
library(ggforestplot)
library(ggbeeswarm)
library(patchwork)
library(lme4)
library(ggstance)
library(DESeq2)
library(knitr)

library(reticulate)
use_python("/projects/home/nealpsmith/.conda/envs/updated_pegasus/bin/python")

setwd('/projects/home/ikernin/github_code/myocarditis/functions')
source('masc.R')
source('plot_masc.R')
source('de.R')
source('stacked_bar.R')
```

Load Python packages

``` python
import pegasus as pg
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions
import scanpy as sc
```

Read in single-cell data

``` python
tissue_nonimmune = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune.zarr')
```

    ## 2024-01-10 14:29:16,199 - pegasusio.readwrite - INFO - zarr file '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune.zarr' is loaded.
    ## 2024-01-10 14:29:16,200 - pegasusio.readwrite - INFO - Function 'read_input' finished in 0.65s.

``` python
tissue_global = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global.zarr')
```

    ## 2024-01-10 14:29:17,070 - pegasusio.readwrite - INFO - zarr file '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global.zarr' is loaded.
    ## 2024-01-10 14:29:17,070 - pegasusio.readwrite - INFO - Function 'read_input' finished in 0.87s.

## Supplemental Figure 8A

``` python
python_functions.plot_umap(tissue_nonimmune, 'Tissue: Non-Immune', python_functions.tissue_nonimmune_pal, fig_size=(6, 5), marker_multiplier=15)
```

<img src="supp_8_files/figure-gfm/fig_s8a-1.png" width="1152" />

## Supplemental Figure 8B

``` python
# split data into endothelial and non-endothelial cells
endothelial_clusters = ['1. Capillary EC: RGCC CA4',
                        '9. Arterial EC: RBP7 BTNL9',
                        '5. Arterial EC: HEY1 SEMA3G',
                        '2. Capillary EC: CXCL2 JUN',
                        '11. Venous EC: ACKR1 PLVAP',
                        '14. Inflammatory EC: CXCL9 IDO1',
                        '13. Endocardial: NPR3 POSTN',
                        '3. Capillary EC: VWF CD36',
                        '15. Lymphatic EC: PDPN CCL21']
non_endothelial_clusters = ['8. Smooth muscle: TAGLN MYH11',
                            '7. Pericytes: KCNJ8 ABCC9',
                            '4. Fibroblasts: DCN LUM',
                            '16. Myofibroblasts: ACTA2 COL14A1',
                            '6. Pericytes: HLA-DRB1 S1PR1',
                            '12. Cardiomyocytes: TNNT2 MB',
                            '17. Neuronal: PLP1 NRXN1',
                            '10. Pericytes: RGS5 KCNJ8-low']

endothelial = tissue_nonimmune[tissue_nonimmune.obs['umap_name'].isin(endothelial_clusters)].copy()
non_endothelial = tissue_nonimmune[tissue_nonimmune.obs['umap_name'].isin(non_endothelial_clusters)].copy()

# plot dotplots

python_functions.make_gene_dotplot(endothelial.to_anndata(),
             cluster_order=endothelial_clusters,
             gene_order=['RGCC', 'CA4',  # Capillary EC1
                         'RBP7', 'BTNL9',  # Arterial EC2
                         'HEY1', 'SEMA3G',  # Arterial EC1
                         'CXCL2', 'JUN',  # Capillary EC2
                         'ACKR1', 'PLVAP',  # Venous EC
                         'CXCL9', 'IDO1',  # Inflammatory EC
                         'NPR3', 'POSTN',  # Endocardia
                         'VWF', 'CD36',  # Capillary EC3
                         'PDPN', 'CCL21'  # Lymphatic ECs
                         ],
             title='Endothelial')
python_functions.make_gene_dotplot(non_endothelial.to_anndata(),
             cluster_order=non_endothelial_clusters,
             gene_order=['TAGLN', 'MYH11',  # Smooth Muscle
                         'KCNJ8', 'ABCC9',  # Pericytes2
                         'DCN', 'LUM',  # Fibroblasts
                         'ACTA2', 'COL14A1', 'ID4',  # Myofibroblasts
                         'HLA-DRB1', 'S1PR1',  # Pericytes1
                         'TNNT2', 'MB',  # Cardiomyoctyes
                         'PLP1', 'NRXN1',  # Neuronal
                         'RGS5'  # Pericytes3
                         ],
             title='Non-Endothelial')
```

<img src="supp_8_files/figure-gfm/fig_s8b-3.png" width="1152" /><img src="supp_8_files/figure-gfm/fig_s8b-4.png" width="1152" />

## Supplemental Figure 8C

``` r
# read in masc tissue results (see figure_2.rmd fig_2c)
tissue_global_obs = read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global_obs.csv')
masc_filtered_df  <- masc_filter(tissue_global_obs)
cluster_masc_res <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/cluster_masc_res.csv')
plot_masc_by_cell_type(cluster_masc_res, masc_filtered_df, lineage='Nonimmune')
```

![](supp_8_files/figure-gfm/fig_s8c-7.png)<!-- -->

## Supplemental Figure 8D

``` python

supp_fig8d_genes = ['VWF', 'CA4', 'ACKR1', "HEY1", "ABCC9", "TAGLN", "MYH11", "ACTA2", "DCN", "PLP1"]
python_functions.multi_hex_featureplot(tissue_nonimmune, supp_fig8d_genes, ncol=3, cmap=python_functions.blues_cmap, gridsize=200)
```

    ##   0%|                                                                                                                                                                                                                | 0/10 [00:00<?, ?it/s] 10%|####################                                                                                                                                                                                    | 1/10 [00:00<00:06,  1.34it/s] 20%|########################################                                                                                                                                                                | 2/10 [00:01<00:05,  1.40it/s] 30%|############################################################                                                                                                                                            | 3/10 [00:02<00:04,  1.40it/s] 40%|################################################################################                                                                                                                        | 4/10 [00:02<00:04,  1.45it/s] 50%|####################################################################################################                                                                                                    | 5/10 [00:03<00:03,  1.43it/s] 60%|########################################################################################################################                                                                                | 6/10 [00:04<00:02,  1.47it/s] 70%|############################################################################################################################################                                                            | 7/10 [00:04<00:02,  1.45it/s] 80%|################################################################################################################################################################                                        | 8/10 [00:05<00:01,  1.48it/s] 90%|####################################################################################################################################################################################                    | 9/10 [00:06<00:00,  1.45it/s]100%|#######################################################################################################################################################################################################| 10/10 [00:06<00:00,  1.49it/s]

<img src="supp_8_files/figure-gfm/fig_s8d-1.png" width="1728" />

## Supplemental Figure 8E

``` python

python_functions.hex_plot(tissue_nonimmune, "% Mito", n_genes=False, cmap=python_functions.blues_cmap)
```

<img src="supp_8_files/figure-gfm/fig_s8e-3.png" width="576" />

``` python
python_functions.hex_plot(tissue_nonimmune, "# Genes", n_genes=True, cmap=python_functions.blues_cmap)
```

<img src="supp_8_files/figure-gfm/fig_s8e-4.png" width="576" />

## Supplemental Figure 8F

``` python

tissue_nonimmune.obs['Condition'] = [x.capitalize() for x in tissue_nonimmune.obs['condition']]
tissue_nonimmune = tissue_nonimmune.to_anndata()
sc.tl.embedding_density(tissue_nonimmune, groupby='Condition')
sc.pl.embedding_density(tissue_nonimmune, basis='umap', key=f'umap_density_Condition')
```

<img src="supp_8_files/figure-gfm/fig_s8f-7.png" width="1514" />

## Supplemental Figure 8G

``` python

tissue_nonimmune = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune.zarr')
```

    ## 2024-01-10 14:30:11,273 - pegasusio.readwrite - INFO - zarr file '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune.zarr' is loaded.
    ## 2024-01-10 14:30:11,273 - pegasusio.readwrite - INFO - Function 'read_input' finished in 0.62s.

``` python
stacked_bar_df = python_functions.get_stacked_bar_df(tissue_nonimmune, 'nonimmune')
```

    ## Getting stacked bar info for: nonimmune

``` python
stacked_bar_order = tissue_nonimmune.obs['umap_name'].cat.categories.values
```

``` r
stacked_bar_order = py$stacked_bar_order[!str_detect(py$stacked_bar_order, 'Doublets')]
plot_clust_perc_by_donor(py$stacked_bar_df, 'nonimmune', cluster_order = stacked_bar_order)
```

    ## Warning in py_to_r.pandas.core.frame.DataFrame(x): index contains duplicated
    ## values: row names not set

![](supp_8_files/figure-gfm/supp_8g_plot-9.png)<!-- -->

## Supplemental Figure 8h

``` python
python_functions.get_pseudobulk_info(tissue_nonimmune, 'tissue_nonimmune')
```

``` r
nonimmune_deres <- run_de_by_condition(counts_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune_pseudocounts.csv',
                    meta_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_nonimmune_metainfo.csv',
                    save_name = 'tissue_nonimmune')
```

``` r
# plot number of up/down-expressed genes
nonimmune_cluster_map <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/nonimmune_cluster_map.csv')
nonimmune_deres <- nonimmune_deres %>%
  left_join(nonimmune_cluster_map, by = c('cluster' = 'cluster_number'))


get_up_down <- function(data){
  n_up <- data %>% filter(log2FoldChange > 0) %>% nrow()
  n_down <- data %>% filter(log2FoldChange < 0) %>% nrow()

  return(list('n_up' = n_up, 'n_down' = n_down))
}

plot_df <- nonimmune_deres %>%
  filter(abs(log2FoldChange) >= 0, padj <= 0.1) %>%
  group_by(cluster_name) %>%
  nest() %>%
  filter(cluster_name != 'Doublets/RBCs') %>%
  mutate(deg = map(data, get_up_down)) %>%
  unnest_wider(deg) %>%
  select(cluster_name, n_up, n_down) %>%
  mutate(n_total = n_up + n_down,
         n_down = -1 * n_down) %>%
  pivot_longer(c(n_up, n_down), names_to = 'direction')

breaks <- pretty(plot_df$value) # get plot break values

ggplot(plot_df, aes(x = reorder(cluster_name, n_total), y = value, fill = direction)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(breaks = breaks,
                     labels = abs(breaks)) +
  scale_fill_manual(values = c('#708090', '#8B3626')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16)) +
  labs(x = element_blank(),
       y = element_blank(),
       title = 'Tissue Non-immune DE Results',
       subtitle = 'Padj < 0.01')
```

![](supp_8_files/figure-gfm/fig_s8h_plot_sig_n_de-1.png)<!-- -->

``` r
nonimmune_clusters <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/nonimmune_cluster_map.csv')
heatmap_genes <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/nonimmune_heatmap_genes.csv')
nonimmune_heatmap_df <- get_heatmap_data(nonimmune_deres, heatmap_genes, nonimmune_clusters)

heatmap_df <- nonimmune_heatmap_df %>%
  mutate(cluster = cluster_name) %>%
  select(!cluster_name)


# Format main body --------------------------------------------------------

# set category order
category_levels <- sort(unique(heatmap_df$category))
category_levels <- c(category_levels[category_levels != 'Other'], 'Other')

# reformat from long to wide
heatmap_df <- heatmap_df %>%
  distinct() %>%
  pivot_wider(names_from = cluster, values_from = c(log2FoldChange, padj)) %>%
  filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = category_levels)) %>%
  arrange(category)

# get information for the main body's cells
heatmap_mtx <- heatmap_df %>%
  select(starts_with("log2FoldChange")) %>%
  replace(is.na(.), 0) %>%
  rename_with(~str_remove(., "log2FoldChange_")) %>%
  select(order(colnames(.))) %>%
  as.matrix()
rownames(heatmap_mtx) <- heatmap_df$gene_symbol
colnames(heatmap_mtx) <- str_remove(colnames(heatmap_mtx), regex(":.*"))

# move CXCL9 before CXCL10
which(rownames(heatmap_mtx) == 'CXCL9')
which(rownames(heatmap_mtx) == 'CXCL10')
row_order <- c(rownames(heatmap_mtx)[1:(which(rownames(heatmap_mtx) == 'CXCL10')-1)],
               'CXCL9',
               rownames(heatmap_mtx)[which(rownames(heatmap_mtx) == 'CXCL10'):(which(rownames(heatmap_mtx) == 'CXCL9')-1)],
               rownames(heatmap_mtx)[(which(rownames(heatmap_mtx) == 'CXCL9')+1): nrow(heatmap_mtx)])
heatmap_mtx <- heatmap_mtx[row_order, ]


# define cell color range
heatmap_col_fun <- colorRamp2(c(floor(min(heatmap_mtx)), 0, ceiling(max(heatmap_mtx))),
                              c("blue", "white", "red"))

# Main body annotation (FDR) ----------------------------------------------

# get fdr values
fdr_mtx <- heatmap_df %>%
  select(starts_with('padj')) %>%
  replace(is.na(.), Inf) %>%
  rename_with(~str_remove(., "padj_")) %>%
  select(order(colnames(.))) %>%
  as.matrix()
colnames(fdr_mtx) <- str_remove(colnames(fdr_mtx), regex(":.*"))
rownames(fdr_mtx) <- heatmap_df$gene_symbol

fdr_mtx <- fdr_mtx[row_order, ]

# make sure same cols and rows
stopifnot(colnames(fdr_mtx) == colnames(heatmap_mtx))
stopifnot(rownames(fdr_mtx) == rownames(heatmap_mtx))

# make function for plotting fdr value
fdr_func <- function(j, i, x, y, width, height, fill){
  if (fdr_mtx[i,j] < 0.1){
    grid.circle(x = x, y = y, r = unit(1.5, 'mm'),
                gp = gpar(fill = 'black', col = NA))
  }
}

# create legend for fdr
lgd_fdr = Legend(pch = 16, type = "points", labels = "FDR < 0.1")


# put together lineage heatmaps -------------------------------------------

# define column splits
lineage_list <- list("endothelial" = c("5. Arterial EC",
                                       "9. Arterial EC",
                                       "1. Capillary EC",
                                       "2. Capillary EC",
                                       "3. Capillary EC",
                                       "13. Endocardial",
                                       "14. Inflammatory EC",
                                       "15. Lymphatic EC",
                                       "11. Venous EC"),
                     "mural" = c("6. Pericytes",
                                 "7. Pericytes",
                                 "10. Pericytes",
                                 "8. Smooth muscle"),
                     "fibroblasts" = c("4. Fibroblasts",
                                       "16. Myofibroblasts"),
                     "cardiomyocytes" = "12. Cardiomyocytes",
                     "nueronal" = "17. Neuronal")


make_lineage_heatmap <- function(lineage, n_lineage){

  # make sure same cols and rows
  stopifnot(colnames(fdr_mtx) == colnames(heatmap_mtx))
  stopifnot(rownames(fdr_mtx) == rownames(heatmap_mtx))

  # subset heatmap for given lineage
  lineage_ht <- heatmap_mtx[, lineage_list[lineage][[1]], drop=F]
  lineage_fdr <- fdr_mtx[, lineage_list[lineage][[1]], drop=F]

  # make function for plotting fdr value
  fdr_func <- function(j, i, x, y, width, height, fill){
    if (lineage_fdr[i,j] < 0.1){
      grid.circle(x = x, y = y, r = unit(1.5, 'mm'),
                  gp = gpar(fill = 'black', col = NA))
    }
  }

  # set column annotation colors
  lineage_cols <- list("endothelial" = c('#fb8072',
                                         '#ffed6f',
                                        '#a65628',
                                        '#f781bf',
                                        '#8dd3c7',
                                        '#4ba93b',
                                        '#5779bb',
                                        '#927acc',
                                        '#d95f02'),
                      "mural" = c('#fdb462',
                                  '#fccde5',
                                  '#c4eaff',
                                  '#bc80bd'),
                      "fibroblasts" = c('#bebada',
                                        '#bf3947'),
                      "cardiomyocytes" = '#737373',
                      "nueronal" = '#f48758')


  clust_col_fun <- lineage_cols[lineage][[1]]
  names(clust_col_fun) <- seq(1, ncol(lineage_ht))

  # create column heatmap annotation
  clust_ha <- HeatmapAnnotation(clust_colors = names(clust_col_fun),
                                col = list(clust_colors = clust_col_fun),
                                show_legend = FALSE,
                                show_annotation_name = FALSE,
                                simple_anno_size = unit(3, "mm"))

  # split rows by gene category
  categories_df <- rownames(lineage_ht) %>%
          as_tibble() %>%
          mutate(gene_symbol = value) %>%
          left_join(heatmap_df %>% select(gene_symbol,category))
  row_split <- str_replace_all(categories_df$category, "_", " ")
  row_split <- factor(row_split, levels = unique(row_split))

  # format for multiple heatmaps
  heatmap_name = if (n_lineage == 1) "Log2FC" else as.character(n_lineage)
  draw_legend = if(n_lineage == 1) TRUE else FALSE

  # create heatmap
  ht <- Heatmap(lineage_ht, name = heatmap_name,
                col = heatmap_col_fun,
                row_split = row_split,
                cell_fun = fdr_func,
                top_annotation = clust_ha,
                show_heatmap_legend = draw_legend,
                cluster_columns = FALSE,  column_names_side = "top",
                show_column_names = T, column_names_rot = 45,
                cluster_rows = FALSE, row_names_side = "left",
                row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                row_gap = unit(2, "mm"), border = TRUE,
                width = ncol(lineage_ht)*unit(6, "mm"),
                height = nrow(lineage_ht)*unit(6, "mm"))
  return(ht)
}

ht_endothelial <- make_lineage_heatmap('endothelial', n_lineage = 1)
ht_mural <- make_lineage_heatmap("mural", 2)
ht_fibroblasts <- make_lineage_heatmap("fibroblasts" , 3)
ht_cardiomyocytes <- make_lineage_heatmap("cardiomyocytes" , 4)
ht_nueronal <- make_lineage_heatmap("nueronal" , 5)

draw(ht_endothelial + ht_mural  + ht_fibroblasts + ht_cardiomyocytes + ht_nueronal,
     annotation_legend_list = list(lgd_fdr),
     merge_legends = TRUE)
```

![](supp_8_files/figure-gfm/fig_s8h_heatmap-1.png)<!-- -->

    ## [1] 15
    ## [1] 11