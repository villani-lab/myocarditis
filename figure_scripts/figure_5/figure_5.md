Figure 5
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
source('tissue_troponin_abundance.R')
source('blood_fatal_abundance.R')
source('de.R')
```

Load Python packages

``` python
import pegasus as pg
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions
```

Read in single-cell data

``` python
tissue_myeloid = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_myeloid.zarr')
```

``` python
blood_myeloid = pg.read_input('/projects/home/ikernin/projects/myocarditis/github_datasets/blood_myeloid.zarr')
```

## Figure 5A

``` python
python_functions.plot_umap(tissue_myeloid, 'Tissue: Myeloid', python_functions.tissue_mnp_pal, wspace=0.9, marker_multiplier=6)
```

<img src="figure_5_files/figure-gfm/fig_5a-1.png" width="1008" />

## Figure 5B

``` python
python_functions.make_gene_dotplot(tissue_myeloid.to_anndata(),
             cluster_order=['7. h-pDC: LILRA4 IRF8',
                            '3. h-MNP: FCGR3A LILRB2',
                            '1. h-MNP: S100A8-low C1QA-low',
                            '6. h-cDC: CLEC9A CD1C',
                            '2. h-MNP: LYVE1 C1QA',
                            '4. h-MNP: S100A12 VCAN',
                            '5. h-MNP: TREM2 APOC1'],
             gene_order=['LILRA4', 'IRF8',
                         'CD14', 'FCGR3A', 'LILRB2',
                         'S100A8',
                         'HLA-DQA1', 'CLEC9A', 'CD1C',
                         'LYVE1', 'C1QA',
                         'S100A12', 'VCAN',
                         'TREM2', 'APOC1'
                         ],
             title='Heart Myeloid')
```

<img src="figure_5_files/figure-gfm/fig_5b-3.png" width="1152" />

## Figure 5c

``` r
# read in masc tissue results (see figure_2.rmd fig_2c)
tissue_global_obs = read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global_obs.csv')
masc_filtered_df  <- masc_filter(tissue_global_obs)
cluster_masc_res <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/cluster_masc_res.csv')
plot_masc_by_cell_type(cluster_masc_res, masc_filtered_df, lineage='Myeloid')
```

![](figure_5_files/figure-gfm/fig_5c-5.png)<!-- -->

## Figure 5D

``` r
tissue_troponin_metadata <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_troponin_metadata.csv')
troponin_filtered_df <- troponin_filter_tissue(tissue_global_obs, tissue_troponin_metadata)

# fit linear model by troponin for DE clusters
select_clusters <- c("h-NK: KLRF1 FCER1G",
                     "h-CD4T: IL7R LTB",
                    "h-CD8T: CD27 LAG3",
                    "h-CD8T: CCL5 NKG7",
                    "h-CD8T: cycling",
                    "h-MNP: S100A8-low C1QA-low",
                    "h-MNP: FCGR3A LILRB2",
                    "h-cDC: CLEC9A CD1C",
                    "Fibroblasts: DCN LUM")
troponin_cluster_percs <- troponin_get_percents_per_level(troponin_filtered_df, level='cluster')
select_cluster_percs <- troponin_cluster_percs %>%
        filter(cluster_names %in% select_clusters)
select_cluster_model <- troponin_fit_model(select_cluster_percs, level='cluster')
kable(select_cluster_model %>%
              select(!c(data, model)) %>%
              unnest(cols = c(trop_coef, trop_se, trop_pval)))

troponin_plot_model(select_cluster_model %>% filter(cluster_names =="h-cDC: CLEC9A CD1C"),
                    select_cluster_percs %>% filter(cluster_names =="h-cDC: CLEC9A CD1C"),
                    "h-cDC: CLEC9A CD1C", level='cluster', point_size = 2.2, type='simple')
```

![](figure_5_files/figure-gfm/fig_5D-1.png)<!-- -->

| cluster\_names             |  trop\_coef |  trop\_se | trop\_pval |      padj |
| :------------------------- | ----------: | --------: | ---------: | --------: |
| Fibroblasts: DCN LUM       |   0.0161216 | 0.0162256 |  0.3438558 | 0.6189405 |
| h-CD4T: IL7R LTB           |   0.0001064 | 0.0061510 |  0.9865381 | 0.9865381 |
| h-CD8T: CCL5 NKG7          |   0.0039344 | 0.0104176 |  0.7135667 | 0.9050230 |
| h-CD8T: CD27 LAG3          |   0.0049633 | 0.0063341 |  0.4514475 | 0.6771712 |
| h-CD8T: cycling            |   0.0058375 | 0.0020689 |  0.0181105 | 0.0814970 |
| h-cDC: CLEC9A CD1C         |   0.0019211 | 0.0006323 |  0.0125033 | 0.0814970 |
| h-MNP: FCGR3A LILRB2       |   0.0021302 | 0.0083788 |  0.8044649 | 0.9050230 |
| h-MNP: S100A8-low C1QA-low |   0.0112636 | 0.0082623 |  0.2027089 | 0.4560950 |
| h-NK: KLRF1 FCER1G         | \-0.0054666 | 0.0030911 |  0.1074143 | 0.3222428 |

## Figure 5F

``` python
# get pseudobulk counts and metadata by donor for all myeloid clusters
# python_functions.get_pseudobulk_info(tissue_myeloid, 'tissue_myeloid')

# get pseudobulk counts and metadata by donor for all mnp clusters
tissue_myeloid.obs['mnp_cell'] = tissue_myeloid.obs['umap_name'].isin(['1. h-MNP: S100A8-low C1QA-low',
                                                               '2. h-MNP: LYVE1 C1QA',
                                                               '3. h-MNP: FCGR3A LILRB2',
                                                               '4. h-MNP: S100A12 VCAN',
                                                               '5. h-MNP: TREM2 APOC1'])
tissue_myeloid.obs['mnp_cell'] = tissue_myeloid.obs['mnp_cell'].replace({True: 'all_mnp', False: 'other'})
# python_functions.get_pseudobulk_info(tissue_myeloid, 'tissue_mnp_grouped', cluster_col='mnp_cell')
```

``` r
# run DE analysis by condition
myeloid_deres <- run_de_by_condition(counts_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_myeloid_pseudocounts.csv',
                               meta_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_myeloid_metainfo.csv',
                               save_name = 'tissue_myeloid')

myeloid_grouped_deres <- run_de_by_condition(counts_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_mnp_grouped_pseudocounts.csv',
                                       meta_filepath = '/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_mnp_grouped_metainfo.csv',
                                       save_name = 'tissue_mnp_grouped')
```

    ## [1] "Cluster 1"
    ## [1] "Cluster 2"
    ## [1] "Cluster 3"
    ## [1] "Cluster 4"
    ## [1] "Cluster 5"
    ## [1] "Cluster 6"
    ## [1] "Cluster 7"
    ## [1] "Cluster 8"
    ## [1] "saving results..."
    ## [1] "Cluster all_mnp"
    ## [1] "Cluster other"
    ## [1] "saving results..."

``` r
# combine de results and meta data for heatmap
myeloid_full_deres <- bind_rows(myeloid_deres %>%
                               mutate(cluster = as.character(cluster)),
                             myeloid_grouped_deres)
myeloid_clusters <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/myeloid_cluster_map.csv')
myeloid_genes <- read_csv('/projects/home/ikernin/projects/myocarditis/github_datasets/myeloid_heatmap_genes.csv')
myeloid_heatmap_df <- get_heatmap_data(myeloid_full_deres, myeloid_genes, myeloid_clusters)
heatmap_df <- myeloid_heatmap_df %>%
  mutate(cluster = cluster_name,
    cluster = case_when(
    cluster == 'all_mnp' ~ 'All MNP',
    TRUE ~ cluster
  )) %>%
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

# define cell color range
heatmap_col_fun <- colorRamp2(c(floor(min(heatmap_mtx)), 0, ceiling(max(heatmap_mtx))),
                              c("blue", "white", "red"))

# split into subclusters and lineage
heatmap_mtx_subcluster <- heatmap_mtx[, !str_detect(colnames(heatmap_mtx), 'All')]
colnames(heatmap_mtx_subcluster) <- str_remove(colnames(heatmap_mtx_subcluster), regex(":.*"))
heatmap_mtx_lineage <- heatmap_mtx[, str_detect(colnames(heatmap_mtx), 'All'), drop=FALSE]



# Main body annotation (FDR) ----------------------------------------------

# get fdr values
fdr_mtx <- heatmap_df %>%
  select(starts_with('padj')) %>%
  replace(is.na(.), Inf) %>%
  rename_with(~str_remove(., "padj_")) %>%
  select(order(colnames(.))) %>%
  as.matrix()

# make sure columns the same
stopifnot(colnames(fdr_mtx) == colnames(heatmap_mtx))

# split into subcluster and lineage
fdr_mtx_subcluster <- fdr_mtx[, !str_detect(colnames(fdr_mtx), 'All')]
fdr_mtx_lineage <- fdr_mtx[, str_detect(colnames(fdr_mtx), 'All'), drop=FALSE]

# make function for plotting fdr value
fdr_func_subcluster <- function(j, i, x, y, width, height, fill){
  if (fdr_mtx_subcluster[i,j] < 0.1){
    grid.circle(x = x, y = y, r = unit(1.5, 'mm'),
                gp = gpar(fill = 'black', col = NA))
  }
}
fdr_func_lineage <- function(j, i, x, y, width, height, fill){
  if (fdr_mtx_lineage[i,j] < 0.1){
    grid.circle(x = x, y = y, r = unit(1.5, "mm"),
                gp = gpar(fill = 'black', col = NA))
  }
}

# create legend for fdr
lgd_fdr = Legend(pch = 16, type = "points", labels = "FDR < 0.1")



# Column annotation (cluster names) ---------------------------------------

# define colors
clust_col_fun <- c('#1b9e77', '#e7298a', '#a6761d', '#252525', '#f43600', '#356F83')
names(clust_col_fun) <- seq(1, 6)
clust_ha <- HeatmapAnnotation(clust_colors = names(clust_col_fun),
                              col = list(clust_colors = clust_col_fun),
                              show_legend = FALSE,
                              show_annotation_name = FALSE,
                              simple_anno_size = unit(3, "mm"))

lineage_col_fun <- c('white')
names(lineage_col_fun) <- 'All'
lineage_ha <- HeatmapAnnotation(lineage = names(lineage_col_fun),
                                col = list(lineage = lineage_col_fun),
                                show_legend = FALSE,
                                show_annotation_name = FALSE,
                                simple_anno_size = unit(3, "mm"),
                                border = T)

# Row annotation (gene names) ---------------------------------------------

# split rows by gene category
stopifnot(rownames(heatmap_mtx_subcluster) == rownames(heatmap_mtx_lineage))
row_split <- str_replace_all(heatmap_df$category, "_", " ")
row_split <- factor(row_split, levels = unique(row_split))


# Plot --------------------------------------------------------------------

ht_subcluster <- Heatmap(heatmap_mtx_subcluster,
                         col = heatmap_col_fun,
                         row_split = row_split,
                         cell_fun = fdr_func_subcluster,
                         top_annotation = clust_ha,
                         name = 'Log2FC',
                         cluster_columns = FALSE,  column_names_side = "top",
                         show_column_names = T, column_names_rot = 45,
                         cluster_rows = FALSE, row_names_side = "left",
                         row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                         row_gap = unit(2, "mm"), border = TRUE,
                         width = ncol(heatmap_mtx_subcluster)*unit(6, "mm"),
                         height = nrow(heatmap_mtx_subcluster)*unit(6, "mm"))

ht_lineage <- Heatmap(heatmap_mtx_lineage,
                      col = heatmap_col_fun,
                      row_split = row_split,
                      cell_fun = fdr_func_lineage,
                      top_annotation = lineage_ha,
                      name = 'Lineage', show_heatmap_legend = FALSE,
                      column_names_gp = gpar(fontface='bold'),
                      cluster_columns = FALSE,  column_names_side = "top",
                      show_column_names = T, column_names_rot = 45,
                      cluster_rows = FALSE, row_names_side = "left",
                      row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                      row_gap = unit(2, "mm"), border = TRUE,
                      width = ncol(heatmap_mtx_lineage)*unit(6, "mm"),
                      height = nrow(heatmap_mtx_lineage)*unit(6, "mm"))

draw(ht_lineage + ht_subcluster,
     annotation_legend_list = list(lgd_fdr),
     merge_legends = TRUE)
```

![](figure_5_files/figure-gfm/fig_5f_heatmap-1.png)<!-- -->

## Figure 5G

``` python
fig_5g_genes = ['IFITM1', 'CXCL10', 'HLA-DQB2', 'STAT1']
python_functions.multi_hexfp_by_condition(tissue_myeloid, fig_5g_genes, cmap = python_functions.blues_cmap, gridsize=200)
```

    ##   0%|                                                                                                                                                                 | 0/4 [00:00<?, ?it/s] 25%|######################################2                                                                                                                  | 1/4 [00:00<00:01,  2.48it/s] 50%|############################################################################5                                                                            | 2/4 [00:00<00:00,  2.56it/s] 75%|##################################################################################################################7                                      | 3/4 [00:01<00:00,  2.61it/s]100%|#########################################################################################################################################################| 4/4 [00:01<00:00,  2.47it/s]

<img src="figure_5_files/figure-gfm/fig_5G-1.png" width="1152" />