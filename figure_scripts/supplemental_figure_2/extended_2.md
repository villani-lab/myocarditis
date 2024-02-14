Supplementary Figure 2
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
library(circlize)
library(Matrix)
library(glue)
library(ggforestplot)
library(ggbeeswarm)
library(ggrepel)
library(patchwork)
library(lme4)
library(ggstance)
library(knitr)
library(ggpubr)

library(reticulate)
use_python("/projects/home/nealpsmith/.conda/envs/updated_pegasus/bin/python")

setwd('/projects/home/ikernin/github_code/myocarditis/functions')
source('masc.R')
source('blood_abundance.R')
source('tissue_plot_masc.R')
source('de.R')
```

Load Python packages

``` python
import pegasus as pg
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions
```

Read in single-cell data

``` python
adt_df = pg.read_input('/projects/home/ikernin/projects/myocarditis/updated_datasets/adt_dataset.zarr')
```

## Figure 2A

``` r
wd = '/projects/home/sramesh/github/myocarditis'

#### first run de
if (!file.exists(glue('{wd}/output/lineage_de_by_deg_case_control_all_results.csv'))) {
  counts <- read_counts(glue("{wd}/output/pb_counts_by_sample_id_and_lineage.csv"))
  meta <- read_meta(glue("{wd}/output/pb_meta_by_sample_id_and_lineage.csv"))
  meta <- meta %>%
    filter(deg_case_control != "NA") %>%
    filter(!str_detect(lineage, "Doublet"))
  case_control_contrast_vec <- c('deg_case_control', 'case', 'control')
  case_control_lineage_results <- run_de_by_comp_var(counts, meta, glue('{wd}/output/lineage'), case_control_contrast_vec,
                                                     deseq_formula = formula("~ deg_case_control + sex + ici_type"))
} else {
    case_control_lineage_results <- read_csv(glue('{wd}/output/lineage_de_by_deg_case_control_all_results.csv'))
}

#### now plot barplot
colors <- c('tomato4', 'slategray')
cluster_order <- c('pDCs', 'cDCs', 'B and plasma', 'MNP', 'CD4', 'CD8 and NK')
res <- case_control_lineage_results %>%
  mutate(cluster = str_replace_all(cluster, '_', ' '))
res <- res[res$padj < 0.1, ]
ups <- downs <- numeric(length(cluster_order))

for (i in seq_along(cluster_order)) {
  clust <- cluster_order[i]
  up_genes <- res[res$cluster == clust & res$log2FoldChange > 0, ]
  down_genes <- res[res$cluster == clust & res$log2FoldChange < 0, ]
  ups[i] <- nrow(up_genes)
  downs[i] <- -nrow(down_genes)
}

df <- data.frame(cluster = cluster_order, up = ups, down = downs)
df <- df %>% mutate(cluster = factor(cluster, levels = rev(cluster_order)))

ggplot(df, aes(x = up, y = cluster)) +
  geom_bar(stat = "identity", fill = colors[1]) +
  geom_bar(stat = "identity", aes(x = down, y = cluster), fill = colors[2]) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        axis.text = element_text(color = 'black')) +
  labs(x = element_blank(),
       y = element_blank())
```

![](extended_2_files/figure-gfm/fig_2A-1.png)<!-- -->

## Figure 2B

Note: please run “adt\_gating.py” for plotting code

``` python
def filter_df(df, x_gene, y_gene, x_line=None, y_line=None, x_gr=True, y_gr=True):
    # make dataframe with gene counts
    x = np.squeeze(np.asarray(df[:, x_gene].copy().get_matrix('X').todense()))
    y = np.squeeze(np.asarray(df[:, y_gene].copy().get_matrix('X').todense()))
    plot_df = pd.DataFrame({'x': x + 1, 'y': y + 1})

    # get info whether passed threshold
    x_mask = plot_df['x'] > x_line if x_gr else plot_df['x'] < x_line
    y_mask = plot_df['y'] > y_line if y_gr else plot_df['y'] < y_line

    return df[x_mask & y_mask].copy()

# get all filtered dataframes to plot all at once
b_df = filter_df(adt_df,
                 x_gene='cite_CD3E',
                 y_gene='cite_CD19',
                 x_line=1.8,
                 y_line=1.7,
                 x_gr=False,
                 y_gr=True)
t_df = filter_df(adt_df,
                 x_gene='cite_CD3E',
                 y_gene='cite_CD19',
                 x_line=2,
                 y_line=1.7,
                 x_gr=True,
                 y_gr=False)
cd8_df = filter_df(t_df,
                   x_gene='cite_CD4',
                   y_gene='cite_CD8A',
                   x_line=2,
                   y_line=1.9,
                   x_gr=False,
                   y_gr=True)
cd4_df = filter_df(t_df,
                   x_gene='cite_CD4',
                   y_gene='cite_CD8A',
                   x_line=2,
                   y_line=1.7,
                   x_gr=True,
                   y_gr=False)
gate_b_df = filter_df(adt_df,
                      x_gene='cite_CD3E',
                      y_gene='cite_CD19',
                      x_line=1.4,
                      y_line=1.5,
                      x_gr=False,
                      y_gr=False)
nk_df = filter_df(gate_b_df,
                  x_gene='cite_CD14',
                  y_gene='cite_NCAM1',
                  x_line=1.4,
                  y_line=1.6,
                  x_gr=False,
                  y_gr=True)
gate_d_df = filter_df(gate_b_df,
                      x_gene='cite_CD14',
                      y_gene='cite_NCAM1',
                      x_line=2.8,
                      y_line=1.5,
                      x_gr=False,
                      y_gr=False)
pdc_df = filter_df(gate_d_df,
                   x_gene='cite_HLA-DRB1',
                   y_gene='cite_IL3RA',
                   x_line=1.5,
                   y_line=2.8,
                   x_gr=True,
                   y_gr=True)
gate_e_df = filter_df(gate_d_df,
                      x_gene='cite_HLA-DRB1',
                      y_gene='cite_IL3RA',
                      x_line=5,
                      x_gr=False,
                      y_line=2.0,
                      y_gr=False
                      )
e_df = filter_df(gate_e_df,
                 x_gene='cite_ITGAX',
                 y_gene='cite_CD1C',
                 x_line=7,  # one very large outlier distorting scale
                 x_gr=False,
                 y_line=0,
                 y_gr=True
                 )
dc2_df = filter_df(gate_e_df,
                   x_gene='cite_ITGAX',
                   y_gene='cite_CD1C',
                   x_line=2,
                   y_line=1.7,
                   x_gr=True,
                   y_gr=True)
gate_f_df = filter_df(gate_e_df,
                      x_gene='cite_ITGAX',
                      y_gene='cite_CD1C',
                      x_line=5,
                      x_gr=False,
                      y_line=1.6,
                      y_gr=False
                      )
mnp_df = filter_df(gate_f_df,
                   x_gene='cite_ITGAX',
                   y_gene='cite_CD94',
                   x_line=1.7,
                   y_line=1.6,
                   x_gr=True,
                   y_gr=False)

# add gate to obs data
b_df.obs['gate'] = 'b_gate'
cd8_df.obs['gate'] = 'cd8_gate'
cd4_df.obs['gate'] = 'cd4_gate'
nk_df.obs['gate'] = 'nk_gate'
pdc_df.obs['gate'] = 'pdc_gate'
dc2_df.obs['gate'] = 'dc2_gate'
mnp_df.obs['gate'] = 'mnp_gate'

# get cell barcodes assigned to each gate
obs_to_get = ['lineage', 'condition', 'sample_id', 'steroid_treatment', 'gate']
df_dict = {"b_df": b_df.obs[obs_to_get],
           "cd4_df": cd4_df.obs[obs_to_get],
           "cd8_df": cd8_df.obs[obs_to_get],
           "nk_df": nk_df.obs[obs_to_get],
           "pdc_df": pdc_df.obs[obs_to_get],
           "dc2_df": dc2_df.obs[obs_to_get],
           "mnp_df": mnp_df.obs[obs_to_get]}
gate_obs = pd.concat(df_dict.values())

# get obs from blood dataset
gate_obs = gate_obs[['gate']]
blood_obs = pd.read_csv('/projects/home/sramesh/myo_final/blood/final/myo_blood_global_obs.csv',
                        usecols=['barcodekey', 'deg_case_control', 'lineage', 'sample_id'])
blood_obs = blood_obs.set_index('barcodekey')
gate_masc_obs = pd.merge(gate_obs, blood_obs, how='left', on='barcodekey')
gate_masc_obs.to_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/gate_masc_obs.csv')

# get number of control and myocarditis cells
print(f'All adt cells {adt_df.shape[0]}')
print(f'Control adt cells {adt_df[adt_df.obs["condition"] == "control"].shape[0]}')
print(f'Myocarditis adt cells {adt_df[adt_df.obs["condition"] == "myocarditis"].shape[0]}')
```

    ## All adt cells 294573
    ## Control adt cells 111676
    ## Myocarditis adt cells 182897

``` r
knitr::include_graphics("/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/adt_scatterplots.png")
```

<img src="../../../../projects/myocarditis/updated_datasets/figures/adt_scatterplots.png" width="800px" height="460px" />

## Figure 2C

Note: please run “adt\_gating.py” for plotting code

``` r
knitr::include_graphics("/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/adt_hexbins.pdf")
```

<embed src="../../../../projects/myocarditis/updated_datasets/figures/adt_hexbins.pdf" width="800px" height="460px" type="application/pdf" />

## Figure 2D

``` r
# numbers from previous python chunk
n_adt_cells <- 294573
n_ctrl <- 111676
n_myo <- 182897
df <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/gate_masc_obs.csv') %>% dplyr::select(-1)

# plot %gex lineages per gate
df %>%
  group_by(gate) %>%
  mutate(n_cells_gate = n()) %>%
  group_by(gate,
           lineage) %>%
  mutate(n_cells_lineage = n(),
         perc_cells_gate = n_cells_lineage/n_cells_gate * 100) %>%
  dplyr::select(gate,
         lineage,
         n_cells_gate,
         n_cells_lineage,
         perc_cells_gate) %>%
  distinct() %>%
  ggplot(aes(x = perc_cells_gate,
             y = lineage)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=sprintf("%2.1f", perc_cells_gate)),
            hjust = -0.1) +
  scale_x_continuous(limits = c(0,100)) +
  facet_wrap(~gate) +
  theme_bw() +
  labs(
    title = 'Gate accuracy',
    x = '% of cells in gate',
    y = 'GEX lineage'
  )
```

![](extended_2_files/figure-gfm/fig_2d-1.png)<!-- -->

## Figure 2E

``` r
df <- df %>%
  mutate(merged_gate = case_when(
    gate %in% c("cd8_gate", "nk_gate") ~ 'cd8_nk_gate',
    TRUE ~ gate
  ))

# for ungrouped gates
ungrouped_masc_res <- run_masc(df, c('deg_case_control', 'control', 'case'), 'gate')
ungrouped_abundance <- create_abundance(df, 'deg_case_control', 'gate')
tissue_odds_ratio_and_box_plot(ungrouped_abundance %>%
                                 dplyr::rename('condition' = 'deg_case_control') %>%
                                 mutate(condition = factor(condition, levels = c('control', 'case'))),
                               ungrouped_masc_res %>%
                                 dplyr::rename("conditionmyocarditis.OR" = "deg_case_controlcase.OR",
                                               "conditionmyocarditis.OR.95pct.ci.lower" = "deg_case_controlcase.OR.95pct.ci.lower",
                                               "conditionmyocarditis.OR.95pct.ci.upper" = "deg_case_controlcase.OR.95pct.ci.upper"),
                               'gate',
                               'Ungrouped Gates')
```

![](extended_2_files/figure-gfm/fig_2e-1.png)<!-- -->

``` r
# for grouped gates
grouped_masc_res <- run_masc(df, c('deg_case_control', 'control', 'case'), 'merged_gate')
grouped_abundance <- create_abundance(df, 'deg_case_control', 'merged_gate')
tissue_odds_ratio_and_box_plot(grouped_abundance %>%
                                 dplyr::rename('condition' = 'deg_case_control') %>%
                                 mutate(condition = factor(condition, levels = c('control', 'case'))),
                               grouped_masc_res %>%
                                 dplyr::rename("conditionmyocarditis.OR" = "deg_case_controlcase.OR",
                                               "conditionmyocarditis.OR.95pct.ci.lower" = "deg_case_controlcase.OR.95pct.ci.lower",
                                               "conditionmyocarditis.OR.95pct.ci.upper" = "deg_case_controlcase.OR.95pct.ci.upper"),
                               'merged_gate',
                               'Grouped Gates')
```

![](extended_2_files/figure-gfm/fig_2e-2.png)<!-- -->