Figure 2
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
library(ggrepel)
library(patchwork)
library(lme4)
library(ggstance)
library(DESeq2)
library(knitr)
library(fgsea)
library(ggpubr)

library(reticulate)
use_python("/projects/home/nealpsmith/.conda/envs/updated_pegasus/bin/python")

setwd('/projects/home/ikernin/github_code/myocarditis/functions')
source('masc.R')
source('de.R')
source('tissue_plot_masc.R')
source('tissue_gsea.R')
source('tissue_troponin_abundance.R')
```

Load Python packages

``` python
import pegasus as pg
import os
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append("/projects/home/ikernin/github_code/myocarditis/functions")
import python_functions
```

Read in single-cell data

``` python
tissue_t = pg.read_input('/projects/home/ikernin/projects/myocarditis/updated_datasets/tissue_t.zarr')
```

``` r
tissue_obs <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_full_obs.csv')
```

## Figure 2A

``` python
python_functions.plot_umap(tissue_t, 'Tissue: T and NK', python_functions.tissue_t_pal, marker_multiplier=6)
```

<img src="figure_2_files/figure-gfm/fig_2a-1.png" width="960" />

## Figure 2B

``` python
python_functions.make_gene_dotplot(tissue_t.to_anndata(),
             cluster_order=['CD8 T 2: CCL5, NKG7',
                            'CD4 T: IL7R, LTB',
                            'CD8 T 1: CD27, LAG3',
                            'CD8 T 4: STMN1, TOP2A',
                            'NK: KLRF1, FCER1G',
                            'CD8 T 3: KLRG1, CX3CR1'
                            ],
             gene_order=['CD8A', 'CCL5', 'NKG7',  # CD8_2 markers
                         'CD4', 'IL7R', 'LTB',  # CD4 markers
                         'CD27', 'LAG3', 'PDCD1',  # CD8_1 markers
                         'STMN1', 'TOP2A',  # CD8_4 markers
                         'KLRF1', 'FCER1G',  # NK markers
                         'KLRG1', 'CX3CR1', 'GZMH'  # CD8_3 markers
                         ],
             title='T/NK')
```

<img src="figure_2_files/figure-gfm/fig_2b-3.png" width="1152" />

## Figure 2C

``` r
# filter
masc_df <- masc_filter(tissue_obs)

# get global cluster numbers
global_number_map <- masc_df %>%
  group_by(lineage_subcluster_name, umap_name) %>%
  summarize(n_cells_clust = n()) %>%
  arrange(desc(n_cells_clust)) %>%
  ungroup() %>%
  mutate(global_subcluster_number = rank(n_cells_clust))

# add global cluster numbers to obs
masc_df <- masc_df %>%
  left_join(global_number_map %>% dplyr::select(lineage_subcluster_name, global_subcluster_number))

## run MASC for subclusters
cluster_masc_res <- MASC(masc_df,
                         cluster = masc_df$global_subcluster_number,
                         contrast = "condition",
                         random_effects = "donor",
                         fixed_effects = "",
                         verbose = TRUE, save_models = FALSE)

## add cluster names to results
cluster_masc_formatted <- cluster_masc_res %>%
  as_tibble() %>%
  mutate(cluster_number = unlist(map(str_split(cluster, 'cluster'), 2)),
         cluster_number = as.numeric(cluster_number)) %>%
  left_join(global_number_map, by = c('cluster_number' = "global_subcluster_number"))

## FDR adjust p-values
cluster_masc_formatted['p.adj'] <- p.adjust(cluster_masc_formatted$model.pvalue, method = 'fdr')
write_csv(cluster_masc_formatted, '/projects/home/ikernin/projects/myocarditis/updated_datasets/masc/cluster_masc_res.csv')

# plot masc results
plot_masc_by_cell_type(cluster_masc_formatted, masc_df, lineage='T and NK')
```

![](figure_2_files/figure-gfm/fig_2c-5.png)<!-- -->

## Figure 2D

``` python
os.chdir('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk')

# get pseudobulk values for DE analysis
python_functions.get_pseudobulk_info(tissue_t, 'tissue_t')
tissue_t.obs['t_cell'] = tissue_t.obs['umap_name'].isin(['CD8 T 1: CD27, LAG3',
                                                         'CD4 T: IL7R, LTB',
                                                         'CD8 T 2: CCL5, NKG7',
                                                         'CD8 T 3: KLRG1, CX3CR1',
                                                         'CD8 T 4: STMN1, TOP2A'])
tissue_t.obs['t_cell'] = tissue_t.obs['t_cell'].replace({True: 'all_t', False: 'other'})
python_functions.get_pseudobulk_info(tissue_t, 'tissue_t_grouped', cluster_col='t_cell')
```

``` r
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/de_analysis')

# run DE analysis by condition
tissue_t_cts <- read_counts('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_t_pseudocounts.csv')
tissue_t_meta <- read_meta('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_t_metainfo.csv')
tissue_t_deres <- run_de_by_comp_var(counts = tissue_t_cts,
                               meta = tissue_t_meta,
                               save_name = 'tissue_t',
                               comp_var_contrast_vec = c('condition', "myocarditis", "control"))

tissue_t_grouped_cts <- read_counts('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_t_grouped_pseudocounts.csv')
tissue_t_grouped_meta <- read_meta('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_t_grouped_metainfo.csv')
tissue_t_grouped_deres <- run_de_by_comp_var(counts = tissue_t_grouped_cts,
                               meta = tissue_t_grouped_meta,
                               save_name = 'tissue_t_grouped',
                               comp_var_contrast_vec = c('condition', "myocarditis", "control"))
```

``` r
tissue_obs <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_full_obs.csv')

# combine de results and meta data for heatmap
t_full_deres <- bind_rows(tissue_t_deres %>%
                               mutate(cluster = as.character(cluster)),
                             tissue_t_grouped_deres)
t_clusters <- tissue_obs %>%
  filter(umap_name == 'T and NK cells') %>%
  dplyr::select(lineage_subcluster_name, lineage_subcluster_number) %>%
  distinct() %>%
  dplyr::rename('cluster_name' = 'lineage_subcluster_name') %>%
  dplyr::rename('cluster_number' = 'lineage_subcluster_number') %>%
  add_row(cluster_name = 'all_t', cluster_number = 'all_t') %>%
    add_row(cluster_name = 'other', cluster_number = 'other') %>%
  add_row(cluster_name = 'Doublets/RBCs', cluster_number = '7')
t_genes <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/de_analysis/t_heatmap_genes.csv') # genes to include in heatmap
t_heatmap_df <- get_heatmap_data(t_full_deres, t_genes, t_clusters)

heatmap_df <- t_heatmap_df %>%
  mutate(cluster = cluster_name,
    cluster = case_when(
    cluster == 'all_t' ~ 'All T',
    TRUE ~ cluster
  )) %>%
  dplyr::select(!cluster_name)

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
  dplyr::select(starts_with("log2FoldChange")) %>%
  replace(is.na(.), 0) %>%
  rename_with(~str_remove(., "log2FoldChange_")) %>%
  dplyr::select(order(colnames(.))) %>%
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
  dplyr::select(starts_with('padj')) %>%
  replace(is.na(.), Inf) %>%
  rename_with(~str_remove(., "padj_")) %>%
  dplyr::select(order(colnames(.))) %>%
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
clust_col_fun <- c('#cf8c00', '#ff4040', '#0097ff', '#00d067', '#bdbdbd', '#8a2be2')
names(clust_col_fun) <- seq(1,6)

lineage_col_fun <- c('white')
names(lineage_col_fun) <- 'All'


clust_ha <- HeatmapAnnotation(clust_colors = names(clust_col_fun),
                              col = list(clust_colors = clust_col_fun),
                              show_legend = FALSE,
                              show_annotation_name = FALSE,
                              simple_anno_size = unit(3, "mm"))
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

![](figure_2_files/figure-gfm/fig_3e_de_heatmap-1.png)<!-- -->

``` r
# filter out doublets
filtered_t_deres <- t_full_deres %>%
        left_join(t_clusters,
                  by = c('cluster' = 'cluster_number')) %>%
        filter(!str_detect(str_to_lower(cluster_name), 'doublets')) %>%
        dplyr::rename('cluster_names' = 'cluster_name') %>%
  mutate(cluster_names =
           case_when(
           str_detect(cluster_names, ':') ~ str_c(cluster, cluster_names, sep = '. '),
           TRUE ~ cluster_names)
  )

# read in gsea pathways
pathways <- gmtPathways("/projects/home/ikernin/projects/myocarditis/updated_datasets/msigdb_symbols.gmt")

# run gsea
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/gsea')
run_gsea_by_cluster(filtered_t_deres, 'tissue_t_gsea')
t_gsea_res <- gsea_combine_xlsx('/projects/home/ikernin/projects/myocarditis/updated_datasets/gsea/tissue_t_gsea_all_gsea.xlsx')

# plot gsea
t_pathways <- c(
        "HALLMARK:INTERFERON_GAMMA_RESPONSE",
        "KEGG:ALLOGRAFT_REJECTION",
        "KEGG:CELL_ADHESION_MOLECULES_CAMS",
        "KEGG:DNA_REPLICATION",
        "KEGG:T_CELL_RECEPTOR_SIGNALING_PATHWAY",
        "KEGG:VIRAL_MYOCARDITIS"
                )
t_cluster_order <- c("all",
                     "3. h-CD4T",
                     "2. h-CD8T",
                     "4. h-CD8T",
                     "5. h-CD8T",
                     "6. h-CD8T",
                     "1. h-NK")
# pre-process data
t_plot_df <- t_gsea_res %>%
  mutate(pathway_name = str_c(geneset, pathway_name, sep=':')) %>%
  filter(pathway_name %in% t_pathways,
         cluster_name %in% t_cluster_order) %>%
  mutate(cluster_name = factor(cluster_name),
         pathway_name = factor(pathway_name)) %>%
  complete(cluster_name, pathway_name)

# make heatmap
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/figures')
plot_heatmap(t_plot_df,
             cluster_order = t_cluster_order,
             col_order = t_pathways,
             't_gsea.pdf',
             split = T)
```

``` r
knitr::include_graphics("/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/t_gsea.pdf")
```

![](../../../../../ikernin/projects/myocarditis/updated_datasets/figures/t_gsea.pdf)<!-- -->

``` r
tissue_troponin_metadata <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_troponin_metadata.csv')

mtx <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global_lineage_pseudocounts.csv",
                row.names = 1)

lin_assign <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/global_lineage_number_to_name_map.csv")

lin_assign$clust <- paste("c", lin_assign$umap_number, sep = "")
# trop_values <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_troponin_metadata.csv")

meta_data <- data.frame(row.names = colnames(mtx))

meta_data$clust <- sapply(rownames(meta_data), function(x) strsplit(x, "_")[[1]][3])
meta_data$donor <- sub("_c10|_c[1-9]", "", rownames(meta_data))

meta_data %<>%
  rownames_to_column() %>%
  dplyr::left_join(tissue_troponin_metadata, by = "donor") %>%
  dplyr::left_join(lin_assign, by = "clust") %>%
  column_to_rownames()


if (!file.exists("/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/data/tissue_lin_model_by_troponin.csv")){
  all_res <- data.frame()
  gset_res <- data.frame()
  for (cl in unique(meta_data$umap_name)){
    pdf(glue("/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/{cl}_results.pdf"))
    meta_temp <- meta_data %>%
      dplyr::filter(umap_name == cl) %>%
      na.omit() %>%
      mutate(log_trop = log(nearest_troponin))

    count_temp <- mtx[,rownames(meta_temp)]

    n_samp <- rowSums(count_temp != 0)
    count_temp <- count_temp[n_samp > round(nrow(meta_temp) / 2),]
    # Okay now we can run DESeq
    dds<- DESeqDataSetFromMatrix(countData = count_temp,
                                colData = meta_temp,
                                design = ~log_trop)
    dds<- DESeq(dds)
    res <- as.data.frame(results(dds))
    res<- res[!is.na(res$padj),]
    res$gene <- rownames(res)
    res$cluster <- cl

    if (nrow(res[res$padj < 0.1,]) > 20 ){
        up_label <- res[res$padj < 0.1,] %>%
          filter(log2FoldChange > 0) %>%
          arrange(pvalue) %>%
          top_n(-20, pvalue) %>%
          .$gene
        down_label <- res[res$padj < 0.1,] %>%
          filter(log2FoldChange < 0) %>%
          arrange(pvalue) %>%
          top_n(-20, pvalue) %>%
          .$gene
        label_genes <- c(up_label, down_label)
      } else if(nrow(res[res$padj < 0.1,]) > 0 ) {
        label_genes <- res[res$padj < 0.1,]$gene
      } else {
        label_genes = c()
      }
    print(
        ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
          geom_point(data = res[res$padj > 0.1,], color = "grey") +
          geom_point(data = res[res$log2FoldChange > 0 & res$padj < 0.1,], color = "red") +
          geom_point(data = res[res$log2FoldChange < 0 & res$padj < 0.1,], color = "blue") +
          geom_text_repel(data = res[res$gene %in% label_genes,], aes(label = gene)) +
          ggtitle("")+
          theme_classic(base_size = 20)
      )

    ## Run GSEA ##
    res2 <- res %>%
      dplyr::select(gene, stat) %>%
      na.omit() %>%
      distinct() %>%
      group_by(gene) %>%
      summarize(stat=mean(stat)) %>%
      dplyr::select(gene, stat) %>%
      na.omit()

    ranks <- deframe(res2)

    gene_sets <- gmtPathways("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/msigdb_symbols.gmt")

    sets <- c("KEGG", "HALLMARK", "BIOCARTA")
    for (s in sets){
      gsets <- gene_sets[grep(s, names(gene_sets))]

      fgseaRes <- fgsea(pathways=gsets, stats=ranks, nperm=1000)
      fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
      fgseaResTidy$cluster <- cl
      gset_res <- rbind(gset_res, fgseaResTidy)

      # Make the GSEA plots
     print(
       ggplot(fgseaResTidy[fgseaResTidy$padj < 0.1,], aes(reorder(pathway, NES), NES)) +
         coord_flip() +
         geom_col() +
         labs(x="Pathway", y="Normalized Enrichment Score",
              title=paste("significant", s,"pathways", sep = " ")) +
         theme_minimal()
          )
      }
    dev.off()
    all_res <- rbind(all_res, res)
  }
  write.csv(all_res, "/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/data/tissue_lin_model_by_troponin.csv",
            row.names = FALSE)
  gset_res$leadingEdge <- as.character(gset_res$leadingEdge)
  write.csv(gset_res, "/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/data/tissue_lin_model_by_troponin_gsea_results.csv",
            row.names = FALSE)

} else {
  gene_sets <- gmtPathways("/projects/home/nealpsmith/projects/kupper/all_data_analysis/data/msigdb_symbols.gmt")
  all_res <- read.csv("/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/data/tissue_lin_model_by_troponin.csv")
  gset_res <- read.csv("/projects/home/nealpsmith/projects/myocarditis/tissue_troponin_gene_modeling/data/tissue_lin_model_by_troponin_gsea_results.csv")
}
```

## Figure 2E

``` r
plot_data <- all_res %>%
  dplyr::filter(cluster == "3. T and NK cells")

label_up <- c("MKI67", "TOP2A", "BIRC5", "LAG3", "HLA-DRA", "HLA-DRB1", "KIR2DL4")
label_down <- c("CX3CR1", "S1PR5", "CCL4", "KLF3", "KLRG1")
label_genes <- c(label_up, label_down)
ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
          geom_point(data = plot_data[plot_data$padj > 0.1,], color = "grey") +
          geom_point(data = plot_data[plot_data$log2FoldChange > 0 & plot_data$padj < 0.1,], pch = 21, fill = "red") +
          geom_point(data = plot_data[plot_data$log2FoldChange < 0 & plot_data$padj < 0.1,], pch = 21, fill = "blue") +
          geom_label_repel(data = plot_data[plot_data$gene %in% label_genes,], aes(label = gene)) +
          ggtitle("")+
          theme_classic(base_size = 20)
```

![](figure_2_files/figure-gfm/fig_2e-1.png)<!-- -->

## Figure 2F

``` r
tcell_data <- all_res %>%
  dplyr::filter(cluster == "3. T and NK cells")
plot_gsets <- c("KEGG_CELL_CYCLE",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                "HALLMARK_G2M_CHECKPOINT",
                "HALLMARK_MTORC1_SIGNALING")

plot_list <- list()
for (gset in plot_gsets){
  ranks <- tcell_data %>%
    dplyr::select(gene, stat) %>%
    na.omit() %>%
    arrange(desc(stat))  %>%
    deframe(.)

  gset_list <- gene_sets[gset]

  fgsea_res <- fgsea(gset_list, stats = ranks, nperm = 10000)

  nes <- round(fgsea_res$NES[fgsea_res$pathway == gset], 3)
  pval <- round(fgsea_res$pval[fgsea_res$pathway == gset], 3)
  if (pval == 0) {
    pval <- "< 0.001"
  }
  n_genes <- fgsea_res$size[fgsea_res$pathway == gset]

  rnk <- rank(-ranks)
  ord <- order(rnk)

  statsAdj <- ranks[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ 1)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(gset_list[[gset]], names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  x=y=NULL

  p <- ggplot(toPlot, aes(x = x, y = y)) +
    # geom_point(color="blue", size=0.1) +
    geom_line(color="blue") +
    geom_hline(yintercept=0, colour="black") +
    geom_segment(data=data.frame(x=pathway),
                     mapping=aes(x=x, y=-0.15,
                                 xend=x, yend=-0.25),
                     size=0.4) +
    scale_y_continuous(limits = c(-0.5, 1), expand = c(0.05,0.05)) +
    xlab("Rank") + ylab("Enrichment score") +
    geom_text(aes(label = "")) +
    annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 1000, y  =0.9) +
    annotate("text", label = glue("p-value : {pval}"), x = length(ranks) - 1000, y = 0.8) +
    annotate("text", label = glue("# genes : {n_genes}"), x = length(ranks) - 1000, y = 0.7) +
    ggtitle(glue("{gset}")) +
    theme_classic(base_size = 12)
  plot_list <- c(plot_list, list(p))

}

plots <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
plots
```

![](figure_2_files/figure-gfm/fig_2f-1.png)<!-- -->