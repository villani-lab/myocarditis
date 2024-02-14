Figure 6
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
tissue_fibroblast = pg.read_input('/projects/home/ikernin/projects/myocarditis/updated_datasets/tissue_fibroblast.zarr')
```

``` r
tissue_obs <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_full_obs.csv')
```

    ## Rows: 84576 Columns: 19
    ## ── Column specification ──────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (13): barcodekey, Channel, fatal, on_steroids, condition, sex, donor, so...
    ## dbl  (6): n_genes, n_counts, percent_mito, scale, leiden_labels, umap_numbers
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

## Figure 6a

``` python
python_functions.plot_umap(tissue_fibroblast, 'Tissue: Fibroblast', python_functions.tissue_fibroblast_pal, marker_multiplier=6)
```

<img src="figure_6_files/figure-gfm/fig_6a-1.png" width="960" />

## Figure 6b

``` python
python_functions.make_gene_dotplot(tissue_fibroblast.to_anndata(),
             cluster_order=[
                 'Fibroblast: CXCL9, HLA-DRA',
                 'Fibroblast: DCN, EGR1low',
                 'Myofibroblast: ACTA2, ID4',
                 'Fibroblast: GPX3, EGR1high',
                 'Fibroblast: POSTN, F2R',
                 'Fibroblast: PCOLCE2, IGFBP6',

             ],
             gene_order=[
                 "CXCL9", "HLA-DRA", "CCL5",
                 "DCN",
                 "ACTA2", "ID4",
                 "GPX3", "EGR1",
                 "POSTN", "F2R",
                 "PCOLCE2", "IGFBP6"
             ],
             title='Tissue: Fibroblast')
```

<img src="figure_6_files/figure-gfm/fig_6b-3.png" width="1152" />

## Figure 6C

``` r
# filter
masc_df <- masc_filter(tissue_obs)

# read in masc res (from figure 3 code)
cluster_masc_res <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/masc/cluster_masc_res.csv')

# plot masc results
plot_masc_by_cell_type(cluster_masc_res, masc_df, lineage='Fibroblast')
```

![](figure_6_files/figure-gfm/fig_6c-5.png)<!-- -->

## Figure 6D

``` python
os.chdir('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk')

# get pseudobulk counts and metadata by donor
python_functions.get_pseudobulk_info(tissue_fibroblast, 'tissue_fibroblast')
tissue_fibroblast.obs['fibroblast_cell'] = tissue_fibroblast.obs['umap_name'].isin(['Fibroblast: CXCL9, HLA-DRA',
                                                                                    'Fibroblast: DCN, EGR1low',
                                                                                    'Fibroblast: GPX3, EGR1high',
                                                                                    'Fibroblast: PCOLCE2, IGFBP6',
                                                                                    'Fibroblast: POSTN, F2R'])
tissue_fibroblast.obs['fibroblast_cell'] = tissue_fibroblast.obs['fibroblast_cell'].replace({True: 'all_fibroblast', False: 'other'})
python_functions.get_pseudobulk_info(tissue_fibroblast, 'tissue_fibroblast_grouped', cluster_col='fibroblast_cell')
```

``` r
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/de_analysis')

# run DE analysis by condition
tissue_fibroblast_cts <- read_counts('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_fibroblast_pseudocounts.csv')
tissue_fibroblast_meta <- read_meta('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_fibroblast_metainfo.csv')
tissue_fibroblast_deres <- run_de_by_comp_var(counts = tissue_fibroblast_cts,
                               meta = tissue_fibroblast_meta,
                               save_name = 'tissue_fibroblast',
                               comp_var_contrast_vec = c('condition', "myocarditis", "control"))

tissue_fibroblast_grouped_cts <- read_counts('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_fibroblast_grouped_pseudocounts.csv')
tissue_fibroblast_grouped_meta <- read_meta('/projects/home/ikernin/projects/myocarditis/updated_datasets/pseudobulk/tissue_fibroblast_grouped_metainfo.csv')
tissue_fibroblast_grouped_deres <- run_de_by_comp_var(counts = tissue_fibroblast_grouped_cts,
                               meta = tissue_fibroblast_grouped_meta,
                               save_name = 'tissue_fibroblast_grouped',
                               comp_var_contrast_vec = c('condition', "myocarditis", "control"))
```

``` r
# combine de results and meta data for heatmap
fibroblast_full_deres <- bind_rows(tissue_fibroblast_deres %>%
                               mutate(cluster = as.character(cluster)),
                             tissue_fibroblast_grouped_deres)
fibroblast_clusters <- tissue_obs %>%
  filter(umap_name == 'Fibroblasts') %>%
  dplyr::select(lineage_subcluster_name, lineage_subcluster_number) %>%
  distinct() %>%
  dplyr::rename('cluster_name' = 'lineage_subcluster_name') %>%
  dplyr::rename('cluster_number' = 'lineage_subcluster_number') %>%
  add_row(cluster_name = 'all_fibroblast', cluster_number = 'all_fibroblast') %>%
    add_row(cluster_name = 'other', cluster_number = 'other') %>%
  add_row(cluster_name = 'Doublets/RBCs', cluster_number = '7') %>%
  mutate(cluster_name =
           case_when(
           str_detect(cluster_name, ':') ~ str_c(cluster_number, cluster_name, sep = '. '),
           TRUE ~ cluster_name)
  )
fibroblast_genes <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/de_analysis/fibroblast_genes.csv') # genes to include in heatmap
fibroblast_heatmap_df <- get_heatmap_data(fibroblast_full_deres, fibroblast_genes, fibroblast_clusters)


heatmap_df <- fibroblast_heatmap_df %>%
  mutate(cluster = cluster_name,
    cluster = case_when(
    cluster == 'all_fibroblast' ~ 'All Fibroblast',
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
clust_col_fun <- c('#A85529', '#F67FBD', '#8DD3C8', '#BEBBD9', '#FC7F72', '#80B1D2')
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

![](figure_6_files/figure-gfm/fig_5f_heatmap-1.png)<!-- -->

``` r
# filter out doublets
filtered_fibroblast_deres <- fibroblast_full_deres %>%
        left_join(fibroblast_clusters,
                  by = c('cluster' = 'cluster_number')) %>%
        filter(!str_detect(str_to_lower(cluster_name), 'doublets')) %>%
        dplyr::rename('cluster_names' = 'cluster_name')

# read in gsea pathways
pathways <- gmtPathways("/projects/home/ikernin/projects/myocarditis/updated_datasets/msigdb_symbols.gmt")

# run gsea
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/gsea')
run_gsea_by_cluster(filtered_fibroblast_deres, 'tissue_fibroblast')
fibroblast_gsea_res <- gsea_combine_xlsx('/projects/home/ikernin/projects/myocarditis/updated_datasets/gsea/tissue_fibroblast_all_gsea.xlsx')

# plot gsea
fibroblast_pathways <- c("HALLMARK:INTERFERON_GAMMA_RESPONSE",
                         "KEGG:ALLOGRAFT_REJECTION",
                         "KEGG:ANTIGEN_PROCESSING_AND_PRESENTATION",
                         "KEGG:CELL_ADHESION_MOLECULES_CAMS",
                         "KEGG:DILATED_CARDIOMYOPATHY",
                         "KEGG:DNA_REPLICATION",
                         "KEGG:HYPERTROPHIC_CARDIOMYOPATHY_HCM",
                         "KEGG:VIRAL_MYOCARDITIS")

fibroblast_order <- c('all',
                      '1. Fibroblast',
                      '2. Myofibroblast',
                      '3. Fibroblast',
                      '4. Fibroblast',
                      '5. Fibroblast',
                      '6. Fibroblast')

fibroblast_plot_df <- fibroblast_gsea_res %>%
  mutate(pathway_name = str_c(geneset, pathway_name, sep=':')) %>%
  filter(pathway_name %in% fibroblast_pathways,
         cluster_name %in% fibroblast_order) %>%
  mutate(cluster_name = factor(cluster_name),
         pathway_name = factor(pathway_name)) %>%
  complete(cluster_name, pathway_name)


# make heatmap
setwd('/projects/home/ikernin/projects/myocarditis/updated_datasets/figures')
plot_heatmap(fibroblast_plot_df,
             cluster_order = fibroblast_order,
             col_order = fibroblast_pathways,
             'fibroblast_gsea.pdf',
             split = T)
```

``` r
knitr::include_graphics("/projects/home/ikernin/projects/myocarditis/updated_datasets/figures/fibroblast_gsea.pdf")
```

<embed src="../../../../projects/myocarditis/updated_datasets/figures/fibroblast_gsea.pdf" width="800px" height="460px" type="application/pdf" />

## Figure 6E

``` r
mtx <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_global_lineage_pseudocounts.csv",
                row.names = 1)

lin_assign <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/global_lineage_number_to_name_map.csv")

lin_assign$clust <- paste("c", lin_assign$umap_number, sep = "")
trop_values <- read.csv("/projects/home/ikernin/projects/myocarditis/github_datasets/tissue_troponin_metadata.csv")

meta_data <- data.frame(row.names = colnames(mtx))

meta_data$clust <- sapply(rownames(meta_data), function(x) strsplit(x, "_")[[1]][3])
meta_data$donor <- sub("_c10|_c[1-9]", "", rownames(meta_data))

meta_data %<>%
  rownames_to_column() %>%
  dplyr::left_join(trop_values, by = "donor") %>%
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

plot_gsets <- c("KEGG_CHEMOKINE_SIGNALING_PATHWAY", "KEGG_ALLOGRAFT_REJECTION",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE", "KEGG_VIRAL_MYOCARDITIS")

fib_data <- all_res %>%
  dplyr::filter(cluster == "5. Fibroblasts")

plot_list <- list()
for (gset in plot_gsets){
  print(gset)
  ranks <- fib_data %>%
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
    scale_y_continuous(expand = c(0.05,0.05)) +
    xlab("Rank") + ylab("Enrichment score") +
    geom_text(aes(label = "")) +
    annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 1000, y  =0.9) +
    annotate("text", label = glue("p-value : {pval}"), x = length(ranks) - 1000, y = 0.8) +
    annotate("text", label = glue("# genes : {n_genes}"), x = length(ranks) - 1000, y = 0.7) +
    ggtitle(glue("{gset}")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size=8))
  plot_list <- c(plot_list, list(p))

}

plots <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
plots
```

![](figure_6_files/figure-gfm/fig_6e-1.png)<!-- -->

    ## [1] "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
    ## [1] "KEGG_ALLOGRAFT_REJECTION"
    ## [1] "HALLMARK_INTERFERON_GAMMA_RESPONSE"
    ## [1] "KEGG_VIRAL_MYOCARDITIS"

## Figure 6H

``` r
tissue_troponin_metadata <- read_csv('/projects/home/ikernin/projects/myocarditis/updated_datasets/metadata/tissue_troponin_metadata.csv')
troponin_filtered_df <- troponin_filter_tissue(tissue_obs, tissue_troponin_metadata)

# fit linear model by troponin for DE clusters
select_clusters <- c("h-NK: KLRF1 FCER1G",
                     "h-CD4T: IL7R LTB",
                    "h-CD8T: CD27 LAG3",
                    "h-CD8T: CCL5 NKG7",
                    "h-CD8T: cycling",
                    "h-MNP: S100A8-low C1QA-low",
                    "h-MNP: FCGR3A LILRB2",
                    "h-cDC: CLEC9A CD1C",
                    "Fibroblast: CXCL9, HLA-DRA")
troponin_cluster_percs <- troponin_get_percents_per_level(troponin_filtered_df, level='cluster')
select_cluster_percs <- troponin_cluster_percs %>%
        filter(cluster_names %in% select_clusters)
select_cluster_model <- troponin_fit_model(select_cluster_percs, level='cluster')
kable(select_cluster_model %>%
              dplyr::select(!c(data, model)) %>%
              unnest(cols = c(trop_coef, trop_se, trop_pval)))

troponin_plot_model(select_cluster_model %>% filter(cluster_names =="Fibroblast: CXCL9, HLA-DRA"),
                    select_cluster_percs %>% filter(cluster_names =="Fibroblast: CXCL9, HLA-DRA"),
                   "Fibroblast: CXCL9, HLA-DRA", level='cluster', point_size = 2.2, type='simple')
```

![](figure_6_files/figure-gfm/fig_6h-1.png)<!-- -->

| cluster\_names             |  trop\_coef |  trop\_se | trop\_pval |      padj |
| :------------------------- | ----------: | --------: | ---------: | --------: |
| Fibroblast: CXCL9, HLA-DRA |   0.0079106 | 0.0024072 |  0.0082023 | 0.0532563 |
| h-CD4T: IL7R LTB           |   0.0001508 | 0.0061707 |  0.9809880 | 0.9809880 |
| h-CD8T: CCL5 NKG7          |   0.0041760 | 0.0104397 |  0.6975596 | 0.8932075 |
| h-CD8T: CD27 LAG3          |   0.0049871 | 0.0063537 |  0.4507020 | 0.6760529 |
| h-CD8T: cycling            |   0.0059728 | 0.0021203 |  0.0182542 | 0.0547626 |
| h-cDC: CLEC9A CD1C         |   0.0019781 | 0.0006443 |  0.0118347 | 0.0532563 |
| h-MNP: FCGR3A LILRB2       |   0.0022458 | 0.0083721 |  0.7939622 | 0.8932075 |
| h-MNP: S100A8-low C1QA-low |   0.0115797 | 0.0083111 |  0.1937265 | 0.3487078 |
| h-NK: KLRF1 FCER1G         | \-0.0054702 | 0.0031077 |  0.1088655 | 0.2449473 |