library(tidyverse)
library(glue)
library(openxlsx)
library(biomaRt)
library(fgsea)
library(ComplexHeatmap)
library(circlize)


run_gsea <- function(cluster_res, name, wb_all){
  # get stat per gene_symbol
  cluster_df <- cluster_res %>%
    dplyr::select(gene_symbol, stat) %>%
    distinct() %>%
    na.omit() %>%
    group_by(gene_symbol) %>%
    summarize(stat = mean(stat))
  
  # format ranked dataset for GSEA
  ranks <- deframe(cluster_df)
  
  # look through well defined gene_symbol sets
  plot_list <- list()
  at_least_one <- FALSE
  for (geneset in c("BIOCARTA", "HALLMARK", "KEGG")){
    # print geneset
    print(geneset)
    
    # get pathways in geneset
    gset_pathways <- pathways[grep(geneset, names(pathways))]
    
    # run gsea
    fgseaRes <- fgsea(pathways=gset_pathways, 
                      stats=ranks,
                      minSize=10,
                      nPermSimple=10000) %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    # format pathway names before plotting
    fgseaRes$pathway_name <- sapply(str_split(fgseaRes$pathway, "_", n=2), tail, 1)
    
    # save all non-empty results
    if (nrow(fgseaRes) > 0){
      df_all <- fgseaRes %>% 
        as_tibble() %>%
        relocate(pathway_name) %>%
        dplyr::select(-pathway)
      shortened_name <- unlist(map(str_split(name, ":"), 1))
      addWorksheet(wb_all, glue('{shortened_name}_{geneset}'))
      writeData(wb_all, glue('{shortened_name}_{geneset}'), x=df_all)
    }
  }
}


run_gsea_by_cluster <- function(df, lineage){
  
  # create workbook to save multiple cluster results
  wb_all <- createWorkbook()
  
  # for each cluster run gsea
  for (clust in unique(df$cluster_names)){
    # print cluster
    print(glue('Current cluster: {clust}'))
    
    # filter for cluster
    gsea_df <- df %>%
      filter(cluster_names == clust)
    
    # run for cluster
    run_gsea(gsea_df, clust, wb_all)
  }
  
  saveWorkbook(wb_all, paste0(lineage, "_all_gsea.xlsx"), overwrite = TRUE)
  dev.off()
}

gsea_combine_xlsx <- function(file_path){

  # get sheet names from xlsx file
  sheetnames <- openxlsx::getSheetNames(file_path)

  # store all cluster results in one dataframe
  all_df <- tibble()
  for (sheet in sheetnames){
    temp_df <- read.xlsx(file_path, sheet)
    temp_df['cluster_name'] <- str_split(sheet, '_')[[1]][1]
    temp_df['geneset'] <- str_split_i(sheet, '_', i = -1)
    all_df <- bind_rows(all_df, temp_df)
  }

  return(all_df)
}

plot_heatmap <- function(plot_df, cluster_order, col_order, save_name, split = F){
  
  # format for easy matrix conversion
  heatmap_df <- plot_df %>% 
    dplyr::select(pathway_name, cluster_name, NES, padj) %>%
    pivot_wider(names_from = cluster_name, values_from = c(NES, padj)) 
  
  # get information for the main body's cells
  heatmap_mtx <- heatmap_df %>% 
    dplyr::select(starts_with("NES")) %>%
    rename_with(~str_remove(., "NES_")) %>%
    dplyr::select(cluster_order) %>%
    as.matrix()
  rownames(heatmap_mtx) <- heatmap_df$pathway_name
  
  # define color map
  heatmap_col_fun <- colorRamp2(c(floor(min(heatmap_mtx, na.rm = T)), 
                                  0, 
                                  ceiling(max(heatmap_mtx, na.rm = T))),
                                brewer.pal(5, 'PiYG')[c(5,3,1)])
  
  # get fdr values for main body annotation
  fdr_mtx <- heatmap_df %>% 
    dplyr::select(starts_with('padj')) %>% 
    rename_with(~str_remove(., "padj_")) %>%
    dplyr::select(cluster_order) %>%
    as.matrix()
  rownames(fdr_mtx) <- heatmap_df$pathway_name

  # make function for plotting fdr value
  fdr_func <- function(mtx) {
    function(j, i, x, y, width, height, fill) {
      if (!is.na(mtx[i, j]) & (mtx[i,j] < 0.1)) {
        grid.circle(x = x, y = y, r = unit(1.5, 'mm'),
                    gp = gpar(fill = 'black', col = NA))
      }
    }
  }
  
  # flip rows and columns in matrices
  heatmap_mtx <- t(heatmap_mtx)
  fdr_mtx <- t(fdr_mtx)
  
  # make sure rows and columns the same
  stopifnot(colnames(fdr_mtx) == colnames(heatmap_mtx))
  stopifnot(rownames(fdr_mtx) == rownames(heatmap_mtx))
  
  # make legends
  gsea_lgd <- Legend(col_fun = heatmap_col_fun, title = 'NES', direction = 'horizontal')
  fdr_lgd <- Legend(pch = 16, type = "points", labels = "FDR < 0.1")
  na_lgd <- Legend(labels = 'N/A', legend_gp = gpar(fill = 'grey'))
  pd <- packLegend(gsea_lgd, fdr_lgd, na_lgd, direction = 'horizontal')
  
  
  if (split){
    # split into individual and grouped clusters
    heatmap_mtx_subcluster <- heatmap_mtx[!str_detect(str_to_lower(rownames(heatmap_mtx)), 'all'), ]
    heatmap_mtx_lineage <- heatmap_mtx[str_detect(str_to_lower(rownames(heatmap_mtx)), 'all'), ,drop=FALSE]
    
    fdr_mtx_subcluster <- fdr_mtx[!str_detect(str_to_lower(rownames(fdr_mtx)), 'all'), ]
    fdr_mtx_lineage <- fdr_mtx[str_detect(str_to_lower(rownames(fdr_mtx)), 'all'), ,drop=FALSE]
    
    # plot heatmaps
    ht_subcluster <- Heatmap(heatmap_mtx_subcluster,
                             col = heatmap_col_fun,
                             cell_fun = fdr_func(fdr_mtx_subcluster),
                             column_order = col_order,
                             top_annotation = NULL, column_title = NULL,
                             name = 'NES', show_heatmap_legend = FALSE,
                             cluster_columns = FALSE, column_names_side = "top",
                             show_column_names = T, column_names_rot = 45,
                             cluster_rows = FALSE, row_names_side = "left",
                             row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                             row_gap = unit(2, "mm"), border = TRUE,
                             width = ncol(heatmap_mtx_subcluster)*unit(6, "mm"),
                             height = nrow(heatmap_mtx_subcluster)*unit(6, "mm"))
    ht_lineage <- Heatmap(heatmap_mtx_lineage,
                          col = heatmap_col_fun,
                          cell_fun = fdr_func(fdr_mtx_lineage),
                          column_order = col_order,
                          top_annotation = NULL, column_title = NULL,
                          name = 'NES', show_heatmap_legend = FALSE,
                          cluster_columns = FALSE, column_names_side = "top",
                          show_column_names = T, column_names_rot = 45,
                          cluster_rows = FALSE, row_names_side = "left",
                          row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                          row_gap = unit(2, "mm"), border = TRUE,
                          width = ncol(heatmap_mtx_lineage)*unit(6, "mm"),
                          height = nrow(heatmap_mtx_lineage)*unit(6, "mm"))
    
    pdf(save_name, width = 12, height = 15)
    draw(ht_lineage %v% ht_subcluster,
         heatmap_legend_list = pd,
         heatmap_legend_side = 'bottom',
         merge_legends = TRUE)
    dev.off()
    
  } else {
    
    # plot heatmap
    ht_gsea <- Heatmap(heatmap_mtx,
                       col = heatmap_col_fun,
                       cell_fun = fdr_func(fdr_mtx),
                       column_order = col_order,
                       top_annotation = NULL, column_title = NULL,
                       name = 'NES', show_heatmap_legend = FALSE,
                       cluster_columns = FALSE, column_names_side = "top",
                       show_column_names = T, column_names_rot = 45,
                       cluster_rows = FALSE, row_names_side = "left",
                       row_title_rot = 0, row_title_gp=gpar(fontface='bold'),
                       row_gap = unit(2, "mm"), border = TRUE,
                       width = ncol(heatmap_mtx)*unit(6, "mm"),
                       height = nrow(heatmap_mtx)*unit(6, "mm"))
    
    pdf(save_name, width = 12, height = 15)
    draw(ht_gsea,
         heatmap_legend_list = pd,
         heatmap_legend_side = 'bottom')
    dev.off()
  }
}
