masc_filter <- function(df){
  df %>%
    filter(!str_detect(umap_name, 'Doublets')) %>%
    filter(!str_detect(lineage_subcluster_name, 'Doublets')) %>%
    filter(on_steroids != 'True') %>%
    filter(source != 'CD45+') %>%
    filter(Channel != 'SIC_182_561_1_heart') %>%  # only want first time point of each patient
    mutate(condition = factor(condition, levels = c('control', 'myocarditis')))
}


tissue_odds_ratio_and_box_plot <- function(abundance, masc_res, group, plot_title) {
  
  # arrange results by p-value
  masc_res <- masc_res %>% 
    arrange(desc(model.pvalue)) %>% 
    mutate(group = !!sym(group),
           group = factor(group, levels = group))
  abundance <- abundance %>%
    mutate(group = !!sym(group),
           group = factor(group, levels = levels(masc_res$group)))
  
  
  # format alternating row colors
  if (nrow(masc_res) %% 2 == 0) {
    fill_vec <- rev(rep(c("#d6d6d6", "#ffffff"), length.out = nrow(masc_res)))
  } else {
    fill_vec <- rev(rep(c("#ffffff", "#d6d6d6"), length.out = nrow(masc_res)))
  }
  
  # create odds ratio plot
  orp <- ggplot(masc_res, aes(x = conditionmyocarditis.OR, y = group)) +
    geom_stripes(odd = "#ffffff00", even = "#33333333") +
    ggplot2::geom_errorbarh(aes(xmax = conditionmyocarditis.OR.95pct.ci.upper, height = 0,
                                xmin = conditionmyocarditis.OR.95pct.ci.lower)) +
    geom_point(size = 3) +
    geom_label(data = masc_res,
               mapping = aes(x =  0,
                             hjust =  0, y = group,
                             label = if_else(model.pvalue > 0.001, 
                                             format(round(model.pvalue, 3), nsmall=2),
                                             str_replace(format(signif(model.pvalue, 1), nsmall=1), '-0', '-'))),
               size = 5, color = 'black', alpha = 1,
               fill = fill_vec,
               label.size = 0, label.padding = unit(0.2, 'lines')) +
    geom_vline(xintercept = 1, color = "red") +
    scale_x_log10(labels = function(x) signif(x, 1)) +
    annotation_logticks(side = 'b') +
    theme_forest(base_size = 18) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), legend.position = 'top',
          legend.margin=margin(0,0,0,0), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
    labs(x = "Odds Ratio", y = "", color = "p-value < 0.05")

  
  # create abundance boxplot
  bp <- ggplot(abundance, aes(x = percent, y = group)) +
    geom_stripes(odd = "#ffffff00", even = "#33333333") +
    scale_x_continuous(trans = 'log10', limits = c(1, 100)) +
    geom_boxploth(outlier.shape = NA, aes(fill = condition)) +
    scale_fill_manual(values = c('slategray', 'tomato4')) +
    geom_point(aes(fill = condition,), size = 2, alpha = 0.9, stroke = 0.3,
               position = position_quasirandom(dodge.width = -0.8, groupOnX = F),
               shape = 21) +
    annotation_logticks(side = 'b') +
    theme_forest(base_size = 18) +
    theme(axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
          legend.margin=margin(0,0,0,0), legend.position = 'top',
          legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
    labs(y = '', x = 'Percent + 1', fill = 'Condition')
  
  
  # put plots together and format
  p <- orp +  bp + 
    plot_annotation(title = plot_title,
                    theme = theme(plot.title.position = 'plot',
                                  plot.title = element_text(size = 18, face = 'bold', hjust = 0.5))) +
    plot_layout(widths = unit(c(4, 4), 'in'), heights = unit(length(masc_res$group) * 0.3 + 1, 'in'))
  print(p)
}


plot_masc_by_cell_type <- function(cluster_masc_res, filtered_df, lineage, filter_pdc=T) {
  
  cluster_masc_res <- cluster_masc_res %>% dplyr::rename('cluster_names' = 'lineage_subcluster_name')
  filtered_df <- filtered_df %>% dplyr::rename('cluster_names' = 'lineage_subcluster_name')
  abundance <- filtered_df %>%
    group_by(donor, cluster_names, .drop = FALSE) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    complete(donor, cluster_names, fill = list(count = 0)) %>%
    group_by(donor) %>%
    mutate(percent = count/sum(count) * 100 + 1) %>%  # adding pseudocount
    mutate(log10_percent = log10(percent)) %>%
    filter(is.finite(log10_percent)) %>%
    left_join(filtered_df %>%
                dplyr::select(c(donor, condition)) %>%
                distinct(),
              by='donor')
  
  # pdc CI is very large due to small sample size 
  if (filter_pdc) {
    cluster_masc_res <- cluster_masc_res %>% filter(cluster_names != 'h-pDC: LILRA4 IRF8')
    abundance <- abundance %>% filter(cluster_names != 'h-pDC: LILRA4 IRF8')
  }
  
  # subset to clusters in lineage
  if (lineage == 'Endothelial'){
    clusters <- cluster_masc_res %>%
    filter(umap_name %in% c('Endothelial cells')) %>%
    pull(cluster_names)
  } else if (lineage == 'Non-Endothelial'){
    clusters <- cluster_masc_res %>%
      filter(umap_name %in% c('Mural cells', 'Fibroblasts', 'Neuronal cells', 'Cardiomyocytes')) %>%
      pull(cluster_names)
  }
  else if (lineage == 'MNP'){
    clusters <- cluster_masc_res %>%
      filter(umap_name %in% c('Myeloid cells', 
                               'pDC',
                               'cDC')) %>%
      pull(cluster_names)
  }
  else {
    clusters <- cluster_masc_res %>%
    filter(str_detect(umap_name, lineage)) %>%
    pull(cluster_names)
  }
  
  tissue_odds_ratio_and_box_plot(abundance %>% filter(cluster_names %in% clusters),
                         cluster_masc_res %>% filter(cluster_names %in% clusters),
                         group = 'cluster_names',
                         plot_title = lineage)
}


# by lineage
plot_masc_by_lineage <- function(lineage_masc_res, filtered_df) {
  
  abundance <- filtered_df %>%
    group_by(donor, umap_name, .drop = FALSE) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    complete(donor, umap_name, fill = list(count = 0)) %>%
    group_by(donor) %>%
    mutate(percent = count/sum(count) * 100 + 1) %>%  # adding pseudocount
    mutate(log10_percent = log10(percent)) %>%
    filter(is.finite(log10_percent)) %>%
    left_join(filtered_df %>%
                dplyr::select(c(donor, condition)) %>%
                distinct(),
              by='donor')
  
  tissue_odds_ratio_and_box_plot(abundance %>% filter(umap_name != 'pDC'), 
                          lineage_masc_res %>% filter(umap_name != 'pDC'), 
                          group = 'umap_name', plot_title = 'Tissue Global')
}
