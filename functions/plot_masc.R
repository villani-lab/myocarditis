masc_filter <- function(df){
  df %>%
    filter(!str_detect(umap_name, 'Doublets')) %>%
    filter(on_steroids != 'True') %>%
    filter(source != 'CD45+') %>%
    filter(Channel != 'SIC_182_561_1_heart') %>%  # only want first time point of each patient
    mutate(condition = factor(condition, levels = c('control', 'myocarditis')))
}


odds_ratio_and_box_plot <- function(abundance, masc_res, comp_var, group, plot_title, colors, row_order=F) {

  # arrange results by p-value
  if (!row_order) {
    masc_res <- masc_res %>%
      arrange(desc(model.pvalue))
  }
  # rows already arranged
  masc_res <- masc_res %>%
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

  # Get the OR variable
  or_var <- paste(comp_var, levels(abundance[[comp_var]])[2], ".OR", sep = "")

  # create odds ratio plot
  orp <- ggplot(masc_res, aes(x = !!sym(or_var), y = group)) +
    geom_stripes(odd = "#ffffff00", even = "#33333333") +
    ggplot2::geom_errorbarh(aes(xmax = !!sym(glue("{or_var}.95pct.ci.upper")), height = 0,
                                xmin = !!sym(glue("{or_var}.95pct.ci.lower")))) +
    geom_point(size = 3) +
    geom_label(data = masc_res,
               mapping = aes(x =  0,
                             hjust =  0, y = group,
                             label = signif(model.pvalue, 3)),
               size = 5, color = 'black', alpha = 1,
               fill = fill_vec,
               label.size = 0, label.padding = unit(0.2, 'lines')) +
    geom_vline(xintercept = 1, color = "red") +
    scale_x_log10(labels = function(x) signif(x, 1)) +
    annotation_logticks(side = 'b') +
    theme_forest(base_size = 18) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), legend.position = 'top',
          legend.margin=margin(0,0,0,0), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
    labs(x = "Odds Ratio", y = "", color = "p-value < 0.05")

  # create abundance boxplot
  bp <- ggplot(abundance, aes(x = percent, y = group)) +
    geom_stripes(odd = "#ffffff00", even = "#33333333") +
    scale_x_continuous(trans = 'log10', limits = c(1, 100)) +
    geom_boxploth(outlier.shape = NA, aes(fill = !!as.symbol(comp_var))) +
    scale_fill_manual(values = colors) +
    geom_point(aes(fill = !!as.symbol(comp_var),), size = 2, alpha = 0.9, stroke = 0.3,
               position = position_quasirandom(dodge.width = -0.8, groupOnX = F),
               shape = 21) +
    annotation_logticks(side = 'b') +
    theme_forest(base_size = 18) +
    theme(axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          legend.margin=margin(0,0,0,0), legend.position = 'top',
          legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
    labs(y = '', x = 'Percent + 1', fill = '')

  # put plots together and format
  p <- orp +
    bp +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title.position = 'plot',
                                  plot.title = element_text(size = 18, face = 'bold', hjust = 0.5))) +
    plot_layout(widths = unit(c(4, 4), 'in'), heights = unit(length(masc_res$group) * 0.3 + 1, 'in'))
  print(p)

}


plot_masc_by_cell_type <- function(cluster_masc_res, filtered_df, lineage, comp_var = "condition", filter_pdc=T) {

  abundance <- filtered_df %>%
    group_by(donor, global_subcluster_name, .drop = FALSE) %>%
    summarize(count = n()) %>%
    group_by(donor) %>%
    mutate(percent = count/sum(count) * 100 + 1) %>%  # adding pseudocount
    mutate(log10_percent = log10(percent)) %>%
    filter(is.finite(log10_percent)) %>%
    left_join(filtered_df %>%
                dplyr::select(c(donor, !!sym(comp_var))) %>%
                distinct(),
              by='donor') %>%
    mutate(cluster_names = unlist(map(str_split(global_subcluster_name, '\\. '), 2)))

  # pdc CI is very large due to small sample size
  if (filter_pdc) {
    cluster_masc_res <- cluster_masc_res %>% filter(cluster_names != 'h-pDC: LILRA4 IRF8')
    abundance <- abundance %>% filter(cluster_names != 'h-pDC: LILRA4 IRF8')
  }

  # subset to clusters in lineage
  if (lineage == 'Nonimmune'){
    clusters <- cluster_masc_res %>%
    filter(!umap_name %in% c('T and NK cells', 'Myeloid cells', 'B and plasma cells', 'Dendritic cells')) %>%
    pull(cluster_names)
  } else {
    clusters <- cluster_masc_res %>%
    filter(str_detect(umap_name, lineage)) %>%
    pull(cluster_names)
  }

  odds_ratio_and_box_plot(abundance %>% filter(cluster_names %in% clusters),
                          cluster_masc_res %>% filter(cluster_names %in% clusters),
                          group = 'cluster_names',
                          plot_title = lineage,
                          comp_var = comp_var, colors = c('slategray', 'tomato4'))

}


# by lineage
plot_masc_by_lineage <- function(lineage_masc_res, filtered_df) {
  
  abundance <- filtered_df %>%
    group_by(donor, umap_name, .drop = FALSE) %>%
    summarize(count = n()) %>%
    group_by(donor) %>%
    mutate(percent = count/sum(count) * 100 + 1) %>%  # adding pseudocount
    mutate(log10_percent = log10(percent)) %>%
    filter(is.finite(log10_percent)) %>%
    left_join(filtered_df %>%
                dplyr::select(c(donor, condition)) %>%
                distinct(),
              by='donor') %>%
    mutate(lineage_names = unlist(map(str_split(umap_name, '\\. '), 2)))
  
  odds_ratio_and_box_plot(abundance, lineage_masc_res, group = 'lineage_names', plot_title = 'Tissue Global')
  
}