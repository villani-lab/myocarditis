plot_clust_perc_by_donor <- function(df, lineage, cluster_order){
  # set umap color palettes for each lineage
  tissue_global_pal = c('#ff0029', '#377eb8', '#66a61e', '#984ea3', '#00d2d5',
                        '#af8d00', '#7f80cd', '#b3e900', '#c42e60', '#ff7f00')
  tissue_mnp_pal = c('#1b9e77', '#e7298a', '#a6761d', '#252525', '#f43600', '#356F83', '#eff26e', '#80b1d3')
  tissue_t_pal = c('#cf8c00', '#ff4040', '#0097ff', '#00d067', '#bdbdbd', '#8a2be2', '#80b1d3')
  tissue_nonimmune_pal = c('#a65628', '#f781bf', '#8dd3c7', '#bebada', '#fb8072', '#fdb462',
                           '#fccde5', '#bc80bd', '#ffed6f', '#c4eaff', '#d95f02', '#737373',
                           '#4ba93b', '#5779bb', '#927acc', '#bf3947', '#f48758', '#80b1d3')
  tissue_b_pal = c('#ff7373', '#f2b94f', '#80b1d3')
  tissue_palettes <- list('mnp' = tissue_mnp_pal,
                          't' = tissue_t_pal,
                          'b' = tissue_b_pal,
                          'nonimmune' = tissue_nonimmune_pal)
  
  # order donors by name and condition
  donor_df <- tibble(donor = c("SIC_48", "SIC_153","SIC_164","SIC_171",
                               "SIC_197", "SIC_199","SIC_217", "SIC_232",
                               "SIC_258","SIC_264", "SIC_317", "SIC_319",
                               "SIC_182","SIC_333",
                               "Sanger_3", "Sanger_4", "Sanger_5",
                               "Sanger_6", "Sanger_7", "Sanger_11"),
                     case = c(rep('Myocarditis', 12), rep('Control', 8))) %>%
    mutate(donor = factor(donor, levels = unique(donor)))
  
  # plot stacked bar of cluster percents by donor
  df %>%
    left_join(donor_df, by = 'donor') %>%
    mutate(donor = factor(donor, levels = levels(donor_df$donor)),
           case = factor(case, levels = c('Myocarditis', 'Control')),
           cluster_names = factor(cluster_names, levels = cluster_order)) %>%
    arrange(cluster_names) %>%
    ggplot(aes(x = donor, y = perc, group = cluster_names, fill = cluster_names)) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = tissue_palettes[lineage][[1]]) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    facet_grid(~case, scales = 'free_x', space = 'free') +
    labs(x = 'Donor', 
         y = 'Percentage of Donor Cells', 
         fill = 'Cluster',
         title = glue('{str_to_title(lineage)}: % Donor Cells per Cluster'))
}

