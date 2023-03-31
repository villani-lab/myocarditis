condition_filter_df <- function(df){
  # lood at CD45 enriched channels and filter out doublets
  df %>%
    select(donor, Channel, timepoint_cat, umap_name, global_subcluster_name) %>%
    filter(!str_detect(umap_name, 'Doublets')) %>%
    filter(str_detect(Channel, 'CD45')) %>%
    filter(timepoint_cat == 'on_ici' | timepoint_cat == 'pre_steroid') %>%
    mutate(lineage_names = unlist(map(str_split(umap_name, '\\. '), 2)),
           cluster_names = unlist(map(str_split(global_subcluster_name, '\\. '), 2)),
           condition = str_replace(timepoint_cat, 'pre_steroid', 'myocarditis'),
           condition = factor(condition, levels = c('on_ici', 'myocarditis'))) %>%
    select(donor, condition, lineage_names, cluster_names)
}


condition_get_percent_per_level <- function(df, level = 'cluster'){
  
  # get total number of donor cells
  ncells_donor <- df %>%
    group_by(donor, condition) %>%
    summarize(ncells_donor_timepoint = n())
  
  # get total number of donor,level cells per timepoint
  onici_timepoints <- df %>%
    filter(condition == 'on_ici') %>%
    group_by(donor, !!sym(glue('{level}_names'))) %>%
    summarize(!!sym(glue('ncells_{level}')) := n()) %>%
    ungroup() %>%
    complete(donor, !!sym(glue('{level}_names')),
             fill = list2(!!sym(glue("ncells_{level}")) := 0)) %>%
    mutate(condition = 'on_ici')
  
  myocarditis_timepoints <- df %>%
    filter(condition == 'myocarditis') %>%
    group_by(donor, !!sym(glue('{level}_names'))) %>%
    summarize(!!sym(glue('ncells_{level}')) := n()) %>%
    ungroup() %>%
    complete(donor, !!sym(glue('{level}_names')),
             fill = list2(!!sym(glue("ncells_{level}")) := 0)) %>%
    mutate(condition = 'myocarditis')
  
  ncells_level <- bind_rows(onici_timepoints, myocarditis_timepoints) %>%
    arrange(donor, condition)
  
  # put data together and get level percentages
  model_df <- ncells_level %>%
    left_join(ncells_donor) %>%
    mutate(!!sym(glue('{level}_prop_by_all')) := !!sym(glue('ncells_{level}'))/ncells_donor_timepoint,
           perc = !!sym(glue('{level}_prop_by_all')) * 100,
           log2_percent = log2(perc + 1)) 
  
  return(model_df)
}


set_factor_order <- function(df, col_name, order){
  df %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), levels = order))
}


# plot the sample percentage of cells in each cluster
condition_plot_sample_perc <- function(perc_df, title, level = 'cluster'){
  
  perc_df['names'] <- perc_df[glue('{level}_names')]
  perc_df %>%
    mutate(condition = factor(condition, levels = c('myocarditis', 'on_ici')),
           names = factor(names, levels = rev(levels(names)))) %>%
    ggplot(aes(x = names, 
               y = perc + 1)) +
    geom_boxplot(aes(fill = condition),
                 outlier.shape = NA) +
    geom_point(aes(fill = condition),
               pch=21,
               size = 2.5,
               position = position_jitterdodge(jitter.width = 0.1)) +
    geom_hline(yintercept = c(1,10,100),
               linetype='dotted',
               lwd= 0.8,
               alpha = 0.5) +
    scale_fill_manual(values = c("#8A3625", "#6F808F")) +
    scale_color_manual(values = c("#8A3625", "#6F808F")) +
    scale_y_log10() +
    annotation_logticks(side = "b",
                        size=1) +
    coord_flip(clip = "off") +
    annotate(
      geom = "rect",
      ymin = 0, ymax = 100,
      xmin = as.integer(seq(unique(perc_df$names))) - 0.5,
      xmax = as.integer(seq(unique(perc_df$names))) + 0.5,
      fill = rep(
        c(NA, "black"),
        length.out = length(unique(perc_df$names))
      ),
      alpha = 0.1
    ) +
    labs(x = element_blank(),
         y = 'Percent + 1',
         title = paste0(title, ': Percentages per Sample'),
         fill = 'Condition',
         color = 'Condition') +
    theme_classic(base_size = 20) +
    theme(legend.position = 'top',
          legend.direction = 'horizontal',
          legend.justification = c(-0.02, 1),
          legend.title = element_text(size=18),
          legend.text = element_text(size=18)) +
    guides(fill = guide_legend(reverse=TRUE))
}


condition_fit_model <- function(df, level = 'cluster'){
  # define model to fit
  fit_lm_by_level <- function(df){
    lm(log2_percent ~ condition, data=df)
  }
  
  # fit model
  model_df <- df %>%
    mutate(condition = factor(condition, levels = c('on_ici', 'myocarditis'))) %>%
    group_by_at(glue("{level}_names")) %>%
    nest() %>%
    mutate(model = map(data, fit_lm_by_level)) %>%
    mutate(condition_coef = map(model, function(model){summary(model)$coefficients[2,1]}),
           condition_se = map(model, function(model){summary(model)$coefficients[2,2]}),
           condition_pval = map(model, function(model){summary(model)$coefficients[2,4]}))
  
  # adjust pvalues using FDR
  pvals <- p.adjust(model_df$condition_pval, method='fdr')
  model_df$padj <- pvals
  
  # add CI info
  get_ci_from_model <- function(df){
    parameters(df) %>% 
      filter(Parameter != "(Intercept)") %>%
      select(CI_low, CI_high)
  }
  
  model_df <- model_df %>%
    mutate(CIs = map(model, get_ci_from_model)) %>%
    unnest(CIs)
  
  return(model_df)
}


condition_plot_ci_interval <- function(model_df, title, level='cluster'){
  
  model_df['names'] <- model_df[glue('{level}_names')]  
  model_df %>%
    mutate(names = factor(names, levels = rev(levels(names)))) %>%
    ggplot(aes(x = as.numeric(condition_coef), y = names)) +
    geom_stripes(odd = "#ffffff00", even = "#33333333") +
    ggplot2::geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, height = 0)) +
    geom_point(size = 3) +
    geom_label(aes(x = min(CI_low) - 0.25 * abs(min(CI_low)),
                   y = names,
                   hjust = 0,
                   label = if_else(condition_pval > 0.001, 
                                   format(round(as.numeric(condition_pval), 3), nsmall=2),
                                   str_replace(format(signif(as.numeric(condition_pval), 1), nsmall=1), '-0', '-'))
    ),
    size = 5,
    color = 'black',
    label.size = 0, 
    label.padding = unit(0.2, 'lines')) +
    geom_vline(xintercept = 0, color = "red") +
    annotation_logticks(side = 'b') +
    xlim(min(model_df$CI_low) -  0.25 * abs(min(model_df$CI_low)), NA) +
    theme_forest(base_size = 18) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
    labs(x = "Log2FC",
         y = "", 
         color = "p-value < 0.05",
         title = paste0(title, ': CI of Log2FC'))
}
