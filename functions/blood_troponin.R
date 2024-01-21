troponin_get_percent_per_level <- function(df, level = 'cluster'){
  # get donor troponin info
  troponin_donor <- df %>%
    select(c(donor, nearest_troponin)) %>%
    distinct()

  # get total number of donor cells
  ncells_donor <- df %>%
    group_by(donor) %>%
    summarize(ncells_donor = n())

  # get total number of donor,level cells
  ncells_level <- df %>%
    group_by(donor, !!sym(glue('{level}_names'))) %>%
    summarize(!!sym(glue('ncells_{level}')) := n()) %>%
    ungroup() %>%
    complete(donor, !!sym(glue('{level}_names')),
             fill = list2(!!sym(glue("ncells_{level}")) := 0))

  # put data together and get level percentages
  model_df <- ncells_level %>%
    left_join(ncells_donor) %>%
    left_join(troponin_donor) %>%
    mutate(!!sym(glue('{level}_prop_by_all')) := !!sym(glue('ncells_{level}'))/ncells_donor)

  return(model_df)
}


troponin_fit_model <- function(df, level = 'cluster'){
  # define model to fit
  fit_lm_by_level <- function(df){
    f <- paste(glue("log({level}_prop_by_all + 1)"), " ~ log(nearest_troponin)")
    lm(f, data=df)
  }

  # fit model
  model_df <- df %>%
    group_by_at(glue("{level}_names")) %>%
    nest() %>%
    mutate(model = map(data, fit_lm_by_level)) %>%
    mutate(trop_coef = map(model, function(model){summary(model)$coefficients[2,1]}),
           trop_se = map(model, function(model){summary(model)$coefficients[2,2]}),
           trop_pval = map(model, function(model){summary(model)$coefficients[2,4]}))

  # adjust pvalues using FDR
  pvals <- p.adjust(model_df$trop_pval, method='fdr')
  model_df$padj <- pvals

  return(model_df)
}


set_factor_order <- function(df, col_name, order){
  df %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), levels = order))
}


troponin_plot_model <- function(model_res, model_data, title,
                       level='cluster', type='simple', point_size=2){

  # format columns for printing
  model_annot <- model_res %>%
    select(c(!!sym(glue("{level}_names")), trop_coef, trop_se, trop_pval, padj)) %>%
    mutate(
      trop_coef = round(as.numeric(trop_coef), 2),
      trop_se = round(as.numeric(trop_se), 2),
      trop_pval = signif(as.numeric(trop_pval), 2),
      padj = signif(as.numeric(padj), 2),
      trop_label = paste0("Coef: ", trop_coef, " +/- ", trop_se),
      pval_label = paste0("pval = ", trop_pval, "; padj = ", padj)
    )

  # plot
  if (type == 'simple'){
    p <- model_data %>%
      ggplot(aes(x = nearest_troponin,
                 y = !!sym(glue("{level}_prop_by_all")))) +
      geom_point(size=point_size) +
      geom_smooth(method = lm, se = FALSE, color='black', linetype="dashed") +
      theme_bw() +
      theme(text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_log10() +
      facet_wrap(as.formula(glue("~ {level}_names")),
                 scales = "free_y") +
      labs(x = "Serum troponin T (ng/L)",
           y = paste("Proportion of all cells"),
           title = title,
           subtitle = glue("Log({str_to_title(level)} Proportion + 1) ~ Log(Serum troponin T (ng/L))"))
  } else if (type == 'detailed'){
    p <- model_data %>%
      ggplot(aes(x = nearest_troponin,
                 y = !!sym(glue("{level}_prop_by_all")))) +
      geom_point(size=point_size) +
      geom_text(data = model_annot,aes(label = trop_label),
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -2 , size=4.5) +
      geom_text(data = model_annot,aes(label = pval_label),
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -1, size=4.5) +
      geom_smooth(method = lm, se = FALSE) +
      theme_bw() +
      theme(text = element_text(size = 16)) +
      scale_x_log10() +
      facet_wrap(as.formula(glue("~ {level}_names")),
                 scales = "free_y") +
      labs(x = "Serum troponin T (ng/L)",
           y = paste("Proportion of all cells"),
           title = title,
           subtitle = glue("Log({str_to_title(level)} Proportion + 1) ~ Log(Serum troponin T (ng/L))"))
  }

  print(p)
}