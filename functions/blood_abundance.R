create_abundance <- function(obs, comp_var, group_var) {
  abundance <- obs %>%
    filter(!!as.symbol(comp_var) != "NA") %>%
    filter(!str_detect(lineage, "Doublet")) %>%
    group_by(sample_id, !!as.symbol(group_var)) %>%
    summarize(n_cells = n(), .groups = 'drop') %>%
    complete(sample_id, !!as.symbol(group_var), fill = list(n_cells = 0)) %>%
    group_by(sample_id) %>%
    mutate(percent = n_cells/sum(n_cells) * 100 + 1) %>%  # adding pseudocount
    mutate(log10_percent = log10(percent)) %>%
    filter(is.finite(log10_percent)) %>%
    left_join(obs %>%
                select(c(sample_id, !!as.symbol(comp_var))) %>%
                distinct(),
              by='sample_id') %>%
    mutate(cluster = !!as.symbol(group_var))
  return(abundance)
}

run_masc <- function(obs, comp_var_contrast_vec, group_var, fixed_effects="") {
  comp_var <- comp_var_contrast_vec[1]
  second_level <- comp_var_contrast_vec[3]
  comp_var_levels <- comp_var_contrast_vec[2:3]
  masc_data <- obs %>%
    mutate(across(!!as.symbol(group_var), \(x) str_replace_all(x, ' ', '_'))) %>%
    filter(!str_detect(lineage, 'Doublet')) %>%
    filter(!!as.symbol(comp_var) != 'NA')

  masc_data[,comp_var] <- factor(masc_data %>% pull(comp_var), levels=comp_var_levels)

  masc_res <- MASC(masc_data,
                   cluster = masc_data %>% pull(group_var),
                   contrast = comp_var,
                   random_effects = "sample_id",
                   fixed_effects = fixed_effects,
                   verbose = TRUE, save_models = FALSE)
  masc_res[,group_var] <- masc_res$cluster %>% str_remove_all('cluster')

  return(masc_res)
}


masc_helper <- function(obs,
                        comp_var_contrast_vec,
                        group_var,
                        group_subset=NULL,
                        masc_res=NULL,
                        fixed_effects="",
                        colors=c('slategray', 'tomato4'),
                        row_order=F) {
  comp_var <- comp_var_contrast_vec[1]

  if (is.null(masc_res)) {  # if i just need to plot masc results
      masc_res <- run_masc(obs, comp_var_contrast_vec, group_var, fixed_effects = fixed_effects)
  }

  abundance <- create_abundance(obs, comp_var, group_var)
  abundance[,group_var] <- str_replace_all(abundance %>% pull(group_var), '_', ' ')
  abundance[,comp_var] <- factor(abundance %>% pull(comp_var), levels=comp_var_contrast_vec[2:3])

  if (!is.null(group_subset)) {
    subset_masc_res <- masc_res %>% filter(!!as.symbol(group_var) %in% group_subset)
    subset_abundance <- abundance %>% filter(!!as.symbol(group_var) %in% group_subset)
  } else {
    subset_masc_res <- masc_res
    subset_abundance <- abundance
  }

  odds_ratio_and_box_plot(subset_abundance,
                          subset_masc_res,
                          comp_var,
                          group_var,
                          '',
                          colors=colors,
                          row_order=row_order)

}
