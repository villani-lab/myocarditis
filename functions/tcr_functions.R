
tcr_assoc_func <- function(dataset, cluster, contrast){
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

   # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]

    message(paste("Creating logistic mixed models for", test_cluster))

    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1")), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + ", contrast)), collapse = ""))

    # Run null and full mixed-effects models
    null_model <- glm(formula = null_fm, data = dataset,
                              family = binomial)
    full_model <- glm(formula = full_fm, data = dataset,
                              family = binomial)
    model_lrt <- anova(null_model, full_model, test = "LRT")
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint(full_model)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
          }
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                               size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chi)"]][2])
  output$model.padj <- p.adjust(output$model.pvalue, method = "fdr")
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(coef(x$full_model)[[paste(contrast_lvl2, "True", sep = "")]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[paste(contrast_lvl2, "True", sep = ""), "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[paste(contrast_lvl2, "True", sep = ""), "97.5 %"]))
  return(output)
}


