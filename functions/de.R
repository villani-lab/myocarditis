# read in gene names to remove
tcr_genes <- read_tsv('tcr_genes.tsv')$'Approved symbol'
bcr_genes <- read_tsv('bcr_genes.tsv')$'Approved symbol'


# Read in pseudobulk counts as a tibble and convert to data frame 
read_counts <- function(counts_filepath){
  counts <- read_csv(counts_filepath)
  counts_rownames <- counts$featurekey
  counts_colnames <- colnames(counts)[-1]
  counts <- counts %>%
    dplyr::select(-featurekey) %>%
    round() %>%
    data.frame()
  rownames(counts) <- counts_rownames
  colnames(counts) <- counts_colnames
  return(counts)
}

# Read in metadata 
read_meta <- function(meta_filepath){
  meta <- read.csv(meta_filepath, row.names=1)
  meta <- meta %>% 
    filter(condition == 'myocarditis' | condition == 'control') %>%
    mutate(condition = factor(condition, levels=c('control', 'myocarditis')))
  return(meta)
}


# Run DE analysis
run_de_by_comp_var <- function(counts_filepath, meta_filepath, save_name, comp_var_contrast_vec,
                                deseq_formula=NULL, cell_cutoff=5){
  # get datasets
  counts <- read_counts(counts_filepath)
  meta <- read_meta(meta_filepath)
  
  # filter for samples in metadata
  counts <- dplyr::select(counts, rownames(meta))

  # get comparison variable
  comp_var <- comp_var_contrast_vec[1]

  # create dataframe to hold all results for all clusters
  all_res <- data.frame(baseMean=numeric(), log2FoldChange=numeric(), lfcSE=numeric(), 
                        stat=numeric(), pvalue=numeric(), padj=numeric(), 
                        gene_symbol=character(), cluster=numeric())
  
  # run DE and make plot for each cluster
  for (cluster in sort(unique(meta$cluster))) {
    print(paste0('Cluster ', cluster))
    
    # filter for cluster
    meta_cluster <- meta[meta$cluster == cluster,]
    counts_cluster <- counts[,rownames(meta)]
    
    # skip clusters with too few cells 
    if (sum(meta_cluster$n_cells) <= 100) {next}
    
    # filtering for sample/clusters that have more than cell_cutoff cells
    meta_cluster <- meta_cluster[meta_cluster$n_cells >= cell_cutoff, ]
    tot_cells <- sum(meta_cluster$n_cells)

    # skip if design variable has only one level
    if (length(unique(meta_cluster[,comp_var])) == 1) {next}

    # filter the count matrix with sample/clusters that have more than cell_cutoff cells
    counts_cluster <- counts_cluster[,rownames(meta_cluster)]
    
    # for each gene, record how many samples have a non-zero count 
    n_samp <- rowSums(counts_cluster != 0)
    
    # keep genes where least half of the samples have a non-zero count for that gene
    counts_cluster <- counts_cluster[n_samp > (nrow(meta_cluster) / 2),]
    
    # filter out tcr and bcr genes
    counts_cluster <- counts_cluster[!(rownames(counts_cluster) %in% c(tcr_genes, bcr_genes)),]
    
    # make sure count column names match metadata rownames for DESeq2
    stopifnot(colnames(counts_cluster) == rownames(meta_cluster))

    # get formula
    if (is.null(deseq_formula)) {
      deseq_formula <- formula(glue("~ {comp_var}"))
    }

    # run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts_cluster,
                                  colData = meta_cluster,
                                  design = deseq_formula)
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds, contrast=comp_var_contrast_vec))
    res <- res[!is.na(res$padj),]
    
    # just for the volcano plot
    res$gene_symbol <- rownames(res)
    
    # concatenate cluster results to all_res
    res$cluster <- cluster
    all_res <- rbind(all_res, res)
  }
  print('saving results...')
  write_csv(all_res, glue('{save_name}_de_by_{comp_var}_all_results.csv'))
  
  return(all_res)
}

get_heatmap_data <- function(de_results, heatmap_genes, cluster_map){
  # filter de results for select heatmap genes and add umap names
  de_results %>%
    mutate(cluster = as.character(cluster)) %>%
    select(gene_symbol, cluster, log2FoldChange, padj) %>%
    inner_join(heatmap_genes, by = 'gene_symbol') %>%
    left_join(cluster_map %>% 
                mutate(cluster_number = as.character(cluster_number)),
              by = c('cluster' = 'cluster_number')) %>%
    filter(!cluster_name %in% c('Doublets/RBCs', 'other')) %>%
    complete(gene_symbol, cluster_name) %>%
    arrange(gene_symbol)
}