# all gene sets from MSigDB
gene_sets <-  gmtPathways('/projects/home/sramesh/gsea_db/msigdb_symbols.gmt')

run_gsea <- function(de_results) {

  ret <- data.frame(pathway=character(), pval=numeric(), padj=numeric(), ES=numeric(), NES=numeric(),
                    nMoreExtreme=numeric(), size=integer(), leadingEdge=list(), cluster=integer(), geneset=character())

  clusters <- de_results$cluster %>% unique() %>% sort()
  for (clust in clusters) {
    res <- de_results %>% filter(cluster == clust)
    # Rank the genes for GSEA
    res2 <- res %>%
      dplyr::select(gene_symbol, stat) %>%
      na.omit() %>%
      distinct() %>%
      group_by(gene_symbol) %>%
      summarize(stat=mean(stat)) # taking mean does nothing here, res2 is just a sorted table of the gene and it's stat (from res)
    ranks <- deframe(res2) # convert res2 from a dataframe to a vector/list

    # Look at 3 different groups of gene sets, plot all gene sets with adjusted p value < 0.1
    for (geneset in c("BIOCARTA", "HALLMARK", "KEGG")){
      pathways <- gene_sets[grep(geneset, names(gene_sets))]
      fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=10000)
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      fgseaResTidy$cluster <- clust
      fgseaResTidy$geneset <- geneset
      fgseaResTidy$leadingEdge <- as.character(fgseaResTidy$leadingEdge)
      ret <- rbind(ret, fgseaResTidy)

    }
  }

  return(ret)

}



