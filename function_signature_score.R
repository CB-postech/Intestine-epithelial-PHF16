somi_select_genes <- function(seurat, Genes){
  
  if(sum(Genes %in% colnames(seurat@meta.data)) == 1){
    selected_geneids = Genes
  }else if(sum(starts_with("ENS", vars=rownames(seurat)[1])) == 1){
    if(sum(starts_with("ENS", vars=Genes[1])) == 1){
      selected_geneids = Genes
    }else{
      selected_geneids = ensemblGenes[match(Genes, ensemblGenes$external_gene_name), ]$ensembl_gene_id
    }
    selected_geneids = selected_geneids[selected_geneids %in% rownames(seurat)]
  }else{
    if(sum(starts_with("ENS", vars=Genes[1])) == 1){
      selected_geneids = ensemblGenes[Genes, ]$external_gene_name
    }else{
      selected_geneids = Genes
    }
    selected_geneids = selected_geneids[selected_geneids %in% rownames(seurat)]
  }
  return(selected_geneids)
}

avg_signature <- function(seurat, genes, group.by, comparison_list=TRUE, col=NULL, y.mul=1.2, symnum=list(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns"))){
  
  use.genes = somi_select_genes(seurat, genes)
  
  signature <- seurat@assays$RNA@data[use.genes,] %>% as.matrix() %>% colMeans()
  
  df <- data.frame(signature = signature,
                   group = seurat@meta.data[, group.by])
  
  y.max = y.mul*max(df$signature)
  y.min = min(df$signature)

  g <- ggplot(df, aes(x=group, y=signature, fill=group)) +
    geom_violin(scale="width", alpha=0.5) +
    geom_boxplot(width=0.1, alpha=0.5, outlier.size=0.5) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="lightgray", colour=NA),
          axis.title.y = element_text(margin = margin(0,10,0,0)),
          axis.text.x = element_text(angle=45, hjust=1),
          text = element_text(size=15),
          plot.margin = margin(10,10,10,10)) +
    ylim(y.min, y.max)
  
  if(sum(is.null(col)) == 0){
    g <- g + scale_fill_manual(values = col)
  }
  
  if(sum(is.list(comparison_list)) == 1){
    g <- g + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum,step.increase = rep(0, length(comparison_list)))
  }
  if(sum(is.list(comparison_list)) == 0){
    if(sum(comparison_list) == 1){
      # if(sum(is.null(split.by)) == 1){
      #   comparison_list = as.list(as.data.frame(combn(as.vector(unique(df$group)),2)))
      # }else if(sum(is.null(split.by)) == 0){
      #   comparison_list = as.list(as.data.frame(combn(as.vector(unique(df$group)),2)))
      # }
      comparison_list = as.list(as.data.frame(combn(as.vector(unique(df$group)),2)))
      g <- g + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum, step.increase = rep(0, length(comparison_list)))
    }
  }
  
  return(g)
}
