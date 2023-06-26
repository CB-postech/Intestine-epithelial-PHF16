library(ggpubr)
library(tidyselect)

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

somi_violinplot.fig.test <- function(seurat, feature, group.by, split.by=NULL, col=NULL, comparison_list=FALSE, assay="RNA", nrow=1, y.mul=1.1, stat=FALSE){
  
  feature_id = somi_select_genes(seurat, feature)
  print(feature_id)
  
  if(sum(sum(feature_id %in% rownames(seurat[[assay]])), sum(feature_id %in% colnames(seurat@meta.data))) == 0){
    print(paste0("no ", feature, " in seurat obj"))
    #next
  }else{
    if(sum(feature_id %in% rownames(seurat[[assay]])) == 1){
      expr = seurat[[assay]][feature_id, ] %>% as.vector()
    }else{
      feature_id = feature
      expr = seurat@meta.data[, feature_id]
      max.value = max(expr)
    }
    y.max = y.mul*max(expr)
    y.min = min(expr)
    
    df <- data.frame(expr = expr,
                     group1 = seurat[[group.by]])
    colnames(df) = c("expr", "group1")
    
    if(sum(is.null(split.by)) == 0){
      df = data.frame(df, 
                      group2 = seurat[[split.by]])
      colnames(df) = c("expr", "group1", "group2")
    }
    
    df.melted = reshape2::melt(df)
    
    v <- ggplot(df.melted, aes(x=group1, y=value, fill=group1)) +
      geom_violin(scale="width", alpha=0.5) +
      geom_boxplot(width=0.1, alpha=0.7, outlier.size=0.5) +
      geom_jitter(height = 0, size = 0.5, show.legend = FALSE) +
      theme_classic() +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            strip.background = element_rect(fill="lightgray", colour=NA),
            axis.title.y = element_text(margin = margin(0,10,0,0)),
            axis.text.x = element_text(angle=45, hjust=1),
            text = element_text(size=15),
            plot.margin = margin(10,10,10,10)) +
      ylim(y.min, y.max) +
      ylab(feature)
    
    if(sum(is.null(col)) == 0){
      v <- v + scale_fill_manual(values = col)
    }
    
    if(sum(is.null(split.by)) == 0){
      v <- v + facet_wrap(~group2, nrow=nrow)
    }
    if(sum(stat) == 1){
      comparisons = as.list(as.data.frame(combn(as.vector(unique(df.melted$group1)),2), stringsAsFactors=FALSE))
      v <- v + stat_compare_means(label="p.signif", method="t.test", comparisons = comparisons, tip.length=0, symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")))
    }
    # if(sum(is.list(comparison_list)) == 1){
    #   v <- v + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum, step.increase = rep(0, length(comparison_list)))
    # }
    # if(sum(is.list(comparison_list)) == 0){
    #   if(sum(comparison_list) == 1){
    #     if(sum(is.null(split.by)) == 1){
    #       comparison_list = as.list(as.data.frame(combn(as.vector(unique(df.melted$group1)),2), stringsAsFactors=FALSE))
    #     }else if(sum(is.null(split.by)) == 0){
    #       comparison_list = as.list(as.data.frame(combn(as.vector(unique(df.melted$group1)),2), stringsAsFactors=FALSE))
    #     }
    #     v <- v + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum, step.increase = rep(0, length(comparison_list)))
    #   }
    # }
    return(v)
  }
}
