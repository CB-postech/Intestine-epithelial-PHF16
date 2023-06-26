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

somi_violinplot_split <- function(seurat, feature, group.by, split.by=NULL, split.2=NULL, col=NULL, assay="RNA", nrow=1, y.mul=1.1, stat=FALSE){
  
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
    if(sum(stat) == 1){
      y.max = y.mul*max(expr)
    }else{
      y.max = max(expr)
    }
    y.min = min(expr)
    
    df <- data.frame(expr = expr,
                     group1 = seurat[[group.by]])
    colnames(df) = c("expr", "group1")
    
    if(sum(is.null(split.by)) == 0){
      df = data.frame(df, 
                      group2 = seurat[[split.by]])
      colnames(df) = c("expr", "group1", "group2")
    }
    if(sum(is.null(split.2)) == 0){
      df = data.frame(df, 
                      group3 = seurat[[split.2]])
      colnames(df) = c("expr", "group1", "group2", "group3")
    }
    
    df.melted = reshape2::melt(df)
    
    v <- ggplot(df.melted, aes(x=group1, y=value, fill=group2)) +
      geom_split_violin(scale="width", alpha=0.5) +
      geom_boxplot(width=0.15, alpha=0.7, outlier.shape=NA) +
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

    if(sum(is.null(split.2)) == 0){
      v <- v + facet_wrap(~group3, nrow=nrow)
    }

    if(sum(stat) == 1){
      v <- v + stat_compare_means(label="p.signif", method="t.test", tip.length=0, symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf), symbols = c("**", "*", "ns")))
    }
    # if(sum(is.list(comparison_list)) == 1){
    #   v <- v + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum, step.increase = rep(0, length(comparison_list)))
    # }
    # if(sum(is.list(comparison_list)) == 0){
    #   if(sum(comparison_list) == 1){
    #     if(sum(is.null(split.by)) == 1){
    #       comparison_list = as.list(as.data.frame(combn(as.vector(unique(df.melted$group2)),2), stringsAsFactors=FALSE))
    #     }else if(sum(is.null(split.by)) == 0){
    #       comparison_list = as.list(as.data.frame(combn(as.vector(unique(df.melted$group2)),2), stringsAsFactors=FALSE))
    #     }
    #     v <- v + stat_compare_means(comparisons = comparison_list, label="p.signif", method="t.test", tip.length=0, symnum.args = symnum, step.increase = rep(0, length(comparison_list)))
    #   }
    # }
    return(v)
  }
}

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(
      trim = trim, scale = scale, draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}

GeomSplitViolin <- ggplot2:::ggproto("GeomSplitViolin",
                                     ggplot2::GeomViolin,
                                     draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                       data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                       grp <- data[1, "group"]
                                       newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                       newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                       newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                                       if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                         stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                                   1))
                                         quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                         aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                         aesthetics$alpha <- rep(1, nrow(quantiles))
                                         both <- cbind(quantiles, aesthetics)
                                         quantile_grob <- GeomPath$draw_panel(both, ...)
                                         ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                       }
                                       else {
                                         ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                       }
                                     }
)
