somi_dimplot <- function(seurat, reduction="umap", group.by, size=1, col=NULL, label=TRUE, legend=TRUE){
  
  d <- DimPlot(seurat,
               reduction = reduction,
               group.by = group.by,
               pt.size = size,
               order=TRUE,
               label = label,
               label.size = 6) +
    labs(title = paste0("by ", group.by)) +
    theme(plot.title = element_text(size=20, face="italic", hjust=0.5, margin=margin(0,0,15,0,"mm")),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()
    )
  if(sum(is.null(col)) == 0){
    d <- d + scale_color_manual(values = col, limits=names(col), breaks=names(col))
  }
  if(legend == FALSE){
    d <- d + theme(legend.position = "none")
  }

  print(d)
}
dimplot_function <- function(seurat, groups, reduction=c("pca", "tsne", "umap"), seuratname, size=1, col=NULL, label=TRUE, legend=TRUE, dir=plotdir, width=NA, height=NA){
  
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  path = paste0(dir, "dimplot/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  for(red in reduction)
    for(i in groups){
      d <- somi_dimplot(seurat, reduction=red, group.by = i, size=size, col=col, label=label, legend=legend)
      ggsave(paste0(path, seuratname, "_", red, "_", i, ".png"), plot=d, width=width, height=height)
    }
}
celldata_somi_featureplot <- function(seurat, feature, reduction="umap", size=1){
  
  FeaturePlot(seurat,
              features = feature,
              reduction = reduction,
              order=TRUE, pt.size=size) +
    labs(title = feature) +
    theme(plot.title = element_text(size=20, face="bold.italic"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), 
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          #legend.position = "none",
          legend.text=element_text(size=15), 
          legend.title=element_blank(),
          legend.key=element_blank(),
          axis.text = element_blank())
  
}
celldataplot_function <- function(seurat, features, seuratname, size=1, dir=plotdir){
  
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  path = paste0(dir, "dimplot/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  for(i in features){
    reduction="umap"
    if(sum(i %in% colnames(seurat@meta.data)) == 0){
      print(paste0("No ", i, ", in metadata slot"))
      next
    }
    c <- celldata_somi_featureplot(seurat, i, size=size)
    ggsave(paste0(path, seuratname, "_", reduction, "_", i, ".png"), plot=c)
  }
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

seurat_highlightplot <- function(seurat, discrete_val, seuratname, size=1, cols=NULL, dir=plotdir){
  
  objectname = deparse(substitute(seurat))
  
  Idents(seurat) = seurat@meta.data[,discrete_val]
  if(sum(is.null(cols)) == 1){
    Idents(seurat) = as.factor(Idents(seurat))
  }else{
    Idents(seurat) = factor(Idents(seurat), levels = names(cols))
  }
  print(levels(Idents(seurat)))
  
  n = length(levels(Idents(seurat)))
  print(n)
  
  if(sum(is.null(cols)) == 1){
    cols = gg_color_hue(n)
    names(cols) = levels(Idents(seurat))
  } 
  
  for(j in 1:n){
    u <- UMAPPlot(seurat,
             pt.size = size, 
             sizes.highlight = size,
             cells.highlight = WhichCells(seurat, idents = names(cols)[j])) +
      #labs(title = paste0("cluster ", j)) +
      scale_color_manual(values = c("#D3D3D3", cols[[j]]),
                         labels = c("unselected", names(cols)[j])) +
      theme_classic() +
      theme(#plot.title = element_text(size=25, face="bold.italic", hjust=0.5),
        #plot.subtitle = element_text(size=15, face="italic", hjust=0.5),
        plot.title = element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none",
        legend.text=element_text(size=15), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text = element_blank())
    ggsave(paste0(dir, seuratname, "_umap_by_", discrete_val, "_", levels(Idents(seurat))[j], ".png"), plot=u)
  }
}
highlight_function <- function(seurat, discrete_val_list, seuratname, dir, size=1, cols=NULL){
  
  objectname = deparse(substitute(seurat))
  
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  path = paste0(dir, "dimplot/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  
  for(k in 1:length(discrete_val_list)){
    
    val = discrete_val_list[k]
    
    path2 = paste0(path, val, "/")
    if(sum(dir.exists(path2)) == 0){
      dir.create(path2)
    }
    
    seurat_highlightplot(seurat, val, seuratname, size=size, dir=path2, cols=cols)
    
  }
}
somi_featureplot <- function(seurat, feature, reduction="umap", col=NULL, size=0.1){
  
  if(sum(startsWith(rownames(seurat)[1], "ENS")) == 0){
    feature_id = feature
  }else if(sum(feature %in% ensemblGenes$external_gene_name) == 1){
    DefaultAssay(seurat) = "RNA"
    feature_id = ensemblGenes[ensemblGenes$external_gene_name == feature,]$ensembl_gene_id[1]
    print(c(feature, feature_id))
  }else{
    print("No feature")
  }
  
  
  if(is.null(col)){
    col = c("#EEEEEE", brewer.pal(9, "Reds"))
  }
  if(sum(feature_id %in% rownames(seurat)) == 0 && sum(feature %in% colnames(seurat@meta.data)) == 0){
    print(paste0("no ", feature))
  }else{
    
    if(sum(feature %in% ensemblGenes$external_gene_name) == 1){
      expr = GetAssayData(seurat, slot="data", assay="RNA")
      max = max(expr[feature_id,])
      min = min(expr[feature_id,])
    }else{
      max = max(seurat@meta.data[, feature_id])
      min = min(seurat@meta.data[, feature_id])
   }
    
    FeaturePlot(seurat,
                features = feature_id,
                reduction = reduction,
                order=TRUE, 
                pt.size=size) +
      labs(title = feature) +
      scale_color_gradientn(colours = col, 
                            breaks=c(min, (max-min)/2 + min, 
                                     max), 
                            labels=c(round(min, digit=1), round((max-min)/2 + min, digit=1), 
                                     round(max, digit=1))) +
      theme(plot.title = element_text(size=20, face="bold.italic"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            #legend.position = "none",
            legend.text=element_text(size=15), 
            legend.title=element_blank(),
            legend.key=element_blank(),
            axis.text = element_blank())
  }
}
feature_iteration <- function(seurat, genes, seuratname, size=1, dir=plotdir, width=6, height=6){
  
  for(i in genes){
    reduction="umap"
    feature_id = ensemblGenes[ensemblGenes$external_gene_name %in% i,]$ensembl_gene_id[1]
    
    if(sum(feature_id %in% rownames(seurat)) == 0){
      print(paste0("No ", i, ", so passed"))
      next
    }else if(sum(seurat@assays$RNA@data[feature_id,]) == 0){
      print(paste0("no ", i))
      next
    }else{
      f <- somi_featureplot(seurat, i, size=size)
      ggsave(paste0(dir, seuratname, "_", reduction, "_", i, ".png"), plot=f, width=width, height=height)
    }
  }
}
featureplot_function <- function(seurat, genes, seuratname, size=1, path=plotdir, width=6, height=6){
  
  if(sum(is.list(genes)) == 1){
    for(l in 1:length(genes)){
      groupname = gsub(" ", "_", names(genes)[l])
      path2 = paste0(path, groupname, "/")
      if(sum(dir.exists(path2)) == 0){
        dir.create(path2)
      }
      
      genes_in_list = genes[[l]]
      feature_iteration(seurat, genes_in_list, seuratname, size=size, dir=path2, width=width, height=height)
    }
    
  }else{
    feature_iteration(seurat, genes, seuratname, size=size, dir=path, width=width, height=height)
  }
}
run_feature <- function(seurat, genes, seuratname=NULL, size=1, pseudodir=NULL, width=6, height=6){
  
  if(sum(is.null(seuratname)) == 1){
    seuratname = deparse(substitute(seurat))
    seuratname_forpath = gsub("oskm_", "", seuratname)
  }
  
  if(sum(is.null(pseudodir)) == 1){
    pseudodir = paste0(plotdir, seuratname_forpath)
  }
  if(sum(dir.exists(pseudodir)) == 0){
    dir.create(pseudodir)
  }
  dir = paste0(pseudodir, "featureplot/")
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  if(sum(is.list(genes)) == 1){
    genelistname = deparse(substitute(genes))
    plotsavedir = paste0(dir, genelistname, "/")
    if(sum(dir.exists(plotsavedir)) == 0){
      dir.create(plotsavedir)
    }            
  }else{
    plotsavedir = dir
  }
  
  featureplot_function(seurat, genes, seuratname, size=size, path=plotsavedir, width=width, height=height)
}
somi_violinplot <- function(seurat, feature, group.by=NULL, split.by=NULL, col=NULL, size=0){
  
  DefaultAssay(seurat) = "RNA"
  # feature_id = ensemblGenes[ensemblGenes$external_gene_name == feature,]$ensembl_gene_id[1]
  # print(c(feature, feature_id))
  
  if(sum(startsWith(rownames(seurat)[1], "ENS")) == 0){
    feature_id = feature
  }else if(sum(feature %in% ensemblGenes$external_gene_name) == 1){
    DefaultAssay(seurat) = "RNA"
    feature_id = ensemblGenes[ensemblGenes$external_gene_name == feature,]$ensembl_gene_id[1]
    print(c(feature, feature_id))
  }else{
    print("No feature")
  }
  
  if(sum(feature_id %in% rownames(seurat)) == 0){
    #print(paste0("no ", feature))
  }else if(sum(seurat@assays$RNA@data[feature_id,]) == 0){
    print(paste0("no ", feature))
  }else{
    
  VlnPlot(seurat, sort=FALSE, cols=col,
          features = feature_id,
          group.by = group.by, 
          split.by = split.by, 
          pt.size=size) +
    geom_boxplot(width=0.05, outlier.shape=NA) +
    #stat_compare_means(comparisons = comparison_list, label="p.signif", test="t.test") +
    labs(title = feature) +
    theme(plot.title = element_text(size=40, face="bold.italic"),
          axis.title = element_text(size=20),
          axis.text = element_text(size=20),
          axis.title.x = element_blank(),
          legend.text = element_text(size=20))
    
  }
}
violin_iteration <- function(seurat, genes, seuratname, group.by=NULL, split.by=NULL, split.plot=FALSE, size=0, col=NULL, dir=plotdir){
  
  for(i in genes){
    reduction="umap"
    feature_id = ensemblGenes[ensemblGenes$external_gene_name %in% i,]$ensembl_gene_id[1]
    
    if(sum(feature_id %in% rownames(seurat)) == 0){
      print(paste0("No ", i, ", so passed"))
      next
    }else if(sum(seurat@assays$RNA@data[feature_id,]) == 0){
      print(paste0("no ", i))
      next
    }else{
      v <- somi_violinplot(seurat, i, group.by=group.by, split.by=split.by, split.plot=split.plot, col=col, size=size)
      ggsave(paste0(dir, seuratname, "_violin_", i, ".png"), plot=v)
    }
  }
}
violinplot_function <- function(seurat, genes, seuratname, group.by=NULL, split.by=NULL, split.plot=FALSE, size=0, col=NULL, path=plotdir){
  
  if(sum(is.list(genes)) == 1){
    for(l in 1:length(genes)){
      groupname = gsub(" ", "_", names(genes)[l])
      path2 = paste0(path, groupname, "/")
      if(sum(dir.exists(path2)) == 0){
        dir.create(path2)
      }
      
      genes_in_list = genes[[l]]
      violin_iteration(seurat, genes_in_list, seuratname, group.by=group.by, split.by=split.by, split.plot=split.plot, size=size, col=col, dir=path2)
    }
    
  }else{
    violin_iteration(seurat, genes, seuratname, group.by=group.by, split.by=split.by, split.plot=split.plot, size=size, col=col, dir=path)
  }
}
run_violin <- function(seurat, genes, group.by=NULL, split.by=NULL, split.plot=FALSE, size=0, col=NULL, pseudodir=NULL){
  
  seuratname = deparse(substitute(seurat))
  seuratname_forpath = gsub("oskm_", "", seuratname)
  if(sum(is.null(pseudodir)) == 1){
    pseudodir = paste0(plotdir, seuratname_forpath)
  }
  if(sum(dir.exists(pseudodir)) == 0){
    dir.create(pseudodir)
  }
  dir = paste0(pseudodir, "/vlnplot/")
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  if(sum(is.list(genes)) == 1){
    genelistname = deparse(substitute(genes))
    plotsavedir = paste0(dir, genelistname, "/")
    if(sum(dir.exists(plotsavedir)) == 0){
      dir.create(plotsavedir)
    }            
  }else{
    plotsavedir = dir
  }
  
  violinplot_function(seurat, genes, seuratname, group.by=group.by, split.by=split.by, split.plot=split.plot, size=size, col=col, path=plotsavedir)
}
comparison_list = list( c("WT", "WTTAA"), c("KO", "KOTAA"), c("WT", "KO"), c("WTTAA", "KOTAA"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
somi_violinplot_2 <- function(seurat, symbol, x, split.by=NULL, comparison_list=NULL, test="t.test", assay="RNA", ncol=3, col=NULL, ylim=NULL){
  
  if(assay == "RNA"){
    geneid = ensemblGenes[ensemblGenes$external_gene_name %in% symbol, ]$ensembl_gene_id[1]
    df <- data.frame(gene = seurat@assays$RNA@data[geneid,])
    print(c(symbol, geneid))
  }else if(assay == "dorothea"){
    df <- data.frame(gene = seurat@assays$dorothea@data[symbol,])
  }
  
  if(sum(is.null(split.by)) == 0){
    vector = unique(c(x, split.by))
  }else{
    vector = x
  }
  
  meta.data = as.data.frame(seurat@meta.data[, vector])
  colnames(meta.data) = vector
  
  df <- data.frame(df, meta.data)
  
  if(sum(is.null(col)) == 1){
    n = length(unique(df[,x]))
    col = gg_color_hue(n)
  }
  
  if(sum(is.null(split.by)) == 0){
    g <- ggplot(df, aes(x=eval(parse(text=x)), y=gene, fill=eval(parse(text=x)))) +
      geom_violin() +
      #stat_summary(fun.y=median, geom="point", size=2, color="black") +
      scale_fill_manual(values = col) +
      facet_wrap(~eval(parse(text=split.by)), ncol=ncol) +
      theme_classic() +
      ylab(symbol) +
      xlab(x) +
      labs(fill=x) +
      stat_compare_means(comparisons = comparison_list, label="p.signif", test=test) +
      theme(strip.background = element_rect(linetype = "blank", fill="lightgray"))
    
  }else{
    g <- ggplot(df, aes(x=eval(parse(text=x)), y=gene, fill=eval(parse(text=x)))) +
      geom_violin() +
      #stat_summary(fun.y=median, geom="point", size=2, color="black") +
      scale_fill_manual(values = col) +
      theme_classic() +
      ylab(symbol) +
      xlab(x) +
      labs(fill=x) +
      stat_compare_means(comparisons = comparison_list, label="p.signif", test="t.test") +
      theme(strip.background = element_rect(linetype = "blank", fill="lightgray"))
  }
  if(sum(is.null(ylim)) == 0){
    g <- g + ylim(0, ylim)
  }
  print(g)
}
violin_iteration_2 <- function(seurat, symbol, x, split.by=NULL, comparison_list, test, assay="RNA", ncol=3, col=NULL, dir=plotdir, width=NA, height=NA){
  
  for(i in symbol){
    feature_id = ensemblGenes[ensemblGenes$external_gene_name %in% i,]$ensembl_gene_id[1]
    
    if(sum(feature_id %in% rownames(seurat)) == 0){
      print(paste0("No ", i, ", so passed"))
      next
    }else if(sum(seurat@assays$RNA@data[feature_id,]) == 0){
      print(paste0("no ", i))
      next
    }else{
      if(sum(is.null(split.by)) == 1){
        v <- somi_violinplot_2(seurat, i, x, split.by, comparison_list, test, assay, ncol, col)
        ggsave(paste0(dir, i, ".png"), plot=v, width=width, height=height)
      }else{
        v <- somi_violinplot_2(seurat, i, x, split.by, comparison_list, test, assay, ncol, col)
        ggsave(paste0(dir, split.by, "_", i, ".png"), plot=v, width=width, height=height)
      }
      
    }
  }
}
violinplot_function_2 <- function(seurat, symbol, x, split.by=NULL, comparison_list, test, assay="RNA", ncol=3, col=NULL, path=plotdir, width=NA, height=NA){
  
  if(sum(is.list(symbol)) == 1){
    for(l in 1:length(symbol)){
      groupname = gsub(" ", "_", names(symbol)[l])
      path2 = paste0(path, groupname, "/")
      if(sum(dir.exists(path2)) == 0){
        dir.create(path2)
      }
      
      genes_in_list = symbol[[l]]
      violin_iteration_2(seurat, genes_in_list, x, split.by, comparison_list, test, assay, ncol, col, dir=path2, width=width, height=height)
    }
    
  }else{
    violin_iteration_2(seurat, symbol, x, split.by, comparison_list, test, assay, ncol, col, dir=path, width=width, height=height)
  }
}
run_violin_2 <- function(seurat, symbol, x, split.by=NULL, comparison_list, test, assay="RNA", ncol=3, col=NULL, pseudodir=NULL, width=NA, height=NA){
  
  seuratname = deparse(substitute(seurat))
  seuratname_forpath = gsub("oskm_", "", seuratname)
  if(sum(is.null(pseudodir)) == 1){
    pseudodir = paste0(plotdir, seuratname_forpath)
  }
  if(sum(dir.exists(pseudodir)) == 0){
    dir.create(pseudodir)
  }
  dir = paste0(pseudodir, "/vlnplot/")
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  if(sum(is.list(symbol)) == 1){
    genelistname = deparse(substitute(symbol))
    plotsavedir = paste0(dir, genelistname, "/")
    if(sum(dir.exists(plotsavedir)) == 0){
      dir.create(plotsavedir)
    }            
  }else{
    genelistname = deparse(substitute(symbol))
    plotsavedir = paste0(dir, genelistname, "/")
    if(sum(dir.exists(plotsavedir)) == 0){
      dir.create(plotsavedir)
    }
  }
  
  violinplot_function_2(seurat, symbol, x, split.by, comparison_list, test, assay, ncol, col, path=plotsavedir, width=width, height=height)
}
somi_dotplot <- function(seurat, genesymbol, group.by=NULL, cluster_order=NULL, cols=rev(brewer.pal(9, "RdBu"))){
  
  if(sum(is.list(genesymbol)) == 1){
    genesymbol = unlist(genesymbol, use.names=FALSE)
  }
  if(sum(is.null(cluster_order)) == 1){
    cluster_order = levels(as.factor(seurat@meta.data[,group.by]))
    print(cluster_order)
  }
  cluster_order = rev(cluster_order)
  print(cluster_order)
  
  df = ensemblGenes[match(genesymbol, ensemblGenes$external_gene_name), ]
  genes = unique(df$ensembl_gene_id)
  genesymbol = unique(df$external_gene_name)
  
  DotPlot(seurat, assay = "RNA", features = genes, group.by = group.by, col.min = -1.5, col.max = 1.5) +
    scale_x_discrete(labels = genesymbol) + 
    scale_y_discrete(limits = as.character(cluster_order)) +
    scale_color_gradientn(colors = cols) +
    theme(plot.title = element_text(size=50, face="bold.italic"),
          #panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank(), 
          panel.border = element_blank(),
          #legend.position = "none",
          legend.text=element_text(size=15), 
          legend.title=element_text(size=18),
          #legend.key=element_blank(),
          axis.text.x = element_text(size=20, angle = 90, hjust=1, vjust=0.5),
          axis.text.y = element_text(size=20),
          #legend.direction = "vertical", legend.box = "horizontal",
          axis.title.y = element_blank())
}
dotplot_function <- function(seurat, genesymbol, group.by=NULL, cluster_order=NULL, cols=rev(brewer.pal(9, "RdBu")), seuratname="", dir, width=NA, height=NA){
  
  if(sum(is.list(genesymbol)) == 1){
    genelistname = deparse(substitute(genesymbol))
    genelistname = paste0(genelistname, "_")
  }else{
    genelistname = ""
  }
  
  dir = paste0(dir, "dotplot/")
  if(sum(dir.exists(dir)) == 0){
    dir.create(dir)
  }
  
  d <- somi_dotplot(seurat, genesymbol, group.by, cluster_order, cols)
  ggsave(paste0(dir, seuratname, "_", genelistname, group.by, "_dotplot.png"), plot=d, width=width, height=height)
}