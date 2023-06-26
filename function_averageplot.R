# cell type annot - dotplot
# somi_select_genes <- function(seurat, Genes){
#   
#   if(sum(starts_with("ENS", vars=rownames(seurat)[1])) == 1){
#     if(sum(starts_with("ENS", vars=Genes[1])) == 1){
#       selected_geneids = Genes
#     }else{
#       selected_geneids = ensemblGenes[match(Genes, ensemblGenes$external_gene_name), ]$ensembl_gene_id
#     }
#     selected_geneids = selected_geneids[selected_geneids %in% rownames(seurat)]
#   }else{
#     if(sum(starts_with("ENS", vars=Genes[1])) == 1){
#       selected_geneids = ensemblGenes[Genes, ]$external_gene_name
#     }else{
#       selected_geneids = Genes
#     }
#     selected_geneids = selected_geneids[selected_geneids %in% rownames(seurat)]
#   }
#   return(selected_geneids)
# }
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
somi_averageMat <- function(seurat, Genes, group.by = "celltype", scale=TRUE){
  
  selected_geneids = somi_select_genes(seurat, Genes)
  
  data <- as.matrix(seurat@assays$RNA@data[selected_geneids, ])
  if(sum(scale) == 1){scale.data <- t(scale(t(data)))}else{scale.data <- data}
  
  df = data.frame(row.names = colnames(seurat), group=seurat@meta.data[,group.by], t(scale.data))
  df = df %>%
    group_by(group) %>%
    summarise_at(vars(selected_geneids), list(mean))
  df = tibble::column_to_rownames(df, var="group")
  colnames(df) = Genes
  mat = as.matrix(df)
  
  if(sum(scale) == 1){
    mat[mat > 1.5] = 1.5
    mat[mat < -1.5] = -1.5
  }
  
  return(mat)
}
somi_propMat <- function(seurat, Genes, group.by = "celltype"){
  
  selected_geneids = somi_select_genes(seurat, Genes)
  
  data <- as.matrix(seurat@assays$RNA@data[selected_geneids, ])
  propdata <- data > 0
  propdata <- propdata * 1
  
  df = data.frame(row.names = colnames(seurat), group=seurat@meta.data[,group.by], t(propdata))
  df = df %>%
    group_by(group) %>%
    summarise_at(vars(selected_geneids), list(mean))
  df = tibble::column_to_rownames(df, var="group")
  colnames(df) = Genes
  prop.mat = as.matrix(df)
  
  return(prop.mat)
}
somi_dotplot <- function(seurat, Genes, group.by = "celltype", cl.order=NULL, cells=NULL, range=c(0,6)){
  
  if(sum(is.null(cells)) == 0){
    seurat = subset(seurat, cells=cells)
  }
  
  if(sum(is.null(cl.order)) == 0){
    cl.level = cl.order
    seurat@meta.data[, group.by] = factor(seurat@meta.data[, group.by], levels = cl.order)
  }else{
    if(sum(is.null(levels(seurat@meta.data[, group.by]))) == 0){
      cl.level = levels(seurat@meta.data[, group.by])
    }else{
      cl.level = unique(seurat@meta.data[, group.by])
    }
  }
  
  mat <- somi_averageMat(seurat, Genes, group.by)
  prop.mat <- somi_propMat(seurat, Genes, group.by)
  
  dot.df <- data.frame(row.names = cl.level)
  dot.df <- cbind.data.frame(dot.df, mat)
  
  data.df <- data.frame(genes = rep(Genes, length(cl.level)))
  data.df$celltype = rep(cl.level, each=length(Genes))
  data.df$expr = as.vector(t(mat))
  data.df$prop = as.vector(t(prop.mat))
  
  min.value = min(data.df$expr)
  max.value = max(data.df$expr)
  abs.value = max(abs(min.value), max.value)
  
  data.df$celltype = as.character(data.df$celltype)
  
  gd <- ggplot(data.df) + 
    geom_point(aes(x=celltype, y=genes, color=expr, size=prop)) +
    scale_y_discrete(limits = rev(Genes)) + 
    scale_x_discrete(limits = as.character(cl.level)) +
    scale_color_gradientn(colors = rev(c(brewer.pal(11, "RdBu")[1:5], rep(brewer.pal(11, "RdBu")[6], 2), brewer.pal(11, "RdBu")[7:11])), limits = c(-abs.value, abs.value)) + 
    scale_radius(range = range) +
    theme_classic() +
    theme(text = element_text(size=20)) +
    ylab("") + xlab("")
  
  return(gd)
}

