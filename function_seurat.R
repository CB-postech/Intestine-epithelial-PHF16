library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(scran)
library(harmony)

make_seurat_function <- function(sce, hvg, npcs=50){
  
  if(is.null(reducedDims(sce)) == FALSE){
    reducedDims(sce) = NULL
  }
  
  Seurat <- as.Seurat(sce,
                      counts = "counts",
                      data = "logcounts")
  VariableFeatures(Seurat) <- rownames(hvg)
  
  Seurat <- ScaleData(Seurat)
  
  Seurat <- RunPCA(Seurat, npcs=npcs,
                   features = VariableFeatures(Seurat),
                   reduction.key = "pca_",
                   verbose = FALSE)
  print(plot(Seurat@reductions$pca@stdev,
             xlab = "PC",
             ylab = "Eigenvalue"))
  
  return(Seurat)
}
seurat_cluster_function <- function(Seurat, PCA, reduction="pca", assay="RNA", resolution=0.8){
  
  objectname = deparse(substitute(Seurat))
  
  if(is.null(names(Seurat@active.ident))){
    names(Seurat@active.ident) = colnames(Seurat)
  }
  else{
    Seurat <- FindNeighbors(Seurat, dims=1:PCA, reduction = reduction, assay=assay,
                            features = VariableFeatures(Seurat))
    Seurat <- FindClusters(Seurat,
                           resolution = resolution)
  }
  
  Seurat <- RunTSNE(Seurat,
                    reduction = reduction,
                    assay = assay,
                    dims = 1:PCA,
                    features = VariableFeatures(Seurat),
                    check_duplicates = FALSE)
  Seurat <- RunUMAP(Seurat,
                    assay = assay,
                    reduction = reduction,
                    dims = 1:PCA)
}
norm_somi_function <- function(sce.qc){
  
  objectname = deparse(substitute(sce.qc))
  
  clusters <- quickCluster(sce.qc)
  sce.qc <- computeSumFactors(sce.qc, clusters=clusters)
  
  print(paste0("summary of size factors of ", objectname))
  print(summary(sizeFactors(sce.qc)))
  
  sce.qc <- normalize(sce.qc)
  return(sce.qc)
}
hvg_function <- function(sce, FDRcut=0.05, biocut=0.01){
  var.fit <- trendVar(sce,
                      method="spline", parametric=TRUE, use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  
  hvg <- var.out[which(var.out$FDR < FDRcut & var.out$bio > biocut),]
  print(nrow(hvg))
  return(hvg)
}
harmony_function <- function(seurat, batch, assay="RNA"){
  
  set.seed(10)
  seurat <- RunHarmony(seurat, batch, plot_convergence=TRUE, assay.use = assay)
  
  #set.seed(10)
  #seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:PCA)
  #seurat <- FindClusters(seurat, resolution = resolution)
  #seurat <- RunTSNE(seurat, reduction = "harmony", dims = 1:PCA, check_duplicates = FALSE)    
  #seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:PCA)
  
  return(seurat)
}
run_seurat_somi_function <- function(seurat, idents=NULL, cells=NULL, invert=FALSE, npcs=50,
                                     batchcorrection = FALSE, batch=NULL,
                                     PCA=15, resolution=0.8, biocut=0.01){
  
  if(sum(!is.null(idents)) + sum(!is.null(cells)) != 0){
    #Idents(seurat) = seurat$seurat_clusters
    seurat <- subset(seurat, idents = idents, cells=cells, invert = invert)
    print("subsetting is done")   
  }
  
  sce <- as.SingleCellExperiment(seurat)
  sce <- norm_somi_function(sce)
  hvg <- hvg_function(sce, biocut=biocut)
  print("highly variable genes done")
  
  Seurat <- make_seurat_function(sce, hvg, npcs=npcs)
  
  if(sum(batchcorrection) == 1){
    Seurat <- harmony_function(Seurat, batch)
    
    Seurat <- seurat_cluster_function(Seurat, PCA=PCA, resolution=resolution, assay = "RNA", reduction = "harmony")
  }else{
    Seurat <- seurat_cluster_function(Seurat, PCA=PCA, resolution=resolution, assay = "RNA")
  }
  print("clustering done")
  
  return(Seurat)
}
