.libPaths(Sys.getenv("R_LIBS_INTESTINE"))
.libPaths()

setwd("D:/Project/INTESTINE")

library(SingleCellExperiment)
library(scater)
library(RColorBrewer)
library(stringr)
library(DropletUtils)
library(Seurat)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

plotdir = "plots_cellbender/"
rdatadir = "data_cellbender/"

filedir = "D:/Project/INTESTINE/cellbender/"

for(i in 1:12){
  h5_data <- hdf5r::H5File$new(paste0(filedir, "E", i, "_cellbender_filtered.h5"), mode = 'r')
  
  mat <- Matrix::sparseMatrix(
    i = h5_data[['matrix/indices']][],
    p = h5_data[['matrix/indptr']][],
    x = h5_data[['matrix/data']][],
    dimnames = list(
      h5_data[['matrix/features/name']][],
      h5_data[['matrix/barcodes']][]
    ),
    dims = h5_data[['matrix/shape']][],
    index1 = FALSE
  )
  
  sce <- SingleCellExperiment(assays = list(counts = mat))
  
  assign(paste0("cb_sce_e", i), sce)
}

rm(h5_data)
rm(sce)
rm(mat)

qc_hist_function <- function(sce, v=0, column, name, path){
  
  hist_column = sce[[column]]
  
  png(paste0(path, "hist_", column, ".png"))
  hist(
    hist_column,
    breaks = 100,
    xlab = name,
    main="")
  if(v != 0){
    abline(v = v, col="red")
  }
  dev.off()
}
qc_dimplot_function <- function(df, column, name, path){
  
  df = df[order(df[,column]),]
  g <- ggplot(df) +
    geom_point(aes(x=PC1, y=PC2, color = eval(parse(text = column)))) +
    scale_colour_gradientn(colours = brewer.pal(9, "Reds"),
                           limits = c(min(df[,column]), max(df[,column]))) +
    ggtitle(name) +
    theme_bw() +
    theme(plot.title = element_text(size = 15, face="italic", hjust = 0.5, margin = margin(0,0,4,0,"mm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size=15, face="italic"),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(path, "dimplot_", column, ".png"), plot=g, width=6, height=5)
  
}
qc_somi_function <- function(sce, total_counts_cutoff=0, total_features_cutoff=0, mt_percent_cutoff=0,
                             save = FALSE, objectname=NULL){
  
  if(sum(is.null(objectname)) == 1){
    objectname = deparse(substitute(sce))
    objectname = gsub("raw", "", objectname)
  }
  
  path = paste0(plotdir, "QC/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  path = paste0(path, objectname, "/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  
  mtGenes = ensemblGenes[ensemblGenes[,3] == "MT",]
  is.mito = rownames(sce) %in% mtGenes$external_gene_name
  print(paste0("There is ", length(is.mito[is.mito == TRUE]), " mitochondrial genes"))
  
  per.cell <- perCellQCMetrics(
    sce,
    subsets = list(MT=is.mito),
    percent_top = c(50, 100, 200, 500)
  )
  
  colData(sce) <- cbind(colData(sce), per.cell)
  
  print(table(sce$sum == 0))
  sce = sce[, sce$sum != 0]
  
  sce$log10_sum = log10(sce$sum + 1)
  sce$log10_detected = log10(sce$detected + 1)
  
  ### plot function
  
  coldata = c("log10_sum", "log10_detected", "subsets_MT_percent")
  columnName = c("log10 total umi counts", "log10 total feature counts", "Mitochondrial percentage (%)")
  cutoffs = c(total_counts_cutoff, total_features_cutoff, mt_percent_cutoff)
  
  # run pca
  vector <-  c(unique(colnames(sce@colData)))[-c(1,2)]
  #vector <- c("sum", "detected", "subsets_MT_percent")
  sce <- runColDataPCA(sce, ncomponents=5, variables = vector)
  
  # hist
  for(i in 1:length(coldata)){
    qc_hist_function(sce, v=cutoffs[i], coldata[i], columnName[i], path=path)
  }
  
  df <- as.data.frame(sce@int_colData$reducedDims$PCA_coldata)
  df$log10_sum <- sce$log10_sum
  df$log10_detected <- sce$log10_detected
  df$subsets_MT_percent <- sce$subsets_MT_percent
  
  # pca plot
  for(i in 1:length(coldata)){
    qc_dimplot_function(df, coldata[i], columnName[i], path=path)
  }
  
  filter_by_total_counts = sce$log10_sum > total_counts_cutoff
  filter_by_feature_counts = sce$log10_detected > total_features_cutoff
  filter_by_mt_percent = sce$subsets_MT_percent < mt_percent_cutoff
  
  sce$use <- (
    filter_by_feature_counts &
      filter_by_mt_percent &
      filter_by_total_counts
  )
  print(table(sce$use))
  
  ### save function
  if(sum(save) == 1){
    
    sce.qc <- sce[, sce$use]
    
    p <- plotReducedDim(sce, "PCA_coldata",
                        colour_by = "use",
                        size_by = "detected") +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            legend.text=element_text(size=13),
            legend.key=element_blank(),
            axis.text = element_blank())
    ggsave(paste0(path, "dimplot_use.png"), plot=p, width=6, height=5)
    
    print(paste0("sce.qc is generated. There are ", ncol(sce.qc), " cells finally."))
    return(sce.qc)
    
  }else{
    return(sce)
  }
}

for(i in 1:12){
  objectname = paste0("cb_sce_e", i)
  sce = eval(parse(text=objectname))
  sce1 <- qc_somi_function(sce, objectname = objectname,
                           total_counts_cutoff = 2.5, mt_percent_cutoff = 20)
  sce.qc <- qc_somi_function(sce, objectname = objectname,
                             total_counts_cutoff = 2.5, mt_percent_cutoff = 20, save = TRUE)
  assign(paste0("cb_sce_e", i), sce1)
  assign(paste0("cb_sce_e", i, ".qc"), sce.qc)
}

rm(sce)
rm(sce1)
rm(sce.qc)

for(i in 1:12){
  objectname = paste0("cb_sce_e", i)
  save(list = objectname, file=paste0(rdatadir, "sce/", objectname, ".RData"))
  objectname = paste0("cb_sce_e", i, ".qc")
  save(list = objectname, file=paste0(rdatadir, "sce.qc/", objectname, ".RData"))
}

rm(list = paste0("cb_sce_e", 1:12))
