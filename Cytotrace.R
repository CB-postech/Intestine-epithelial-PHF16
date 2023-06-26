.libPaths(Sys.getenv("R_LIBS_INTESTINE_4"))

library(CytoTRACE)
library(Seurat)
library(dplyr)
library(ggplot2)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "data2/"
plotdir = "plots2/"

load("data2/Seurat/seurat_7dpi_5.RData")

# # row genes column cells
# expr <- as.matrix(seurat_7dpi_5@assays$RNA@data)

# run cytotrace
results <- CytoTRACE(as.matrix(seurat_7dpi_5@assays$RNA@data), enableFast = TRUE, ncores=1, subsamplesize = 1000)
cytotrace_res <- results
rm(results)

cytotraceScore <- cytotrace_res$CytoTRACE

save(cytotrace_res, file=paste0(rdatadir, "Seurat/cytotrace_res.RData"))
save(cytotraceScore, file=paste0(rdatadir, "Seurat/cytotraceScore.RData"))
