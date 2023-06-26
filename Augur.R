.libPaths(Sys.getenv("R_LIBS_INTESTINE_4"))
.libPaths()

library(Augur)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(cartography)
library(rcartocolor)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "data2/"
plotdir = "plots2/"

load(paste0(rdatadir, "Seurat/seurat_7dpi_5.RData"))

unique(seurat_7dpi_5$day)

Idents(seurat_7dpi_5) = seurat_7dpi_5$day
seurat_D0 = subset(seurat_7dpi_5, idents = "D0")
seurat_D3 = subset(seurat_7dpi_5, idents = "D3")
seurat_D7 = subset(seurat_7dpi_5, idents = "D7")

run_augur <- function(seurat, group="celltype"){
  
  print(dim(seurat))
  expr = as.matrix(seurat@assays$RNA@data)
  meta = seurat@meta.data
  meta$cell_type = meta[, group]
  meta$label = meta$type
  
  augur = calculate_auc(expr, meta, label_col="label", cell_type_col="cell_type")
  
  return(augur)
}

# augur_D0 = run_augur(seurat_D0, group="celltype4")
# augur_D3 = run_augur(seurat_D3)
# augur_D7 = run_augur(seurat_D7)
# 
# augur_AUC_D0 = augur_D0$AUC
# augur_AUC_D3 = augur_D3$AUC
# augur_AUC_D7 = augur_D7$AUC
# 
# save(augur_D0, file=paste0(rdatadir, "Seurat/augur_D0.RData"))
# save(augur_D3, file=paste0(rdatadir, "Seurat/augur_D3.RData"))
# save(augur_D7, file=paste0(rdatadir, "Seurat/augur_D7.RData"))
# 
# save(augur_AUC_D0, file=paste0(rdatadir, "Seurat/augur_AUC_D0.RData"))
# save(augur_AUC_D3, file=paste0(rdatadir, "Seurat/augur_AUC_D3.RData"))
# save(augur_AUC_D7, file=paste0(rdatadir, "Seurat/augur_AUC_D7.RData"))

#
augur_PClin_D0 = run_augur(seurat_D0, group="celltype4")
augur_PClin_D3 = run_augur(seurat_D3, group="celltype4")
augur_PClin_D7 = run_augur(seurat_D7, group="celltype4")

augur_PClin_AUC_D0 = augur_PClin_D0$AUC
augur_PClin_AUC_D3 = augur_PClin_D3$AUC
augur_PClin_AUC_D7 = augur_PClin_D7$AUC

save(augur_PClin_D0, file=paste0(rdatadir, "Seurat/augur_PClin_D0.RData"))
save(augur_PClin_D3, file=paste0(rdatadir, "Seurat/augur_PClin_D3.RData"))
save(augur_PClin_D7, file=paste0(rdatadir, "Seurat/augur_PClin_D7.RData"))

save(augur_AUC_D0, file=paste0(rdatadir, "Seurat/augur_AUC_D0.RData"))
save(augur_AUC_D3, file=paste0(rdatadir, "Seurat/augur_AUC_D3.RData"))
save(augur_AUC_D7, file=paste0(rdatadir, "Seurat/augur_AUC_D7.RData"))
