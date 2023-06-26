.libPaths(Sys.getenv("R_LIBS_INTESTINE"))
.libPaths()

setwd("D:/Project/INTESTINE")

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

source('source_cellbender/function_visualization.R')

library(Seurat)
library(RColorBrewer)
library(ggplot2)

plotdir = "plots_cellbender/"
rdatadir = "data_cellbender/"

dup_genes = rownames(sce_total.norm)[duplicated(rownames(sce_total.norm))]
sce_total.norm = sce_total.norm[!rownames(sce_total.norm) %in% dup_genes,]

cb_seurat <- make_seurat_somi_function(sce_total.norm, hvg_total)
cb_seurat <- seurat_cluster_somi_function(cb_seurat, 20)

DimPlot(cb_seurat, label=TRUE, label.size=8)
DimPlot(cb_seurat, group.by="day", split.by="day")

VlnPlot(cb_seurat, "log10_sum", pt.size=0)
VlnPlot(cb_seurat, "subsets_MT_percent", pt.size=0)

somi_featureplot(cb_seurat, "Ptprc")
somi_featureplot(cb_seurat, "Epcam")
somi_featureplot(cb_seurat, "Ppbp")
somi_featureplot(cb_seurat, "Pf4")
somi_featureplot(cb_seurat, "Msln")
somi_featureplot(cb_seurat, "Pecam1")
somi_featureplot(cb_seurat, "Cdh5")
somi_featureplot(cb_seurat, "Vwf")
somi_featureplot(cb_seurat, "Hba-a2")
somi_featureplot(cb_seurat, "Hbb-bt")
somi_violinplot(cb_seurat, "Hbb-bt")

somi_featureplot(cb_seurat, "Jade3")
somi_violinplot(cb_seurat, "Jade3")

somi_featureplot(cb_seurat, "Cd3d")
somi_featureplot(cb_seurat, "Ms4a1")
somi_featureplot(cb_seurat, "Nkg7")
somi_featureplot(cb_seurat, "Cd300e")
somi_featureplot(cb_seurat, "Jchain")

somi_featureplot(cb_seurat, "Lgr5")
somi_featureplot(cb_seurat, "Slc12a2")
somi_featureplot(cb_seurat, "Olfm4")
somi_featureplot(cb_seurat, "Ascl2")
somi_featureplot(cb_seurat, "Gkn3")
somi_featureplot(cb_seurat, "Mcm6")
somi_featureplot(cb_seurat, "Mcm5")
somi_featureplot(cb_seurat, "Cdk4")
somi_featureplot(cb_seurat, "Mki67")
somi_featureplot(cb_seurat, "Arg2")
somi_featureplot(cb_seurat, "Il18")
somi_featureplot(cb_seurat, "Car4")
somi_featureplot(cb_seurat, "Ccl25")
somi_featureplot(cb_seurat, "Apoa4")
somi_featureplot(cb_seurat, "Apoa1")
somi_featureplot(cb_seurat, "Alpi")
somi_featureplot(cb_seurat, "Fabp1")
somi_featureplot(cb_seurat, "Dll1")
somi_featureplot(cb_seurat, "Atoh1")
somi_featureplot(cb_seurat, "Neurog3")
somi_featureplot(cb_seurat, "Neurod1")
somi_featureplot(cb_seurat, "Chga")
somi_featureplot(cb_seurat, "Chgb")
somi_featureplot(cb_seurat, "Dclk1")
somi_featureplot(cb_seurat, "Trpm5")
somi_featureplot(cb_seurat, "Gfi1")
somi_featureplot(cb_seurat, "Spdef")
somi_featureplot(cb_seurat, "Tff3")
somi_featureplot(cb_seurat, "Muc2")
somi_featureplot(cb_seurat, "Agr2")
somi_featureplot(cb_seurat, "Clu")
somi_featureplot(cb_seurat, "Ly6a")
somi_featureplot(cb_seurat, "Anxa1")
somi_featureplot(cb_seurat, "Anxa3")
somi_featureplot(cb_seurat, "Cd44")
somi_featureplot(cb_seurat, "Ccnd1")
somi_featureplot(cb_seurat, "Ccnd2")
somi_featureplot(cb_seurat, "Bmi1")
somi_featureplot(cb_seurat, "Lrig1")
somi_featureplot(cb_seurat, "Hopx")
somi_featureplot(cb_seurat, "Tert")
somi_featureplot(cb_seurat, "Lyz1")
somi_featureplot(cb_seurat, "Defa17")
somi_featureplot(cb_seurat, "Defa22")
somi_featureplot(cb_seurat, "Defa24")

somi_featureplot(cb_seurat, "Clu")

save(cb_seurat, file=paste0(rdatadir, "seurat/cb_seurat.RData"))
# load(paste0(rdatadir, "seurat/cb_seurat.RData"))