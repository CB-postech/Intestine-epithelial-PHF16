.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(ggplot2)
library(RColorBrewer)

source('source_cellbender/function_visualization.R')
source('source_cellbender/function_seurat.R')
source('source_cellbender/function_deviance.R')

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

seurat_7dpi_3 <- run_seurat_somi_function(seurat_7dpi_2, idents = c(14,18,25,26), invert = TRUE,
                                          batchcorrection = TRUE, batch = "seq_date")

DimPlot(seurat_7dpi_3, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_3, group.by = "day", split.by = "day")
FeaturePlot(seurat_7dpi_3, "total_counts")
FeaturePlot(seurat_7dpi_3, "log10_total_counts")
FeaturePlot(seurat_7dpi_3, "pct_counts_MT")
VlnPlot(seurat_7dpi_3, "log10_total_counts", pt.size=0)

save(seurat_7dpi_3, file="data2/Seurat/seurat_7dpi_3.RData")
# load("data2/Seurat/seurat_7dpi_3.RData")

somi_violinplot(seurat_7dpi_3, "Hbb-bt")
DimPlot(seurat_7dpi_3, cells.highlight = colnames(seurat_7dpi_3[, seurat_7dpi_3$seurat_clusters %in% c(17)]))

somi_featureplot(seurat_7dpi_3, "Ptprc")
somi_featureplot(seurat_7dpi_3, "Epcam")
somi_featureplot(seurat_7dpi_3, "Ppbp")
somi_featureplot(seurat_7dpi_3, "Pf4")
somi_featureplot(seurat_7dpi_3, "Msln")
somi_featureplot(seurat_7dpi_3, "Pecam1")
somi_featureplot(seurat_7dpi_3, "Vwf")
somi_featureplot(seurat_7dpi_3, "Cdh5")
somi_featureplot(seurat_7dpi_3, "Hba-a2")
somi_featureplot(seurat_7dpi_3, "Hbb-bt")

somi_featureplot(seurat_7dpi_3, "Lgr5")
somi_featureplot(seurat_7dpi_3, "Olfm4")
somi_featureplot(seurat_7dpi_3, "Mki67")
somi_featureplot(seurat_7dpi_3, "Apoa4")
somi_featureplot(seurat_7dpi_3, "Neurod1")
somi_featureplot(seurat_7dpi_3, "Chga")
somi_featureplot(seurat_7dpi_3, "Dclk1")
somi_featureplot(seurat_7dpi_3, "Clu")
somi_featureplot(seurat_7dpi_3, "Anxa1")
somi_featureplot(seurat_7dpi_3, "Defa17")
somi_featureplot(seurat_7dpi_3, "Muc2")
somi_featureplot(seurat_7dpi_3, "Agr2")
somi_featureplot(seurat_7dpi_3, "Dll1")
somi_featureplot(seurat_7dpi_3, "Atoh1")

marker <- FindMarkers(seurat_7dpi_3, ident.1 = c(2), logfc.threshold = 1.5)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC > 0), n=40)

somi_featureplot(seurat_7dpi_3, "Gzma")

#
load("D:/Project/INTESTINE/data/Seurat.3.dev/intestine_seurat.3.dev.RData")
head(colnames(intestine_seurat.3.dev))

seurat_7dpi_3$orig.celltype = "Undefined"
seurat_7dpi_3$orig.celltype[colnames(intestine_seurat.3.dev)] = intestine_seurat.3.dev$celltype
seurat_7dpi_3$orig.celltype = as.factor(seurat_7dpi_3$orig.celltype)
seurat_7dpi_3$orig.celltype = relevel(seurat_7dpi_3$orig.celltype, ref="Undefined")
DimPlot(seurat_7dpi_3[,colnames(intestine_seurat.3.dev)], group.by = "orig.celltype")

#
seurat_sub <- run_seurat_somi_function(seurat_7dpi_3, idents = c(2,12,20), invert = FALSE,
                                       batchcorrection = TRUE, batch = "seq_date", biocut = 0.1)
DimPlot(seurat_sub, label=TRUE, label.size=8)
DimPlot(seurat_sub, group.by = "day", split.by = "day")
FeaturePlot(seurat_sub, "log10_total_counts")
FeaturePlot(seurat_sub, "pct_counts_MT")
VlnPlot(seurat_sub, "log10_total_counts", pt.size=0)

somi_featureplot(seurat_sub, "Ptprc")
somi_featureplot(seurat_sub, "Epcam")
somi_featureplot(seurat_sub, "Ppbp")
somi_featureplot(seurat_sub, "Pf4")
somi_featureplot(seurat_sub, "Msln")
somi_featureplot(seurat_sub, "Pecam1")
somi_featureplot(seurat_sub, "Vwf")
somi_featureplot(seurat_sub, "Cdh5")
somi_featureplot(seurat_sub, "Hba-a2")
somi_featureplot(seurat_sub, "Hbb-bt")
somi_featureplot(seurat_sub, "Gzma")

somi_featureplot(seurat_sub, "Lgr5")
somi_featureplot(seurat_sub, "Olfm4")
somi_featureplot(seurat_sub, "Mki67")
somi_featureplot(seurat_sub, "Apoa4")
somi_featureplot(seurat_sub, "Neurod1")
somi_featureplot(seurat_sub, "Chga")
somi_featureplot(seurat_sub, "Dclk1")
somi_featureplot(seurat_sub, "Clu")
somi_featureplot(seurat_sub, "Anxa1")
somi_featureplot(seurat_sub, "Defa17")
somi_featureplot(seurat_sub, "Muc2")
somi_featureplot(seurat_sub, "Agr2")
somi_featureplot(seurat_sub, "Dll1")
somi_featureplot(seurat_sub, "Atoh1")

remove.cells = colnames(seurat_sub[, seurat_sub$seurat_clusters == "7"])
use.cells = colnames(seurat_sub)[!colnames(seurat_sub) %in% remove.cells]
length(use.cells)
length(remove.cells)

save(seurat_sub, file="data2/Seurat/seurat_sub.RData")
