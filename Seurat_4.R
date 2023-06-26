.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(ggplot2)
library(RColorBrewer)

source('source_cellbender/function_visualization.R')
source('source_cellbender/function_seurat.R')
source('source_cellbender/function_deviance.R')

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

seurat_7dpi_4 <- run_seurat_somi_function(seurat_7dpi_3, cells = remove.cells, invert = TRUE,
                                          batchcorrection = TRUE, batch = "seq_date")

DimPlot(seurat_7dpi_4, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_4, group.by = "seurat_clusters", split.by = "day")
FeaturePlot(seurat_7dpi_4, "total_counts")
FeaturePlot(seurat_7dpi_4, "log10_total_counts")
FeaturePlot(seurat_7dpi_4, "pct_counts_MT")
VlnPlot(seurat_7dpi_4, "log10_total_counts", pt.size=0)
VlnPlot(seurat_7dpi_4, "log10_total_features_by_counts", pt.size=0)
VlnPlot(seurat_7dpi_4, "pct_counts_MT", pt.size=0)

save(seurat_7dpi_4, file="data2/Seurat/seurat_7dpi_4.RData")
# load("data2/Seurat/seurat_7dpi_4.RData")

somi_violinplot(seurat_7dpi_4, "Hbb-bt")
DimPlot(seurat_7dpi_4, cells.highlight = colnames(seurat_7dpi_4[, seurat_7dpi_4$seurat_clusters %in% c(17)]))

somi_featureplot(seurat_7dpi_4, "Ptprc")
somi_featureplot(seurat_7dpi_4, "Epcam")
somi_featureplot(seurat_7dpi_4, "Ppbp")
somi_featureplot(seurat_7dpi_4, "Pf4")
somi_featureplot(seurat_7dpi_4, "Msln")
somi_featureplot(seurat_7dpi_4, "Pecam1")
somi_featureplot(seurat_7dpi_4, "Vwf")
somi_featureplot(seurat_7dpi_4, "Cdh5")
somi_featureplot(seurat_7dpi_4, "Hba-a2")
somi_featureplot(seurat_7dpi_4, "Hbb-bt")

somi_featureplot(seurat_7dpi_4, "Lgr5")
somi_featureplot(seurat_7dpi_4, "Slc12a2")
somi_featureplot(seurat_7dpi_4, "Olfm4")
somi_featureplot(seurat_7dpi_4, "Ascl2")
somi_featureplot(seurat_7dpi_4, "Gkn3")
somi_featureplot(seurat_7dpi_4, "Mcm6")
somi_featureplot(seurat_7dpi_4, "Mcm5")
somi_featureplot(seurat_7dpi_4, "Cdk4")
somi_featureplot(seurat_7dpi_4, "Mki67")
somi_featureplot(seurat_7dpi_4, "Arg2")
somi_featureplot(seurat_7dpi_4, "Il18")
somi_featureplot(seurat_7dpi_4, "Car4")
somi_featureplot(seurat_7dpi_4, "Ccl25")
somi_featureplot(seurat_7dpi_4, "Apoa4")
somi_featureplot(seurat_7dpi_4, "Apoa1")
somi_featureplot(seurat_7dpi_4, "Alpi")
somi_featureplot(seurat_7dpi_4, "Fabp1")
somi_featureplot(seurat_7dpi_4, "Dll1")
somi_featureplot(seurat_7dpi_4, "Atoh1")
somi_featureplot(seurat_7dpi_4, "Neurog3")
somi_featureplot(seurat_7dpi_4, "Neurod1")
somi_featureplot(seurat_7dpi_4, "Chga")
somi_featureplot(seurat_7dpi_4, "Chgb")
somi_featureplot(seurat_7dpi_4, "Dclk1")
somi_featureplot(seurat_7dpi_4, "Trpm5")
somi_featureplot(seurat_7dpi_4, "Gfi1")
somi_featureplot(seurat_7dpi_4, "Spdef")
somi_featureplot(seurat_7dpi_4, "Tff3")
somi_featureplot(seurat_7dpi_4, "Muc2")
somi_featureplot(seurat_7dpi_4, "Agr2")
somi_featureplot(seurat_7dpi_4, "Clu")
somi_featureplot(seurat_7dpi_4, "Ly6a")
somi_featureplot(seurat_7dpi_4, "Anxa1")
somi_featureplot(seurat_7dpi_4, "Anxa3")
somi_featureplot(seurat_7dpi_4, "Cd44")
somi_featureplot(seurat_7dpi_4, "Ccnd1")
somi_featureplot(seurat_7dpi_4, "Ccnd2")
somi_featureplot(seurat_7dpi_4, "Bmi1")
somi_featureplot(seurat_7dpi_4, "Lrig1")
somi_featureplot(seurat_7dpi_4, "Hopx")
somi_featureplot(seurat_7dpi_4, "Tert")
somi_featureplot(seurat_7dpi_4, "Lyz1")
somi_featureplot(seurat_7dpi_4, "Defa17")
somi_featureplot(seurat_7dpi_4, "Defa22")
somi_featureplot(seurat_7dpi_4, "Defa24")

marker <- FindMarkers(seurat_7dpi_4, ident.1 = c(22), logfc.threshold = 1)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC > 0), n=40)

somi_featureplot(seurat_7dpi_4, "Fabp2")

#
load("D:/Project/INTESTINE/data/Seurat.3.dev/intestine_seurat.3.dev.RData")
head(colnames(intestine_seurat.3.dev))

seurat_7dpi_4$orig.celltype = "Undefined"
common_cells = intersect(colnames(seurat_7dpi_4), colnames(intestine_seurat.3.dev))
seurat_7dpi_4$orig.celltype[common_cells] = intestine_seurat.3.dev$celltype[common_cells]
seurat_7dpi_4$orig.celltype = as.factor(seurat_7dpi_4$orig.celltype)
seurat_7dpi_4$orig.celltype = relevel(seurat_7dpi_4$orig.celltype, ref="Undefined")
DimPlot(seurat_7dpi_4[,colnames(intestine_seurat.3.dev)], group.by = "orig.celltype")

#
DimPlot(seurat_7dpi_4, label=TRUE, label.size=8)
VlnPlot(seurat_7dpi_4, "log10_total_counts", pt.size=0)
VlnPlot(seurat_7dpi_4, "pct_counts_MT", pt.size=0)

remove.cells = colnames(seurat_7dpi_4[, seurat_7dpi_4$seurat_clusters %in% c(4,22)])
use.cells = colnames(seurat_7dpi_4)[!colnames(seurat_7dpi_4) %in% remove.cells]
length(use.cells)
length(remove.cells)

# deviance residuals for correcting UMI bias
seurat_7dpi_4_dev <- run_deviance_residual(seurat_7dpi_4, batch="seq_date")
# seurat_7dpi_4_dev <- FindNeighbors(seurat_7dpi_4_dev, reduction="harmony", assay="deviance", dims=1:10)
# seurat_7dpi_4_dev <- FindClusters(seurat_7dpi_4_dev, resolution=0.5)

save(seurat_7dpi_4_dev, file="data2/Seurat/seurat_7dpi_4_dev.RData")
# load("data2/Seurat/seurat_7dpi_4_dev.RData")

DimPlot(seurat_7dpi_4_dev, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_4_dev, group.by = "day")
DimPlot(seurat_7dpi_4_dev, group.by = "seq_date")
DimPlot(seurat_7dpi_4_dev, group.by = "sample")
DimPlot(seurat_7dpi_4_dev, group.by = "condition")

FeaturePlot(seurat_7dpi_4_dev, "log10_total_counts")
FeaturePlot(seurat_7dpi_4_dev, "pct_counts_MT")
VlnPlot(seurat_7dpi_4_dev, "log10_total_counts")
VlnPlot(seurat_7dpi_4_dev, "pct_counts_MT", pt.size = 0)

somi_featureplot(seurat_7dpi_4_dev, "Ptprc")
somi_featureplot(seurat_7dpi_4_dev, "Epcam")
somi_featureplot(seurat_7dpi_4_dev, "Ppbp")
somi_featureplot(seurat_7dpi_4_dev, "Pf4")
somi_featureplot(seurat_7dpi_4_dev, "Msln")
somi_featureplot(seurat_7dpi_4_dev, "Pecam1")
somi_featureplot(seurat_7dpi_4_dev, "Vwf")
somi_featureplot(seurat_7dpi_4_dev, "Cdh5")
somi_featureplot(seurat_7dpi_4_dev, "Hba-a2")
somi_featureplot(seurat_7dpi_4_dev, "Hbb-bt")

somi_featureplot(seurat_7dpi_4_dev, "Lgr5")
somi_featureplot(seurat_7dpi_4_dev, "Olfm4")
somi_featureplot(seurat_7dpi_4_dev, "Mki67")
somi_featureplot(seurat_7dpi_4_dev, "Apoa4")
somi_featureplot(seurat_7dpi_4_dev, "Neurod1")
somi_featureplot(seurat_7dpi_4_dev, "Chga")
somi_featureplot(seurat_7dpi_4_dev, "Dclk1")
somi_featureplot(seurat_7dpi_4_dev, "Clu")
somi_featureplot(seurat_7dpi_4_dev, "Anxa1")
somi_featureplot(seurat_7dpi_4_dev, "Defa17")
somi_featureplot(seurat_7dpi_4_dev, "Muc2")
somi_featureplot(seurat_7dpi_4_dev, "Agr2")
somi_featureplot(seurat_7dpi_4_dev, "Dll1")
somi_featureplot(seurat_7dpi_4_dev, "Atoh1")

# highMT = c(7,16,23,25)
# Immune = c(0,10,11,18,19)
# 
# marker <- FindMarkers(seurat_7dpi_4_dev, ident.1 = c(), logfc.threshold = 1)
# marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
# head(subset(marker, avg_logFC > 0), n=40)
# 
# somi_featureplot(seurat_7dpi_4_dev, "Jchain")