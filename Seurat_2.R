.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(ggplot2)
library(RColorBrewer)

source('source_cellbender/function_visualization.R')
source('source_cellbender/function_seurat.R')
source('source_cellbender/function_deviance.R')

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

seurat_7dpi_2 <- run_seurat_somi_function(seurat_7dpi_bc, idents = c(0,9,10,19,25,28,29), invert = TRUE,
                                          batchcorrection = TRUE, batch = "seq_date")

DimPlot(seurat_7dpi_2, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_2, group.by = "day", split.by = "day")
FeaturePlot(seurat_7dpi_2, "log10_total_counts")
FeaturePlot(seurat_7dpi_2, "pct_counts_MT")

save(seurat_7dpi_2, file="data2/Seurat/seurat_7dpi_2.RData")
# load("data2/Seurat/seurat_7dpi_2.RData")

somi_violinplot(seurat_7dpi_2, "Hbb-bt")
DimPlot(seurat_7dpi_2, cells.highlight = colnames(seurat_7dpi_2[, seurat_7dpi_2$seurat_clusters %in% c(17)]))

somi_featureplot(seurat_7dpi_2, "Ptprc")
somi_featureplot(seurat_7dpi_2, "Epcam")
somi_featureplot(seurat_7dpi_2, "Ppbp")
somi_featureplot(seurat_7dpi_2, "Pf4")
somi_featureplot(seurat_7dpi_2, "Msln")
somi_featureplot(seurat_7dpi_2, "Pecam1")
somi_featureplot(seurat_7dpi_2, "Vwf")
somi_featureplot(seurat_7dpi_2, "Cdh5")
somi_featureplot(seurat_7dpi_2, "Hba-a2")
somi_featureplot(seurat_7dpi_2, "Hbb-bt")
somi_violinplot(seurat_7dpi_2, "Hbb-bt")

somi_featureplot(seurat_7dpi_2, "Lgr5")
somi_featureplot(seurat_7dpi_2, "Olfm4")
somi_featureplot(seurat_7dpi_2, "Mki67")
somi_featureplot(seurat_7dpi_2, "Apoa4")
somi_featureplot(seurat_7dpi_2, "Neurod1")
somi_featureplot(seurat_7dpi_2, "Chga")
somi_featureplot(seurat_7dpi_2, "Dclk1")
somi_featureplot(seurat_7dpi_2, "Clu")
somi_featureplot(seurat_7dpi_2, "Anxa1")
somi_featureplot(seurat_7dpi_2, "Defa17")
somi_featureplot(seurat_7dpi_2, "Muc2")
somi_featureplot(seurat_7dpi_2, "Agr2")
somi_featureplot(seurat_7dpi_2, "Dll1")
somi_featureplot(seurat_7dpi_2, "Atoh1")

marker <- FindMarkers(seurat_7dpi_2, ident.1 = c(18), logfc.threshold = 1.5)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC < 0), n=40)

somi_featureplot(seurat_7dpi_2, "Bpgm")
