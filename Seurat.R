.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(ggplot2)
library(RColorBrewer)

source('source_cellbender/function_visualization.R')
source('source_cellbender/function_seurat.R')
source('source_cellbender/function_deviance.R')

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

seurat_7dpi <- make_seurat_somi_function(sce_7dpi.norm, hvg_7dpi)
seurat_7dpi <- seurat_cluster_somi_function(seurat_7dpi, 20)
UMAPPlot(seurat_7dpi, label=TRUE, label.size=6)

DimPlot(seurat_7dpi, group.by = "sample")
DimPlot(seurat_7dpi, group.by = "day")

somi_featureplot(seurat_7dpi, "Ptprc")
somi_featureplot(seurat_7dpi, "Epcam")
somi_featureplot(seurat_7dpi, "Ppbp")
somi_featureplot(seurat_7dpi, "Pf4")
somi_featureplot(seurat_7dpi, "Msln")
somi_featureplot(seurat_7dpi, "Hba-a2")
somi_featureplot(seurat_7dpi, "Hbb-bt")
somi_violinplot(seurat_7dpi, "Hbb-bt")

somi_featureplot(seurat_7dpi, "Jade3")
somi_violinplot(seurat_7dpi, "Jade3", "seurat_clusters")

somi_featureplot(seurat_7dpi, "Cd3d")
somi_featureplot(seurat_7dpi, "Ms4a1")
somi_featureplot(seurat_7dpi, "Nkg7")
somi_featureplot(seurat_7dpi, "Cd300e")
somi_featureplot(seurat_7dpi, "Jchain")

somi_featureplot(seurat_7dpi, "Lgr5")
somi_featureplot(seurat_7dpi, "Olfm4")
somi_featureplot(seurat_7dpi, "Mki67")
somi_featureplot(seurat_7dpi, "Apoa4")
somi_featureplot(seurat_7dpi, "Neurod1")
somi_featureplot(seurat_7dpi, "Chga")
somi_featureplot(seurat_7dpi, "Dclk1")
somi_featureplot(seurat_7dpi, "Clu")
somi_featureplot(seurat_7dpi, "Anxa1")
somi_featureplot(seurat_7dpi, "Defa17")

seurat_7dpi$seq_date = seurat_7dpi$day
seurat_7dpi$seq_date = gsub("D0", "batch1", seurat_7dpi$seq_date)
seurat_7dpi$seq_date = gsub("D3", "batch1", seurat_7dpi$seq_date)
seurat_7dpi$seq_date = gsub("D7", "batch2", seurat_7dpi$seq_date)

UMAPPlot(seurat_7dpi, group.by = "seq_date")

save(seurat_7dpi, file=paste0("data2/Seurat/seurat_7dpi.RData"))
# load(paste0("data2/Seurat/seurat_7dpi.RData"))

# batch correction with harmony
library(harmony)
set.seed(10)
seurat_7dpi_bc <- RunHarmony(seurat_7dpi, "seq_date", plot_convergence=TRUE, assay.use = "RNA")
seurat_7dpi_bc <- seurat_cluster_somi_function(seurat_7dpi_bc, 20, "harmony")

DimPlot(seurat_7dpi_bc, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_bc, group.by = "day", split.by = "day")
FeaturePlot(seurat_7dpi_bc, "log10_total_counts")
FeaturePlot(seurat_7dpi_bc, "pct_counts_MT")

save(seurat_7dpi_bc, file="data2/Seurat/seurat_7dpi_bc.RData")
# load("data2/Seurat/seurat_7dpi_bc.RData")

somi_featureplot(seurat_7dpi_bc, "Ptprc")
somi_featureplot(seurat_7dpi_bc, "Epcam")
somi_featureplot(seurat_7dpi_bc, "Ppbp")
somi_featureplot(seurat_7dpi_bc, "Pf4")
somi_featureplot(seurat_7dpi_bc, "Msln")
somi_featureplot(seurat_7dpi_bc, "Pecam1")
somi_featureplot(seurat_7dpi_bc, "Vwf")
somi_featureplot(seurat_7dpi_bc, "Cdh5")
somi_featureplot(seurat_7dpi_bc, "Hba-a2")
somi_featureplot(seurat_7dpi_bc, "Hbb-bt")
somi_violinplot(seurat_7dpi_bc, "Hbb-bt")

somi_featureplot(seurat_7dpi_bc, "Lgr5")
somi_featureplot(seurat_7dpi_bc, "Olfm4")
somi_featureplot(seurat_7dpi_bc, "Mki67")
somi_featureplot(seurat_7dpi_bc, "Apoa4")
somi_featureplot(seurat_7dpi_bc, "Neurod1")
somi_featureplot(seurat_7dpi_bc, "Chga")
somi_featureplot(seurat_7dpi_bc, "Dclk1")
somi_featureplot(seurat_7dpi_bc, "Clu")
somi_featureplot(seurat_7dpi_bc, "Anxa1")
somi_featureplot(seurat_7dpi_bc, "Defa17")
somi_featureplot(seurat_7dpi_bc, "Muc2")
somi_featureplot(seurat_7dpi_bc, "Agr2")
somi_featureplot(seurat_7dpi_bc, "Dll1")
somi_featureplot(seurat_7dpi_bc, "Atoh1")

FeaturePlot(seurat_7dpi_bc, "log10_total_counts")

marker <- FindMarkers(seurat_7dpi_bc, ident.1 = c(19), logfc.threshold = 1.5)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC > 0), n=40)

somi_featureplot(seurat_7dpi_bc, "Malat1")

# deviance residuals for correcting UMI bias
seurat_7dpi_dev <- run_deviance_residual(seurat_7dpi, batch="seq_date")
seurat_7dpi_dev <- FindNeighbors(seurat_7dpi_dev, reduction="harmony", assay="deviance", dims=1:10)
seurat_7dpi_dev <- FindClusters(seurat_7dpi_dev, resolution=0.5)

save(seurat_7dpi_dev, file="data2/Seurat/seurat_7dpi_dev.RData")
# load("data2/Seurat/seurat_7dpi_dev.RData")

DimPlot(seurat_7dpi_dev, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_dev, group.by = "day")
DimPlot(seurat_7dpi_dev, group.by = "seq_date")
DimPlot(seurat_7dpi_dev, group.by = "sample")
DimPlot(seurat_7dpi_dev, group.by = "condition")

FeaturePlot(seurat_7dpi_dev, "log10_total_counts")
FeaturePlot(seurat_7dpi_dev, "pct_counts_MT")
VlnPlot(seurat_7dpi_dev, "log10_total_counts")
VlnPlot(seurat_7dpi_dev, "pct_counts_MT", pt.size = 0)

somi_featureplot(seurat_7dpi_dev, "Ptprc")
somi_featureplot(seurat_7dpi_dev, "Epcam")
somi_featureplot(seurat_7dpi_dev, "Ppbp")
somi_featureplot(seurat_7dpi_dev, "Pf4")
somi_featureplot(seurat_7dpi_dev, "Msln")
somi_featureplot(seurat_7dpi_dev, "Pecam1")
somi_featureplot(seurat_7dpi_dev, "Vwf")
somi_featureplot(seurat_7dpi_dev, "Cdh5")
somi_featureplot(seurat_7dpi_dev, "Hba-a2")
somi_featureplot(seurat_7dpi_dev, "Hbb-bt")

somi_featureplot(seurat_7dpi_dev, "Lgr5")
somi_featureplot(seurat_7dpi_dev, "Olfm4")
somi_featureplot(seurat_7dpi_dev, "Mki67")
somi_featureplot(seurat_7dpi_dev, "Apoa4")
somi_featureplot(seurat_7dpi_dev, "Neurod1")
somi_featureplot(seurat_7dpi_dev, "Chga")
somi_featureplot(seurat_7dpi_dev, "Dclk1")
somi_featureplot(seurat_7dpi_dev, "Clu")
somi_featureplot(seurat_7dpi_dev, "Anxa1")
somi_featureplot(seurat_7dpi_dev, "Defa17")
somi_featureplot(seurat_7dpi_dev, "Muc2")
somi_featureplot(seurat_7dpi_dev, "Agr2")
somi_featureplot(seurat_7dpi_dev, "Dll1")
somi_featureplot(seurat_7dpi_dev, "Atoh1")

# highMT = c(7,16,23,25)
# Immune = c(0,10,11,18,19)
# 
# marker <- FindMarkers(seurat_7dpi_dev, ident.1 = c(), logfc.threshold = 1)
# marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
# head(subset(marker, avg_logFC > 0), n=40)
# 
# somi_featureplot(seurat_7dpi_dev, "Jchain")