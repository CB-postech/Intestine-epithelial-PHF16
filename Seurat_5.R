.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(ggplot2)
library(RColorBrewer)

source('source_cellbender/function_visualization.R')
source('source_cellbender/function_violin_split.R')
source('source_cellbender/function_violin_fig.R')
source('source_cellbender/function_seurat.R')
source('source_cellbender/function_deviance.R')

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "data2/"
plotdir = "plots2/"

seurat_7dpi_5 <- run_seurat_somi_function(seurat_7dpi_4, cells = remove.cells, invert = TRUE,
                                          batchcorrection = TRUE, batch = "seq_date")

DimPlot(seurat_7dpi_5, label=TRUE, label.size=8)
DimPlot(seurat_7dpi_5, group.by = "day", split.by = "day")
FeaturePlot(seurat_7dpi_5, "total_counts")
FeaturePlot(seurat_7dpi_5, "log10_total_counts")
FeaturePlot(seurat_7dpi_5, "pct_counts_MT")
VlnPlot(seurat_7dpi_5, "log10_total_counts", pt.size=0)
VlnPlot(seurat_7dpi_5, "log10_total_features_by_counts", pt.size=0)
VlnPlot(seurat_7dpi_5, "pct_counts_MT", pt.size=0)

save(seurat_7dpi_5, file="data2/Seurat/seurat_7dpi_5.RData")
# load("data2/Seurat/seurat_7dpi_5.RData")

somi_featureplot(seurat_7dpi_5, "Ptprc")
somi_featureplot(seurat_7dpi_5, "Epcam")
somi_featureplot(seurat_7dpi_5, "Ppbp")
somi_featureplot(seurat_7dpi_5, "Pf4")
somi_featureplot(seurat_7dpi_5, "Msln")
somi_featureplot(seurat_7dpi_5, "Pecam1")
somi_featureplot(seurat_7dpi_5, "Vwf")
somi_featureplot(seurat_7dpi_5, "Cdh5")
somi_featureplot(seurat_7dpi_5, "Hba-a2")
somi_featureplot(seurat_7dpi_5, "Hbb-bt")

somi_featureplot(seurat_7dpi_5, "Lgr5")
somi_featureplot(seurat_7dpi_5, "Slc12a2")
somi_featureplot(seurat_7dpi_5, "Olfm4")
somi_featureplot(seurat_7dpi_5, "Ascl2")
somi_featureplot(seurat_7dpi_5, "Gkn3")
somi_featureplot(seurat_7dpi_5, "Mcm6")
somi_featureplot(seurat_7dpi_5, "Mcm5")
somi_featureplot(seurat_7dpi_5, "Cdk4")
somi_featureplot(seurat_7dpi_5, "Mki67")
somi_featureplot(seurat_7dpi_5, "Arg2")
somi_featureplot(seurat_7dpi_5, "Il18")
somi_featureplot(seurat_7dpi_5, "Car4")
somi_featureplot(seurat_7dpi_5, "Ccl25")
somi_featureplot(seurat_7dpi_5, "Apoa4")
somi_featureplot(seurat_7dpi_5, "Apoa1")
somi_featureplot(seurat_7dpi_5, "Alpi")
somi_featureplot(seurat_7dpi_5, "Fabp1")
somi_featureplot(seurat_7dpi_5, "Dll1")
somi_featureplot(seurat_7dpi_5, "Atoh1")
somi_featureplot(seurat_7dpi_5, "Neurog3")
somi_featureplot(seurat_7dpi_5, "Neurod1")
somi_featureplot(seurat_7dpi_5, "Chga")
somi_featureplot(seurat_7dpi_5, "Chgb")
somi_featureplot(seurat_7dpi_5, "Dclk1")
somi_featureplot(seurat_7dpi_5, "Trpm5")
somi_featureplot(seurat_7dpi_5, "Gfi1")
somi_featureplot(seurat_7dpi_5, "Spdef")
somi_featureplot(seurat_7dpi_5, "Tff3")
somi_featureplot(seurat_7dpi_5, "Muc2")
somi_featureplot(seurat_7dpi_5, "Agr2")
somi_featureplot(seurat_7dpi_5, "Clu")
somi_featureplot(seurat_7dpi_5, "Ly6a")
somi_featureplot(seurat_7dpi_5, "Anxa1")
somi_featureplot(seurat_7dpi_5, "Anxa3")
somi_featureplot(seurat_7dpi_5, "Cd44")
somi_featureplot(seurat_7dpi_5, "Ccnd1")
somi_featureplot(seurat_7dpi_5, "Ccnd2")
somi_featureplot(seurat_7dpi_5, "Bmi1")
somi_featureplot(seurat_7dpi_5, "Lrig1")
somi_featureplot(seurat_7dpi_5, "Hopx")
somi_featureplot(seurat_7dpi_5, "Tert")
somi_featureplot(seurat_7dpi_5, "Lyz1")
somi_featureplot(seurat_7dpi_5, "Defa17")
somi_featureplot(seurat_7dpi_5, "Defa22")
somi_featureplot(seurat_7dpi_5, "Defa24")

marker <- FindMarkers(seurat_7dpi_5, ident.1 = c(13), logfc.threshold = 0.5)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC > 0), n=40)

somi_featureplot(seurat_7dpi_5, "Ly6a")

marker <- FindMarkers(seurat_7dpi_5, ident.1 = c(13), ident.2 = c(9,25), logfc.threshold = 0.5)
marker$symbol = ensemblGenes[rownames(marker), "external_gene_name"]
head(subset(marker, avg_logFC > 0), n=40)

somi_featureplot(seurat_7dpi_5, "Msln")

#
load("D:/Project/INTESTINE/data/Seurat.3.dev/intestine_seurat.3.dev.RData")
head(colnames(intestine_seurat.3.dev))

seurat_7dpi_5$orig.celltype = "Undefined"
common_cells = intersect(colnames(seurat_7dpi_5), colnames(intestine_seurat.3.dev))
seurat_7dpi_5$orig.celltype[common_cells] = intestine_seurat.3.dev$celltype[common_cells]
seurat_7dpi_5$orig.celltype = as.factor(seurat_7dpi_5$orig.celltype)
seurat_7dpi_5$orig.celltype = relevel(seurat_7dpi_5$orig.celltype, ref="Undefined")
DimPlot(seurat_7dpi_5[,colnames(intestine_seurat.3.dev)], group.by = "orig.celltype")

#
CBC = c(8,10)
TA = c(1,3,12,14,23)
TC = c(17)
EE = c(18,20,24)
EC = c(4,5,6,11,15,21,22,26)
GC = c(0,2,16,19)
PC = c(7)
RevSC = c(9,13,25)

seurat_7dpi_5$celltype = NA

types = c("CBC", "TA", "TC", "EE", "EC", "GC", "PC", "RevSC")
for(t in types){
  cl = eval(parse(text = t))
  cells = seurat_7dpi_5[, seurat_7dpi_5$seurat_clusters %in% cl] %>% colnames()
  seurat_7dpi_5$celltype[cells] = rep(t, length(cells))
}

UMAPPlot(seurat_7dpi_5, group.by = "celltype", label=TRUE, label.size=8)

#
DimPlot(seurat_7dpi_5, group.by = "celltype")
DimPlot(seurat_7dpi_5, group.by = "seurat_clusters", label=TRUE, label.size=6)

seurat_7dpi_5$celltype2 = seurat_7dpi_5$celltype
seurat_7dpi_5$celltype2[colnames(seurat_7dpi_5[,seurat_7dpi_5$seurat_clusters %in% c(9,25)])] = "RevSC_1"
seurat_7dpi_5$celltype2[colnames(seurat_7dpi_5[,seurat_7dpi_5$seurat_clusters %in% c(13)])] = "RevSC_2"

DimPlot(seurat_7dpi_5, group.by = "celltype2")

#
seurat_7dpi_5$celltype3 = as.character(seurat_7dpi_5$celltype2)
cells = colnames(seurat_7dpi_5[, seurat_7dpi_5$celltype == "RevSC"])
seurat_7dpi_5$celltype3[cells] = seurat_7dpi_5$seurat_clusters[cells]
seurat_7dpi_5$celltype3 = gsub("26", "RevSC_1", seurat_7dpi_5$celltype3)
seurat_7dpi_5$celltype3 = gsub("10", "RevSC_2", seurat_7dpi_5$celltype3)
seurat_7dpi_5$celltype3 = gsub("14", "RevSC_3", seurat_7dpi_5$celltype3)

UMAPPlot(seurat_7dpi_5, group.by = "celltype3")

seurat_D0$celltype3 = seurat_7dpi_5$celltype3[colnames(seurat_D0)]
seurat_D3$celltype3 = seurat_7dpi_5$celltype3[colnames(seurat_D3)]
seurat_D7$celltype3 = seurat_7dpi_5$celltype3[colnames(seurat_D7)]

# save(seurat_7dpi_5, file="data2/Seurat/seurat_7dpi_5.RData")
# save(seurat_D0, file="data2/Seurat/seurat_D0.RData")
# save(seurat_D3, file="data2/Seurat/seurat_D3.RData")
# save(seurat_D7, file="data2/Seurat/seurat_D7.RData")
# load("data2/Seurat/seurat_7dpi_5.RData")

#
Idents(seurat_7dpi_5) = seurat_7dpi_5$celltype2
RevSC_1_marker <- FindMarkers(seurat_7dpi_5, ident.1 = "RevSC_1", logfc.threshold = 0, min.pct = 0)
RevSC_1_marker$symbol = ensemblGenes[rownames(RevSC_1_marker), "external_gene_name"]
save(RevSC_1_marker, file=paste0(rdatadir, "Seurat/RevSC_1_marker.RData"))

RevSC_compare <- FindMarkers(seurat_7dpi_5, ident.1 = "RevSC_1", ident.2 = "RevSC_2", logfc.threshold = 0, min.pct = 0)
RevSC_compare$symbol = ensemblGenes[rownames(RevSC_compare), "external_gene_name"]
save(RevSC_compare, file=paste0(rdatadir, "Seurat/RevSC_compare.RData"))

#
Idents(seurat_7dpi_5) = seurat_7dpi_5$celltype
celltype_markers <- FindAllMarkers(seurat_7dpi_5, logfc.threshold = 0, min.pct = 0)
celltype_markers$symbol = ensemblGenes[celltype_markers$gene, "external_gene_name"]
save(celltype_markers, file=paste0(rdatadir, "Seurat/celltype_markers.RData"))

#
Idents(seurat_7dpi_5) = seurat_7dpi_5$seurat_clusters
cluster_markers <- FindAllMarkers(seurat_7dpi_5, logfc.threshold = 0, min.pct = 0)
cluster_markers$symbol = ensemblGenes[cluster_markers$gene, "external_gene_name"]
save(cluster_markers, file=paste0(rdatadir, "Seurat/cluster_markers.RData"))

#
Idents(seurat_7dpi_5) = seurat_7dpi_5$seurat_clusters
cluster9 <- FindMarkers(seurat_7dpi_5, ident.1 = c(9), ident.2 = c(13, 25), logfc.threshold = 0, min.pct = 0)
cluster13 <- FindMarkers(seurat_7dpi_5, ident.1 = c(13), ident.2 = c(9, 25), logfc.threshold = 0, min.pct = 0)
cluster25 <- FindMarkers(seurat_7dpi_5, ident.1 = c(25), ident.2 = c(13, 9), logfc.threshold = 0, min.pct = 0)
cluster9$symbol = ensemblGenes[rownames(cluster9), "external_gene_name"]
cluster13$symbol = ensemblGenes[rownames(cluster13), "external_gene_name"]
cluster25$symbol = ensemblGenes[rownames(cluster25), "external_gene_name"]

head(subset(cluster25, avg_logFC > 0), n=20)

#
somi_violinplot.fig.test(seurat_7dpi_5, "log10_total_counts", group.by = "celltype", col = celltype_color, y.mul=1)
somi_violinplot.fig.test(seurat_7dpi_5, "log10_total_features_by_counts", group.by = "celltype", col = celltype_color, y.mul=1)

#
UMAPPlot(seurat_7dpi_5, group.by = "seurat_clusters", label=TRUE) # revSC_1 = 9,25 / revSC_2 = 13
load(paste0(rdatadir, "Seurat/RevSC_1_marker.RData"))
load(paste0(rdatadir, "Seurat/cluster_markers.RData"))
load(paste0(rdatadir, "Seurat/RevSC_compare.RData"))

head(RevSC_compare)
head(RevSC_1_marker)
head(subset(cluster_markers, cluster == 25), n=40)

somi_featureplot(seurat_7dpi_5, "Ly6c1")
UMAPPlot(seurat_7dpi_5, group.by = "seurat_clusters", label=TRUE)

#
# head(RevSC_compare)
# 
# comparedf <- RevSC_compare
# comparedf$log10P = -log10(comparedf$p_val_adj)
# comparedf$log10P[comparedf$log10P > 300] = 300
# 
# somi_volcanoplot(comparedf, top_n=0, color_up=, color_down=)

# barcodes
# for velocity
unique(seurat_7dpi_5$sample)
for(s in unique(seurat_7dpi_5$sample)){
  barcodes <- colnames(seurat_7dpi_5[, seurat_7dpi_5$sample %in% s])
  barcodes = gsub("[A-Za-z0-9]+_", "", barcodes)
  barcodes = paste0(barcodes, "-1") %>% as.data.frame()
  write.table(barcodes, file=paste0(rdatadir, "velocyto/barcodes_", s, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

#
markers <- FindMarkers(seurat_7dpi_5, ident.1=c(13), ident.2=c(9,25), logfc.threshold = 0, min.pct = 0)

Idents(seurat_7dpi_5) = seurat_7dpi_5$type
markers_type <- FindMarkers(seurat_7dpi_5, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, min.pct = 0)
markers_type$symbol = ensemblGenes[rownames(markers_type), "external_gene_name"]
save(markers_type, file=paste0(rdatadir, "Seurat/markers_type.RData"))

#
umap <- seurat_7dpi_5@reductions$umap@cell.embeddings
write.csv(umap, file=paste0(rdatadir, "Seurat/umap_embed_Total.csv"), row.names = TRUE)
harmony <- seurat_7dpi_5@reductions$harmony@cell.embeddings
write.csv(harmony, file=paste0(rdatadir, "Seurat/harmony_embed.csv"), row.names = TRUE)
meta <- seurat_7dpi_5@meta.data
write.csv(meta, file=paste0(rdatadir, "Seurat/meta_data_Total.csv"), row.names = TRUE)

# cell cycle scoring
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes = convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes = convertHumanGeneList(cc.genes$g2m.genes)

s.genes.id = ensemblGenes[ensemblGenes$external_gene_name %in% m.s.genes,]$ensembl_gene_id
g2m.genes.id = ensemblGenes[ensemblGenes$external_gene_name %in% m.g2m.genes,]$ensembl_gene_id

seurat_7dpi_5 <- CellCycleScoring(seurat_7dpi_5,
                           s.features = s.genes.id,
                           g2m.features = g2m.genes.id)
DimPlot(seurat_7dpi_5, group.by = "Phase")
seurat_7dpi_5$CC.Difference <- seurat_7dpi_5$S.Score - seurat_7dpi_5$G2M.Score

df <- seurat_7dpi_5@meta.data
df$sample = factor(df$sample, levels = rev(paste0("Sample", 1:12)))
df$condition = factor(df$condition, levels = rev(names(condition_color)))
df$celltype = factor(df$celltype, levels = rev(names(celltype_color)))
df$celltype = factor(df$celltype, levels = rev(c("CBC", "TA", "EC", "EE", "GC", "TC", "RevSC", "PC")))
# df$celltype2 = factor(df$celltype2, levels = rev(names(celltype2_color)))
# df$celltype3 = factor(df$celltype3, levels = rev(names(celltype3_color)))
df$Phase = factor(df$Phase, levels = c("G1", "G2M", "S"))

df$rep = df$sample
df$rep = gsub("^Sample1$|^Sample3$|^Sample5$|^Sample7$|^Sample9$|^Sample11$", "Rep1", df$rep)
df$rep = gsub("^Sample2$|^Sample4$|^Sample6$|^Sample8$|^Sample10$|^Sample12$", "Rep2", df$rep)

library(dplyr)
df2 <- df %>%
  group_by(day, type, celltype, Phase, .drop=FALSE) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n))
g_bar <- ggplot(df2, aes(freq, type, fill=Phase)) +
  geom_bar(position='fill', stat='identity') +
  # scale_fill_manual(values = unlist(celltype_color)) +
  facet_wrap(~day) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")
g_bar

g_bar <- ggplot(df2, aes(n, celltype, fill=Phase)) +
  geom_bar(stat='identity') +
  # scale_fill_manual(values = unlist(celltype_color)) +
  facet_wrap(~day+type, nrow=3) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")
g_bar

print_plot_to_ppt(g_bar, path, width=8, height=10)

#
df2 = df[df$day == "D7",]
df2 <- df2 %>%
  group_by(celltype, type, Phase, .drop=FALSE) %>%
  summarise(n=n()) %>%
  mutate(freq = n/sum(n))

g_bar <- ggplot(df2, aes(freq, celltype, fill=Phase)) +
  geom_bar(position='fill', stat='identity') +
  # scale_fill_manual(values = unlist(celltype_color)) +
  facet_wrap(~type) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")
g_bar

g_bar <- ggplot(df2, aes(n, celltype, fill=Phase)) +
  geom_bar(stat='identity') +
  # scale_fill_manual(values = unlist(celltype_color)) +
  facet_wrap(~type) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")
g_bar

table(seurat_7dpi_5$type, seurat_7dpi_5$Phase, seurat_7dpi_5$day, seurat_7dpi_5$celltype)

ggplot(df2, aes(n, celltype, fill=Phase)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = ) +
  facet_wrap(~type) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none")

#
Idents(seurat_7dpi_5) = seurat_7dpi_5$day
subseurat <- subset(seurat_7dpi_5, idents = "D3", invert=TRUE)

v <- somi_violinplot_split(subseurat, "Mki67", "day", "type", "celltype", col=type_color, y.mul = 1, stat=TRUE)
print_plot_to_ppt(v, path, width=8, height=3)

# TEAD 1-4
somi_violinplot_split(seurat_7dpi_5, "Tead1", "day", "type", "celltype", col = type_color, y.mul=1, stat=TRUE)

genes = c("Tead1", "Tead2", "Tead3", "Tead4")
for(gene in genes){
  v <- somi_violinplot_split(seurat_7dpi_5, gene, "day", "type", "celltype", col = type_color, y.mul=1, stat=TRUE)
  ggsave(paste0("plots2/violin_", gene, ".png"), plot=v, width=15, height=5)
}

#
wt_color = condition_color[grep("WT", names(condition_color))]
names(wt_color) = levels(seurat_7dpi_5$day)
somi_violinplot.fig.test(seurat_7dpi_5[, seurat_7dpi_5$type == "WT"], "Jade3", group.by = "day", split.by = "celltype", stat=TRUE, y.mul = 1.2, col = wt_color)

somi_featureplot(seurat_7dpi_5[, seurat_7dpi_5$type == "WT" & seurat_7dpi_5$day == "D7"], "Jade3")
somi_featureplot(seurat_7dpi_5[, seurat_7dpi_5$type == "WT"], "Neurog3")
