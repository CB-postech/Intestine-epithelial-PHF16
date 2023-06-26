.libPaths(Sys.getenv("R_LIBS_INTESTINE"))
.libPaths()

setwd("D:/Project/INTESTINE")

library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

plotdir = "plots_cellbender/"
rdatadir = "data_cellbender/"

for(i in 1:12){
  load(paste0(rdatadir, "sce.qc/", "cb_sce_e", i, ".qc.RData"))
}

for(i in 1:12){
  sce <- eval(parse(text = paste0("cb_sce_e", i, ".qc")))
  colnames(sce) = paste0("e", i, "_", colnames(sce))
  
  assign(paste0("cb_sce_e", i, ".qc"), sce)
}
rm(sce)

merged_counts <- cbind(counts(cb_sce_e1.qc), counts(cb_sce_e2.qc), counts(cb_sce_e3.qc), counts(cb_sce_e4.qc),
                       counts(cb_sce_e5.qc), counts(cb_sce_e6.qc), counts(cb_sce_e7.qc), counts(cb_sce_e8.qc),
                       counts(cb_sce_e9.qc), counts(cb_sce_e10.qc), counts(cb_sce_e11.qc), counts(cb_sce_e12.qc))
# logmc <- log2(merged_counts+1)

col_data = rbind(colData(cb_sce_e1.qc), colData(cb_sce_e2.qc), colData(cb_sce_e3.qc), colData(cb_sce_e4.qc),
                 colData(cb_sce_e5.qc), colData(cb_sce_e6.qc), colData(cb_sce_e7.qc), colData(cb_sce_e8.qc),
                 colData(cb_sce_e9.qc), colData(cb_sce_e10.qc), colData(cb_sce_e11.qc), colData(cb_sce_e12.qc))
head(col_data)
cell_info = data.frame(cell = colnames(merged_counts),
                       sample = c(rep("Sample1", ncol(cb_sce_e1.qc)), rep("Sample2", ncol(cb_sce_e2.qc)),
                                  rep("Sample3", ncol(cb_sce_e3.qc)), rep("Sample4", ncol(cb_sce_e4.qc)),
                                  rep("Sample5", ncol(cb_sce_e5.qc)), rep("Sample6", ncol(cb_sce_e6.qc)),
                                  rep("Sample7", ncol(cb_sce_e7.qc)), rep("Sample8", ncol(cb_sce_e8.qc)),
                                  rep("Sample9", ncol(cb_sce_e9.qc)), rep("Sample10", ncol(cb_sce_e10.qc)),
                                  rep("Sample11", ncol(cb_sce_e11.qc)), rep("Sample12", ncol(cb_sce_e12.qc))),
                       type = c(rep("WT", ncol(cb_sce_e1.qc) + ncol(cb_sce_e2.qc)),
                                rep("KO", ncol(cb_sce_e3.qc) + ncol(cb_sce_e4.qc)),
                                rep("WT", ncol(cb_sce_e5.qc) + ncol(cb_sce_e6.qc)),
                                rep("KO", ncol(cb_sce_e7.qc) + ncol(cb_sce_e8.qc)),
                                rep("WT", ncol(cb_sce_e9.qc) + ncol(cb_sce_e10.qc)),
                                rep("KO", ncol(cb_sce_e11.qc) + ncol(cb_sce_e12.qc))),
                       irr = c(rep("Nor", ncol(cb_sce_e1.qc) + ncol(cb_sce_e2.qc) + ncol(cb_sce_e3.qc) + ncol(cb_sce_e4.qc)),
                               rep("Irr", ncol(cb_sce_e5.qc) + ncol(cb_sce_e6.qc) + ncol(cb_sce_e7.qc) + ncol(cb_sce_e8.qc) +
                                     ncol(cb_sce_e9.qc) + ncol(cb_sce_e10.qc) + ncol(cb_sce_e11.qc) + ncol(cb_sce_e12.qc))),
                       day = c(rep("D0", ncol(cb_sce_e1.qc) + ncol(cb_sce_e2.qc) + ncol(cb_sce_e3.qc) + ncol(cb_sce_e4.qc)),
                               rep("D3", ncol(cb_sce_e5.qc) + ncol(cb_sce_e6.qc) + ncol(cb_sce_e7.qc) + ncol(cb_sce_e8.qc)),
                               rep("D7", ncol(cb_sce_e9.qc) + ncol(cb_sce_e10.qc) + ncol(cb_sce_e11.qc) + ncol(cb_sce_e12.qc))),
                       seq_date = c(rep("batch1", ncol(cb_sce_e1.qc) + ncol(cb_sce_e2.qc) + ncol(cb_sce_e3.qc) + ncol(cb_sce_e4.qc) +
                                          ncol(cb_sce_e5.qc) + ncol(cb_sce_e6.qc) + ncol(cb_sce_e7.qc) + ncol(cb_sce_e8.qc)),
                                    rep("batch2", ncol(cb_sce_e9.qc) + ncol(cb_sce_e10.qc) + ncol(cb_sce_e11.qc) + ncol(cb_sce_e12.qc)))
)
cell_info$condition = paste(cell_info$day, cell_info$type, cell_info$irr, sep="_")

colData = cbind(cell_info, col_data)
rm(cell_info)
rm(col_data)

sce_total <- SingleCellExperiment(assays = list(counts = merged_counts),
                                 colData = colData)

rm(merged_counts)
rm(colData)

keep_feature = rowSums(counts(sce_total) > 0) > 0
table(keep_feature) # TRUE 35379 FALSE 20042
sce_total = sce_total[keep_feature, ]
dim(sce_total) # 35379 genes, 60875 cells
rm(keep_feature)

#
norm_somi_function <- function(sce.qc){
  
  objectname = deparse(substitute(sce.qc))
  
  clusters <- quickCluster(sce.qc)
  sce.qc <- computeSumFactors(sce.qc, clusters=clusters)
  
  print(paste0("summary of size factors of ", objectname))
  print(summary(sizeFactors(sce.qc)))
  
  sce.qc <- normalize(sce.qc)
  return(sce.qc)
}

sce_total.norm <- norm_somi_function(sce_total)

rownames(sce_total.norm@assays@data$logcounts) = rownames(sce_total.norm@assays@data$counts)
colnames(sce_total.norm@assays@data$logcounts) = colnames(sce_total.norm@assays@data$counts)

save(sce_total, file=paste0(rdatadir, "sce/sce_total.RData"))
save(sce_total.norm, file=paste0(rdatadir, "sce/sce_total.norm.RData"))

hvg_function <- function(sce, FDRcut=0.05, biocut=0.01){
  var.fit <- trendVar(sce,
                      method="spline", parametric=TRUE, use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)

  hvg <- var.out[which(var.out$FDR < FDRcut & var.out$bio > biocut),]
  print(nrow(hvg)) # 1467 genes
  return(hvg)
}
hvg_total <- hvg_function(sce_total.norm, biocut = 0.05) # 2273 genes

save(hvg_total, file=paste0(rdatadir, "sce/hvg_total.RData"))
