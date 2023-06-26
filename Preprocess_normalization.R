.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

# load packages
library(scater)
library(scran)

for(i in 1:8){load(paste0("D:/Project/BOWEL/data/QC_Sce/bowel_sce", i, ".qc.RData"))}
for(i in 9:12){load(paste0("D:/Project/INTESTINE/data2/Sce/bowel_sce", i, ".qc.RData"))}

colnames(bowel_sce9.qc) = paste("e9", colnames(bowel_sce9.qc), sep="_")
colnames(bowel_sce10.qc) = paste("e10", colnames(bowel_sce10.qc), sep="_")
colnames(bowel_sce11.qc) = paste("e11", colnames(bowel_sce11.qc), sep="_")
colnames(bowel_sce12.qc) = paste("e12", colnames(bowel_sce12.qc), sep="_")

merged_counts <- cbind(counts(bowel_sce1.qc), counts(bowel_sce2.qc), counts(bowel_sce3.qc), counts(bowel_sce4.qc),
                       counts(bowel_sce5.qc), counts(bowel_sce6.qc), counts(bowel_sce7.qc), counts(bowel_sce8.qc),
                       counts(bowel_sce9.qc), counts(bowel_sce10.qc), counts(bowel_sce11.qc), counts(bowel_sce12.qc))
# logmc <- log2(merged_counts+1)

col_data = rbind(colData(bowel_sce1.qc), colData(bowel_sce2.qc), colData(bowel_sce3.qc), colData(bowel_sce4.qc),
                 colData(bowel_sce5.qc), colData(bowel_sce6.qc), colData(bowel_sce7.qc), colData(bowel_sce8.qc),
                 colData(bowel_sce9.qc), colData(bowel_sce10.qc), colData(bowel_sce11.qc), colData(bowel_sce12.qc))
col_data = col_data[,-1]
cell_info = data.frame(cell = colnames(merged_counts),
                       sample = c(rep("Sample1", ncol(bowel_sce1.qc)), rep("Sample2", ncol(bowel_sce2.qc)),
                                  rep("Sample3", ncol(bowel_sce3.qc)), rep("Sample4", ncol(bowel_sce4.qc)),
                                  rep("Sample5", ncol(bowel_sce5.qc)), rep("Sample6", ncol(bowel_sce6.qc)),
                                  rep("Sample7", ncol(bowel_sce7.qc)), rep("Sample8", ncol(bowel_sce8.qc)),
                                  rep("Sample9", ncol(bowel_sce9.qc)), rep("Sample10", ncol(bowel_sce10.qc)),
                                  rep("Sample11", ncol(bowel_sce11.qc)), rep("Sample12", ncol(bowel_sce12.qc))),
                       type = c(rep("WT", ncol(bowel_sce1.qc) + ncol(bowel_sce2.qc)),
                                rep("KO", ncol(bowel_sce3.qc) + ncol(bowel_sce4.qc)),
                                rep("WT", ncol(bowel_sce5.qc) + ncol(bowel_sce6.qc)),
                                rep("KO", ncol(bowel_sce7.qc) + ncol(bowel_sce8.qc)),
                                rep("WT", ncol(bowel_sce9.qc) + ncol(bowel_sce10.qc)),
                                rep("KO", ncol(bowel_sce11.qc) + ncol(bowel_sce12.qc))),
                       irr = c(rep("Nor", ncol(bowel_sce1.qc) + ncol(bowel_sce2.qc) + ncol(bowel_sce3.qc) + ncol(bowel_sce4.qc)),
                               rep("Irr", ncol(bowel_sce5.qc) + ncol(bowel_sce6.qc) + ncol(bowel_sce7.qc) + ncol(bowel_sce8.qc) +
                                     ncol(bowel_sce9.qc) + ncol(bowel_sce10.qc) + ncol(bowel_sce11.qc) + ncol(bowel_sce12.qc))),
                       day = c(rep("D0", ncol(bowel_sce1.qc) + ncol(bowel_sce2.qc) + ncol(bowel_sce3.qc) + ncol(bowel_sce4.qc)),
                               rep("D3", ncol(bowel_sce5.qc) + ncol(bowel_sce6.qc) + ncol(bowel_sce7.qc) + ncol(bowel_sce8.qc)),
                               rep("D7", ncol(bowel_sce9.qc) + ncol(bowel_sce10.qc) + ncol(bowel_sce11.qc) + ncol(bowel_sce12.qc)))
)
cell_info$condition = paste(cell_info$day, cell_info$type, cell_info$irr, sep="_")
colData = Matrix::rBind(colData(bowel_sce1.qc), colData(bowel_sce2.qc), colData(bowel_sce3.qc), colData(bowel_sce4.qc),
                        colData(bowel_sce5.qc), colData(bowel_sce6.qc), colData(bowel_sce7.qc), colData(bowel_sce8.qc),
                        colData(bowel_sce9.qc), colData(bowel_sce10.qc), colData(bowel_sce11.qc), colData(bowel_sce12.qc))
colData = cbind(cell_info, colData)

sce_7dpi <- SingleCellExperiment(assays = list(counts = merged_counts),
                                colData = colData)

rm(merged_counts)
rm(cell_info)

keep_feature = rowSums(counts(sce_7dpi) > 0) > 0
table(keep_feature) # TRUE 35533 FALSE 19888
sce_7dpi = sce_7dpi[keep_feature, ]
dim(sce_7dpi) # 35533 genes, 59652 cells
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

sce_7dpi.norm <- norm_somi_function(sce_7dpi)

save(sce_7dpi, file="data2/Sce/sce_7dpi.RData")
save(sce_7dpi.norm, file="data2/Sce/sce_7dpi.norm.RData")
