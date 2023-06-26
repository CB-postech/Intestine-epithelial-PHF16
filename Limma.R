.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(Seurat)
library(limma)
library(dplyr)
library(scran)
library(scater)
library(SingleCellExperiment)

rdatadir = "data2/"
plotdir = "plots2/"

load(paste0(rdatadir, "Seurat/seurat_7dpi_5.RData"))

load("D:/somi_function/ensemblGenes_mmusculus_2021-01-06.RData")

celltypes = levels(as.factor(seurat_7dpi_5$celltype))
for(i in celltypes){
  cells = seurat_7dpi_5[,seurat_7dpi_5$celltype == i] %>% colnames
  print(i)
  tab = table(seurat_7dpi_5[,cells]$condition)
  print(tab)
  print(length(tab))
}

DefaultAssay(seurat_7dpi_5) = "RNA"

norm_somi_function <- function(sce.qc){
  
  objectname = deparse(substitute(sce.qc))
  
  clusters <- quickCluster(sce.qc)
  sce.qc <- computeSumFactors(sce.qc, clusters=clusters)
  
  print(paste0("summary of size factors of ", objectname))
  print(summary(sizeFactors(sce.qc)))
  
  sce.qc <- normalize(sce.qc)
  return(sce.qc)
}
run_limma <- function(seurat, subset, group){
  
  use.cells = seurat[,seurat[[subset]] == group] %>% colnames()
  seurat = subset(seurat, cells = use.cells)
  
  sce <- as.SingleCellExperiment(seurat)
  
  min.size=100
  if(length(use.cells) < 100){
    min.size=length(use.cells)
  }
  sce <- norm_somi_function(sce)
  print(paste0("dimension is ", nrow(counts(sce)), " X ", ncol(counts(sce)), "."))
  
  expr <- assay(sce, "logcounts")
  
  keep <- rowSums(as.matrix(expr)) > 0
  expr = expr[keep,]
  
  condition <- as.factor(seurat$condition)
  sample <- as.factor(seurat$sample)
  
  print("base data generation is done.")
  
  design <- model.matrix(~0+condition)
  colnames(design) = levels(condition)
  
  print(length(levels(condition)))
  if(sum(grep("D0", levels(condition))) == 0){
    cont.matrix <- makeContrasts(
      "KOeffect_D7" = D7_KO_Irr - D7_WT_Irr,
      "KOeffect_D3" = D3_KO_Irr - D3_WT_Irr,
      "Ireffect_KO_D7_D3" = D7_KO_Irr - D3_KO_Irr,
      "Ireffect_WT_D7_D3" = D7_WT_Irr - D3_WT_Irr,
      "Interact_D7_D3" = (D7_KO_Irr - D7_WT_Irr) - (D3_KO_Irr - D3_WT_Irr),
      levels = design
    )
  }else if(sum(grep("D3", levels(condition))) == 0){
    cont.matrix <- makeContrasts(
      "KOeffect_D7" = D7_KO_Irr - D7_WT_Irr,
      "KOeffect_D0" = D0_KO_Nor - D0_WT_Nor,
      "Ireffect_KO_D7_D0" = D7_KO_Irr - D0_KO_Nor,
      "Ireffect_WT_D7_D0" = D7_WT_Irr - D0_WT_Nor,
      "Interact_D7_D0" = (D7_KO_Irr - D7_WT_Irr) - (D0_KO_Nor - D0_WT_Nor),
      levels = design
    )
  }else{
    cont.matrix <- makeContrasts(
      "KOeffect_D7" = D7_KO_Irr - D7_WT_Irr,
      "KOeffect_D3" = D3_KO_Irr - D3_WT_Irr,
      "KOeffect_D0" = D0_KO_Nor - D0_WT_Nor,
      "Ireffect_KO_D3_D0" = D3_KO_Irr - D0_KO_Nor,
      "Ireffect_WT_D3_D0" = D3_WT_Irr - D0_WT_Nor,
      "Ireffect_KO_D7_D0" = D7_KO_Irr - D0_KO_Nor,
      "Ireffect_WT_D7_D0" = D7_WT_Irr - D0_WT_Nor,
      "Ireffect_KO_D7_D3" = D7_KO_Irr - D3_KO_Irr,
      "Ireffect_WT_D7_D3" = D7_WT_Irr - D3_WT_Irr,
      "Interact_D7_D0" = (D7_KO_Irr - D7_WT_Irr) - (D0_KO_Nor - D0_WT_Nor),
      "Interact_D7_D3" = (D7_KO_Irr - D7_WT_Irr) - (D3_KO_Irr - D3_WT_Irr),
      "Interact_D3_D0" = (D3_KO_Irr - D3_WT_Irr) - (D0_KO_Nor - D0_WT_Nor),
      levels = design
    )
  }
  
  corfit <- duplicateCorrelation(expr,
                                 design = design,
                                 block = sample)
  cor = corfit$consensus.correlation
  
  print("Batch effect correlation is calculated. Start fitting.")
  
  fit <- lmFit(expr,
               design = design,
               block = sample,
               correlation = cor)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
  
  return(fit2)
}

for(i in celltypes){
  fit <- run_limma(seurat_7dpi_5, 'celltype', i)
  assign(paste0("fit_", i), fit)
}
for(i in celltypes){ 
  save(list = paste0("fit_", i), file=paste0(rdatadir, "Seurat/fit_", i, ".RData")) 
}
for(i in celltypes){load(paste0(rdatadir, "Seurat/fit_", i, ".RData"))}

fit_PC_lin <- run_limma(seurat_7dpi_5, 'celltype4', "PC_lineage")

for(i in celltypes){
  fit <- eval(parse(text = paste0("fit_", i)))
  
  deglist_full <- lapply(colnames(fit$coefficients), function(x) topTable(fit, coef=x, number=Inf, p.value=1, lfc=0))
  deglist <- lapply(colnames(fit$coefficients), function(x) topTable(fit, coef=x, number=Inf, p.value=0.05, lfc=0.1))
  
  deglist_full <- lapply(deglist_full, function(x) transform(x, symbol = ensemblGenes[rownames(x), "external_gene_name"]))
  deglist <- lapply(deglist, function(x) transform(x, symbol = ensemblGenes[rownames(x), "external_gene_name"]))
  
  names(deglist_full) = colnames(fit$coefficients)
  names(deglist) = colnames(fit$coefficients)
  
  assign(paste0("deglist_full_", i), deglist_full)
  assign(paste0("deglist_", i), deglist)
  
  rm(fit)
  rm(deglist_full)
  rm(deglist)
}

coefs = colnames(fit_EC$coefficients)
degMat <- matrix(nrow = length(celltypes), ncol = length(coefs), dimnames = list(celltypes, coefs))

for(i in celltypes){
  deglist <- eval(parse(text = paste0("deglist_", i)))
  degMat[i, names(deglist)] = sapply(deglist, function(x) nrow(x))
}

celltype_num <- table(seurat_7dpi_5$celltype)

degRatio <- apply(degMat, 2, function(x) x / celltype_num)

save(degMat, file=paste0(rdatadir, "Seurat/degMat.RData"))
save(degRatio, file=paste0(rdatadir, "Seurat/degRatio.RData"))

#
degMat2 <- degMat[,grep("KOeffect", colnames(degMat))]
degMat2

cellNum <- matrix(data = c(table(seurat_D7$celltype)[rownames(degMat2)],
                           table(seurat_D3$celltype)[rownames(degMat2)],
                           table(seurat_D0$celltype)[rownames(degMat2)]), ncol=3, dimnames = list(rownames(degMat2), colnames(degMat2)))

degMat2 / cellNum

degRatio2 <- degMat2 / cellNum
degRatio2["TA", "KOeffect_D3"] = NA
degRatio2["TC", "KOeffect_D3"] = NA

colors = c("lightgray", colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(255))
colors = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(255)

use_terms = colnames(degRatio)[grep("KOeffect", colnames(degRatio))]
dg <- pheatmap(degRatio2[c("RevSC", "EC", "CBC", "TA", "TC", "PC", "GC", "EE"),sort(colnames(degRatio2))], cluster_rows = FALSE, cluster_cols = FALSE, color = colors)
pheatmap(degMat2[c("RevSC", "EC", "CBC", "TA", "TC", "PC", "GC", "EE"),sort(colnames(degRatio2))], cluster_rows = FALSE, cluster_cols = FALSE, color = colors)

t(round(degRatio[,rev(grep("KOeffect", colnames(degRatio)))], digit=3))

#
for(i in celltypes){ 
  save(list = paste0("deglist_", i), file=paste0(rdatadir, "Seurat/deglist_", i, ".RData")) 
  save(list = paste0("deglist_full_", i), file=paste0(rdatadir, "Seurat/deglist_full_", i, ".RData")) 
}

library(openxlsx)
make_xlsx <- function(list, listname, option=NULL){
  
  workbook <- createWorkbook(listname)
  
  for(i in 1:length(list)){
    if(nrow(list[[i]]) == 0){
      print("pass")
      next
    }
    
    sheetname = names(list)[i]
    
    addWorksheet(workbook, sheetname)
    
    if(sum(is.null(option)) == 1){
      writeDataTable(workbook, sheetname, list[[i]], keepNA=TRUE, tableStyle = "none", rowNames = TRUE)
    }else if(option == "up"){
      degDF = subset(list[[i]], logFC > 0)
      degDF = degDF[, c("logFC", "AveExpr", "adj.P.Val", "symbol")]
      writeDataTable(workbook, sheetname, degDF, keepNA=TRUE, tableStyle = "none", rowNames = TRUE)
    }else if(option == "down"){
      degDF = subset(list[[i]], logFC < 0)
      degDF = degDF[, c("logFC", "AveExpr", "adj.P.Val", "symbol")]
      writeDataTable(workbook, sheetname, degDF, keepNA=TRUE, tableStyle = "none", rowNames = TRUE)
    }else{print("option should be null or up or down")}
  }
  
  if(sum(is.null(option)) == 0){
    if(option == "up"){ listname = paste0(listname, "_up") }
    if(option == "down"){ listname = paste0(listname, "_down") }
  }
  
  saveWorkbook(workbook, file=paste0(rdatadir, "Seurat/deg_csv/", listname, ".xlsx"), overwrite = TRUE)
}
make_multiple_xlsx <- function(..., option=NULL){
  names <- as.character(substitute(list(...)))[-1L]
  List = list(...)
  names(List) = names
  
  for(j in 1:length(List)){
    make_xlsx(List[[j]], names(List)[j], option=option)
  }
}
make_multiple_xlsx(deglist_CBC, deglist_RevSC, deglist_EC, deglist_EE, deglist_GC, deglist_PC, deglist_TA, deglist_TC)
make_multiple_xlsx(deglist_full_CBC, deglist_full_RevSC, deglist_full_EC, deglist_full_EE, deglist_full_GC, deglist_full_PC, deglist_full_TA, deglist_full_TC)
make_multiple_xlsx(deglist_CBC, deglist_RevSC, deglist_EC, deglist_EE, deglist_GC, deglist_PC, deglist_TA, deglist_TC, option="up")
make_multiple_xlsx(deglist_CBC, deglist_RevSC, deglist_EC, deglist_EE, deglist_GC, deglist_PC, deglist_TA, deglist_TC, option="down")
