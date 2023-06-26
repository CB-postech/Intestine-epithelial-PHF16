.libPaths(Sys.getenv("R_LIBS_INTESTINE_4"))

rdatadir = "data2/"
plotdir = "plots2/"

library(oligo)

celfiles <- list.files("Q:/Public/mouse/intestine/E-MTAB-5249/", full.names = TRUE)
celfiles

rawData <- read.celfiles(celfiles)

geneCore <- rma(rawData, target="core")
featureData(geneCore) <- getNetAffx(geneCore, "probeset")

#
library(limma)
library(dplyr)

eset <- exprs(geneCore)

colnames(geneCore)
targets = geneCore@phenoData@data
targets$level = c(rep("Homeostasis", 3), rep("Repair", 3))
targets

TS = targets$level
TS <- factor(TS, levels = c("Homeostasis", "Repair"))
design <- model.matrix(~0+TS)
colnames(design) = levels(TS)

fit <- lmFit(eset, design)
fit

cont.matrix <- makeContrasts(
  Repair_vs_Homeo = Repair - Homeostasis,
  levels=design
)
fit_2 <- contrasts.fit(fit, cont.matrix)
fit_2 <- eBayes(fit_2, trend=TRUE)
fit_2
fit_2$coefficients %>% colnames()

#
plotSA(fit_2)

#
coefficients(fit_2) %>% colnames()
output_RepairVSHomeo <- topTable(fit_2, coef=1, number=Inf)
output_RepairVSHomeo$probe = rownames(output_RepairVSHomeo)
head(output_RepairVSHomeo)

#
library(biomaRt)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset="mmusculus_gene_ensembl", 
                   host="www.ensembl.org")
View(listAttributes(ensembl))
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'affy_mogene_2_1_st_v1'), 
                      mart=ensembl)
head(ensemblGenes)

#
symbols = sapply(rownames(output_RepairVSHomeo), function(x) {ensemblGenes[ensemblGenes$affy_mogene_2_1_st_v1 %in% x,]$external_gene_name[1] } )
output_RepairVSHomeo$symbols = symbols
#output_RepairVSHomeo$symbols = gsub("character.*", NA, output_RepairVSHomeo$symbols)

Repair_output <- output_RepairVSHomeo[!is.na(output_RepairVSHomeo$symbols),]
Repair_output <- Repair_output[-which(duplicated(Repair_output$symbols)),]
rownames(Repair_output) = Repair_output$symbols
head(Repair_output)
save(Repair_output, file=paste0(rdatadir, "Microarray/Repair_output.RData"))

Repair_up = subset(Repair_output, logFC > 0)
Repair_down = subset(Repair_output, logFC < 0)

head(Repair_up)

nrow(subset(Repair_up, logFC > 1 & adj.P.Val < 0.05))
