.libPaths(Sys.getenv("R_LIBS_INTESTINE_4"))

rdatadir = "data2/"
plotdir = "plots2/"

library(oligo)

celfiles <- list.files("Q:/Public/mouse/intestine/GSE45155/", full.names = TRUE)
celfiles = celfiles[-grep("sh", celfiles)]
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
targets$level = rep(c("Control", "YAPTAZ_KD"), 3)
targets

TS = targets$level
TS <- factor(TS, levels = c("Control", "YAPTAZ_KD"))
design <- model.matrix(~0+TS)
colnames(design) = levels(TS)

fit <- lmFit(eset, design)
fit

cont.matrix <- makeContrasts(
  YAPTAZ_KD_vs_Control = YAPTAZ_KD - Control,
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
output_YAPTAZ_KDVSControl <- topTable(fit_2, coef=1, number=Inf)
output_YAPTAZ_KDVSControl$probe = rownames(output_YAPTAZ_KDVSControl)
head(output_YAPTAZ_KDVSControl)

#
library(biomaRt)

# ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#                    dataset="mmusculus_gene_ensembl", 
#                    host="www.ensembl.org")
View(listAttributes(ensembl))
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'affy_mogene_1_0_st_v1'), 
                      mart=ensembl)
head(ensemblGenes)

#
symbols = sapply(rownames(output_YAPTAZ_KDVSControl), function(x) {ensemblGenes[ensemblGenes$affy_mogene_1_0_st_v1 %in% x,]$external_gene_name[1] } )
output_YAPTAZ_KDVSControl$symbols = symbols

KD_output <- output_YAPTAZ_KDVSControl[!is.na(output_YAPTAZ_KDVSControl$symbols),]
KD_output <- KD_output[-which(duplicated(KD_output$symbols)),]
rownames(KD_output) = KD_output$symbols
head(KD_output)
save(KD_output, file=paste0(rdatadir, "Microarray/KD_output.RData"))

KD_up = subset(KD_output, logFC > 0)
KD_down = subset(KD_output, logFC < 0)

head(KD_up)

nrow(subset(Repair_up, logFC > 1 & adj.P.Val < 0.05))