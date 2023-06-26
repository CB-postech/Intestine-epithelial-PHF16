.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(SCENT)
library(Seurat)
library(scran)
library(scater)
library(dplyr)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "data2/"
plotdir = "plots2/"

load("data2/Seurat/seurat_7dpi_5.RData")

# load network data
data(net17Jan16)
print(dim(net17Jan16.m))

#
rawmat <- seurat_7dpi_5@assays$RNA@counts
logmat <- log2(rawmat + 1.1)

rm(seurat_7dpi_5)
rm(rawmat)

#
library(biomaRt)

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset="mmusculus_gene_ensembl", 
                   host="www.ensembl.org")
ensembl_ref <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'gene_biotype'), 
                     mart=ensembl, useCache = FALSE)
head(ensembl_ref)

table(rownames(logmat) %in% ensembl_ref$ensembl_gene_id)

use.genes = intersect(rownames(logmat), ensembl_ref$ensembl_gene_id)
ensembl_ref = ensembl_ref[match(use.genes, ensembl_ref$ensembl_gene_id),]

ensembl_ref = ensembl_ref[!is.na(ensembl_ref$entrezgene_id),]

use.genes = ensembl_ref$ensembl_gene_id
table(use.genes %in% rownames(logmat))

logmat = logmat[use.genes, ]

logcd <- as.matrix(logmat)
identical(rownames(logcd), ensembl_ref$ensembl_gene_id)
rownames(logcd) = ensembl_ref$entrezgene_id
rm(logmat)

table(is.na(ensembl_ref$entrezgene_id))
table(is.na(rownames(logcd)))

#
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = x, mart = mouse, attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
  return(genesV2)
}
mouse_to_human_entrezIDs <- convertHumanGeneList(rownames(logcd))
head(mouse_to_human_entrezIDs)

table(rownames(logcd) %in% mouse_to_human_entrezIDs$NCBI.gene..formerly.Entrezgene..ID)

entrez.df <- mouse_to_human_entrezIDs[!is.na(mouse_to_human_entrezIDs$NCBI.gene..formerly.Entrezgene..ID) & 
                                        !is.na(mouse_to_human_entrezIDs$NCBI.gene..formerly.Entrezgene..ID.1),]
use.ids = intersect(rownames(logcd), entrez.df$NCBI.gene..formerly.Entrezgene..ID)
table(use.ids %in% rownames(logcd))

logcd = logcd[rownames(logcd) %in% use.ids,]
entrez.df = entrez.df[match(rownames(logcd), use.ids),]
head(entrez.df)

rownames(logcd) <- entrez.df[,2]
table(!is.na(rownames(logcd)))
dim(logcd)

#
integ.1 <- DoIntegPPI(exp.m = logcd, ppiA.m = net17Jan16.m)
str(integ.1)

sr.o <- CompSRana(integ.1, local=FALSE)

scentScore <- sr.o$SR

save(integ.1, file=paste0(rdatadir, "Seurat/integ.1.RData"))
save(sr.o, file=paste0(rdatadir, "Seurat/sr.o.RData"))
save(scentScore, file=paste0(rdatadir, "Seurat/scentScore.RData"))
