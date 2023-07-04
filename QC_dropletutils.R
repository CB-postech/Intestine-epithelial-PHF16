.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

library(DropletUtils)

dropletutils_somi_function <- function(sample, sampleinfo, lower=100, FDR_cutoff=0.05){
  
  sampleinfo = as.character(sampleinfo)
  sample = deparse(substitute(sample))
  
  dir.name <- paste0("D:/Project/INTESTINE/rawdata/", sample, "/")
  list.files(dir.name)
  
  sce <- read10xCounts(dir.name)
  
  br.out <- barcodeRanks(counts(sce))
  
  set.seed(2019)
  e.out <- emptyDrops(counts(sce), lower = lower)  ## Cells that have UMI counts lower than 100 are empty cells.
  table(Sig=e.out$FDR <= FDR_cutoff, Limited=e.out$Limited)
  is.cell <- e.out$FDR <= FDR_cutoff
  
  print(sum(is.cell, na.rm=TRUE))
  
  print(table(br.out$rank == sum(is.cell, na.rm=TRUE)))
  
  png(paste0("plots2/dropletutils/dropletutils_", sample, ".png"))
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  # abline(h=br.out[br.out$rank == sum(is.cell, na.rm=TRUE),]$total, col="red", lty=2)
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  # legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"),
  #        legend=c("knee", "inflection", "FDR_0.05"))
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))
  dev.off()
  
  colnames(sce) = colData(sce)$Barcode
  sce <- sce[,which(e.out$FDR <= FDR_cutoff)]
  return(sce)
}

sce9 <- dropletutils_somi_function(WT1, "WT 1 mouse crypt cell w irradiation DAY7", lower=500)
sce10 <- dropletutils_somi_function(WT2, "WT 2 mouse crypt cell w irradiation DAY7", lower=500)
sce11 <- dropletutils_somi_function(KO1, "Phf16-/y 1 mouse crypt cell w irradiation DAY7", lower=500)
sce12 <- dropletutils_somi_function(KO2, "Phf16-/y 2 mouse crypt cell w irradiation DAY7", lower=500)

for(i in 9:12){
  save(list=paste("sce", i, sep=""),
       file=paste0("data2/Raw_Sce/sce", i, ".RData"))
}
