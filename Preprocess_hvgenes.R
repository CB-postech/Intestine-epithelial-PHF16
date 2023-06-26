.libPaths(Sys.getenv("R_LIBS_INTESTINE"))

# load packages
library(scater)
library(scran)

var.fit <- trendVar(sce_7dpi.norm,
                    method="spline", parametric=TRUE, use.spikes=FALSE)
var.out <- decomposeVar(sce_7dpi.norm, var.fit)
FDRcut = 0.05
biocut = 0.1
hvg <- var.out[which(var.out$FDR < FDRcut & var.out$bio > biocut),]
nrow(hvg) # 1467 genes

rm(var.fit)
rm(var.out)
rm(o)
rm(biocut)
rm(FDRcut)

hvg_7dpi <- hvg
rm(hvg)

save(hvg_7dpi, file="data2/Sce/hvg_7dpi.RData")
