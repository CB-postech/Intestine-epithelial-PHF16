.libPaths(Sys.getenv("R_LIBS_INTESTINE_4"))

rdatadir = "data2/"
plotdir = "plots2/"

library(oligo)
library(limma)
library(dplyr)

idatfiles <- list.files("Q:/Public/mouse/intestine/E-MTAB-5246/", full.names = TRUE)
idatfiles

bgxfiles <- list.files("Q:/Public/mouse/intestine/MouseWG-6_V2_0_R3_11278593_A/", full.names = TRUE)
bgxfiles

rawData <- read.idat(idatfiles, bgxfiles)

rawData$other$Detection <- detectionPValues(rawData)
propexpr(rawData)
y <- neqc(rawData)
y
dim(y)

expressed <- rowSums(y$other$Detection < 0.05) >= 3
length(expressed)

#
# TS = c("Fetal", "Adult", "Fetal", "Adult", "Adult", "Fetal")
TS = c("Adult", "Fetal", "Adult", "Fetal", "Fetal", "Adult")
TS <- factor(TS, levels = c("Fetal", "Adult"))
design <- model.matrix(~0+TS)
colnames(design) = levels(TS)

#
fit <- lmFit(y, design)
fit

cont.matrix <- makeContrasts(
  Fetal_vs_Adult = Fetal - Adult,
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
output_FetalVSAdult <- topTable(fit_2, coef=1, number=Inf)
output_FetalVSAdult$probe = rownames(output_FetalVSAdult)
head(output_FetalVSAdult)

#
Fetal_output = output_FetalVSAdult[!is.na(output_FetalVSAdult$Symbol),]
save(Fetal_output, file=paste0(rdatadir, "Microarray/Fetal_output.RData"))

Fetal_up = subset(Fetal_output, logFC > 0)
Fetal_down = subset(Fetal_output, logFC < 0)

head(Fetal_up)

nrow(subset(Fetal_up, logFC > 2 & adj.P.Val < 0.05))
