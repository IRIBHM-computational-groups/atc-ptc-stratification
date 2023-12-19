## VARIABLES
set.seed(42)
md.setting <- 0.01
dimensionality <- 1:30

## DIRECTORIES
article.dir <- "/mnt/iribhm/ngs/ST/article/"
output.dir <- paste0(article.dir, "output/")
input.dir <- paste0(article.dir, "input/")
fig.dir <- paste0(output.dir, "figures/")
table.dir <- paste0(output.dir, "tables/")
tcga.dir <- paste0(input.dir, "TCGA/")
sn.dir <- paste0(input.dir, "snRNA-seq/")
sp.dir <- paste0(input.dir, "spRNA-seq/")
misc.dir <- paste0(input.dir, "misc/")

## COLOR PALLETTES
samples.vec <- c("ATC1","ATC2","ATC3A","ATC3B","ATC4A","ATC4B",paste0("PTC",1:9))
samples.pal <- glasbey(1:length(samples.vec) + 1)
names(samples.pal) <- samples.vec
exp.vec <- c(paste0("ATC",1:6,"SN1"), paste0("PTC",1:9,"SN1"))
exp.pal <- glasbey(1:length(samples.vec) + 1)
names(exp.pal) <- exp.vec
sequentialPal <- function(n1 = 0, n2 = 1){colorRamp2(c(n1,n2), c("white","red"))}
divergingPal <- function(n1 = -1, n2 = 0, n3 = 1){colorRamp2(c(n1,n2,n3), c("blue","white","red"))}
divergingPalBin <- function(n1 = 0, n2 = 1){colorRamp2(c(n1,n2), c("blue","red"))}

## ANNOTATION COLOR PALETTTES
ct.low.vec <- c("Epith TSHR-","Myeloid","T-cell","B-cell","Fibroblast","Dendritic","Endothelial","High mitochdr","Epith TSHR+ TPO-","Pericyte","Epith TSHR+ TPO+")
ct.low.pal <- glasbey(1:length(ct.low.vec) + 1)
names(ct.low.pal) <- ct.low.vec
ct.medium.vec <- c("B-cell","B-naive","Dendritic","Endothelial","Epith TSHR low","Epith TSHR+ TPO low","Epith TSHR+ TPO+","Epith TSHR+ TPO+ NIS+","Epith TSHR+ TPO-","Epith TSHR-","Fibroblast","High mitochdr","Myeloid","Pericyte","T-cell","T-regs")
ct.medium.pal <- glasbey(1:length(ct.medium.vec) + 1)
names(ct.medium.pal) <- ct.medium.vec
ct.high.vec <- c("B IgA","B IgG PRDM1 high","B IgG PRDM1 low","B naive","Dblt T cell","Dblt epith","Dblt macroph","Endoth arterial","Endoth lymphatic","Endoth vascular","Endoth vascular capillary","Endoth vascular tip","Endoth venous","Epith TSHR low","Epith TSHR+ TPO low","Epith TSHR+ TPO+","Epith TSHR+ TPO+ NIS+","Epith TSHR+ TPO-","Epith TSHR+ TPO- FN1","Epith TSHR+ TPO- PTC2","Epith TSHR+ TPO- PTC4","Epith TSHR+ TPO- PTC9","Epith TSHR-","Fibrobl ATC","Fibrobl ATC1 2","Fibrobl PTC","High mitochr EMT1","High mitochr EMT2","High mitochr EMT2 ATC5","High mitochr EMT2 ATC6","High mitochr T cell","High mitochr epithelial ATC3","High mitochr epithelial ATC5","High mitochr epithelial PTC","High mitochr macrophage","Macroph ATC","Macroph ATC1","Macroph ATC2","Macroph PTC","Macroph prolif","Mast cell","Mature cDC","Pericyte","T CD4","T CD4 PTC2","T CD4 exhausted","T CD8","T CD8 ATC","T CD8 ATC1","T CD8 prolif","T reg","cDC1","moDC")
ct.high.pal <- glasbey(1:length(ct.high.vec) + 1)
names(ct.high.pal) <- ct.high.vec