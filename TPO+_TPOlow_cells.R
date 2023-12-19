setwd("/mnt/iribhm/ngs/ST/article/R/")
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")

merged.3.6 <- readRDS(paste0(sn.dir, "merged.3.6.Rds"))
samples.vec.norm <- c("PTC3N", "PTC6N")
samples.pal.norm <- glasbey(c(21,20))
names(samples.pal.norm) <- samples.vec.norm
samples.pal <- c(samples.pal, samples.pal.norm)
DimPlotSave(merged.3.6, group.var = "Sample")
lapply(c("TSHR","TPO"), FeaturePlotSave, merged.3.6)

s.genes <- c("TPO","TSHR")
PTC6S1 <- readRDS(paste0(sp.dir, "PTC6S1_seurat.Rds"))
for (s.gene in s.genes){
    SpFtPlotSave(PTC6S1, s.gene = s.gene)
}
