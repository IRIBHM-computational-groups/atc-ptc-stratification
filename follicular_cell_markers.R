setwd(script.dir) ## Set working directory to the folder containing the scripts of this repository
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))

## EPITHELIAL PLOTS
lapply(c("SLC5A5","TPO","TG","TSHR"), FeaturePlotSaveAll, merged.experiments)
keep <- merged.experiments$cell.types.medium.resolution %in% c("Epith TSHR low","Epith TSHR+ TPO low","Epith TSHR+ TPO+",
                                                               "Epith TSHR+ TPO+ NIS+","Epith TSHR+ TPO-","Epith TSHR-")
keep <- colnames(merged.experiments)[keep]
merged.experiments@project.name <- "all_epith"
DimPlotSave(merged.experiments, group.var = "cell.types.medium.resolution", cells = keep)
table(merged.experiments$Sample[merged.experiments$cell.types.medium.resolution == "Epith TSHR+ TPO+ NIS+"]) / sum(merged.experiments$cell.types.medium.resolution == "Epith TSHR+ TPO+ NIS+")
table(merged.experiments$Sample[merged.experiments$cell.types.medium.resolution == "Epith TSHR+ TPO+"]) / sum(merged.experiments$cell.types.medium.resolution == "Epith TSHR+ TPO+")
table(merged.experiments$cell.types.medium.resolution[merged.experiments$Sample == "PTC6"]) / sum(merged.experiments$Sample == "PTC6")

FeaturePlotSaveAll(merged.experiments, ft.var = "SLC5A5", order = T)

# ORGANOID PLOTS
orga <- readRDS(paste0(sn.dir, "HS_ORGANOIDSC1_seurat_mt26.Rds"))
orga@project.name <- "organoids"
orga$Cell.type <- orga$Cell.Types
DimPlotSave(orga, group.var = "Cell.type")
lapply(c("SLC5A5","TPO","TG","TSHR"), FeaturePlotSaveOrganoids, orga)

origin = merged.experiments$Sample
origin[grepl("ATC", origin)] = "ATC"
origin[grepl("PTC", origin)] = "PTC"
table(origin[merged.experiments$cell.types.medium.resolution == "High mitochdr"]) / sum(merged.experiments$cell.types.medium.resolution == "High mitochdr")
