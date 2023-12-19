setwd("/mnt/iribhm/ngs/ST/article/R/")
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))
## FIG2BD - SUPFIG2ABC - HEATMAPS FOR OUR DATASET
HeatmapSave(merged.experiments, group.var = "cell.types.low.resolution")
HeatmapSave(merged.experiments, group.var = "cell.types.medium.resolution")
HeatmapSave(merged.experiments, group.var = "cell.types.high.resolution")
HeatmapSave(merged.experiments, group.var = "cell.types.low.resolution", norm.by = "cluster")
HeatmapSave(merged.experiments, group.var = "cell.types.medium.resolution", norm.by = "cluster")
HeatmapSave(merged.experiments, group.var = "cell.types.high.resolution", norm.by = "cluster")
HeatmapSave(merged.experiments, group.var = "cell.types.low.resolution", norm.by = "cluster", do.cor = TRUE)
HeatmapSave(merged.experiments, group.var = "cell.types.medium.resolution", norm.by = "cluster", do.cor = TRUE)

## FIG2C - PERCENTAGE CELLS OF EPITHELIAL ORIGIN
order.vec <- rev(c("PTC1","PTC4","PTC9","PTC2","PTC5","ATC3A","ATC3B","ATC2","PTC8","ATC1","ATC4A","ATC4B","PTC3","PTC6","PTC7"))
ct.med <- merged.experiments$cell.types.medium.resolution
ct.high <- merged.experiments$cell.types.high.resolution
ct.unique <- unique(ct.med)
ct.high.unique <- unique(ct.high)
is.epi.orig <- c(ct.med %in% grep("Epith", ct.unique, value = T))
is.epi.orig[ct.high %in% grep("High",ct.high.unique,value=T)[-c(2,5)]] <- TRUE
Sample <- merged.experiments$Sample
prop <- sapply(unique(Sample), function(z){sum(is.epi.orig[Sample == z]) / sum(Sample == z)})
prop <- prop[order.vec]
pdf(paste0(fig.dir, "barplot_epithelial_origin_all.pdf"), width = 2)
barplot(prop, las = 2, horiz = TRUE, xlab = "% epithelial")
dev.off()

## FIG2E - CORRELOGRAM FOR TCGA
tcga <- readRDS(paste0(tcga.dir, "percentages_16clusters_tcga.Rds"))
rownames(tcga)[c(1,7)] <- c("Epith TSHR+ TPO+", "Epith TSHR+ TPO-")
annot.tcga <- readRDS(paste0(tcga.dir, "tcga_annotations.Rds"))
annot.tcga <- annot.tcga[rownames(annot.tcga) %in% colnames(tcga),]
tcga <- tcga[,annot.tcga$BRAFV600E_RAS %in% "BRAF_V600E"]
annot.tcga <- annot.tcga[rownames(annot.tcga) %in% colnames(tcga),]
mcor <- cor(t(tcga), method = "spearman")
pdf(paste0(fig.dir, paste0("heatmap_corTRUE_", gsub("\\.","_","cell.types.low.resolution"), "_clusterNorm_tcga"), ".pdf"))
Heatmap(mcor, col = divergingPal(), show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()
pdf(paste0(fig.dir, paste0("heatmap_props_cell.types.low.resolution_tcga"), ".pdf"), width = 20)
ciberHeatmap()
dev.off()
