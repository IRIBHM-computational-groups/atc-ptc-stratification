setwd("/mnt/iribhm/ngs/ST/article/R/")
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))

##Fig 7a, UMAP with CDH11, FAP, PRRX1, ZEB1
mes.genes <- c("CDH11","FAP","ZEB1","PRRX1")
lapply(mes.genes, FeaturePlotSaveAll, merged.experiments, order = T)

##Fig. 7c, barplot of % cells expressing TOP2A, use medium res. annotation
top2a.counts <- tapply(merged.experiments@assays$SCT@counts["TOP2A",], merged.experiments$cell.types.medium.resolution, function(z){sum(z > 0)})
cluster.counts <- table(merged.experiments$cell.types.medium.resolution)
top2a.percentages <- top2a.counts / cluster.counts
pdf(paste0(fig.dir, "barplot_top2apos_all.pdf"), height = 4)
par(mar=c(11,4,1,1)+.1)
barplot(sort(top2a.percentages), ylab = "% TOP2A positive", las = 2)
dev.off()

##Fig. 7b, heatmap showing % expressing cells
Idents(merged.experiments) <- merged.experiments$cell.types.medium.resolution
low.genes <- FindMarkers(merged.experiments, ident.1 = "Epith TSHR low", ident.2 = "Epith TSHR+ TPO+", only.pos = T)
top10 <- rownames(low.genes)[order(low.genes$avg_log2FC, decreasing = T)][1:10]
other.markers <- c("FN1", "SLC5A5", "TPO", "TG", "TSHR")
markers <- c(top10, other.markers)
ctypes.keep <- c("Epith TSHR+ TPO+ NIS+","Epith TSHR+ TPO+","Epith TSHR+ TPO low","Epith TSHR+ TPO-","Epith TSHR low","Epith TSHR-","Fibroblast")
percentage.matrix <- matrix(data = 0, nrow = length(ctypes.keep), ncol = length(markers))
rownames(percentage.matrix) <- ctypes.keep
colnames(percentage.matrix) <- markers
for (r in ctypes.keep){
    for (c in markers){
        counts <- merged.experiments@assays$RNA@counts[c, merged.experiments$cell.types.medium.resolution %in% r]
        number.expressing <- sum(counts > 0)
        percentage.matrix[r, c] <- number.expressing / length(counts)
    }
}
pdf(paste0(fig.dir, "heatmap_percentage_pos_nuclei_epith.pdf"), height = 3)
Heatmap(percentage.matrix, name = "% nuclei expressing markers", col = divergingPalBin(), 
                    show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

m <- FindAllMarkers(merged.experiments, only.pos = T)
saveRDS(m, file = paste0(metatc,"markers_article_annotation.Rds"))
m <- readRDS(paste0(metatc,"markers_article_annotation.Rds"))
m[m$cluster == "EMT1", c(2,5,6,7)][1:50,]
m12 <- FindMarkers(merged.experiments, ident.1 = "EMT1", ident.2 = "EMT2")
m12 <- m12[order(m12$avg_log2FC),]
m12[grep("^PAX", rownames(m12)),]
FeaturePlot(merged.experiments, features = "CDH1")
epi.genes <- c("CDH1","")
PrintEntrez("MEG8")
sum(m.samples[c("Myeloid","T-prolif","T-cell","B-naive","T-regs","Fibroblast","Dendritic",
                "Endothelial","Pericyte","B-cell"),"ATC5"])
sum(m.samples[c("Myeloid","T-prolif","T-cell","B-naive","T-regs","Fibroblast","Dendritic",
                "Endothelial","Pericyte","B-cell"),"ATC6"])
sum(m.samples[c("Myeloid","T-prolif","T-cell","B-naive","T-regs","Fibroblast","Dendritic",
                "Endothelial","Pericyte","B-cell"),"ATC3"])
sum(m.samples[c("Myeloid","T-prolif","T-cell","B-naive","T-regs","Fibroblast","Dendritic",
                "Endothelial","Pericyte","B-cell"),"ATC4"])

##Fig. 7e, spRNA-seq, mean MT genes, CD68
s.genes <- c("CD68","percent.mt")
sp.filenames <- list.files(sp.dir, pattern = "S[0-9]_seurat.Rds")
sp.sectors <- gsub("_seurat.Rds", "", sp.filenames)
sp.bundle <- lapply(sp.filenames, function(z){readRDS(paste0(sp.dir, z))})
names(sp.bundle) <- sp.sectors
for (s.g in s.genes){
    for (i in sp.sectors){
        SpFtPlotSave(s.object = sp.bundle[[i]], s.gene = s.g)
    }
}
sp.cors <- lapply(sp.bundle, function(z){cor.test(x = z@assays$SCT@counts["CD68",], y = z$percent.mt, method = "spearman")})
sp.rhos <- unlist(lapply(sp.cors, function(z){return(z$estimate)}))
sp.ptc <- grepl("PTC",names(sp.rhos))*1
pdf(paste0(fig.dir, "barplot_cor_cd68_mt_ST.pdf"), height = 6)
par(mar=c(6,4,1,1)+.1)
barplot(sp.rhos, las = 2, ylab = "Spearman's rho between CD68 and %MT genes across ST spots", col = sp.ptc)
dev.off()

genes1 <- VariableFeatures(sp.bundle$PTC2S5)
genes2 <- genes1
for (s in sp.bundle){
    genes2 <- genes2[genes2 %in% rownames(s)]
}
genes2 <- c(genes2, "CD68")
sp.all.rhos <- lapply(sp.bundle, function(sector){
    sp.cors <- lapply(genes2, function(gene){cor.test(x = sector@assays$SCT@counts[gene,], y = sector$percent.mt, method = "spearman")})
    sp.rhos <- unlist(lapply(sp.cors, function(z){return(z$estimate)}))
    return(sp.rhos)
})

sp.ptc <- grepl("PTC",names(sp.rhos))*1 + 2
pdf(paste0(fig.dir, "boxlot_cor_var_mt_ST.pdf"), height = 9, width = 9)
par(mar=c(6,4,1,1)+.1)
boxplot(sp.all.rhos, las = 2, ylab = "Spearman's rho between common variable genes and %MT genes across ST spots", col = sp.ptc, main = "Abline at -0.1")
abline(h = -0.1)
dev.off()
