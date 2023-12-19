setwd("/mnt/iribhm/ngs/ST/article/R/")
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))
gs <- list(REACTOME = gmtPathways(paste0(misc.dir, "c2.cp.reactome.v7.5.1.symbols.gmt")),
          BIOCARTA = gmtPathways(paste0(misc.dir, "c2.cp.biocarta.v7.5.1.symbols.gmt")))

## MYELOID CELL ABUNDANCE
is.myelo.orig <- c(merged.experiments$cell.types.medium.resolution %in% "Myeloid")
is.myelo.orig[merged.experiments$cell.types.high.resolution %in% c("High mitochr macrophage")] <- TRUE
Sample <- merged.experiments$Sample
pp <- sort(sapply(unique(Sample), function(z){sum(is.myelo.orig[Sample == z]) / sum(Sample == z)})) * 100
pdf(paste0(fig.dir, "barplot_myeloid_all.pdf"))
par(mar=c(5,8,4,1)+.1)
pp <- sort(m.samples["Myeloid",] * 100)
colors <- rep("gray", times = length(pp))
colors[grepl("ATC",names(pp))] <- "red"
barplot(pp, horiz = T, las = 2, xlab = "% of myeloid cells", col = colors)
dev.off()

Idents(merged.experiments) <- merged.experiments$cell.types.low.resolution
myeloid <- subset(merged.experiments, idents = c("Myeloid","Dendritic"))
myeloid <- basicSeurat(myeloid)
myeloid@project.name <- "myeloid"
DimPlotSave(myeloid, group.var = "Sample")
DimPlotSave(myeloid, group.var = "cell.types.high.resolution")
FeaturePlotSave(myeloid, ft.var = "HIF1A")

Idents(merged.experiments) <- merged.experiments$Sample
gg <- VlnPlot(merged.experiments, features = "HIF1A", group.by = "cell.types.medium.resolution", split.by = "Sample", idents = c("ATC1","ATC2"), cols = ct.medium.pal)
savePlotsLegend(gg, "vlnplot_HIF1A_myeloid")

## FIG4D FGSEA
myeloid <- SetIdent(object = myeloid, value = "SCT_snn_res.0.1")
m <- FindMarkers(myeloid, ident.1 = 1, ident.2 = 2)
atc.clusters <- myeloid$SCT_snn_res.0.1 %in% c(1,2)
exp.mat <- GetAssayData(object = myeloid, slot = "counts", assay = "SCT")[VariableFeatures(myeloid), atc.clusters]
exp.mat <- as.matrix(exp.mat)
exp.mat <- exp.mat[rowSums(exp.mat) > 0,]
categs <- as.numeric(myeloid$SCT_snn_res.0.1[atc.clusters]) - 1
fgsea.res <- fgseaLabel(pathways = gs$REACTOME, mat = exp.mat, labels = categs, nperm = 1000)

Idents(merged.experiments) <- merged.experiments$Sample
markers <- FindMarkers(merged.experiments, ident.1 = "ATC1", logfc.threshold = 0, only.pos = F)
markers.keep <- markers[markers$p_val_adj < 0.05,]
preranks <- markers.keep$avg_log2FC
names(preranks) <- rownames(markers.keep)
ranks <- rank(preranks)
fgsea.res <- fgsea(pathways = gs$BIOCARTA, stats = preranks)
selected <- "BIOCARTA_P53HYPOXIA_PATHWAY"
in.gset <- gs$BIOCARTA[[selected]][gs$BIOCARTA[[selected]] %in% names(preranks)]
ordered.genes <- in.gset[order(preranks[in.gset],decreasing = T)]
escore <- round(fgsea.res[fgsea.res$pathway == selected, "ES"], digits = 3)
pval <- round(fgsea.res[fgsea.res$pathway == selected, "pval"], digits = 3)
gg <- plotEnrichment(pathway = gs$BIOCARTA[[selected]], stats = preranks) + labs(title = paste(selected, "ES =", escore, "Pval =", pval), 
                                                                                                 y = "Enrichment score", x = "Rank") + 
                                                                                            theme(axis.text = element_text(size=15),
                                                                                                axis.title = element_text(size=15))
ggsave(filename = paste0(fig.dir, "gseaplot_hypoxia_all.pdf"), plot = gg)

sort.markers <- markers.keep[order(markers.keep$avg_log2FC, decreasing = T),]
n <- which(rownames(sort.markers) == "HIF1A") 
n
sort.markers[n,]
fgsea.res[grep("HYPOX",fgsea.res$pathway),]
length(ranks)
