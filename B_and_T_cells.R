setwd("/mnt/iribhm/ngs/ST/article/R/")
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))

## B-cells vertical
b.abund <- computeFractions(merged.experiments, "cell.types.low.resolution")["B-cell",]
pdf(paste0(fig.dir, "barplot_bcell_frac_all.pdf"), width = 2)
par(mar=c(5,8,4,1)+.1)
barplot(sort(b.abund * 100), horiz = T, las = 2, xlab = "% of B cells")
dev.off()
## T-cells vertical
t.abund <- computeFractions(merged.experiments, "cell.types.low.resolution")["T-cell",]
pdf(paste0(fig.dir, "barplot_tcell_frac_all.pdf"), width = 2)
par(mar=c(5,8,4,1)+.1)
barplot(sort(t.abund * 100), horiz = T, las = 2, xlab = "% of T cells")
dev.off()
## B-cells horizontal
b.abund <- computeFractions(merged.experiments, "cell.types.low.resolution")["B-cell",]
pdf(paste0(fig.dir, "barplot_bcell_frac_all.pdf"), height = 3, width = 12)
barplot(sort(b.abund * 100), horiz = F, las = 1, ylab = "% of B cells")
dev.off()
## T-cells horizontal
t.abund <- computeFractions(merged.experiments, "cell.types.low.resolution")["T-cell",]
pdf(paste0(fig.dir, "barplot_tcell_frac_all.pdf"), height = 3, width = 12)
barplot(sort(t.abund * 100), horiz = F, las = 1, ylab = "% of T cells")
dev.off()

## IMPORT TCGA DATA AND FILTER
tcga <- readRDS(paste0(tcga.dir, "percentages_16clusters_tcga.Rds"))
rownames(tcga)[c(1,7)] <- c("Epith TSHR+ TPO+", "Epith TSHR+ TPO-")
annot.tcga <- readRDS(paste0(tcga.dir, "tcga_annotations.Rds"))
annot.tcga <- annot.tcga[rownames(annot.tcga) %in% colnames(tcga),]
tcga <- tcga[,annot.tcga$BRAFV600E_RAS %in% "BRAF_V600E"]
annot.tcga <- annot.tcga[rownames(annot.tcga) %in% colnames(tcga),]

## wilcoxon thyroiditis b-cells
w.res.b <- wilcox.test(x = tcga["B-cell",annot.tcga$thyroiditis %in%
                                            c("thyroiditis","hashimoto","graves")],
                     y = tcga["B-cell",annot.tcga$thyroiditis %in% c("other")],
                     alternative = "greater")
## boxplot thyroiditis b-cells
thyroiditis.or.not <- rep("Thyroiditis", times = ncol(tcga))
thyroiditis.or.not[annot.tcga$thyroiditis %in% c("other")] <- "No Thyroiditis"
pdf(file = paste0(fig.dir, "boxplot_bcell_tcga_brafonly.pdf"))
boxplot(formula = tcga["B-cell",] ~ thyroiditis.or.not)
dev.off()
## wilcoxon thyroiditis T-cells
w.res.t <- wilcox.test(x = tcga["T-cell",annot.tcga$thyroiditis %in%
                                            c("thyroiditis","hashimoto","graves")],
                     y = tcga["T-cell",annot.tcga$thyroiditis %in% c("other")],
                     alternative = "greater")
## boxplot thyroiditis T-cells
thyroiditis.or.not <- rep("Thyroiditis", times = ncol(tcga))
thyroiditis.or.not[annot.tcga$thyroiditis %in% c("other")] <- "No Thyroiditis"
pdf(file = paste0(fig.dir, "boxplot_tcell_tcga_brafonly.pdf"))
boxplot(formula = tcga["T-cell",] ~ thyroiditis.or.not)
dev.off()

## CORRELATION IN ST
s.genes <- c("CD3E","BLK","TSHR")
sp.sectors <- c("PTC3S3","PTC6S1","PTC7S1")
sp.filenames <- paste0(sp.sectors, "_seurat.Rds")
sp.bundle <- lapply(sp.filenames, function(z){readRDS(paste0(sp.dir, z))})
names(sp.bundle) <- sp.sectors
for (s.gene in s.genes){
    for (i in sp.sectors){
        SpFtPlotSave(sp.bundle[[i]], s.gene = s.gene)
    }
}

## B-T CORRELATIONS ACROSS ST SPOTS
samples <- grep("TC[0-9]S[0-9]", list.files(sp.dir), value = T)
stobj <- lapply(samples, function(z){readRDS(paste0(sp.dir, z))})
names(stobj) <- samples
stcors <- lapply(stobj, function(z){em <- z@assays$SCT@data
# IGKC can be used as an alternative to POU2AF1 but less expressed in naive B lymphocytes
    if (all(c("CD3E","POU2AF1") %in% rownames(em))){
        return(cor(em["CD3E",], em["POU2AF1",], method = "spearman"))
    } else { return(NA) }
})
stc <- unlist(stcors)
names(stc) <- gsub("_seurat.Rds","",names(stc))
keep <- !is.na(stc)
stcf <- stc[keep]
stcn <- as.numeric(stc)
stcnf <- stcn[keep]
barplot(stcf,las=3)
mean(stcnf[grepl("PTC",names(stcf))])
median(stcnf)
median(stcnf[grepl("PTC",names(stcf))])
