setwd(script.dir) ## Set working directory to the folder containing the scripts of this repository
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))

merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))
merged.3.6 <- readRDS(paste0(sn.dir, "merged.3.6.Rds"))
qc.metrics <- readRDS(paste0(sn.dir, "merged_qc_metrics.Rds"))
names(qc.metrics) <- gsub("SN1_seurat","", names(qc.metrics))
qc.metrics <- qc.metrics[order(names(qc.metrics))]
sn.stats1 <- data.frame(cell.counts = as.data.frame(table(merged.experiments$Sample))$Freq,
                       median.UMI.per.cell = tapply(X = merged.experiments$nCount_RNA, INDEX = merged.experiments$Sample, FUN = median), 
                       median.feature.per.cell = tapply(X = merged.experiments$nFeature_RNA, INDEX = merged.experiments$Sample, FUN = median), 
                       median.percentage.mt = tapply(X = merged.experiments$percent.mt, INDEX = merged.experiments$Sample, FUN = median))
sn.stats2 <- data.frame(cell.counts = as.data.frame(table(merged.3.6$orig.ident))$Freq,
                       median.UMI.per.cell = tapply(X = merged.3.6$nCount_RNA, INDEX = merged.3.6$orig.ident, FUN = median), 
                       median.feature.per.cell = tapply(X = merged.3.6$nFeature_RNA, INDEX = merged.3.6$orig.ident, FUN = median), 
                       median.percentage.mt = tapply(X = merged.3.6$percent.mt, INDEX = merged.3.6$orig.ident, FUN = median))
sn.stats <- rbind(sn.stats1, sn.stats2)
sn.stats <- sn.stats[order(rownames(sn.stats)),]
sn.stats$percentage.UMI.removed.soupx <- unlist(lapply(qc.metrics, function(z){z$soupx.rho}))
sn.stats$number.doublets <- unlist(lapply(qc.metrics, function(z){z$cells.removed["doublets"]}))
write.table(x = sn.stats, file = paste0(table.dir, "sn_qc_metrics.tsv"), sep = "\t", quote = F)

sp.filenames <- grep("_seurat.Rds", list.files(sp.dir), value = T)
sp.list <- lapply(sp.filenames, function(z){readRDS(paste0(sp.dir, z))})
names(sp.list) <- sp.filenames
sp.stats <- data.frame(row.names = sp.filenames)
for(sp in sp.filenames){
    sp.stats[sp,"spot.counts"] <- ncol(sp.list[[sp]])
    sp.stats[sp,"median.UMI.per.spot"] <- median(sp.list[[sp]]$nCount_RNA)
    sp.stats[sp,"median.feature.per.spot"] <- median(sp.list[[sp]]$nFeature_RNA)
    sp.stats[sp,"median.percentage.mt"] <- median(sp.list[[sp]]$percent.mt)
}
write.table(x = sp.stats, file = paste0(table.dir, "sp_qc_metrics.tsv"), sep = "\t", quote = F)

Idents(merged.experiments) <- merged.experiments$cell.types.low.resolution
lowres <- FindAllMarkers(merged.experiments, only.pos = T)
Idents(merged.experiments) <- merged.experiments$cell.types.medium.resolution
medres <- FindAllMarkers(merged.experiments, only.pos = T)
Idents(merged.experiments) <- merged.experiments$cell.types.high.resolution
highres <- FindAllMarkers(merged.experiments, only.pos = T)
 <- 
highres.short <- 
write.table(lowres, file = paste0(table.dir, "markers_low_all.csv"), quote = F, sep = ",")
write.table(medres, file = paste0(table.dir, "markers_medium_all.csv"), quote = F, sep = ",")
write.table(highres, file = paste0(table.dir, "markers_high_all.csv"), quote = F, sep = ",")
resolution <- list(low = read.csv(paste0(table.dir, "markers_low_all.csv")),
                  medium = read.csv(paste0(table.dir, "markers_medium_all.csv")),
                  high = read.csv(paste0(table.dir, "markers_high_all.csv")))
shortened <- lapply(resolution, function(z){
    short <- data.frame()
    for (c in unique(z$cluster)){
        new <- z[z$cluster %in% c,][1:100,]
        short <- rbind(short, new)
        }
    return(short)
})
lapply(names(shortened), function(z){write.table(shortened[[z]], file = paste0(table.dir, "markers_",z,"_short_all.csv"), quote = F, sep = ",")})

median(c(merged.experiments$nCount_RNA,merged.3.6$nCount_RNA))
median(c(merged.experiments$nFeature_RNA,merged.3.6$nFeature_RNA))
