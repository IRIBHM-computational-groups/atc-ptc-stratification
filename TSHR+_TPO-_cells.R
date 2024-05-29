setwd(script.dir) ## Set working directory to the folder containing the scripts of this repository
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))
gs <- list(WIKIPATHWAYS = gmtPathways(paste0(misc.dir, "c2.cp.wikipathways.v7.5.1.symbols.gmt")))

## PTC2 - FIG 8B
Idents(merged.experiments) <- "Sample"
ptc2 <- subset(merged.experiments, idents = "PTC2")
ptc2 <- basicSeurat(ptc2)
ptc2@project.name <- "PTC2"
tmp <- as.character(ptc2$cell.types.medium.resolution)
tmp[ptc2$SCT_snn_res.0.1 %in% 0] <- "Epith TSHR+ TPO- CNA"
tmp[ptc2$SCT_snn_res.0.1 %in% 1] <- "Epith TSHR+ TPO- No CNA"
ptc2$cell.types.medium.resolution <- as.factor(tmp)
Idents(ptc2) <- "cell.types.medium.resolution"
DimPlotSave(ptc2, group.var = "cell.types.medium.resolution")
gg <- VlnPlot(ptc2, features = "TG", group.by = "cell.types.medium.resolution", idents = c("Epith TSHR+ TPO- CNA", "Epith TSHR+ TPO- No CNA"))
savePlotsLegend(gg, "vlnplot_TG_PTC2")
boxplot(formula = ptc2@assays$SCT@data["TG", ptc2$SCT_snn_res.0.1 %in% c(0,1)] ~ as.factor(ptc2$cell.types.medium.resolution[ptc2$SCT_snn_res.0.1 %in% c(0,1)]))
counts.cna <- ptc2@assays$SCT@counts["FN1", ptc2$cell.types.medium.resolution %in% "Epith TSHR+ TPO- CNA"]
counts.no.cna <- ptc2@assays$SCT@counts["FN1", ptc2$cell.types.medium.resolution %in% "Epith TSHR+ TPO- No CNA"]
counts.cna <- ptc2@assays$SCT@data["TG", ptc2$cell.types.medium.resolution %in% "Epith TSHR+ TPO- CNA"]
counts.no.cna <- ptc2@assays$SCT@data["TG", ptc2$cell.types.medium.resolution %in% "Epith TSHR+ TPO- No CNA"]
wilcox.test(x = counts.cna, y = counts.no.cna, alternative = "less")
mean(counts.cna)
mean(counts.no.cna)
mean(counts.no.cna) / mean(counts.cna)

## PTC2S6 CNV,TG - FIG 8C
ptc2s6 <- readRDS(paste0(sp.dir, "PTC2S6_seurat.Rds"))
load("/mnt/iribhm//software/spatial_transcriptomics/reference_files/gene_order_table.Rda")
g.sets <- list(Chr3 = rownames(gene_order_table[gene_order_table$chr %in% "chr3",]), 
               Chr7 = rownames(gene_order_table[gene_order_table$chr %in% "chr7",]), 
               Chr8 = rownames(gene_order_table[gene_order_table$chr %in% "chr8",]), 
               Chr22 = rownames(gene_order_table[gene_order_table$chr %in% "chr22",]))
ptc2s6 <- AddModuleScore(ptc2s6, features = g.sets)
colnames(ptc2s6@meta.data)[grepl("Cluster", colnames(ptc2s6@meta.data))] <- names(g.sets)
legend.param <- c(0.075, 0.87)
for (s.gene in c("Chr3","Chr7","Chr8","Chr22","TG")){
    SpFtPlotSave(ptc2s6, s.gene = s.gene)
}

## FIG 9A - BOXPLOT FN1 TCGA
tcga <- readRDS(paste0(tcga.dir, "counts_tcga.Rds"))
tcga <- apply(X = tcga, MARGIN = 2, FUN = function(z){z / sum(z) * 10^6})
rownames(tcga)[c(1,7)] <- c("Epith TSHR+ TPO+", "Epith TSHR+ TPO-")
annot.tcga <- readRDS(paste0(tcga.dir, "tcga_annotations.Rds"))
annot.tcga <- annot.tcga[rownames(annot.tcga) %in% colnames(tcga),]
annot.tcga$BRAF_STATUS <- annot.tcga$BRAFV600E_RAS
annot.tcga$BRAF_STATUS[annot.tcga$normal == 1] <- "NORMAL"
annot.tcga$BRAF_STATUS[annot.tcga$BRAF_STATUS == "RAS"] <- "OTHER"
pdf(paste0(fig.dir, "boxplot_FN1_tcga.pdf"), width = 4)
boxplot(tcga["FN1",] ~ annot.tcga$BRAF_STATUS, xlab = "Mutational status", ylab = "FN1 transcripts per million")
dev.off()

## FN1 IN ALL - FIG 9B-C
fn1.counts <- tapply(merged.experiments@assays$SCT@counts["FN1",], merged.experiments$cell.types.medium.resolution, sum)
pdf(paste0(fig.dir, "barplot_FN1_all.pdf"))
par(mar=c(5.1, 13 ,4.1 ,2.1))
barplot(fn1.counts, col = ct.medium.pal, horiz = TRUE, las = 2, main = "Total FN1 UMIs per cell type")
dev.off()

gg <- VlnPlot(merged.experiments, features = "FN1", group.by = "cell.types.medium.resolution", cols = ct.medium.pal)
savePlotsLegend(gg, "vlnplot_FN1_all")

## TPO- subcluster - FIG 9D-E-F
Idents(merged.experiments) <- merged.experiments$cell.types.medium.resolution
tponeg <- subset(x = merged.experiments, idents = "Epith TSHR+ TPO-")
tponeg$tumor.type <- substr(tponeg$Sample, 1, 3)
tponeg <- subset(x = tponeg, subset = tumor.type == "PTC")
tponeg <- basicSeurat(tponeg)
tponeg@active.ident <- tponeg$SCT_snn_res.0.5
pdf(paste0(fig.dir, "scatterplot_FN1_TG_tponeg.pdf"))
plot.scatter(x = tponeg@assays$SCT@data["FN1",], y = tponeg@assays$SCT@data["TG",], col = rgb(red = 0, green = 0, blue = 0, alpha = 0.01))
dev.off()

if (any(grepl("tponeg", ls()))) {saveRDS(tponeg, file = "/tmp/tponeg.Rds")} else {tponeg <- readRDS("/tmp/tponeg.Rds")}

ptcs.names <- paste0("PTC", 1:9)
ptcs <- SplitObject(tponeg, split.by = "Sample")
names(ptcs) <- ptcs.names
ptcs <- ptcs[! names(ptcs) %in% "PTC3"] ## remove PTC3 tponeg because there's only 34 nuclei
ptcs <- lapply(ptcs, basicSeurat)
saveRDS(ptcs, file = "~/tponeg_processed.Rds")
if (any(grepl("ptcs", ls()))) {saveRDS(ptcs, file = "/tmp/ptcs.Rds")} else {ptcs <- readRDS("/tmp/ptcs.Rds")}
print("checkpoint 1")

xtg <- lapply(gs['WIKIPATHWAYS'], function(z) lapply(ptcs, doFgseaMulti, genes='TG', gs=z, topn=1000))
if (any(grepl("xtg", ls()))) {saveRDS(xtg, file = "/tmp/xtg.Rds")} else {xtg <- readRDS("/tmp/xtg.Rds")}
ytg <- getPathways(xtg[['WIKIPATHWAYS']])      
ytg <- lapply(ytg, as.data.frame) 
recurrent.pathways <- getRecurrentPathways(ytg)
p <- recurrent.pathways
colnames(p) <- gsub('WP_', '', colnames(p))
colnames(p) <- sapply(colnames(p), function(z) {if (nchar(z)<=30) return(z) else return(paste0(substr(z,1,27), '...'))})
options(repr.plot.width=10, repr.plot.height=5)
pdf(paste0(fig.dir, "heatmap_gsets_TG_tponeg.pdf"), width = 15)
Heatmap(p, row_order = 1:nrow(p), column_order = 1:ncol(p), 
        column_names_rot = 45, column_names_side = "top",
        col=colorRamp2(c(min(p, na.rm=TRUE), 0, max(p, na.rm=TRUE)), c("blue", "white", "red")), name = 'NES')
dev.off()
                                               
print("checkpoint 2")
                                                
v <- sapply(ptcs, varianceExplained)
x <- lapply(gs, function(z) lapply(ptcs, doFgseaPcaMulti, gs=z, topn=1000))
if (any(grepl("x", ls()))) {saveRDS(xtg, file = "/tmp/x.Rds")} else {x <- readRDS("/tmp/x.Rds")}
y <- getPathways(x[['WIKIPATHWAYS']])
p <- getRecurrentPathwaysPca(y)
rownames(p) <- gsub('.padj', '', rownames(p))
colnames(p) <- gsub('WP_', '', colnames(p))
colnames(p) <- gsub('REACTOME_', '', colnames(p))
colnames(p) <- sapply(colnames(p), function(z) {if (nchar(z)<=30) return(z) else return(paste0(substr(z,1,27), '...'))})
p <- -log10(p)
p <- p[,1:10]
p <- p[, order(apply(p, 2, function(z) sum(is.na(z))))]
ha <- rowAnnotation(Variance = anno_barplot(v*100), gp = gpar(fill = 2:3, col = 2:3), width = unit(1.5, "cm"))
pdf(paste0(fig.dir, "heatmap_gsets_PCA_tponeg.pdf"), width = 20)
Heatmap(p, row_order = 1:nrow(p), column_order = 1:ncol(p), 
        column_names_rot = 45, column_names_side = "top",
        col=colorRamp2(c(min(p, na.rm=TRUE), 0, max(p, na.rm=TRUE)), c("blue", "white", "red")), 
        left_annotation = ha,
#        name = '-log(adj. p)'
       name = 'NES')
dev.off()
                     
print("checkpoint 3")

## FIG. 9G, spRNA-seq PTC7, TG, FN1
ptc7s1 <- readRDS(paste0(sp.dir, "PTC7S1_seurat.Rds"))
SpFtPlotSave(ptc7s1, "TG")
SpFtPlotSave(ptc7s1, "FN1")

## CELLCHAT FN1 - FIG 9H
samples <- unique(merged.experiments$Sample)
object <- merged.experiments
cchdb <- CellChatDB.human
cellchat.results <- list()
for (sample in samples){
    subsetted <- CreateSeuratObject(object@assays$RNA@counts[,object$Sample %in% sample])
    subsetted <- NormalizeData(subsetted)
    data.input <- GetAssayData(subsetted, assay = "RNA", slot = "data")
    identity <- object@meta.data[object$Sample %in% sample,]
    identity$cell.types.medium.resolution = droplevels(identity$cell.types.medium.resolution, 
                                                       exclude = setdiff(levels(identity$cell.types.medium.resolution),unique(identity$cell.types.medium.resolution)))
    cellchat <- createCellChat(object = data.input, meta = identity, group.by = "cell.types.medium.resolution")
    ## cellchat <- addMeta(cellchat, meta = identity, meta.name = "metadata")
    ## cellchat <- setIdent(cellchat, ident.use = "Annotation")
    groupSize <- as.numeric(table(cellchat@idents))
    cellchat@DB <- cchdb
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multicore", workers = 20) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat.results[[sample]] <- cellchat
    print(paste(sample, "done"))
}
samps <- grep("TC[0-9]SN1", list.files("/mnt/iribhm/ngs/ST/experiments"), value = T)
cellchat.results <- lapply(samps, function(z){readRDS(paste0("/mnt/iribhm/ngs/ST/experiments/", z, "/cellchat/cellchat_", gsub("SN1","",z), ".Rds"))})
names(cellchat.results) <- samples

cellchat.results <- lapply(samples, function(z){processCC(object = merged.experiments, sample = z)})
names(cellchat.results) <- samples
for (sample in samples[-c(4)]){
    groupSize <- as.numeric(table(cellchat.results[[sample]]@idents))
    pdf(file = paste0(fig.dir, "/circplot_FN1_", sample, ".pdf"))
    netVisual_aggregate(cellchat.results[[sample]], vertex.weight = groupSize, signaling = "FN1", layout = "circle")
    dev.off()
    print(sample)
}

## LIGANDS-FN1 CORRELATIONS - FIG 9I
receptors <- c("ITGA3_ITGB1","ITGA4_ITGB1","ITGA5_ITGB1","ITGA8_ITGB1","ITGAV_ITGB1","ITGA4_ITGB7","ITGAV_ITGB3","ITGA2B_ITGB3","ITGAV_ITGB6","ITGAV_ITGB8","CD44","SDC1","SDC4")
receptors.single <- c("ITGA3","ITGB1","ITGA4","ITGB1","ITGA5","ITGB1","ITGA8","ITGB1","ITGAV","ITGB1","ITGA4","ITGB7","ITGAV","ITGB3","ITGA2B","ITGB3","ITGAV","ITGB6","ITGAV","ITGB8","CD44","SDC1","SDC4")
seurats <- grep("TC[0-9]S[0-9]", list.files(sp.dir), value = T)
sectors <- gsub("_seurat.Rds", "", seurats)
seurats.list <- lapply(seurats, function(z){readRDS(paste0(sp.dir,z))})
names(seurats.list) <- sectors
mat.list <- lapply(seurats.list, function(z){return(z@assays$SCT@data)})
df <- matrix(nrow = length(receptors.single), ncol = length(seurats.list))
colnames(df) <- names(seurats.list)
rownames(df) <- receptors.single
for (gene in rownames(df)){
    df[gene,] <- sapply(mat.list, function(z){
        if(gene %in% rownames(z)) {cor(z["FN1",], z[gene,])
        }else{return(0)}
    }
    )
}
dff <- df[!apply(as.matrix(df), 1, function(z) is.na(sum(z))),]
pdf(paste0(fig.dir, "heatmap_FN1_integrins_corr_st.pdf"))
Heatmap(dff, col = divergingPal())
dev.off()
