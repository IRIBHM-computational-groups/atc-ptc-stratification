setwd(script.dir) ## Set working directory to the folder containing the scripts of this repository
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")
merged.experiments <- readRDS(paste0(sn.dir, "merged_seurat_basic_SCT_MT_annotated.Rds"))

# PLOTS WHOLE DATASET
DimPlotSave(merged.experiments, group.var = "Sample")
DimPlotSave(merged.experiments, group.var = "Experiment")
DimPlotSave(merged.experiments, group.var = "SCT_snn_res.1")
DimPlotSave(merged.experiments, group.var = "cell.types.medium.resolution")
DimPlotSave(merged.experiments, group.var = "cell.types.low.resolution")
DimPlotSave(merged.experiments, group.var = "cell.types.high.resolution")
# STATS WHOLE DATASET
ssize.table <- table(merged.experiments$Sample)
highest.pct <- list()
cor.table <- list()
for (cluster in levels(merged.experiments$SCT_snn_res.1)){
    full.cluster.table <- rep(0, length(ssize.table))
    names(full.cluster.table) <- names(ssize.table)
    selected <- merged.experiments$SCT_snn_res.1 %in% cluster
    cluster.table <- table(merged.experiments$Sample[selected])
    full.cluster.table[names(cluster.table)] <- cluster.table
    highest.pct[[cluster]] <- max(cluster.table) / sum(selected)
    cor.table[[cluster]] <- cor(full.cluster.table, ssize.table, method = "spearman")
}
highest.pct <- unlist(highest.pct)
sort(highest.pct)
barplot(sort(highest.pct), xlab = "Cluster number", ylab = "Percentage of cells from the same sample")
abline(h = 0.95)
barplot(sort(unlist(cor.table)))
mean(unlist(cor.table))

# HARMONY WHOLE DATASET
harmo <- RunHarmony(merged.experiments, group.by.vars = "Experiment")
harmo <- RunUMAP(harmo, reduction = "harmony", dims = dimensionality)
harmo <- FindNeighbors(harmo, dims = dimensionality, reduction = "harmony")
harmo <- FindClusters(harmo, resolution = 0.1)
harmo <- FindClusters(harmo, resolution = 0.5)
harmo <- FindClusters(harmo, resolution = 1)
harmo@project.name <- "all_harmo"
saveRDS(harmo, file = paste0(sn.dir, "merged_seurat_harmony.Rds"))
DimPlotSave(harmo, group.var = "Sample")
DimPlotSave(harmo, group.var = "SCT_snn_res.1")
DimPlotSave(harmo, group.var = "cell.types.medium.resolution")
DimPlotSave(harmo, group.var = "Experiment")
# STATS HARMONY WHOLE DATASET
ssize.table.harmo <- table(harmo$Sample)
highest.pct.harmo <- list()
cor.table.harmo <- list()
for (cluster in levels(harmo$SCT_snn_res.1)){
    full.cluster.table.harmo <- rep(0, length(ssize.table.harmo))
    names(full.cluster.table.harmo) <- names(ssize.table.harmo)
    selected <- harmo$SCT_snn_res.1 %in% cluster
    cluster.table.harmo <- table(harmo$Sample[selected])
    full.cluster.table.harmo[names(cluster.table.harmo)] <- cluster.table.harmo
    highest.pct.harmo[[cluster]] <- max(cluster.table.harmo) / sum(selected)
    cor.table.harmo[[cluster]] <- cor(full.cluster.table.harmo, ssize.table.harmo, method = "spearman")
}
highest.pct.harmo <- unlist(highest.pct.harmo)
sort(highest.pct.harmo)
barplot(sort(highest.pct.harmo), xlab = "Cluster number", ylab = "Percentage of cells from the same sample")
abline(h = 0.95)
barplot(sort(unlist(cor.table.harmo)))
mean(unlist(cor.table.harmo))
wilcox.test(x = unlist(cor.table), y = unlist(cor.table.harmo))

# SEURAT'S CCA WHOLE DATASET
merged.experiments <- SetIdent(merged.experiments, value = "Experiment")
merged.list <- SplitObject(object = merged.experiments, split.by = "Experiment")
merged.list <- lapply(X = merged.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = merged.list)
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, anchor.features = features)
cca.merged <- IntegrateData(anchorset = merged.anchors)
cca.merged@project.name <- "all_cca"
DefaultAssay(cca.merged) <- "integrated"
cca.merged <- ScaleData(cca.merged, verbose = FALSE)
cca.merged <- RunPCA(cca.merged, npcs = 30, verbose = FALSE)
cca.merged <- RunUMAP(cca.merged, reduction = "pca", dims = dimensionality)
cca.merged <- FindNeighbors(cca.merged, reduction = "pca", dims = dimensionality)
cca.merged <- FindClusters(cca.merged, resolution = 1)
DimPlotSave(cca.merged, group.var = "Experiment")
DimPlotSave(cca.merged, group.var = "Sample")

# SEURAT'S RPCA WHOLE DATASET
merged.experiments <- SetIdent(merged.experiments, value = "Experiment")
merged.list <- SplitObject(merged.experiments, split.by = "Experiment")
merged.list <- lapply(X = merged.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = merged.list)
merged.list <- lapply(X = merged.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, anchor.features = features, reduction = "rpca")
cca.atcs <- IntegrateData(anchorset = merged.anchors)
cca.atcs <- ScaleData(cca.atcs, verbose = FALSE)
cca.atcs <- RunPCA(cca.atcs, npcs = 30, verbose = FALSE)
cca.atcs <- RunUMAP(cca.atcs, reduction = "pca", dims = dimensionality)
cca.atcs <- FindNeighbors(cca.atcs, reduction = "pca", dims = dimensionality)
cca.atcs <- FindClusters(cca.atcs, resolution = 1)
cca.atcs@project.name <- "all_rpca"
DimPlotSave(cca.atcs, group.var = "Experiment")
DimPlotSave(cca.atcs, group.var = "Sample")

# VIREO ATC RESULTS
atcs <- readRDS(paste0(sn.dir, "seurat_merged_ATC1_ATC2.Rds"))
atcs@project.name <- "atcs"
atcs$Cell.type <- atcs$labels
plast <- names(atcs$Experiment[atcs$Experiment %in% "ATC2SN1"])
DimPlotSave(atcs, group.var = "Experiment", order = plast)
DimPlotSave(atcs, group.var = "Sample")
DimPlotSave(atcs, group.var = "Cell.type")
DimPlotSave(atcs, group.var = "vireo")
DimPlotSave(atcs, group.var = "vireo.best")
a <- 0:max(as.numeric(levels(m$cluster))) ; a <- as.character(a)
names(a)[a %in% c(0)] <- "Cancer2" #
names(a)[a %in% c(1)] <- "Cancer1" #
names(a)[a %in% c(2)] <- "T-cell" #CD247,CD2,TTN,THEMIS,PTPRC
names(a)[a %in% c(3)] <- "Cancer2" #TG,PAX8,TPO
names(a)[a %in% c(4)] <- "Myeloid1" #HLAs,CD163
names(a)[a %in% c(5)] <- "Myeloid2" #CD163
names(a)[a %in% c(6)] <- "Endothelial" #VWF, FLT1
names(a)[a %in% c(7)] <- "Myeloid1" #XCR1, HLAs
labels <- as.character(atcs@active.ident)
names(labels) <- names((atcs@active.ident))
for (i in levels(atcs@active.ident)) {labels[labels == i] <- names(a[a == i])}
atcs$labels <- labels
DimPlot(atcs, group.by = "labels", label = T)

# ATC HARMONY
atcs.harmo <- RunHarmony(atcs, group.by.vars = "Experiment")
atcs.harmo <- RunUMAP(atcs.harmo, reduction = "harmony", dims = dimensionality)
atcs.harmo <- FindNeighbors(atcs.harmo, dims = dimensionality, reduction = "harmony")
atcs.harmo <- FindClusters(atcs.harmo, resolution = 1)
atcs.harmo@project.name <- "atcs_harmo"
saveRDS(atcs.harmo, file = paste0(sn.dir, "merged_seurat_atcs_harmony.Rds"))
DimPlotSave(atcs.harmo, group.var = "Experiment")
DimPlotSave(atcs.harmo, group.var = "labels")
DimPlotSave(atcs.harmo, group.var = "vireo")
DimPlotSave(atcs.harmo, group.var = "vireo.best")

# SEURAT'S CCA ATC
atcs <- SetIdent(atcs, value = "Experiment")
atcs.list <- SplitObject(object = atcs, split.by = "Experiment")
atcs.list <- lapply(X = atcs.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = atcs.list)
atcs.anchors <- FindIntegrationAnchors(object.list = atcs.list, anchor.features = features)
cca.atcs <- IntegrateData(anchorset = atcs.anchors)
cca.atcs <- ScaleData(cca.atcs, verbose = FALSE)
cca.atcs <- RunPCA(cca.atcs, npcs = 30, verbose = FALSE)
cca.atcs <- RunUMAP(cca.atcs, reduction = "pca", dims = dimensionality)
cca.atcs <- FindNeighbors(cca.atcs, reduction = "pca", dims = dimensionality)
cca.atcs <- FindClusters(cca.atcs, resolution = 1)
cca.atcs@project.name <- "atcs_cca"
DimPlotSave(cca.atcs, group.var = "Experiment")
DimPlotSave(cca.atcs, group.var = "labels")
DimPlotSave(cca.atcs, group.var = "vireo")
DimPlotSave(cca.atcs, group.var = "vireo.best")

# SEURAT'S RPCA ATC
atcs <- SetIdent(atcs, value = "Experiment")
atcs.list <- SplitObject(atcs, split.by = "Experiment")
atcs.list <- lapply(X = atcs.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = atcs.list)
atcs.list <- lapply(X = atcs.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
atcs.anchors <- FindIntegrationAnchors(object.list = atcs.list, anchor.features = features, reduction = "rpca")
rpca.atcs <- IntegrateData(anchorset = atcs.anchors)
rpca.atcs <- ScaleData(rpca.atcs, verbose = FALSE)
rpca.atcs <- RunPCA(rpca.atcs, npcs = 30, verbose = FALSE)
rpca.atcs <- RunUMAP(rpca.atcs, reduction = "pca", dims = dimensionality)
rpca.atcs <- FindNeighbors(rpca.atcs, reduction = "pca", dims = dimensionality)
rpca.atcs <- FindClusters(rpca.atcs, resolution = 1)
rpca.atcs@project.name <- "atcs_rpca"
DimPlotSave(rpca.atcs, group.var = "Experiment")
DimPlotSave(rpca.atcs, group.var = "labels")
DimPlotSave(rpca.atcs, group.var = "vireo")
DimPlotSave(rpca.atcs, group.var = "vireo.best")
