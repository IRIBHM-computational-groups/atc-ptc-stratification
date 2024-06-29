setwd(script.dir) ## Set working directory to the folder containing the scripts of this repository
source("setup.R")
source("constants.R")
source("utils.R")
source("graphics.R")

setwd(GSE193581_folder) ## Set working directory to the folder containing the expression matrix files for GSE193581
mat.files <- list.files("GSE193581_RAW/")
counts.list <- lapply(mat.files, function(mat.file){read.table(paste0("GSE193581_RAW/", mat.file)) ; print(mat.file)})
sample.names <- gsub("GSM[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_", "", gsub("_UMI.txt", "", mat.files))
names(counts.list) <- sample.names
seurat.list <- lapply(sample.names, function(z){CreateSeuratObject(counts = counts.list[[z]], project = z)})
names(seurat.list) <- sample.names
saveRDS(seurat.list, file = "list_seurat_all.Rds")

seurat.merge <- merge(seurat.list[[1]], seurat.list[-1], add.cell.ids = NULL, project = "JCI")
annotation <- read.table("GSE193581_celltype_annotation.txt", sep = "\t")
rownames(annotation) <- gsub("-", ".", rownames(annotation))
seurat.merge$Cell.Types <- annotation[colnames(seurat.merge),"celltype"]
seurat.merge$Sample <- seurat.merge$orig.ident 
seurat.merge@project.name <- "jci"
seurat.merge[["percent.mt"]] <- PercentageFeatureSet(seurat.merge, pattern = "^MT-")
seurat.merge <- SCTransform(seurat.merge, vars.to.regress = "percent.mt")
seurat.merge <- RunPCA(seurat.merge, seed.use = 1)
seurat.merge <- FindNeighbors(seurat.merge, dims = 1:30)
seurat.merge <- FindClusters(seurat.merge, resolution = 0.1)
seurat.merge <- FindClusters(seurat.merge, resolution = 0.5)
seurat.merge <- FindClusters(seurat.merge, resolution = 1)
seurat.merge <- RunUMAP(seurat.merge, dims = 1:30, seed.use = 1)
saveRDS(seurat.merge, file = "seurat_merged.Rds")

seurat.merge <- readRDS("seurat_merged.Rds")
DimPlotSaveOther(seurat.merge, "Cell.Types")
DimPlotSaveOther(seurat.merge, "Sample")
DimPlotSaveOther(seurat.merge, "Sample")
lapply(c("SLC5A5","TPO","TG","TSHR"), FeaturePlotSave, seurat.merge)

qc.metrics <- seurat.merge@meta.data[,c("nCount_RNA","nFeature_RNA")]
qc.metrics$tumor.samples <- gsub("_..................$","",rownames(qc.metrics))
qc.summary <- data.frame(median_UMI = tapply(X = qc.metrics$nCount_RNA, INDEX = qc.metrics$tumor.samples, FUN = median),
                         median_features = tapply(X = qc.metrics$nFeature_RNA, INDEX = qc.metrics$tumor.samples, FUN = median),
                         cell_count = tapply(X = qc.metrics$nFeature_RNA, INDEX = qc.metrics$tumor.samples, FUN = length),
                         tumor_type = c(rep("ATC", 10), rep("NORM", 6), rep("PTC", 7)))
wilcox.test(qc.summary$median_UMI[qc.summary$tumor_type %in% "ATC"], qc.summary$median_UMI[qc.summary$tumor_type %in% "PTC"]) #pvalue = 0.2295
wilcox.test(qc.summary$median_features[qc.summary$tumor_type %in% "ATC"], qc.summary$median_features[qc.summary$tumor_type %in% "PTC"]) #pvalue = 0.4173
wilcox.test(qc.summary$cell_count[qc.summary$tumor_type %in% "ATC"], qc.summary$cell_count[qc.summary$tumor_type %in% "PTC"]) #pvalue = 0.8868

mutations <- data.frame(sample = rownames(qc.summary), mutation = c("BRAF","BRAF","P53","RAS","RAS","RAS","P53",NA,"RAS",NA, rep("None",6),NA,NA,"BRAF",NA,"BRAF","BRAF",NA))
seurat.merge$Mutations <- mutations[match(seurat.merge$Sample, mutations$sample), "mutation"]

DimPlotSaveOther(seurat.merge, "Mutations")

t.cell.counts <- tapply(seurat.merge$Cell.Types[seurat.merge$Cell.Types %in% c("T cell", "NK cell")], seurat.merge$Sample[seurat.merge$Cell.Types %in% c("T cell", "NK cell")], length)
b.cell.counts <- tapply(seurat.merge$Cell.Types[seurat.merge$Cell.Types %in% "B cell"], seurat.merge$Sample[seurat.merge$Cell.Types %in% "B cell"], length)
t.cell.prop <- t.cell.counts / table(seurat.merge$Sample)[names(t.cell.counts)]
b.cell.prop <- b.cell.counts / table(seurat.merge$Sample)[names(b.cell.counts)]

png(paste0(fig.dir, "barplot_tcell_frac_jci.png"), width = 200)
par(mar=c(5,8,4,1)+.1)
barplot(sort(t.cell.prop * 100), las = 2, horiz = T, xlab = "% of T cells")
dev.off()
png(paste0(fig.dir, "barplot_bcell_frac_jci.png"), width = 200)
par(mar=c(5,8,4,1)+.1)
barplot(sort(b.cell.prop * 100), las = 2, horiz = T, xlab = "% of B cells")
dev.off()