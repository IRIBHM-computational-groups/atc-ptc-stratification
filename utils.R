basicSeurat <- function(object){
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    object <- SCTransform(object, vars.to.regress = "percent.mt")
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000)
    object <- RunPCA(object, features = VariableFeatures(object = object))
    object <- FindNeighbors(object, dims = dimensionality)
    object <- FindClusters(object, resolution = 0.1)
    object <- FindClusters(object, resolution = 0.5)
    object <- FindClusters(object, resolution = 1)
    object <- RunUMAP(object, dims = dimensionality, md = md.setting)
    return(object)
}
printEntrez <- function(gene){
    entrezid <- mapIds(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL')
    entrezid <- entrezid[!is.na(entrezid)]
    summary <- unlist(lapply(entrezid, FUN = function(z) entrez_summary(db = "gene",id = z)[["summary"]]))
    return(summary)
}
processCiberResults <- function(results.dir, write.barplots = FALSE){
    reference <- str_match(results.dir, "cibersortx_(.*?)_")[2]
    norownames <- read.table(file = paste0(results.dir, "/CIBERSORTx_Adjusted.txt"), header = TRUE)
    results <- norownames[,-1]
    rownames(results) <- norownames[,1]
    colnames(results) <- gsub("\\.", "-", colnames(results))
    nc <- ncol(results)
    ncc <- nc - 3
    results.verif <- t(results[,1:ncc])
    ## addendum <- matrix(data = 0, nrow = 2, ncol = ncol(results.verif))
    ## rownames(addendum) <- c("Mesenchymal", "EMT")
    ## results.verif <- rbind(results.verif, addendum)
    ground.truth <- readRDS(paste0("/mnt/iribhm/ngs/TCGA-THCA/celltype_sample_counts_", reference, ".Rds"))
    ground.truth.sel <- ground.truth[rownames(results.verif),]
    ground.truth.sel <- apply(ground.truth.sel, 2, function(z) z / sum(z))
    results.verif <- cbind(results.verif, ground.truth.sel)
    if (write.barplots){
        blue <- rgb(0, 0, 1, alpha=0.5)
        red <- rgb(1, 0, 0, alpha=0.5)
        for (i in colnames(ground.truth)){
            sn.exp <- paste0(i, "SN1")
            cor.sn <- cor(results.verif[,sn.exp], results.verif[,i])
            print(paste("Correlation of", i, "is", cor.sn))
            png(paste0(results.dir, "/barplot_cor_", i, ".png"), width = 1200, height = 800)
            barplot(unlist(results.verif[,i]), col = blue,
                    main = paste("Prediction in red for", i, "with cor =", cor.sn), las = 2)
            barplot(unlist(results.verif[,sn.exp]), col = red, add = TRUE, las = 2)
            dev.off()
        }
    }
    return(results.verif)
}
hcCluster <- function(results.verif, n = 5){
    dist_mat <- dist(t(results.verif), method = 'euclidean')
    hclust_avg <- hclust(dist_mat, method = 'average')
    cut_avg <- cutree(hclust_avg, k = n)
    cut_avg[order(names(cut_avg))]
    return(cut_avg)
}
processCC <- function(object, sample, cchdb = CellChatDB.human){
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
    return(cellchat)
}
SVDNormalization <- function(dataset, n.dimentions = 10) {
    ## Args:
    ## dataset: a Seurat spatial transcriptomics object
    ## n.dimensions: number of dimention of SVD

    ## Returns:
    ## A Seurat spatial transcriptomics object with SVD as default assay

    ## Description:
    ## Smooth and impute read counts with Sigular Value Decomposition

    ## Take log
    DefaultAssay(object = dataset) <- "RNA"
    counts <- GetAssayData(object = dataset, slot = "counts")
    counts <- log10(1 + counts)

    ## Impute data
    svd <- svd(counts)
    imputed <- svd$u[, 1:n.dimentions] %*% diag(svd$d[1:n.dimentions]) %*% t(svd$v[,1:n.dimentions])
    rownames(imputed) <- rownames(counts)
    imputed <- 10^imputed - 1

    ## Include imupted matrix in Seurat object
    colnames(imputed) <- colnames(counts)
    dataset[['SVD']] <- CreateAssayObject(count = imputed)
    DefaultAssay(object = dataset) <- "SVD"
    return(dataset)
}
computeFractions <- function(s.object, group.var, norm.by = "sample"){
    df <- data.frame(orig = s.object$Sample, ident = s.object[[group.var]])
    names(df) <- c("orig","ident")
    m <- matrix(nrow = length(unique(df$ident)), ncol = length(unique(df$orig)))
    rownames(m) <- unique(df$ident)
    colnames(m) <- unique(df$orig)
    for (r in rownames(m)){
        for (c in colnames(m)){
            m[r,c] <- nrow(df[df$ident == r & df$orig == c,])
            }
        }
    m.samples <- apply(m, MARGIN = 2, FUN = function(z){z / sum(z)})
    m.clusters <- apply(m, MARGIN = 1, FUN = function(z){z / sum(z)})
    if (norm.by == "sample"){return(m.samples)}
    else {return(m.clusters)}
}
doFgsea <- function(dataset, gene, gs, topn=100, minSize = 10, maxSize = 500)
{
    DefaultAssay(object = dataset) <- 'SCT'
    expression <- GetAssayData(object = dataset, slot = "counts")
    vg <- dataset@assays$SCT@var.features
    if (!(gene %in% vg)) {
        vg = c(gene, vg)
    }
    head(vg)
    expression <- as.matrix(expression[vg,])
    stats <- apply(expression, 1, function(z) cor(z, expression[gene,], method='sp', use='pa'))
    x <- fgsea(pathways = gs, stats = stats, minSize = minSize, maxSize  = maxSize)
    clps <- collapsePathways(x, gs, stats)
    x <- x[order(x$pval),]
    return(list(FGSEA=x[1:topn,], COLLAPSED=clps))
}
doFgseaMulti <- function(dataset, genes, gs, topn=100)
{
    l <- lapply(genes, function(z) doFgsea(dataset, z, gs, topn=topn))
    return(l)
}
simplifyPathways <- function(z) {return(z$FGSEA[z$FGSEA$pathway %in% z$COLLAPSED$mainPathways,])}
getPathways <- function(z) 
{
    res <- lapply(z, function(z) lapply(z, function(zz) {
        y <- simplifyPathways(zz);
        y <- y[y$padj<0.1, c('pathway', 'padj', 'NES')]}))
    return(res)
}
getRecurrentPathways <- function(y)
{
    # Get pathways expressed in >2 tumors 
    tab <- sort(table(unlist(sapply(y, function(z) z$pathway))), dec=TRUE)
    tp <- names(tab[tab>1])
    sapply(tp, function(z) sapply(y, function(s) if (z %in% s$pathway) s[s$pathway==z, 'NES'] else NA))
}
getPCARank <- function(dataset, pc) 
{
    l <- Loadings(dataset, reduction ='pca')[,pc]
    rnk <- rank(l)
    names(rnk) <- names(l)
    return(rnk)
}
doFgseaPca <- function(dataset, pc, gs, topn = 100, minSize = 10, maxSize = 500)
{
    stats <- getPCARank(dataset, pc)
    x <- fgsea(pathways = gs, stats = stats, minSize = minSize, maxSize  = maxSize)
    clps <- collapsePathways(x, gs, stats)
    x <- x[order(x$pval),]
    return(list(FGSEA = x[1:topn,], COLLAPSED = clps))
}
doFgseaPcaMulti <- function(dataset, gs, topn = 100)
{
    l <- list(PC_1 = doFgseaPca(dataset, 'PC_1', gs, topn = topn),
         PC_2 = doFgseaPca(dataset, 'PC_2', gs, topn = topn),
         PC_3 = doFgseaPca(dataset, 'PC_3', gs, topn = topn))
    return(l)
}
getRecurrentPathwaysPca <- function(y)
{
    # Get pathways expressed in >2 tumors 
    tab <- sort(table(unlist(sapply(y, function(z) unique(unlist(sapply(z, function(zz) zz[, 'pathway'])))))), dec=TRUE)
    tp <- names(tab[tab>2])
    # Merge PCs per samples
    y1 <- lapply(y, function(z) rbind(z$PC_1, z$PC_2, z$PC_3))
    y1 <- lapply(y1, function(z) z[z$pathway %in% tp,])
    tp <- unique(unlist(lapply(y1, function(z) z$pathway)))
    #sapply(tp, function(z) sapply(y1, function(s) {if (!z %in% s$pathway) NA else s[pathway, 'padj']})) 
    final.df <- sapply(tp, function(z) sapply(y1, function(s) min(s[s$pathway==z, 'NES'])))
    final.df[final.df == Inf] <- NA
    return(final.df)
}
varianceExplained = function(SeuratObj, pc = 1:3)
{
    # On Seurat 3:
    SeuratObj <- RunPCA(SeuratObj, verbose = FALSE)
    mat <- Seurat::GetAssayData(SeuratObj, assay = "SCT", slot = "scale.data")
    pca <- SeuratObj[["pca"]]

    # Get the total variance:
    total_variance <- sum(matrixStats::rowVars(mat))
    
    eigValues = (pca@stdev)^2  ## EigenValues
    varExplained = eigValues / total_variance
    return(sum(varExplained[pc]))
}
