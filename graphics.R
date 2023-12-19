savePlotsLegend <- function(ggobject, fname){
    legend <- cowplot::get_legend(ggobject)
    ggobject.nolegend <- ggobject + NoLegend()
    ggsave(filename = paste0(fig.dir, fname, ".pdf"), plot = ggobject.nolegend)
    ggsave(filename = paste0(fig.dir, fname, "_legend.pdf"), plot = legend)
}
plot.scatter <- function(x, y, col=transparent.rgb("black", 85), pch=19, log="",
                         main="", xlab="", ylab="", fit="lm", cor=TRUE, cor.x=0.2, cor.y=0.9,
                         id.line=FALSE, id.x=0.4, id.y=0.5, xat=NULL, yat=NULL, file="", ...){
    if (file != "") {
        do.call(gsub(".*([a-z]+{3})$", "\\1", file), list(file))
        par(mar=c(6.1, 6.5, 4.1, 1.1))
    }
    plot(x, y, col=col, pch=pch, cex=2, axes=FALSE, frame=T, log=log,
         main=main, xlab="", ylab="", cex.main=2.5, ...)
    axis(1, at=xat, cex.axis=2.5, padj=0.25)
    axis(2, at=yat, cex.axis=2.5)
    mtext(side=1, text=xlab, line=4, cex=2.5)
    mtext(side=2, text=ylab, line=4, cex=2.5)
    if (fit != "") {
        if (log=="") {
            k <- y
            l <- x
        } else if (log=="x") {
            l <- log10(x)
            k <- y
        } else if (log=="y") {
            k <- log10(y)
            l <- x
        } else if (log=="xy" || log=="yx") {
            k = log10(y)
            l = log10(x)
        } else {
            cat("Error: incorrect 'log' argument\n")
            return()
        }
        if (fit == "lm") {
            abline(lm(k ~ l), col="red")
        } else {
            abline(lqs(k ~ l), col="red")
        }
    }
    if (id.line) {
        abline(0, 1, col="green", lwd=2)
        pos.x <- min(x, na.rm=TRUE) + id.x*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
        pos.y <- min(y, na.rm=TRUE) + id.y*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE)) 
        text(pos.x, pos.y, "x=y", col="green", cex=4)
    }
    z <- cor.test(x, y, method="spearman", exact = FALSE)
    rho <- signif(z$estimate, d=1)
    if ((p <- z$p.value)==0) {
        p <- paste("(p<2e-16)", sep="")
    } else {
        p <- paste("(p=", signif(p, d=1), ")", sep="")
    }
    if (cor) {
        pos.x <- max(x, na.rm=TRUE) + cor.x*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE)) - 2
        pos.y <- min(y, na.rm=TRUE) + cor.y*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE)) 
        p.pos.y <- pos.y - 0.125*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        text(pos.x, pos.y, cex=2, labels=bquote(rho== .(signif(rho, d=2))))
        text(pos.x, p.pos.y, cex=1, labels=p)
    }
    if (file != "") {
        dev.off()
    }
    return(c(rho=rho, p=p))
}
makeColumnAnnotation <- function(annot.object, categ.vars, cont.vars){
    columnAnnotation(df = lapply(c(categ.vars,cont.vars), function(name){
        l <- list()
        l[[name]] <- annot.object[[name]]
        return(l)
    }),
                     col = c(
                         sapply(categ.vars, function(name){
                         vec <- annot.object[[name]]
                         vec.unique <- length(unique(vec[!is.na(vec)]))
                         result <- glasbey(2:c(vec.unique+1))
                         names(result) <- unique(vec[!is.na(vec)])
                         return(result)}),
                         sapply(cont.vars, function(name){
                         vec <- annot.object[[name]]
                         vec.range <- range(vec[!is.na(vec)])
                         result <- colorRamp2(vec.range,c("white","red"))
                         return(result)})))
}
transparent.rgb <- function(col, alpha=75){
    tmp <- c(col2rgb(col), alpha, 255)
    names(tmp) <- c("red", "green", "blue", "alpha", "maxColorValue")
    do.call("rgb", as.list(tmp))
}
makePal <- function(vec){
    vec.u <- unique(vec)
    length.vec <- 1:length(unique(vec.u)) + 1 
    return(glasbey(length.vec))
}
DimPlotSave <- function(s.object, group.var, ...){
    if (group.var == "cell.types.low.resolution"){cols.pal <- ct.low.pal ; tt <- "Cell type"}
    else if (group.var == "cell.types.medium.resolution"){cols.pal <- ct.medium.pal ; tt <- "Cell type"}
    else if (group.var == "cell.types.high.resolution"){cols.pal <- ct.high.pal ; tt <- "Cell type"}
    else if (group.var == "Sample"){cols.pal <- samples.pal ; tt <- group.var}
    else if (group.var == "Experiment"){cols.pal <- exp.pal ; tt <- group.var}
    else {cols.pal <- makePal(s.object[[group.var]][,1]) ; tt <- gsub("\\."," ",group.var)}
    ggobject <- DimPlot(object = s.object, group.by = group.var, cols = cols.pal, shuffle = T, ...) + labs(title = "", color = tt)
    legend <- cowplot::get_legend(ggobject)
    ggobject.nolegend <- ggobject + NoLegend()
    fname <- paste0("dimplot_", gsub("\\.","_",group.var), "_", s.object@project.name)
    ggsave(filename = paste0(fig.dir, fname, ".pdf"), plot = ggobject.nolegend)
    ggsave(filename = paste0(fig.dir, fname, "_legend.pdf"), plot = legend)
}
FeaturePlotSave <- function(ft.var, s.object, ...){
    ggobject <- FeaturePlot(object = s.object, features = ft.var, ...) + labs(title = "", color = ft.var)
    fname <- paste0("ftplot_", gsub("\\.","_",paste(ft.var, collapse = "_")), "_", s.object@project.name)
    ggsave(filename = paste0(fig.dir, fname, ".pdf"), plot = ggobject)
}
FeaturePlotSaveAll <- function(ft.var, s.object, ...){
    ggobject <- FeaturePlot(object = s.object, features = ft.var, ...) + labs(title = "", color = ft.var) + theme(legend.position = c(0.9, 0.2), 
                                                                           legend.title = element_text(size=20,face = "bold"), 
                                                                           axis.text = element_text(size=15),
                                                                          legend.text = element_text(size=15))
    fname <- paste0("ftplot_", gsub("\\.","_",paste(ft.var, collapse = "_")), "_", s.object@project.name)
    ggsave(filename = paste0(fig.dir, fname, ".pdf"), plot = ggobject)
}
FeaturePlotSaveOrganoids <- function(ft.var, s.object, ...){
    ggobject <- FeaturePlot(object = s.object, features = ft.var, ...) + labs(title = "", color = ft.var) + theme(legend.position = c(0.2, 0.8), 
                                                                           legend.title = element_text(size=20,face = "bold"), 
                                                                           axis.text = element_text(size=15),
                                                                          legend.text = element_text(size=15))
    fname <- paste0("ftplot_", gsub("\\.","_",paste(ft.var, collapse = "_")), "_", s.object@project.name)
    ggsave(filename = paste0(fig.dir, fname, ".pdf"), plot = ggobject)
}
HeatmapSave <- function(s.object, group.var, norm.by = "sample", do.cor = FALSE, ...){
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
    if (do.cor){
        if (norm.by == "sample"){
            h <- Heatmap(cor(m.samples, method = "spearman"), name = "Correlations sample normalized", col = divergingPal(), 
                    show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
        } else if (norm.by == "cluster") {
            h <- Heatmap(cor(m.clusters, method = "spearman"), name = "Correlations cluster normalized", col = divergingPal(), 
                    show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
            }
    } else {
            if (norm.by == "sample"){
            h <- Heatmap(t(m.samples), name = "Sample normalized", col = divergingPalBin(), 
                    show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
        } else if (norm.by == "cluster") {
            h <- Heatmap(t(m.clusters), name = "Cluster normalized", col = divergingPalBin(), 
                    show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
        }
    }
    fname <- paste0("heatmap_cor", do.cor, "_", gsub("\\.","_",group.var), "_", norm.by, "Norm_", s.object@project.name)
    pdf(paste0(fig.dir, fname, ".pdf"))
    draw(h)
    dev.off()
}
ciberHeatmap <- function(prop.mat = tcga, annot = annot.tcga, categ.vars.input = c("histological_type", "thyroiditis",
                                                                                   "Neoplasm_Disease_Stage"),
                                                                cont.vars.input = c("differentiation_score","purity",
                                                                                    "BRAF_RAF_score")){
    made.annot <- makeColumnAnnotation(annot.object = annot, categ.vars = categ.vars.input, cont.vars = cont.vars.input)
    h <- Heatmap(prop.mat, name='Composition', col = divergingPalBin(),
        show_column_names = FALSE, top_annotation = made.annot,
        show_column_dend = FALSE, show_row_dend = FALSE, cluster_rows = TRUE, cluster_columns = TRUE)
    return(h)
}
SpFtPlotSave <- function(s.object, s.gene, d.assay = "SCT"){
        s.name <- s.object@project.name
        legend.param <- c(0.065, 0.135)
        t.size <- c(20,15)
        legend.orient <- "vertical"
        if (DefaultAssay(s.object) != "SVD"){s.object <- SVDNormalization(s.object)}
        if (s.name %in% c("ATC2S1","ATC2S3","ATC3AS2","PTC2S1","PTC7S1","PTC9S1")){legend.param <- c(0.065, 0.135)}
        else if (s.name %in% c("ATC2S2","ATC3BS1","PTC2S5","PTC3S4") | grepl("PTC0S", s.name)){legend.param <- c(0.935, 0.135)}
        else if (s.name %in% c("PTC2S2") | grepl("ATC1S", s.name)){legend.param <- c(0.85, 0.04) ; t.size <- c(13,8) ; legend.orient = "horizontal"}
        else if (s.name %in% c("PTC2S3","PTC2S4","PTC2S6","PTC4S1","PTC6S1")){legend.param <- c(0.065, 0.87)}
        else if (s.name %in% c("PTC3S1","PTC3S2","PTC3S3","PTC8S1","ATC3AS1","ATC3AS3") | grepl("PTC1S", i)){legend.param <- c(0.935, 0.87)}
        else if (s.name %in% "PTC5S1"){legend.param <- c(0.96, 0.87) ; t.size <- c(15,12)}
        else if (s.name %in% "ATC3BS2"){legend.param <- c(0.045, 0.87) ; t.size <- c(15,12)}
        else if (s.name %in% "ATC3BS3"){legend.param <- c(0.96, 0.125) ; t.size <- c(15,12)}
        gg <- SpatialFeaturePlot(s.object, features = s.gene, pt.size.factor = 2.1, crop = F) + theme(legend.position = legend.param, 
                                                                          legend.title = element_text(size=t.size[1],face = "bold"), 
                                                                          axis.text = element_text(size=15),
                                                                          legend.text = element_text(size=t.size[2]),
                                                                          legend.direction = legend.orient)
        ggsave(paste0(fig.dir, "spftplot_", s.gene, "_", s.name, ".pdf"), plot = gg, dpi = 600)
}