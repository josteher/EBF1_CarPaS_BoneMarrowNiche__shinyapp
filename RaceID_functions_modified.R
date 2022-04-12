# Modified functions from the RaceID package


plotsymbolsmap_m <- function (object, types, subset = NULL, samples_col = NULL, cex = 0.25, 
                              fr = FALSE, leg = TRUE, map = TRUE, leg.pos="topleft", leg.cex=0.75, axes=F) 
    {
      if (is.null(subset)) 
        subset <- unique(types)
      h <- sort(unique(types)) %in% subset
      if (!is.null(subset)) {
        fp <- rep(FALSE, length(types))
        fp[types %in% subset] <- TRUE
      }
      if (is.null(samples_col)) {
        samples_col <- rainbow(length(unique(types[fp])))
      }
      else {
        samples_col <- samples_col[h]
      }
      if (fr | dim(object@tsne)[1] == 0) 
        d <- object@fr
      else d <- object@tsne
      if (map) {
        plot(d, xlab = "", ylab = "", axes = axes, cex = cex, 
             pch = 20, col = "grey")
        for (i in 1:length(unique(types[fp]))) {
          f <- types == sort(unique(types[fp]))[i]
          points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
                 cex = cex)
        }
      }
      else {
        plot(d, xlab = "", ylab = "", axes = axes, cex = 0, 
             pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                                 1])), ylim = c(min(d[, 2]), max(d[, 2])))
      }
      if (leg) 
        legend(leg.pos, legend = sort(unique(types[fp])), col = samples_col, 
               pch = 20, cex = leg.cex, bty = "n")
    }





plotmap_m <- function (object, final = TRUE, tp = 1, fr = FALSE, cex = 0.5, my_part=1:max(part)) 
{
  if (length(object@tsne) == 0 & length(object@fr) == 0) 
    stop("run comptsne/compfr before plotmap")
  if (final & length(object@cpart) == 0) 
    stop("run findoutliers before plotmap")
  if (!final & length(object@cluster$kpart) == 0) 
    stop("run clustexp before plotmap")
  if (!is.numeric(tp) | (is.numeric(tp) & tp > 1 | tp < 0)) 
    stop("tp has to be a number between 0 and 1 (transparency)")
  if (!is.logical(fr)) 
    stop("fr has to be TRUE or FALSE")
  part <- if (final) 
    object@cpart
  else object@cluster$kpart
  if (fr | dim(object@tsne)[1] == 0) 
    d <- object@fr
  else d <- object@tsne
  row.names(d) <- names(part)
  plot(d, xlab = "", ylab = "", cex = cex, axes = FALSE, col="lightgrey", pch=20)
  for (i in my_part) {
    if (sum(part == i) > 0) 
      points(d[part == i, 1], d[part == i, 2], col = adjustcolor(object@fcol[i], 
                                                                 tp), pch = 20, cex = cex)
  }
  for (i in my_part) {
    if (sum(part == i) > 0) 
      points(d[object@medoids[i], 1], d[object@medoids[i], 
                                        2], col = adjustcolor(object@fcol[i], tp), pch = 20, 
             cex = 4)
    if (sum(part == i) > 0) 
      points(d[object@medoids[i], 1], d[object@medoids[i], 
                                        2], col = adjustcolor("white", tp), pch = 20, 
             cex = 3)
    if (sum(part == i) > 0) 
      text(d[object@medoids[i], 1], d[object@medoids[i], 
                                      2], i, col = adjustcolor("black", tp), cex = 0.75, 
           font = 4)
  }
}





# If differentially expressed genes betweeen cluster 1 vs cluster 2 are already calculated, one can 
# use the following function in order to adapt the fold changes when cluster 2 against cluster 1 needs to be comapred 
calc_inv_order_diffex <- function(diffex_res){
  diffex_res$foldChange <- 1/diffex_res$foldChange
  diffex_res$log2FoldChange <- -diffex_res$log2FoldChange
  return(diffex_res)
}




# MA Plot with readable labels
maplot <- function(diffexp, 
                   title = NULL, 
                   fc = 1, 
                   size = 1, 
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = rownames(new.table),
                   legend = "top",
                   xlab = "Log2 mean expression",
                   ylab = "Log2 fold change",
                   top = 60,
                   font.label = c("italic", 11),
                   font.legend = "bold",
                   font.main   = "bold",
                   ggtheme     = ggplot2::theme_classic()){
  
  new.table <- diffexp$res[,c("baseMean","log2FoldChange","padj")]
  ggmaplot(new.table, 
           main        = title,
           fc          = fc,
           size        = size,
           palette     = palette,
           genenames   = genenames,
           legend      = legend,
           top         = top,
           font.label  = font.label,
           font.legend = font.legend,
           font.main   = font.main,
           xlab = xlab,
           ylab = ylab,
           ggtheme     = ggtheme)
} 



# Modified plotexpmap function
plotexpmap_m <- function (object, g, n = NULL, logsc = FALSE, imputed = FALSE, 
                          fr = FALSE, um = FALSE, cells = NULL, cex = 1, map = TRUE, 
                          leg = TRUE, noise = FALSE) 
{
  if (length(object@tsne) == 0 )
    stop("run comptsne/compfr/compumap before plotlabelsmap")
  if (!is.logical(fr)) 
    stop("fr has to be TRUE or FALSE")
  if (!is.logical(um)) 
    stop("um has to be TRUE or FALSE")
  if (length(intersect(g, rownames(object@ndata))) < length(unique(g))) 
    stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
  if (!is.numeric(logsc) & !is.logical(logsc)) 
    stop("argument logsc has to be logical (TRUE/FALSE)")
  if (!is.null(cells)) {
    if (sum(!cells %in% colnames(object@ndata)) > 0) 
      stop("cells has to be a subset of cell ids, i.e. column names of slot ndata")
  }
  # if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
  #   if (dim(object@fr)[1] != 0) {
  #     fr <- TRUE
  #   }
  #   else if (dim(object@umap)[1] != 0) {
  #     um <- TRUE
  #   }
  # }
  if (imputed & length(object@imputed) == 0) 
    stop("imputing needs to be done by running compdist with knn > 0")
  if (is.null(n)) 
    n <- g[1]
  if (is.null(cells)) 
    cells <- colnames(object@ndata)
  knn <- object@imputed$knn
  if (!noise) {
    if (length(g) == 1) {
      l <- as.vector(object@ndata[g, ] * min(object@counts) + 
                       0.1)
    }
    else {
      l <- apply(as.data.frame(as.matrix(object@ndata)[g, 
      ]) * min(object@counts), 2, sum) + 0.1
    }
    if (imputed) {
      l <- apply(rbind(object@imputed$nn, object@imputed$probs), 
                 2, function(y) {
                   ind <- y[1:(knn + 1)]
                   p <- y[(knn + 2):(2 * knn + 2)]
                   sum(l[ind] * p)
                 })
    }
  }
  else {
    if (is.null(object@noise)) 
      stop("run noise analysis first!")
    if (length(g) == 1) {
      l <- as.vector(object@noise[g, ] + 0.1)
    }
    else {
      l <- apply(as.data.frame(as.matrix(object@noise)[g, 
      ]), 2, sum) + 0.1
    }
  }
  if (logsc) {
    f <- l == 0
    l <- log2(l)
    l[f] <- NA
  }
  h <- colnames(object@ndata) %in% cells
  mi <- min(l, na.rm = TRUE)
  ma <- max(l, na.rm = TRUE)
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  # if (fr) {
  #   d <- object@fr
  # }
  # else if (um) {
  #   d <- object@umap
  # }
  # else {
    d <- object@tsne
  # }
  pardefault <- par()
  layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
         widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
  par(mar = c(3, 5, 2.5, 2))
  if (!leg) 
    n <- NA
  plot(c(min(d[, 1]), max(d[, 1])), c(min(d[, 2]), max(d[, 
                                                         2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
       col = "lightgrey", axes = FALSE)
  if (map) {
    v <- v[h]
    d <- d[h, ]
    kk <- order(v, decreasing = F)
    points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
           cex = cex)
  }
  if (leg) {
    par(mar = c(10, 2.5, 2.5, 4))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                 nrow = 1), col = ColorRamp, xlab = "", ylab = "", 
          xaxt = "n")
    # layout(1)
    par(mar = pardefault$mar)
  }
}






# Plotmarkergenes

plotmarkergenes_m <- function (object, genes, imputed = FALSE, cthr = 0, cl = NULL, 
                               cells = NULL, order.cells = FALSE, aggr = FALSE, norm = FALSE, 
                               cap = NULL, flo = NULL, samples = NULL, cluster_cols = FALSE, 
                               cluster_rows = TRUE, cluster_set = FALSE, samples_col = NULL, 
                               zsc = FALSE, logscale = TRUE, fontsize=10) 
{
  if (imputed & length(object@imputed) == 0) 
    stop("imputing needs to be done by running compdist with knn > 0")
  if (!is.null(cl)) {
    if (sum(!cl %in% object@cpart) > 0) 
      stop("cl has to be a subset of clusters in slot cpart")
  }
  if (!is.null(cells)) {
    if (sum(!cells %in% names(object@cpart)) > 0) 
      stop("cells has to be a subset of cell ids, i.e. names of slot cpart")
  }
  m <- aggregate(rep(1, length(object@cpart)), by = list(object@cpart), 
                 sum)
  pt <- object@cpart[object@cpart %in% m[m[, 2] > cthr, 1]]
  if (is.null(cl)) 
    cl <- 1:max(object@cpart)
  if (is.null(cells)) 
    cells <- names(object@cpart)
  pt <- pt[pt %in% cl]
  pt <- pt[names(pt) %in% cells]
  x <- as.matrix(object@ndata)[genes, ]
  if (imputed) {
    knn <- object@imputed$knn
    dd <- apply(x, 1, function(x) {
      apply(rbind(object@imputed$nn, object@imputed$probs), 
            2, function(y) {
              ind <- y[1:(knn + 1)]
              p <- y[(knn + 2):(2 * knn + 2)]
              sum(x[ind] * p)
            })
    })
    dd <- t(dd)
    colnames(dd) <- colnames(x)
    rownames(dd) <- rownames(x)
    x <- dd
  }
  x <- x[, names(pt)] * min(object@counts[names(pt)])
  f <- apply(x, 1, var) > 0 & apply(x > 0.1, 1, sum) > 1
  x <- x[f, ]
  if (norm) 
    x <- x/apply(x, 1, sum)
  x <- as.data.frame(as.matrix(x)) + 0.1
  z <- object@ndata[object@cluster$features, names(pt)] * min(object@counts[names(pt)]) + 
    0.1
  z <- as.matrix(z)
  if (aggr) {
    y <- aggregate(t(x), by = list(cl = pt), mean)
    z <- as.data.frame(as.matrix(t(y[, -1])))
    names(z) <- as.character(y[, 1])
    anc <- data.frame(cluster = paste("c", y[, 1], sep = ""))
    rownames(anc) <- names(z)
    v <- object@fcol[sort(unique(y[, 1]))]
    names(v) <- paste("c", sort(unique(y[, 1])), sep = "")
    xl <- log2(z)
    if (!is.null(cap)) {
      for (i in 1:ncol(xl)) xl[xl[, i] > cap, i] <- cap
    }
    if (!is.null(flo)) {
      for (i in 1:ncol(xl)) xl[xl[, i] < flo, i] <- flo
    }
    pheatmap(xl, cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
             border_color = NA, fontsize=fontsize )
  }
  else {
    if (length(unique(pt)) == 1) {
      n <- names(pt)
    }
    else {
      y <- aggregate(t(as.matrix(z)), by = list(cl = pt), 
                     mean)
      k <- hclust(as.dist(dist.gen(y[, -1], method = object@clusterpar$metric)))
      if (cluster_set) 
        set <- y[k$order, 1]
      else set <- cl
      n <- c()
      for (i in set) {
        p <- names(pt)[pt == i]
        if (length(p) >= 2) {
          k <- hclust(as.dist(dist.gen(as.matrix(t(z[apply(z[, 
                                                             p], 1, var) >= 0, p])), method = object@clusterpar$metric)))
          n <- append(n, p[k$order])
        }
        else {
          n <- append(n, p)
        }
      }
    }
    if (!is.null(samples)) {
      names(samples) <- colnames(object@ndata)
      anc <- data.frame(cluster = paste("c", pt[n], sep = ""), 
                        samples = samples[n])
    }
    else {
      anc <- data.frame(cluster = paste("c", pt[n], sep = ""))
    }
    rownames(anc) <- n
    v <- object@fcol[unique(pt[n])]
    names(v) <- paste("c", unique(pt[n]), sep = "")
    if (logscale) 
      xl <- log2(x[, n])
    else xl <- x[, n]
    if (zsc) 
      xl <- zscore(xl)
    if (!is.null(cap)) {
      for (i in 1:ncol(xl)) xl[xl[, i] > cap, i] <- cap
    }
    if (!is.null(flo)) {
      for (i in 1:ncol(xl)) xl[xl[, i] < flo, i] <- flo
    }
    if (order.cells) {
      g <- order(colnames(xl))
      xl <- xl[, g]
    }
    if (!is.null(samples)) {
      f <- object@cpart %in% sort(unique(pt[n]))
      h <- sort(unique(samples)) %in% unique(samples[f])
      samples <- samples[f]
      if (is.null(samples_col)) {
        saCol <- rainbow(length(unique(samples)))
        names(saCol) <- unique(samples)
      }
      else {
        saCol <- samples_col[h]
        names(saCol) <- sort(unique(samples))
      }
      pheatmap(xl, annotation_col = anc, annotation_colors = list(cluster = v, 
                                                                  samples = saCol), cluster_cols = cluster_cols, 
               cluster_rows = cluster_rows, show_colnames = FALSE, 
               border_color = NA, fontsize=fontsize)
    }
    else {
      pheatmap(xl, annotation_col = anc, annotation_colors = list(cluster = v), 
               cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
               show_colnames = FALSE, border_color = NA, fontsize=fontsize)
    }
  }
}




