library(bioDist)
library(gplots)
library(gtools)
library(limma)
library(DT)
library(png)
library(plotrix)
library(plyr)

heatmap.22 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun.row = dist, distfun.col = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
        w), symm = FALSE, scale = c("none", "row", "column"), 
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
    symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
    col = "heat.colors", colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("column", "row", "both", 
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
    colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5, 
    density.info = c("histogram", "density", "none"), denscol = tracecol, 
    symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
    key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, 
    key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && any(duplicated(breaks))) 
        stop("breaks may not contian duplicate values")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (all(Colv == "Rowv")) 
        Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
            (dendrogram %in% c("both", "row"))) {
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
            if (dendrogram == "both") 
                dendrogram <- "column"
            else dendrogram <- "none"
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
            (dendrogram %in% c("both", "column"))) {
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
            if (dendrogram == "both") 
                dendrogram <- "row"
            else dendrogram <- "none"
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
            nr)) 
            stop("Rowv dendrogram doesn't match size of x")
        if (length(rowInd) < nr) 
            nr <- length(rowInd)
    }
    else if (is.integer(Rowv)) {
        distr <- distfun.row(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        distr <- distfun.row(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Rowv)) {
        rowInd <- nr:1
        ddr <- as.dendrogram(hclust(dist(diag(nr))))
    }
    else {
        rowInd <- nr:1
        ddr <- as.dendrogram(Rowv)
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd > 
            nc)) 
            stop("Colv dendrogram doesn't match size of x")
        if (length(colInd) < nc) 
            nc <- length(colInd)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        distc <- distfun.col(if (symm) 
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        distc <- distfun.col(if (symm) 
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Colv)) {
        colInd <- 1:nc
        ddc <- as.dendrogram(hclust(dist(diag(nc))))
    }
    else {
        colInd <- 1:nc
        ddc <- as.dendrogram(Colv)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (!is.null(colRow)) 
        colRow <- colRow[rowInd]
    if (!is.null(colCol)) 
        colCol <- colCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    plot.index <- 1
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        plot.index <- plot.index + 1
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        plot.index <- plot.index + 1
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol) && is.null(colCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.null(srtCol) || is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            if (is.null(srtCol)) 
                srtCol <- 90
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol, col = colCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow) && is.null(colRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.null(srtRow) || is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow, col = colRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in 1:length(colInd)) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in 1:length(rowInd)) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol, 
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    plot.index <- plot.index + 1
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(stats:::plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
            yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(stats:::plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
            leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab)) 
            mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab)) 
            mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title)) 
            mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0) 
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min.breaks
            max.raw <- max.breaks
        }
        z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row") 
                key.xlab <- "Row Z-Score"
            else if (scale == "column") 
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
                cex = par("cex") * par("cex.lab"))
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE, 
                from = min.scale, to = max.scale)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[!omit]
            dens$y <- dens$y[!omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                  0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Density"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                  labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Count"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (is.null(key.title)) 
            title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else {
        par(mar = c(0, 0, 0, 0))
        plot.new()
    }
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
    if (!is.null(extrafun)) 
        extrafun()
    invisible(retval)
}

leftstr <- function(x, n = 1) substr(x, 1, min(n, nchar(x)))
rightstr <- function(x, n = 1) substr(x, nchar(x) - min(n - 1, nchar(x) - 1), nchar(x))
plot.null <- function()
	plot(0, 0, col = "transparent", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
mousecase <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

quarter.circle <- function (x, y, radius, quadrant = 1, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1){ 
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
        col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y
        if (quadrant == 1)
        	polygon(c(0, radius[circle] + x, xv[xv >= x & yv >= y], 0), 
        		c(0, 0, yv[xv >= x & yv >= y], radius[circle] * ymult + y), 
        		border = border, col = col[circle], lty = lty, lwd = lwd)
        if (quadrant == 2)
        	polygon(c(0, 0, xv[xv < x & yv >= y], -radius[circle] + x),
        		c(0, radius[circle] * ymult + y, yv[xv < x & yv >= y], 0), 
        		border = border, col = col[circle], lty = lty, lwd = lwd)
        if (quadrant == 3)
        	polygon(c(0, -radius[circle] + x, xv[xv < x & yv < y], 0),
        		c(0, 0, yv[xv < x & yv < y], -radius[circle] * ymult + y), 
        		border = border, col = col[circle], lty = lty, lwd = lwd)
        if (quadrant == 4)
        	polygon(c(0, 0, xv[xv >= x & yv < y], radius[circle] + x),
        		c(0, -radius[circle] * ymult + y, yv[xv >= x & yv < y], 0), 
        		border = border, col = col[circle], lty = lty, lwd = lwd)}}

ui <- fluidPage(
  titlePanel("Gene Expression Comparison App"),
  sidebarLayout(
  div(class = "col-sm-5",
  tags$form(class = "well",
  	selectInput("dataset", "Dataset", "GSE67985", "GSE67985", width = "30%"),
  		div(style="padding-bottom:3px", "Chai, Yang; Parada, Carolina; Grimaldi, Alexandre; Ho, Thach-Vu; Harunaga, Jill; Samuels, Bridget; Johnson, Daniel."),
   		div(style="padding-bottom:3px", "Distal/proximal mandible/maxilla at embryonic day 10.5/11.5/12.5/13.5/14.5, 3 replicates/condition."),
   		div(style="padding-bottom:3px", "Affymetrix Mouse Genome 430 2.0 Array."),
   		div(style="padding-bottom:7px", "Mouse drawing by Darrin Lunde and Nguyen Truong Son."),
 	conditionalPanel("!output.datloaded",
    		helpText("Data is loading (may take a minute) ...")),
    conditionalPanel("output.datloaded", 
    		helpText("Data is loaded."))),
   tags$form(class = "well",
  	  div(tags$label("Select genes/probesets and colors (optional for global analysis)")),
  	  bootstrapPage(
  		div(style="display:inline-block; width:15%", textInput("gene1", NULL, "name")),
  		div(style="display:inline-block; width:15%", textInput("gene2", NULL, "name")),
  		div(style="display:inline-block; width:15%", textInput("gene3", NULL, "name")),
   		div(style="display:inline-block; width:15%", textInput("gene4", NULL, "name")),
  		div(style="display:inline-block; width:15%", textInput("gene5", NULL, "name")),
  		div(style="display:inline-block; width:15%", textInput("gene6", NULL, "name"))),
      bootstrapPage(
       div(style="display:inline-block; width:15%", textInput("col1", NULL, "blue")),
	   div(style="display:inline-block; width:15%", textInput("col2", NULL, "red")),
	   div(style="display:inline-block; width:15%", textInput("col3", NULL, "forestgreen")),
	   div(style="display:inline-block; width:15%", textInput("col4", NULL, "magenta")),
	   div(style="display:inline-block; width:15%", textInput("col5", NULL, "darkorange")),
	   div(style="display:inline-block; width:15%", textInput("col6", NULL, "gold")))),
  tags$form(class = "well",
  	checkboxGroupInput("sel", "Select samples (each is three replicates)",
      choices = c(
        "E10.5_Mnd_D", "E11.5_Mnd_D", "E12.5_Mnd_D", "E13.5_Mnd_D", "E14.5_Mnd_D",
        "E10.5_Mnd_P", "E11.5_Mnd_P", "E12.5_Mnd_P", "E13.5_Mnd_P", "E14.5_Mnd_P",
        "E10.5_Max_D", "E11.5_Max_D", "E12.5_Max_D", "E13.5_Max_D", "E14.5_Max_D",
        "E10.5_Max_P", "E11.5_Max_P", "E12.5_Max_P", "E13.5_Max_P", "E14.5_Max_P"), inline = T)),
  tags$form(class = "well",
  splitLayout(cellWidths = c("35%", "65%"),
    radioButtons("comp", "Compare condition",
    		choices = c(
    			"distal vs. proximal" = "place",
    			"maxilla vs. mandible" = "bone",
    			"age (earliest vs. latest selected)" = "age"),
    		selected = "place"),
    list(conditionalPanel("input.comp == 'place'", 
	    	radioButtons("invert_place", br(), c("proximal down, distal up" = "normal", 
    				"distal down, proximal up" = "inverted"), "normal", inline = T)),
		conditionalPanel("input.comp == 'bone'",
    		radioButtons("invert_bone", list(br(), br()), c("mandible down, maxilla up" = "normal", 
    				"maxilla down, mandible up" = "inverted"), "normal", inline = T)),
		conditionalPanel("input.comp == 'age'",
    		radioButtons("invert_age", list(br(), br(), br(), br()),
    				c("earliest down, latest up" = "normal", 
    				"latest down, earliest up" = "inverted"), "normal", inline = T)))),
   		conditionalPanel("output.numsel > 6",
   			checkboxInput("heatadjust", "For heatmap, factor out comparisons other than the one selected",
   				value = T))),
  tags$form(class= "well",
      bootstrapPage(
    	div(style="display:inline-block; width:30%", textInput("fc", "Fold change cut-off", 2)),
        div(style="display:inline-block", checkboxInput("log", "log2", T))),    
    textInput("fdr", "False discovery rate", 0.01, width = "30%"),
    textInput("max", "Max. DE genes/probesets", Inf, width = "30%"),
    radioButtons("summary", "For each gene, show",
    		choices = c("average probeset" = "AVE",	"most highly expressed probeset" = "A",
    			"most differentially expressed probeset" = "M", "all probesets" = "Z"),
    			selected = "AVE", inline = T),
    radioButtons("heatcol", "For heatmap, high values are",
    		choices = c("green" = "rg", "red" = "gr"), selected = "rg", inline = T),
    radioButtons("distfun", "For heatmap, cluster genes using",
    		choices = c("absolute correlation" = "AC", "correlation" = "C", "Euclidean distance" = "E"),
    			selected = "AC", inline = T),
    radioButtons("heatscale", "For heatmap, scale genes by",
    		choices = c("mean-centering" = "MC", "Z-score" = "Z", "none" = "N"), selected = "MC", inline = T),
	radioButtons("basecol", "For individual plots, zero point is",
    		choices = c("E10.5", "mean of ages"), selected = "E10.5", inline = T))),
  mainPanel(width = 7,
  		tabsetPanel(
  			tabPanel("Global", 
				downloadButton("download.heatmap", "Download heatmap"),
				plotOutput("heatmap", height = "450px"), br(),
				downloadButton("download.ma.plot", "Download MA plot"),
				plotOutput("ma.plot", height = "500px"), br(),
    			downloadButton("download.table", "Download table"),	br(), br(),
    			dataTableOutput("table", width = "90%")),
    		tabPanel("Individual",
	    		downloadButton("download.bullseyes", "Download individual plot(s)"),
	    		plotOutput("bullseyes")
	    		)))
    		))

server <- function(input, output){
	age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
	bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
	place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
	full.design <- model.matrix(~age + bone + place)
	rownames(full.design) <- paste0(age, bone, place, 1:3)
	load("data/jaws3.R")
	load("data/exprs.R")	
	rownames(dat) <- dat[, "probeset"]
	genes.tab <- xtabs(~unlist(dat$symbol))
	single.genes <- names(genes.tab)[genes.tab == 1]
	output$datloaded <- reactive(is.numeric(nrow(dat)))
	outputOptions(output, "datloaded", suspendWhenHidden = FALSE)
	
	output$table <- renderDataTable({
		sel <- input$sel
		sel <- substr(sel, 2, nchar(sel))
		sel <- gsub("_", "", sel)
		sels <- character()
		for (s in sel) sels <- c(sels, paste0(s, 1:3))
		if (input$comp == "place"){ 
			target.col <- "placeD"
			ones <- grep("P", sels, value = T)
			twos <- grep("D", sels, value = T)
			ones.reduced <- gsub("P", "", ones)
			twos.reduced <- gsub("D", "", twos)
			ones <- ones[ones.reduced %in% twos.reduced]
			twos <- twos[twos.reduced %in% ones.reduced]
			invert <- ifelse(input$invert_place == "inverted", T, F)
			one <- ifelse(invert, "distal", "proximal")
			two <- ifelse(invert, "proximal", "distal")}
		if (input$comp == "bone"){
			target.col <- "boneMnd"
			ones <- grep("Max", sels, value = T)
			twos <- grep("Mnd", sels, value = T)
			ones.reduced <- gsub("Max", "", ones)
			twos.reduced <- gsub("Mnd", "", twos)
			ones <- ones[ones.reduced %in% twos.reduced]
			twos <- twos[twos.reduced %in% ones.reduced]
			invert <- ifelse(input$invert_bone == "inverted", T, F)
			one <- ifelse(invert, "maxilla", "mandible")
			two <- ifelse(invert, "mandible", "maxilla")}
		ages <- unique(substr(sels, 1, 4))
		if (length(ages)) age.1 <- min(ages) else age.1 <- 10.5
		young.col <- paste0("age", age.1)
		if (input$comp == "age"){
			age.2 <- max(ages)
			target.col <- paste0("age", age.2)
			ones <- grep(age.1, sels, value = T)
			twos <- grep(age.2, sels, value = T)
			ones.reduced <- gsub(age.1, "", ones)
			twos.reduced <- gsub(age.2, "", twos)
			ones <- ones[ones.reduced %in% twos.reduced]
			twos <- twos[twos.reduced %in% ones.reduced]
			invert <- ifelse(input$invert_age == "inverted", T, F)
			one <- ifelse(invert, paste0("E", age.2), paste0("E", age.1))
			two <- ifelse(invert, paste0("E", age.1), paste0("E", age.2))}
		full.des <- full.design[c(ones, twos), ]
		constant <- apply(full.des, 2, function(x) length(unique(x)) == 1)
		constant[1] <- F
		constant[names(constant) == young.col] <- T
		design.cols <- names(constant[constant == F])
		design <- data.frame(lapply(design.cols, function(x){col <- list(full.des[, x])}))
		names(design) <- design.cols
		control.cols <- design.cols[design.cols != target.col]
		control <- data.frame(lapply(control.cols, function(x){col <- list(full.des[, x])}))
		names(control) <- control.cols
		if (!nrow(design) | !target.col %in% colnames(design) | !all(sels %in% c(ones, twos))) return({
			output$ma.plot <- renderPlot(plot.null())
			output$heatmap <- renderPlot(plot.null())
			data.frame(NULL)})
		output$numsel <- reactive(length(c(ones, twos)))	
		outputOptions(output, "numsel", suspendWhenHidden = FALSE)
		dat.sel <- dat[, c(ones, twos, "symbol")]
		if (input$summary == "AVE"){
			dat.sel <- ddply(dat.sel, .(symbol), numcolwise(mean))
			dat.sel <- dat.sel[dat.sel$symbol != "NA", ]
			rownames(dat.sel) <- dat$probeset[match(dat.sel$symbol, dat$symbol)]}
		if (input$summary == "A"){
			dat.sel$A <- rowMeans(dat.sel)
			dat.sel$P <- rownames(dat.sel)
			dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
			dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
			dat.agg <- aggregate(A ~ symbol, dat.multiple, max)
			dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
			rownames(dat.multiple) <- dat.multiple$P
			dat.sel <- rbind(dat.single, dat.multiple)
			dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("A", "P")]]}
		if (input$summary == "M"){
			dat.sel$MM <- abs(rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones]))
			dat.sel$P <- rownames(dat.sel)
			dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
			dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
			dat.agg <- aggregate(MM ~ symbol, dat.multiple, max)
			dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
			rownames(dat.multiple) <- dat.multiple$P
			dat.sel <- rbind(dat.single, dat.multiple)
			dat.sel <- dat.sel[, colnames(dat.sel)[!colnames(dat.sel) %in% c("MM", "P")]]}
		fit <- lmFit(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]], design)
		efit <- eBayes(fit)
		dat.sel$A <- rowMeans(dat.sel[, colnames(dat.sel)[colnames(dat.sel) != "symbol"]])
		dat.sel$M <- rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones])
		if (invert) dat.sel$M <- -dat.sel$M
		if (!is.null(dat.sel$symbol))
			dat.sel$symbol[is.na(dat.sel$symbol)] <-  rownames(dat.sel)[is.na(dat.sel$symbol)]
		dat.sel$color <- "black"
		dat.sel$color[dat.sel$symbol == mousecase(input$gene1) | rownames(dat.sel) == mousecase(input$gene1)] <- 
			tolower(input$col1)
		dat.sel$color[dat.sel$symbol == mousecase(input$gene2) | rownames(dat.sel) == mousecase(input$gene2)] <- 
			tolower(input$col2)
		dat.sel$color[dat.sel$symbol == mousecase(input$gene3) | rownames(dat.sel) == mousecase(input$gene3)] <- 
			tolower(input$col3)
		dat.sel$color[dat.sel$symbol == mousecase(input$gene4) | rownames(dat.sel) == mousecase(input$gene4)] <- 
			tolower(input$col4)
		dat.sel$color[dat.sel$symbol == mousecase(input$gene5) | rownames(dat.sel) == mousecase(input$gene5)] <- 
			tolower(input$col5)
		dat.sel$color[dat.sel$symbol == mousecase(input$gene6) | rownames(dat.sel) == mousecase(input$gene6)] <- 
			tolower(input$col6)
		lfc <- ifelse(input$log, abs(as.numeric(input$fc)), log2(abs(as.numeric(input$fc))))
		top <- topTable(efit, coef = which(colnames(design) == target.col), 
			lfc = lfc, p.value = as.numeric(input$fdr), num = Inf)
		top <- merge(top, dat[, c("probeset", "symbol", "desc")], by.x = 0, by.y = "probeset")
		top.colored <- top[rownames(top) %in% dat.sel$probeset[dat.sel$color != "black"], ]
		if (nrow(top.colored)){
			top.black <- top[1:min(nrow(top), as.numeric(input$max)), ]
			top <- merge(top.colored, top.black, all = T)}
		if (nrow(top)){
			if (invert) top$logFC <- -top$logFC			
			top <- top[order(top$logFC, decreasing = T), ]
			top$FC <- 2 ^ top$logFC
			top$FC[top$logFC < 0] <- -1 / top$FC[top$logFC < 0]
			top$AveExpr <- round(top$AveExpr, 3)
			top$logFC <- round(top$logFC, 2)
			top$FC <- round(top$FC, 2)
			top$adj.P.Val <- format(top$adj.P.Val, scientific = T, digits = 2)
			top$desc <- format(top$desc, justify = "left")
			top$symbol <- format(top$symbol, justify = "left")
			top <- top[, c("Row.names", "symbol", "AveExpr", "logFC", "FC", "adj.P.Val", "desc")]
			colnames(top) <- c("Probeset", "Gene", "Ave.Expr", "Log2FC", "Fold.Change", "Adj.P.Val", "Description")}
		else top <- data.frame(NULL)
		dat.top <- dat.sel[rownames(dat.sel) %in% top$Probeset, ]
		
		output$ma.plot <- renderPlot({
			ma.plot <- function(){
				par(mar = c(5, 5, 5, 8))
				lim <- max(abs(range(dat.sel$M)))
				plot(NULL, cex = 0.9,
					xaxt = "n", yaxt = "n",
					xlab = "Average Expression",
					xlim = range(dat.sel$A),
					ylim = c(-lim * 1.05, lim * 1.05), 
					ylab = paste0(one,
						paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""),
						ifelse(input$log, "Log2 Fold Change", "Fold Change"),
						paste(rep(" ", ifelse(input$log, 15, 20)), collapse = ""), two),
					cex.lab = 1.2, font.main = 1, cex.main = 1.2,
					main = paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " "))
				points(dat.sel$A[dat.sel$color == "black"], dat.sel$M[dat.sel$color == "black"], pch = 21, cex = 0.9)
				points(dat.sel$A[dat.sel$color != "black"], dat.sel$M[dat.sel$color != "black"], pch = 21, cex = 0.9,
					col = dat.sel$color[dat.sel$color != "black"])
				if (nrow(dat.top)){
					points(dat.top$A, dat.top$M, pch = 21, cex = 0.9, col = dat.top$color, bg = dat.top$color)
					text(dat.top$A, dat.top$M, dat.top$symbol, col = dat.top$color, 
						pos = ifelse(dat.top$M > 0, 3, 1), offset = 0.4, cex = 0.9)}
				ats.0 <- seq(1, 9, 1)
				ats <- c(-1 * rev(ats.0), 0, ats.0)
				log.labs <- ats
				log.labs[log.labs > 0] <- paste0("+", log.labs[log.labs > 0])
				raw.labs <- 2 ^ ats
				raw.labs[raw.labs < 1] <- -1 / raw.labs[raw.labs < 1]
				raw.labs[raw.labs > 1] <- paste0("+", raw.labs[raw.labs > 1])
				axis(1, seq(0, 20, 1))
				if (input$log) axis(2, ats, log.labs)
				else axis(2, ats, raw.labs, cex.axis = ifelse(lim > 6, 0.75, 1))
				abline(h = lfc * c(-1, 1), lty = 2)}
 			plotPNG(ma.plot, "ma.plot.png", width = 1500, height = 750, res = 300, pointsize = 4)
			ma.plot()
		})
			
		output$download.ma.plot <- downloadHandler("Facebase_Microarray_MA_Plot.png", 
			function(file) file.copy("ma.plot.png", file))	

		output$heatmap <- renderPlot({
			heatmap <- function(){
			withProgress({
				if (nrow(dat.top) < 2) return(plot.null())
				dat.heat <- as.matrix(dat.top[, sels])
				dat.heat <- dat.heat - mean(dat.heat)
				if (input$heatadjust){
					for (col in colnames(control)[colnames(control) != "(Intercept)"]){
						for (v in unique(control[, col])){
							mm <- mean(dat.heat[, rownames(control)[control[, col] == v]])
							dat.heat[, rownames(control)[control[, col] == v]] <- 
								dat.heat[, rownames(control)[control[, col] == v]] - mm}}}
				if (input$heatscale %in% c("MC", "Z")) dat.heat <- sweep(dat.heat, 1, rowMeans(dat.heat))
				if (input$heatscale == "Z") dat.heat <- sweep(dat.heat, 1, apply(dat.heat, 1, sd), "/")
				weights <- (ncol(dat.heat):1) + 10000
				weights[colnames(dat.heat) %in% twos] <- 
					weights[colnames(dat.heat) %in% twos] + ifelse(xor(invert, input$comp == "bone"), -100, 100)
				extreme <- ceiling(10 *  max(dat.heat)) / 10
				heatmap.22(t(dat.heat),
					Rowv = as.integer(weights),
					col = ifelse(input$heatcol == "rg", redgreen, greenred),
					cexCol = min(1.5, 0.2 + 1/log10(nrow(dat.heat))),
					scale = "none", key = T, key.title = NA, lwid = c(1, 4), lhei = c(1.5, 4),
					density.info = "density", densadj = 0.5, denscol = "white",
					key.ylab = NA, key.xlab = ifelse(input$log, "Log2 Fold Change", "Fold Change"),
					key.par = list(cex.lab = 1.25, cex.axis = 1, tcl = -0.35, mgp = c(2, 1, 0)),
					key.ytickfun = function() list(labels = F, tick = F),
					key.xtickfun = function(){
						breaks <- parent.frame()$breaks
						list(at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
							 labels = round(seq(-extreme, extreme, length = 6), 1), 
							 padj = -0.25)},
					trace = "none", margins = c(5, 10),
					distfun.row = dist,
					distfun.col = function(...) {
						if (input$distfun == "E") return(dist(...))
						if (input$distfun == "AC") return(cor.dist(..., abs = T))
						if (input$distfun == "C") return(cor.dist(..., abs = F))},
					labCol = dat.top$symbol, colCol = dat.top$color)
			incProgress(1)}, value = NULL, message = "Calculations in progress...")}
			plotPNG(heatmap, "heatmap.png", width = 1500, height = 750, res = 300, pointsize = 6)
			heatmap()
		})
			
		output$download.heatmap <- downloadHandler("Facebase_Microarray_Heatmap.png",
			function(file) file.copy("heatmap.png", file))
			
		rownames(top) <- NULL
		output$download.table <- downloadHandler("Facebase_Microarray_Table.csv",
			function(file) write.csv(as.data.frame(top), file, row.names = F))
		top$Color <- dat.top$color[match(top$Probeset, rownames(dat.top))]
		top <- formatStyle(datatable(top), "Color", target = "row", 
			color = styleEqual(
				c(0, input$col1, input$col2, input$col3, input$col4, input$col5),
				c("white", input$col1, input$col2, input$col3, input$col4, input$col5)))
	})
		
	output$bullseyes <- renderPlot({
	bullseyes <- function(){
		genes <- c(input$gene1, input$gene2, input$gene3, input$gene4, input$gene5, input$gene6)
		good.genes <- which(genes %in% dat$symbol | genes %in% dat$probeset)
		genes <- genes[good.genes]
		if (length(genes) == 0) return(plot.null)
		colors <- c(input$col1, input$col2, input$col3, input$col4, input$col5, input$col6)
		colors <- colors[good.genes]
		exprs <- dat[dat$symbol %in% genes, ]
		exprs <- data.frame(
			symbol = exprs$symbol,
			probeset = exprs$probeset,
			E14.5MaxD = rowMeans(exprs[, grep("14.5MaxD", colnames(exprs))]),
			E13.5MaxD = rowMeans(exprs[, grep("13.5MaxD", colnames(exprs))]),
			E12.5MaxD = rowMeans(exprs[, grep("12.5MaxD", colnames(exprs))]),
			E11.5MaxD = rowMeans(exprs[, grep("11.5MaxD", colnames(exprs))]),
			E10.5MaxD = rowMeans(exprs[, grep("10.5MaxD", colnames(exprs))]),
			E14.5MaxP = rowMeans(exprs[, grep("14.5MaxP", colnames(exprs))]),
			E13.5MaxP = rowMeans(exprs[, grep("13.5MaxP", colnames(exprs))]),
			E12.5MaxP = rowMeans(exprs[, grep("12.5MaxP", colnames(exprs))]),
			E11.5MaxP = rowMeans(exprs[, grep("11.5MaxP", colnames(exprs))]),
			E10.5MaxP = rowMeans(exprs[, grep("10.5MaxP", colnames(exprs))]),
			E14.5MndD = rowMeans(exprs[, grep("14.5MndD", colnames(exprs))]),
			E13.5MndD = rowMeans(exprs[, grep("13.5MndD", colnames(exprs))]),
			E12.5MndD = rowMeans(exprs[, grep("12.5MndD", colnames(exprs))]),
			E11.5MndD = rowMeans(exprs[, grep("11.5MndD", colnames(exprs))]),
			E10.5MndD = rowMeans(exprs[, grep("10.5MndD", colnames(exprs))]),
			E14.5MndP = rowMeans(exprs[, grep("14.5MndP", colnames(exprs))]),
			E13.5MndP = rowMeans(exprs[, grep("13.5MndP", colnames(exprs))]),
			E12.5MndP = rowMeans(exprs[, grep("12.5MndP", colnames(exprs))]),
			E11.5MndP = rowMeans(exprs[, grep("11.5MndP", colnames(exprs))]),
			E10.5MndP = rowMeans(exprs[, grep("10.5MndP", colnames(exprs))]),
			stringsAsFactors = F)
		if (input$summary %in% c("AVE")){
			exprs.full <- exprs
			exprs <- ddply(exprs, .(symbol), numcolwise(mean))
			for (i in 1:nrow(exprs)) exprs$probeset[i] <- 
				paste(exprs.full$probeset[exprs.full$symbol == exprs$symbol[i]], collapse = "\n")
			exprs <- exprs[, c(1, 22, 2:21)]}
		if (input$summary == "A"){
			exprs$rms <- rowMeans(exprs[, 3:21])
			maxs <- ddply(exprs, .(symbol), summarise, mx = which.max(rms))
			expr <- exprs[0, ]
			for (sym in maxs$symbol)
				expr <- rbind(expr, exprs[exprs$symbol == sym, ][maxs$mx[maxs$symbol == sym], ])
			exprs <- expr}
		if (input$summary == "M"){
			exprs$range <- apply(exprs[, -(1:2)], 1, function(x) diff(range(x)))
			maxs <- ddply(exprs, .(symbol), summarise, mx = which.max(range))
			expr <- exprs[0, ]
			for (sym in maxs$symbol)
				expr <- rbind(expr, exprs[exprs$symbol == sym, ][maxs$mx[maxs$symbol == sym], ])
			exprs <- expr}
		exprs <- exprs[order(match(exprs$symbol, genes), exprs$probeset), ]
		if (input$basecol == "mean of ages")
			exprs[, -(1:2)] <- sweep(exprs[, -(1:2)], 1, rowMeans(exprs[, -(1:2)]))
		if (input$basecol == "E10.5") exprs[, -(1:2)] <-
			sweep(exprs[, -(1:2)], 1, rowMeans(exprs[, grep("E10.5", colnames(exprs))]))
		if (nrow(exprs) == 1) par(mfrow = c(1, 1)) else par(mfrow = c(ceiling(nrow(exprs) / 2), 2))
		par(mar = c(1, 1, 1, 1))					
		mn <- min(unlist(exprs[, -(1:2)]))
		mx <- max(unlist(exprs[, -(1:2)]))
		two <- ifelse(nrow(exprs) == 2, T, F)
		for (i in 1:nrow(exprs)){
			expr <- exprs[i, ]
			pal <- colorRampPalette(c("white", colors[which(genes == expr$symbol)]))(1001)
			cols <- function(x) pal[round((x - mn) / (mx - mn) * 1000) + 1]
			plot(c(-1.25, 1.25), c(-1.875, 1.875), type = "n", xaxt = "n", yaxt = "n", frame.plot = F)
		 	rasterImage(bg, par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
			dmax <- unlist(expr[, grep("MaxD", colnames(expr))])
			pmax <- unlist(expr[, grep("MaxP", colnames(expr))])
			dmnd <- unlist(expr[, grep("MndD", colnames(expr))])
			pmnd <- unlist(expr[, grep("MndP", colnames(expr))])
			all <- c(dmax, pmax, dmnd, pmnd)
 			if (length(dmax) == 0) return(plot.null())
			quarter.circle(0, 0, c(1.025, 0.825, 0.625, 0.425, 0.225),
				quadrant = 1, border = NA, col = cols(dmax))
			quarter.circle(0, 0, c(1.025, 0.825, 0.625, 0.425, 0.225),
				quadrant = 2, border = NA, col = cols(pmax))
			quarter.circle(0, 0, c(1.025, 0.825, 0.625, 0.425, 0.225),
				quadrant = 3, border = NA, col = cols(pmnd))
			quarter.circle(0, 0, c(1.025, 0.825, 0.625, 0.425, 0.225),
				quadrant = 4, border = NA, col = cols(dmnd))
			draw.circle(0, 0, c(1.025, 0.825, 0.625, 0.425, 0.225), border = "black", lwd = 0.5)
			segments(-1.025, c(-0.028, 0.028), 1.025, c(-0.028, 0.028), lwd = 0.5)
			segments(c(-0.03, 0.03), -1.025 * getYmult(), c(-0.03, 0.03), 1.025 * getYmult(), lwd = 0.5)
			rect(-2, -0.025, 2, 0.025, border = NA, col = "white")
			rect(-0.025, -1.1, 0.025, 1.1, border = NA, col = "white")
			text(c(0.125, 0.325, 0.525, 0.725, 0.925), 0, c("E10.5", "E11.5", "E12.5", "E13.5", "E14.5"),
				cex = ifelse(two, 0.75, 1))
			arctext("distal maxilla", c(0, 0), 1.1, middle = pi / 4, cex = ifelse(two, 1, 1.5))
			arctext("proximal maxilla", c(0, 0), 1.1, middle = 3 * pi / 4, cex = ifelse(two, 1, 1.5))
			arctext("proximal mandible", c(0, 0), 1.1, middle = 5 * pi / 4, clockwise = F, cex = ifelse(two, 1, 1.5))
			arctext("distal mandible", c(0, 0), 1.1, middle = 7 * pi / 4, clockwise = F, cex = ifelse(two, 1, 1.5))
			text(1.2, 1.95, expr$symbol, adj = 1, cex = ifelse(two, 2, 3))
			text(1.2, 1.85, expr$probeset, adj = c(1, 1), cex = ifelse(two, 1, 1.5))
			rect(seq(0.575, 1.245, 0.005), -1.9, seq(0.58, 1.25, 0.005), -1.75, border = NA,
				col = colorRampPalette(cols(c(min(all), max(all))))(135))
			segments(0.575 - min(all) / diff(range(all)) * (1.245 - 0.575), -1.9, 
				0.575 - min(all) / diff(range(all)) * (1.245 - 0.575), -1.75, col = "white", lwd = ifelse(two, 3, 4))
			text(c(0.63, 1.19), -1.95,
				format(round(c(min(all), max(all)), 2), nsmall = 2), cex = ifelse(two, 0.75, 1))
			text(0.915, -1.7, "Log2 Relative Expression", cex = ifelse(two, 0.75, 1))
			}}				
	plotPNG(bullseyes, "bullseyes.png", width = 2000, height = 1000, res = 300)
	bullseyes()
	}, height = function(){
				genes <- c(input$gene1, input$gene2, input$gene3, input$gene4, input$gene5, input$gene6)
				good.genes <- which(genes %in% dat$symbol | genes %in% dat$probeset)
				genes <- genes[good.genes]
				if (input$summary == "Z") num.probesets <- length(dat$probeset[dat$symbol %in% genes]) else
					num.probesets <- length(genes)
				if (num.probesets == 0) return(1300)
				if (num.probesets == 1) return(1300)
				if (num.probesets == 2) return(700)
				if (num.probesets %in% 3:4) return(1400)
				if (num.probesets > 4) return(1400 + (num.probesets - 4) * 450)})

	output$download.bullseyes <- downloadHandler("Facebase_Microarray_Bullseyes.png", 
	function(file) file.copy("bullseyes.png", file))	
	
}
		
shinyApp(ui = ui, server = server)