
# Loads R packages:
library(bioDist) 
library(gplots) 
library(limma) 
library(DT)

## these are for heatplot
#library(ggplot2)
#library(made4)

## base on Dan's app.R for shiny server 

processCELdata_f <- function(inputCONFIG) {

leftstr <- function(x, n = 1) substr(x, 1, min(n, nchar(x)))
rightstr <- function(x, n = 1) substr(x, nchar(x) - min(n - 1, nchar(x) - 1), nchar(x))
plot.null <- function()
plot(0, 0, col = "transparent", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
mousecase <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))


# This sets up the contrasts in the data. Five ages, two bones, 
# two regions per bone = 20 conditions.
# Three replicates per condition = 60 samples.

age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
full.design <- model.matrix(~age + bone + place)
rownames(full.design) <- paste0(age, bone, place, 1:3)


# This design matrix is for later when we do regressions on each row
# of the data table, using only certain columns.


# This just opens the file from the web.
repeat{
      con <- url(inputCONFIG$inputFile)
      t <- try(load(con))
      close(con)
      if (!inherits(t, "try-error")) break}
rownames(dat) <- dat[, "probeset"]
genes.tab <- xtabs(~unlist(dat$symbol))
single.genes <- names(genes.tab)[genes.tab == 1]

# This identifies which samples the user has selected.
# The samples called "ones" will be contrasted with those called "twos".
# The rest of the code here tries to make sure that you’re comparing 
# like with like, e.g. comparing two against two, not two against three, etc.

sel <- inputCONFIG$sel
sel <- substr(sel, 2, nchar(sel))
sel <- gsub("_", "", sel)
sels <- character()
for (s in sel) sels <- c(sels, paste0(s, 1:3))
## if the comparison is "place" (distal vs. proximal):
if (inputCONFIG$comp == "place"){ 
target.col <- "placeD"
ones <- grep("P", sels, value = T)
twos <- grep("D", sels, value = T)
ones.reduced <- gsub("P", "", ones)
twos.reduced <- gsub("D", "", twos)
ones <- ones[ones.reduced %in% twos.reduced]
twos <- twos[twos.reduced %in% ones.reduced]
invert <- ifelse(inputCONFIG$invert_place == "inverted", T, F)
one <- ifelse(invert, "distal", "proximal")
two <- ifelse(invert, "proximal", "distal")}

## if the comparison is "bone" (maxilla vs. mandible):

if (inputCONFIG$comp == "bone"){
target.col <- "boneMnd"
ones <- grep("Max", sels, value = T)
twos <- grep("Mnd", sels, value = T)
ones.reduced <- gsub("Max", "", ones)
twos.reduced <- gsub("Mnd", "", twos)
ones <- ones[ones.reduced %in% twos.reduced]
twos <- twos[twos.reduced %in% ones.reduced]
invert <- ifelse(inputCONFIG$invert_bone == "inverted", T, F)
one <- ifelse(invert, "mandible", "maxilla")
two <- ifelse(invert, "maxilla", "mandible")}

## if the comparison is "age" (e.g. E10.5 vs. E11.5, or E12.5 vs. E14.5):
ages <- unique(substr(sels, 1, 4))
if (length(ages)) age.1 <- min(ages) else age.1 <- 10.5
young.col <- paste0("age", age.1)
if (inputCONFIG$comp == "age"){
age.2 <- max(ages)
target.col <- paste0("age", age.2)
ones <- grep(age.1, sels, value = T)
twos <- grep(age.2, sels, value = T)
ones.reduced <- gsub(age.1, "", ones)
twos.reduced <- gsub(age.2, "", twos)
ones <- ones[ones.reduced %in% twos.reduced]
twos <- twos[twos.reduced %in% ones.reduced]
invert <- ifelse(inputCONFIG$invert_age == "inverted", T, F)
one <- ifelse(invert, paste0("E", age.2), paste0("E", age.1))
two <- ifelse(invert, paste0("E", age.1), paste0("E", age.2))}

# All of the above section is just getting "ones" and "twos" to 
# represent the right samples.
# If I’m going to compare E10.5 distal mandible with E10.5 proximal 
# mandible (3 samples vs. 3 samples), this is simple.
# But if I wanted to compare distal mandible vs. proximal mandible 
# using all age stages, I want to control for age.
# The function lmFit() is going to do the regression models, but 
# you have to give it all the necessary information.
# So "constant" are the contrasts that don’t vary within the whole 
# comparison (e.g. Age, if all samples are E10.5).
# "Control" are the columns that do vary, but are not the primary 
# one/two contrast.


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

# This just means don’t show the MA plot or heat map if something 
# has gone wrong:

if (!nrow(design) | !target.col %in% colnames(design) | !all(sels %in% c(ones, twos))) return({
print("BAD BAD BAD")
print(nrow(design))
print(target.col)
print(colnames(design))
print(c(ones,twos))
print(sels)
data.frame(NULL)})


# "dat.sel" is just the columns we will be comparing:

dat.sel <- dat[, c(ones, twos)]
cn <- colnames(dat.sel)


# There is a user option for how to summarize probesets. Many 
# genes have more than one probeset,
# so we have more than one measurement of a gene expression level
# for some genes.
# The default is to show all the probesets separately.
# If summary = "A" then we choose, for each gene, the probeset 
# with the highest expression:


if (inputCONFIG$summary == "A"){
dat.sel$A <- rowMeans(dat.sel)
dat.sel$P <- rownames(dat.sel)
dat.sel$symbol <- unlist(dat$symbol)
dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
dat.agg <- aggregate(A ~ symbol, dat.multiple, max)
dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
rownames(dat.multiple) <- dat.multiple$P
dat.sel <- rbind(dat.single, dat.multiple)
dat.sel <- dat.sel[, cn]}


# If summary = "M" then we choose, for each gene, the probeset 
# with the most DIFFERENTIAL expression
# that is, the biggest difference between the ones and the twos.
# (People don’t usually use this.)

if (inputCONFIG$summary == "M"){
dat.sel$MM <- abs(rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones]))
dat.sel$P <- rownames(dat.sel)
dat.sel$symbol <- unlist(dat$symbol)
dat.single <- dat.sel[dat.sel$symbol %in% single.genes, ]
dat.multiple <- dat.sel[!dat.sel$symbol %in% single.genes, ]
dat.agg <- aggregate(MM ~ symbol, dat.multiple, max)
dat.multiple <- merge(dat.agg, dat.multiple, sort = F)
rownames(dat.multiple) <- dat.multiple$P
dat.sel <- rbind(dat.single, dat.multiple)
dat.sel <- dat.sel[, cn]}

# Now we fit the regression models with lmFit(). The variable A
# refers to the average between all samples in the comparison 
# (ones and twos) and M refers to the difference between the 
# average of ones and the average of twos.

fit <- lmFit(dat.sel, design)

efit <- eBayes(fit)

dat.sel$A <- rowMeans(dat.sel)

dat.sel$M <- rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones])


# This part has to do with coloring and details of the output.

if (invert) dat.sel$M <- -dat.sel$M

dat.sel$symbol <- unlist(sapply(rownames(dat.sel), function(x) dat$symbol[dat$probeset == x]))

dat.sel$symbol[is.na(dat.sel$symbol)] <-  rownames(dat.sel)[is.na(dat.sel$symbol)]

dat.sel$color <- "black"

dat.sel$color[dat.sel$symbol == mousecase(inputCONFIG$gene1) | rownames(dat.sel) == mousecase(inputCONFIG$gene1)] <- tolower(inputCONFIG$col1)

dat.sel$color[dat.sel$symbol == mousecase(inputCONFIG$gene2) | rownames(dat.sel) == mousecase(inputCONFIG$gene2)] <- tolower(inputCONFIG$col2)

dat.sel$color[dat.sel$symbol == mousecase(inputCONFIG$gene3) | rownames(dat.sel) == mousecase(inputCONFIG$gene3)] <- tolower(inputCONFIG$col3)

dat.sel$color[dat.sel$symbol == mousecase(inputCONFIG$gene4) | rownames(dat.sel) == mousecase(inputCONFIG$gene4)] <- tolower(inputCONFIG$col4)

dat.sel$color[dat.sel$symbol == mousecase(inputCONFIG$gene5) | rownames(dat.sel) == mousecase(inputCONFIG$gene5)] <- tolower(inputCONFIG$col5)


# Based on user input, we have "lfc" which is a fold-change cut-off
# for differentially expressed genes.
# Also based on user input, we have "fdr" which is a false discovery 
# rate from the Benjamini-Hochberg correction.
# Notice that the FDR p-values are the only thing that we get out of 
# the regression models.


lfc <- ifelse(inputCONFIG$log, abs(as.numeric(inputCONFIG$fc)), log2(abs(as.numeric(inputCONFIG$fc))))

# "top" refers to the table in the bottom right of the output. 
# This is the primary output.
# The MA plot is the second most important output.
# The heat map is the third most important output.

# This part is just formatting the "top" table, which is generated 
# by topTable() and the fit or efit object from lmFit():

top <- topTable(efit, coef = which(colnames(design) == target.col), 
              lfc = lfc, p.value = as.numeric(inputCONFIG$fdr), num = Inf)

top.colored <- top[rownames(top) %in% dat.sel$probeset[dat.sel$color != "black"], ]

if (nrow(top.colored)){

top.black <- top[1:min(nrow(top), as.numeric(inputCONFIG$max)), ]
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
top$rownames <- row.names(top)
top <- top[, c("rownames", "symbol", "AveExpr", "logFC", "FC", "adj.P.Val", "desc")]
colnames(top) <- c("Probeset", "Gene", "Ave.Expr", "Log2FC", "Fold.Change", "Adj.P.Val", "Description")}
dat.top <- dat.sel[rownames(dat.sel) %in% top$Probeset, ]

### renderPlot, ma.plot function 
{
    par(mar = c(5, 5, 5, 8))
    lim <- max(abs(range(dat.sel$M)))
    
    plot(NULL, cex = 0.9,
        xaxt = "n", yaxt = "n",
        xlab = "Average Expression",
        xlim = range(dat.sel$A),
        ylim = c(-lim * 1.05, lim * 1.05), 
        ylab = paste0(one,
        paste(rep(" ", ifelse(inputCONFIG$log, 15, 20)), collapse = ""),
        ifelse(inputCONFIG$log, "Log2 Fold Change", "Fold Change"),
        paste(rep(" ", ifelse(inputCONFIG$log, 15, 20)), collapse = ""), two),
    
        cex.lab = 1.2, 
        font.main = 1,
        cex.main = 1.2,
        main = paste(unique(gsub("1$|2$|3$", "", c(ones, twos))), collapse = " ")
    )

    points(dat.sel$A[dat.sel$color == "black"],
           dat.sel$M[dat.sel$color == "black"],
           pch = 21, cex = 0.9)
    points(dat.sel$A[dat.sel$color != "black"],
           dat.sel$M[dat.sel$color != "black"],
           pch = 21, cex = 0.9,
           col = dat.sel$color[dat.sel$color != "black"])

    if (nrow(dat.top)){
     points(dat.top$A, dat.top$M, pch = 21, 
            cex = 0.9, col = dat.top$color, bg = dat.top$color)
     text(dat.top$A, dat.top$M, dat.top$symbol, col = dat.top$color, 
      pos = ifelse(dat.top$M > 0, 3, 1), offset = 0.4, cex = 0.9)
    }

#MEI, JSON file
generateMAPLOTjson_f(ones,twos,inputCONFIG,dat.sel,dat.top)

    ats.0 <- seq(1, 9, 1)
    ats <- c(-1 * rev(ats.0), 0, ats.0)

    log.labs <- ats
    log.labs[log.labs > 0] <- paste0("+", log.labs[log.labs > 0])
    raw.labs <- 2 ^ ats
    raw.labs[raw.labs < 1] <- -1 / raw.labs[raw.labs < 1]
    raw.labs[raw.labs > 1] <- paste0("+", raw.labs[raw.labs > 1])
    axis(1, seq(0, 20, 1))
    if (inputCONFIG$log) axis(2, ats, log.labs)
    else axis(2, ats, raw.labs, cex.axis = ifelse(lim > 6, 0.75, 1))
    abline(h = lfc * c(-1, 1), lty = 2)
} 


# The next part is the heat map. Part of this involves normalizing the data.
# Sometimes people standardize each row (gene) so that the ups and downs (M)
# can be visualized apart from the differences in A (which don’t mean much).
# I did this in a different way. I think the user should have the option
# of whether to scale the data for the heat map. Another issue is how to
# cluster the rows and columns, there are several options for the distance
# matrix, like Euclidean, Pearson correlation, and absolute Pearson correlation.
# Again I think this could be a user choice. (I used Pearson correlation,
# which makes sense for samples (columns). For rows, I’d prefer 
# absolute Pearson correlation, but heatmap.2() doesn’t allow a different 
# distance metric for rows and columns…


# The function heatmap.2() is from the gplots package. It draws the 
# heat map, does the clustering, and prints the dendrograms, etc.

dat.heat <- as.matrix(dat.top[, sels])

for (col in colnames(control)){
for (v in unique(control[, col])){
mm <- mean(dat.heat[, rownames(control)[control[, col] == v]])
dat.heat[, rownames(control)[control[, col] == v]] <- 
dat.heat[, rownames(control)[control[, col] == v]] - mm}}
dat.heat <- sweep(dat.heat, 1, rowMeans(dat.heat))

weights <- (ncol(dat.heat):1) + 10000
weights[colnames(dat.heat) %in% twos] <- 
  weights[colnames(dat.heat) %in% twos] + ifelse(invert, -100, 100)

if (inputCONFIG$log) { 
   extreme <- max(abs(dat.top$M)) / 2 
} else { extreme <- 2 ^ max(abs(dat.top$M)) / 2 }

extreme <- ceiling(10 * extreme) / 10

#MEI, make JSON file
generateHEATMAPjson_f(ones,twos,inputCONFIG,dat.top,dat.heat)

heatmap.2(t(dat.heat),
          Rowv = as.integer(weights),
          col = ifelse(inputCONFIG$heatcol == "rg", redgreen, greenred),
          cexCol = min(1.5, 0.2 + 1/log10(nrow(dat.heat))),
          scale = "none", key = T, key.title = NA, lwid = c(1, 4), lhei = c(1.5, 4),
          density.info = "density", densadj = 0.5, denscol = "white",
          key.ylab = NA, key.xlab = ifelse(inputCONFIG$log, "Log2 Fold Change", "Fold Change"),
          key.par = list(cex.lab = 1.25, cex.axis = 1, tcl = -0.35, mgp = c(2, 1, 0)),
          key.ytickfun = function() list(labels = F, tick = F),
          key.xtickfun = function(){
          breaks <- parent.frame()$breaks
          list(at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
          labels = round(seq(-extreme, extreme, length = 6), 1), 
          padj = -0.25)},
          trace = "none", margins = c(5, 10),
          distfun = function(...) cor.dist(..., abs = F),
          labCol = dat.top$symbol, colCol = dat.top$color)

} ## end of processCELdataFunction 
