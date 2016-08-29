library(bioDist)
library(gplots)
library(limma)
library(DT)

leftstr <- function(x, n = 1) substr(x, 1, min(n, nchar(x)))
rightstr <- function(x, n = 1) substr(x, nchar(x) - min(n - 1, nchar(x) - 1), nchar(x))
plot.null <- function()
	plot(0, 0, col = "transparent", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
mousecase <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))

ui <- fluidPage(
  titlePanel("Gene Expression Comparison App"),
  sidebarLayout(
  div(class = "col-sm-5",
  tags$form(class = "well",
  	selectInput("dataset", "Dataset", "GSE67985", "GSE67985", width = "30%"),
  		div(style="padding-bottom:3px", "Chai, Yang; Parada, Carolina; Grimaldi, Alexandre; Ho, Thach-Vu; Harunaga, Jill; Samuels, Bridget."),
   		div(style="padding-bottom:3px", "Distal/proximal mandible/maxilla at embryonic day 10.5/11.5/12.5/13.5/14.5, 3 replicates/condition."),
   		div(style="padding-bottom:7px", "Affymetrix Mouse Genome 430 2.0 Array."),
 	conditionalPanel("!output.datloaded",
    		helpText("Data is loading (may take a minute) ...")),
    conditionalPanel("output.datloaded", 
    		helpText("Data is loaded."))),
  tags$form(class = "well",
  	checkboxGroupInput("sel", "Select",
      choices = c(
        "E10.5_Mnd_D", "E11.5_Mnd_D", "E12.5_Mnd_D", "E13.5_Mnd_D", "E14.5_Mnd_D",
        "E10.5_Mnd_P", "E11.5_Mnd_P", "E12.5_Mnd_P", "E13.5_Mnd_P", "E14.5_Mnd_P",
        "E10.5_Max_D", "E11.5_Max_D", "E12.5_Max_D", "E13.5_Max_D", "E14.5_Max_D",
        "E10.5_Max_P", "E11.5_Max_P", "E12.5_Max_P", "E13.5_Max_P", "E14.5_Max_P"), inline = T)),
  tags$form(class = "well",
  splitLayout(cellWidths = c("35%", "65%"),
    radioButtons("comp", "Compare",
    		choices = c(
    			"proximal vs. distal" = "place",
    			"maxilla vs. mandible" = "bone",
    			"age (earliest vs. latest selected)" = "age"),
    		selected = "place"),
    list(conditionalPanel("input.comp == 'place'", 
	    	radioButtons("invert_place", br(), c("proximal down, distal up" = "normal", 
    				"distal down, proximal up" = "inverted"), "normal", inline = T)),
		conditionalPanel("input.comp == 'bone'",
    		radioButtons("invert_bone", list(br(), br()), c("maxilla down, mandible up" = "normal", 
    				"mandible down, maxilla up" = "inverted"), "normal", inline = T)),
		conditionalPanel("input.comp == 'age'",
    		radioButtons("invert_age", list(br(), br(), br(), br()),
    				c("earliest down, latest up" = "normal", 
    				"latest down, earliest up" = "inverted"), "normal", inline = T))))),
  tags$form(class= "well",
      bootstrapPage(
    	div(style="display:inline-block; width:30%", textInput("fc", "Fold change cut-off", 2)),
        div(style="display:inline-block", checkboxInput("log", "log2", T))),    
    textInput("fdr", "False discovery rate", 0.01, width = "30%"),
    textInput("max", "Max. DE genes/probesets", Inf, width = "30%"),
    radioButtons("summary", "For each gene, show",
    		choices = c("all probesets" = "Z", "most highly expressed probeset" = "A",
    					"most differentially expressed probeset" = "M"), selected = "Z", inline = T),
    radioButtons("heatcol", "For heatmap, high values are",
    		choices = c("green" = "rg", "red" = "gr"), selected = "rg", inline = T)),
  tags$form(class = "well",
  	  div(tags$label("Highlight any gene/probeset in any color")),
  	  bootstrapPage(
  		div(style="display:inline-block; width:18%", textInput("gene1", NULL, "name")),
  		div(style="display:inline-block; width:18%", textInput("gene2", NULL, "name")),
  		div(style="display:inline-block; width:18%", textInput("gene3", NULL, "name")),
   		div(style="display:inline-block; width:18%", textInput("gene4", NULL, "name")),
  		div(style="display:inline-block; width:18%", textInput("gene5", NULL, "name"))),
      bootstrapPage(
       div(style="display:inline-block; width:18%", textInput("col1", NULL, "blue")),
	   div(style="display:inline-block; width:18%", textInput("col2", NULL, "red")),
	   div(style="display:inline-block; width:18%", textInput("col3", NULL, "forestgreen")),
	   div(style="display:inline-block; width:18%", textInput("col4", NULL, "magenta")),
	   div(style="display:inline-block; width:18%", textInput("col5", NULL, "darkorange"))))),
  mainPanel(width = 7, 
		downloadButton("download.heatmap", "Download heatmap"),
		plotOutput("heatmap", height = "450px"), br(),
		downloadButton("download.ma.plot", "Download MA plot"),
		plotOutput("ma.plot", height = "500px"), br(),
    	downloadButton("download.table", "Download table"),	br(), br(),
    	dataTableOutput("table", width = "90%"))))

server <- function(input, output){
	age <- factor(rep(c(10.5, 11.5, 12.5, 13.5, 14.5), each = 12))
	bone <- factor(rep(rep(c("Max", "Mnd"), each = 6), 5))
	place <- relevel(factor(rep(rep(c("D", "P"), each = 3), 10)), "P")
	full.design <- model.matrix(~age + bone + place)
	rownames(full.design) <- paste0(age, bone, place, 1:3)
	repeat{
		con <- url("http://www.danielezrajohnson.com/exprs.R")
		t <- try(load(con))
		close(con)
		if (!inherits(t, "try-error")) break}
	
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
			one <- ifelse(invert, "mandible", "maxilla")
			two <- ifelse(invert, "maxilla", "mandible")}
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
		dat.sel <- dat[, c(ones, twos)]
		cn <- colnames(dat.sel)
		if (input$summary == "A"){
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
		if (input$summary == "M"){
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
		fit <- lmFit(dat.sel, design)
		efit <- eBayes(fit)
		dat.sel$A <- rowMeans(dat.sel)
		dat.sel$M <- rowMeans(dat.sel[, twos]) - rowMeans(dat.sel[, ones])
		if (invert) dat.sel$M <- -dat.sel$M
		dat.sel$symbol <- unlist(sapply(rownames(dat.sel), function(x) dat$symbol[dat$probeset == x]))
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
		lfc <- ifelse(input$log, abs(as.numeric(input$fc)), log2(abs(as.numeric(input$fc))))
		top <- topTable(efit, coef = which(colnames(design) == target.col), 
			lfc = lfc, p.value = as.numeric(input$fdr), num = Inf)
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
			plotPNG(ma.plot, "ma.plot.png", width = 2000, height = 1000, res = 300)
			ma.plot()
		})
			
		output$download.ma.plot <- downloadHandler("Facebase_Microarray_MA_Plot.png", 
			function(file) file.copy("ma.plot.png", file))	

		output$heatmap <- renderPlot({
			heatmap <- function(){
			withProgress({
				if (nrow(dat.top) < 2) return(plot.null())
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
				if (input$log) extreme <- max(abs(dat.top$M)) / 2
				else extreme <- 2 ^ max(abs(dat.top$M)) / 2
				extreme <- ceiling(10 * extreme) / 10
				heatmap.2(t(dat.heat),
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
					distfun = function(...) cor.dist(..., abs = F),
					labCol = dat.top$symbol, colCol = dat.top$color)
			incProgress(1)}, value = NULL, message = "Calculations in progress...")}
			plotPNG(heatmap, "heatmap.png", width = 2000, height = 1000, res = 300)
			heatmap()
		})
			
		output$download.heatmap <- downloadHandler("Facebase_Microarray_Heatmap.png",
			function(file) file.copy("heatmap.png", file))
					
		rownames(top) <- NULL
		output$download.table <- downloadHandler("Facebase_Microarray_Table.csv",
			function(file) write.csv(as.data.frame(top), file, row.names = F))
		top$Color <- unlist(sapply(top$Probeset, function(x) dat.top$color[rownames(dat.top) == x]))
		top <- formatStyle(datatable(top), "Color", target = "row", 
			color = styleEqual(
				c(0, input$col1, input$col2, input$col3, input$col4, input$col5),
				c("white", input$col1, input$col2, input$col3, input$col4, input$col5)))
	})	

	}
		
shinyApp(ui = ui, server = server)