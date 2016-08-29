
setConfig_f <- function() {

localCONFIG<-list()

localCONFIG$inputFile <- "http://localhost/data/CEL/danfile.R"
## "E10.5_Mnd_D", "E11.5_Mnd_D", "E12.5_Mnd_D", "E13.5_Mnd_D", "E14.5_Mnd_D",
## "E10.5_Mnd_P", "E11.5_Mnd_P", "E12.5_Mnd_P", "E13.5_Mnd_P", "E14.5_Mnd_P",
## "E10.5_Max_D", "E11.5_Max_D", "E12.5_Max_D", "E13.5_Max_D", "E14.5_Max_D",
## "E10.5_Max_P", "E11.5_Max_P", "E12.5_Max_P", "E13.5_Max_P", "E14.5_Max_P"),
localCONFIG$sel <- c("E10.5_Mnd_D","E10.5_Mnd_P")
## "place", "age", "bone"
localCONFIG$comp <- "place"
## "normal", "inverted"
localCONFIG$invert_place <- "normal"
localCONFIG$invert_bone <- "normal"
localCONFIG$invert_age <- "normal"
## "all probesets" = "Z",
## "most highly expressed probeset" = "A",
## "most differentially expressed probeset" = "M")
localCONFIG$summary <- "Z"
localCONFIG$fdr <- 0.01
localCONFIG$fc <- 2.0
localCONFIG$max <- "Inf"
## "log"(0), "log2"(1) 
localCONFIG$log <- 1
## "rg", "gr"
localCONFIG$heatcol <- "rg"
# user specified color to particular gene
localCONFIG$col1 <- ""
localCONFIG$col2 <- ""
localCONFIG$col3 <- ""
localCONFIG$col4 <- ""
localCONFIG$col5 <- ""
localCONFIG$gene1 <- ""
localCONFIG$gene2 <- ""
localCONFIG$gene3 <- ""
localCONFIG$gene4 <- ""
localCONFIG$gene5 <- ""

return(localCONFIG)
}

