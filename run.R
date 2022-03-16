args = commandArgs(trailingOnly=TRUE)
circname <- args[1]

source('annotate.R')
source('goterms.R')

library("openxlsx")

resultsToDataFrame <- function(results) {
    names <- slotNames('annotationResult')
    df <- sapply(names, function(name) lapply(results, function(r) { if(is.null(r)) NA else slot(r, name) }))
    # order the data frame
    #df <- df[is.null(df[,"pvalue"]),]
    df <- df[order(unlist(df[,"pvalue"])),]
    merge(df, as.data.frame(goTerms), by="row.names")
}

runOnCirc <- function(circRNA, goLower = 10, goUpper = 1000) {
    results <- annotateCirc(circRNA, goTermNames)
    df <- resultsToDataFrame(results)
    df <- df[order(unlist(df[,"pvalue"])),]
    df <- cbind(df, bonferroni=p.adjust(df$pvalue, method="bonferroni"), fdr=p.adjust(df$pvalue, method="fdr"))
    write.xlsx(df, paste(circRNA, ".xlsx", sep=""))
    
    df <- df[!is.na(df$pvalue) & df$goSize >= goLower & df$goSize <= goUpper,]
    df$bonferroni <- p.adjust(df$pvalue, method="bonferroni")
    df$fdr <- p.adjust(df$pvalue, method="fdr")
    if(nrow(df) >= 1) {
        write.xlsx(df, paste(circRNA, "-short.xlsx", sep=""))
    }
}

#runOnCirc("hsa_circ_0007694")
#runOnCirc("hsa_circ_0000228")
#runOnCirc("hsa_circ_0003793")

runOnCirc(circname)

