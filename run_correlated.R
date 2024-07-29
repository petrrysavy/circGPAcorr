source('annotateCorrelated.R')
source('goterms.R')

source('weightdiscretization.R')

library(openxlsx)

args = commandArgs(trailingOnly=TRUE)
circrna <- args[1]

# Load the sample sheet
sampleSheet <- read.xlsx("../std/sampleSheet20211214.xlsx", sheet=1, rowNames=TRUE)
rownames(sampleSheet) <- sapply(rownames(sampleSheet), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# circRNA expression file
circRNA <- read.table(file='../std/200904_circRNA_fromRNAseq_annot_hg19.tsv', header=TRUE)
circrnaIndex <- which(circRNA$circRNA_ID == circrna)
circexpr <- circRNA[circrnaIndex, 10:ncol(circRNA)]
colnames(circexpr) <- sapply(colnames(circexpr), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])
circexpr[is.na(circexpr)] <- 0
rownames(circexpr) <- c(circrna)

# gene expression file
genes <- read.table(file='../std/200904_allRNA_fromRNAseq_annot_hg38.tsv', header=TRUE)
genesMatrix <- genes[,7:ncol(genes)]
genesMatrix <- mapply(genesMatrix, FUN=as.integer)
rownames(genesMatrix) <- genes[,6]
colnames(genesMatrix) <- sapply(colnames(genesMatrix), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# miRNA expression file
mirnaMatrix <- read.xlsx('../std2/final_all_samples_miRNA_seq.xlsx', sheet=1, rowNames=TRUE, colNames = TRUE)
colnames(mirnaMatrix) <- sapply(colnames(mirnaMatrix), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# filter so that the rows/cols are the same
validNames <- intersect(colnames(mirnaMatrix), intersect(colnames(genesMatrix), intersect(colnames(circexpr), rownames(sampleSheet))))
sampleSheet <- sampleSheet[rownames(sampleSheet) %in% validNames,]

sortAndFilter <- function(matrix) {
    matrix <- matrix[,colnames(matrix) %in% validNames]
    index = sapply(rownames(sampleSheet), function(x) which(colnames(matrix) == x))
    return(matrix[, index])
}

# order the rows/cols so that they are in the same order
mirnaMatrix <- sortAndFilter(mirnaMatrix)
genesMatrix <- sortAndFilter(genesMatrix)
circexpr <- sortAndFilter(circexpr)

# indices to filter samples which are WT vs SPL by 3 mutations column (1 vs 2)
filterInterestingIndex <- sampleSheet[,4] %in% c(1,2)

resultsToDataFrame <- function(results) {
    names <- slotNames('annotationResult')
    df <- sapply(names, function(name) lapply(results, function(r) { if(is.null(r)) NA else slot(r, name) }))
    # order the data frame
    df <- df[order(unlist(df[,"pvalue"])),]
    merge(df, as.data.frame(goTerms), by="row.names")
}

runOnCirc <- function(circRNA, goLower = 10, goUpper = 1000) {
    results <- annotateCircCorrelated(circRNA, goTermNames, goLower, goUpper, circexpr, mirnaMatrix, genesMatrix)
    df <- resultsToDataFrame(results)
    df <- df[order(unlist(df[,"pvalue"])),]
    df <- cbind(df, bonferroni=p.adjust(df$pvalue, method="bonferroni"), fdr=p.adjust(df$pvalue, method="fdr"))
    write.xlsx(df, paste(circRNA, "-correlated.xlsx", sep=""))
    
    df <- df[!is.na(df$pvalue) & df$goSize >= goLower & df$goSize <= goUpper,]
    df$bonferroni <- p.adjust(df$pvalue, method="bonferroni")
    df$fdr <- p.adjust(df$pvalue, method="fdr")
    if(nrow(df) >= 1) {
        write.xlsx(df, paste(circRNA, "-correlated-short.xlsx", sep=""))
    }
    
    results <- annotateCircCorrelated(circRNA, goTermNames, goLower, goUpper, circexpr[,filterInterestingIndex], mirnaMatrix[,filterInterestingIndex], genesMatrix[,filterInterestingIndex])
    df <- resultsToDataFrame(results)
    df <- df[order(unlist(df[,"pvalue"])),]
    df <- cbind(df, bonferroni=p.adjust(df$pvalue, method="bonferroni"), fdr=p.adjust(df$pvalue, method="fdr"))
    write.xlsx(df, paste(circRNA, "-correlated-filtered.xlsx", sep=""))
    
    df <- df[!is.na(df$pvalue) & df$goSize >= goLower & df$goSize <= goUpper,]
    df$bonferroni <- p.adjust(df$pvalue, method="bonferroni")
    df$fdr <- p.adjust(df$pvalue, method="fdr")
    if(nrow(df) >= 1) {
        write.xlsx(df, paste(circRNA, "-correlated-filtered-short.xlsx", sep=""))
    }
}

runOnCirc(circrna)

