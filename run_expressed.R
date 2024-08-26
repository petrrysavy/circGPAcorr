source('annotateExpressed.R')
source('goterms.R')

source('weightdiscretization.R')

library(openxlsx)

args = commandArgs(trailingOnly=TRUE)
circrna <- args[1]

# Load the sample sheet
sampleSheet <- read.xlsx("data/sampleSheet20211214.xlsx", sheet=1, rowNames=TRUE)
rownames(sampleSheet) <- sapply(rownames(sampleSheet), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# circRNA expression file
circRNA <- read.table(file='data/200904_circRNA_fromRNAseq_annot_hg19.tsv', header=TRUE)
# circrnaIndex <- which(circRNA$circRNA_ID == circrna)
# circexpr <- circRNA[circrnaIndex, 10:ncol(circRNA)]
# colnames(circexpr) <- sapply(colnames(circexpr), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])
# circexpr[is.na(circexpr)] <- 0
# rownames(circexpr) <- c(circrna)
usablecircs <- !is.na(circRNA$circRNA_ID)
circMatrix <- circRNA[usablecircs, 10:ncol(circRNA)]
colnames(circMatrix) <- sapply(colnames(circMatrix), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])
circMatrix[is.na(circMatrix)] <- 0
rownames(circMatrix) <- c(circRNA$circRNA_ID[usablecircs])

# gene expression file
genes <- read.table(file='data/200904_allRNA_fromRNAseq_annot_hg38.tsv', header=TRUE)
genesMatrix <- genes[,7:ncol(genes)]
genesMatrix <- mapply(genesMatrix, FUN=as.integer)
rownames(genesMatrix) <- genes[,6]
colnames(genesMatrix) <- sapply(colnames(genesMatrix), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# miRNA expression file
mirnaMatrix <- read.xlsx('data/final_all_samples_miRNA_seq.xlsx', sheet=1, rowNames=TRUE, colNames = TRUE)
colnames(mirnaMatrix) <- sapply(colnames(mirnaMatrix), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])

# filter so that the rows/cols are the same
validNames <- intersect(colnames(mirnaMatrix), intersect(colnames(genesMatrix), intersect(colnames(circMatrix), rownames(sampleSheet))))
sampleSheet <- sampleSheet[rownames(sampleSheet) %in% validNames,]

sortAndFilter <- function(matrix) {
    matrix <- matrix[,colnames(matrix) %in% validNames]
    index = sapply(rownames(sampleSheet), function(x) which(colnames(matrix) == x))
    return(matrix[, index])
}

# order the rows/cols so that they are in the same order
mirnaMatrix <- sortAndFilter(mirnaMatrix)
genesMatrix <- sortAndFilter(genesMatrix)
circMatrix <- sortAndFilter(circMatrix)
# circexpr <- circMatrix[which(rownames(circMatrix) == circrna), ]

# indices to filter samples which are WT vs SPL by 3 mutations column (1 vs 2)
filterInterestingIndex <- sampleSheet[,4] %in% c(1,2)

# NOW Deseq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mirnaMatrix, colData=sampleSheet, design = ~0 + GROUP)
dds <- DESeq(dds)#, betaPrior=TRUE)
resMiRNA <- results(dds, contrast= list(c('GROUPSPL'), c("GROUPWT")))
# mirnaExpressionVector <- abs(resMiRNA[all_mirna_list, "log2FoldChange"])
# mirnaExpressionVector[is.na(mirnaExpressionVector)] <- 0.0

dds <- DESeqDataSetFromMatrix(countData = genesMatrix, colData=sampleSheet, design = ~0 + GROUP)
dds <- DESeq(dds)#, betaPrior=TRUE)
resMRNA <- results(dds, contrast= list(c('GROUPSPL'), c("GROUPWT")))
# mrnaExpressionVector <- abs(resMRNA[all_gene_list, "log2FoldChange"])
# mrnaExpressionVector[is.na(mrnaExpressionVector)] <- 0.0

dds <- DESeqDataSetFromMatrix(countData = circMatrix, colData=sampleSheet, design = ~0 + GROUP)
dds <- DESeq(dds)#, betaPrior=TRUE)
resCirc <- results(dds, contrast= list(c('GROUPSPL'), c("GROUPWT")))

#circrnaIndex <- which(rownames(circMatrix) == circrna)

resultsToDataFrame <- function(results) {
    names <- slotNames('annotationResultCorrelated')
    df <- sapply(names, function(name) lapply(results, function(r) { if(is.null(r)) NA else slot(r, name) }))
    # order the data frame
    df <- df[order(unlist(df[,"pvalueUB"])),]
    merge(df, as.data.frame(goTerms), by="row.names")
}

runOnCirc <- function(circRNA, goLower = 10, goUpper = 1000) {
    results <- annotateCircExpressed(circRNA, goTermNames, resCirc[circrna,"log2FoldChange"], resMiRNA, resMRNA, goLower, goUpper, pvalueThreshold=0.01, trialsMin=100, trialsMax=10000)
    print("annotation finished")
    df <- resultsToDataFrame(results)
    print("results in data frame")
    df <- df[order(unlist(df[,"pvalueUB"])),]
    print("ordered")
    df <- cbind(df, bonferroni=p.adjust(df$pvalueUB, method="bonferroni"), fdr=p.adjust(df$pvalueUB, method="fdr"))
    print("pvalue adjustment")
    df2 <- df
    df2[is.na(df2)]<-'NA'
    write.xlsx(df2, paste(circRNA, "-expressed.xlsx", sep=""), overwrite=TRUE)
    print("output written")
    
    df <- df[!is.na(df$pvalueUB) & df$goSize >= goLower & df$goSize <= goUpper,]
    print("filtered shorter list")
    df$bonferroni <- p.adjust(df$pvalueUB, method="bonferroni")
    print("pvalue adjustment")
    df$fdr <- p.adjust(df$pvalueUB, method="fdr")
    if(nrow(df) >= 1) {
        df[is.na(df)]<-'NA'
        write.xlsx(df, paste(circRNA, "-expressed-short.xlsx", sep=""), overwrite=TRUE)
        print("output written")
    }
}

runOnCirc(circrna)

