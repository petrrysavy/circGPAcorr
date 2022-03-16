# sbatch -p cpulong run.sh hsa_circ_0000228
args = commandArgs(trailingOnly=TRUE)
circname <- args[1]

source('../annotate.R', chdir=TRUE)
source('../goterms.R', chdir=TRUE)

# goTermNames <- c('GO:0045595', 'GO:0014013')
# annotateCirc('hsa_circ_0000228')

library("openxlsx")

library('pracma')
library('Rcpp')
library('ggpubr')
library('readr')
sourceCpp('bootstrapSpearmanExperiment.cpp')

# sampling a log-uniform distribution, then convert to integers
# ... this will be our points of interrest for the plot
samplePoints <- c(1:99, logspace(2, 6, 901))
samplePoints <- as.integer(samplePoints)

annotateLog <- function(muCirc, goMu, goM, Im) {
    ImmuCirc <- Im %*% muCirc
    score <- dot(muCirc, goMu) + dot(ImmuCirc, goM) # ignore the halfs, I love whole numbers
    pvalue <- pValueFull(muCirc, goMu, goM, ImmuCirc, score)
    if(score == 0) { bootstrap <- rep(1.0, length(samplePoints)) }
    else {  bootstrap <- pvalue_bootstrap_logging(muCirc, sum(goMu), ImmuCirc, sum(goM), score, samplePoints) }
    return(list(pvalue=pvalue, bootstrap=bootstrap))
}

annotateCirc <- function(circRNA) {
    muCirc <- buildMuCirc(circRNA)
    sapply(goTermNames, function(goTerm) {
        tryCatch({
            print(paste("Probing GO term:", goTerm))
            goMu <- buildGOMu(goTerm)
            goM <- buildGOM(goTerm)
            return(annotateLog(muCirc, goMu, goM, Im))
        }, error=function(cond) {
            print(cond)
            return(NULL)
        })
    })
}

res <- annotateCirc(circname)
saveRDS(res, file=paste(circname, ".RDS", sep=""))
# res <- readRDS('hsa_circ_0000228.RDS')
# sapply(res["bootstrap",], function(x) unlist(x)[999])
#  unlist(res["pvalue",])

#unlisted <- unlist(sapply(names(res), function(key) unlist(res[key])))


#pval <- unlist(sapply(names(res), function(key) as.numeric(unlisted[paste(key, '.', key, '.pvalue', sep="")])))
#pval <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key])[paste(key, '.pvalue', sep="")])))
pval <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key]$pvalue))))
nona <- !is.na(pval)
#pval <- unlist(res["pvalue",])

#getCorrelation <- function(i) {
#    boot <- sapply(res["bootstrap",], function(x) unlist(x)[i])
#    return(cor(pval, boot, method="spearman"))
#}

getCorrelation <- function(i) {
    print(i)
    boot <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
    return(cor(pval[nona], boot[nona], method="spearman"))
}

spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation(i)})
write_lines(spearmans, paste(circname, "-spearmans.dat", sep=""))
# write_lines(spearmans, paste("fsadfaa", ".dat", sep=""))

pval2 <- pval[pval<0.5]
nona2 <- !is.na(pval2) & !is.na(names(pval2))
pval2 <- pval2[nona2]
getCorrelation2 <- function(i) {
    print(i)
    boot <- unlist(sapply(names(pval2), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
    return(cor(pval2, boot, method="spearman"))
}
spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation2(i)})
write_lines(spearmans, paste(circname, "-spearman2.dat", sep=""))


pval3 <- pval[pval<0.05]
nona3 <- !is.na(pval3) & !is.na(names(pval3))
pval3 <- pval3[nona3]
getCorrelation3 <- function(i) {
    print(i)
    boot <- unlist(sapply(names(pval3), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
    return(cor(pval3, boot, method="spearman"))
}
spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation3(i)})
write_lines(spearmans, paste(circname, "-spearman3.dat", sep=""))

