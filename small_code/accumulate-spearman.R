library(openxlsx)
library(readr)
library(pracma)
library(stringr)

circs <- c("hsa_circ_0005232", "hsa_circ_0124919", "hsa_circ_0000994", "hsa_circ_0010023", "hsa_circ_0008101", "hsa_circ_0002059", "hsa_circ_0003889", "hsa_circ_0001523", "hsa_circ_0002212", "hsa_circ_0078150", "hsa_circ_0000551", "hsa_circ_0046760", "hsa_circ_0002711")
# circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0004624", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")
# circstex <- c("hsa-circ-0000005", "hsa-circ-0002816", "hsa-circ-0000006", "hsa-circ-0001897", "hsa-circ-0004624", "hsa-circ-0024604", "hsa-circ-0000228", "hsa-circ-0001540", "hsa-circ-0044708", "hsa-circ-0003583")
circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")
circstex <- c("hsa-circ-0000005", "hsa-circ-0002816", "hsa-circ-0000006", "hsa-circ-0001897", "hsa-circ-0024604", "hsa-circ-0000228", "hsa-circ-0001540", "hsa-circ-0044708", "hsa-circ-0003583")

samplePoints <- c(1:99, logspace(2, 6, 901))
samplePoints <- as.integer(samplePoints)

emptymatrix <- function() {
    m <- matrix(0.0, nrow=1000, ncol=10)
    m[,1] <- samplePoints
    colnames(m) <- c("points", circstex)
    return(m)
}

spearmans1 <- emptymatrix()
spearmans2 <- emptymatrix()
spearmans3 <- emptymatrix()

for(ci in 1:length(circs)) {
    circ <- circs[ci]
    res <- readRDS(paste(circ, ".RDS", sep=""))
    
    pval <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key]$pvalue))))
    nona <- !is.na(pval)
    
    getCorrelation <- function(i) {
        if( (i %% 10) == 0) { print(i) }
        boot <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
        boot <- (boot * samplePoints[i] + 1) / (samplePoints[i] + 1)
        return(cor(pval[nona], boot[nona], method="spearman"))
    }

    spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation(i)})
    spearmans1[,ci+1] <- spearmans
    write_lines(spearmans, paste(circ, "-spearmans.dat", sep=""))

    pval2 <- pval[pval<0.5]
    nona2 <- !is.na(pval2) & !is.na(names(pval2))
    pval2 <- pval2[nona2]
    getCorrelation2 <- function(i) {
        if( (i %% 10) == 0) { print(i) }
        boot <- unlist(sapply(names(pval2), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
        boot <- (boot * samplePoints[i] + 1) / (samplePoints[i] + 1)
        return(cor(pval2, boot, method="spearman"))
    }
    spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation2(i)})
    spearmans2[,ci+1] <- spearmans
    write_lines(spearmans, paste(circ, "-spearman2.dat", sep=""))


    pval3 <- pval[pval<0.05]
    nona3 <- !is.na(pval3) & !is.na(names(pval3))
    pval3 <- pval3[nona3]
    getCorrelation3 <- function(i) {
        if( (i %% 10) == 0) { print(i) }
        boot <- unlist(sapply(names(pval3), function(key) as.numeric(unlist(res[,key]$bootstrap[i]))))
        boot <- (boot * samplePoints[i] + 1) / (samplePoints[i] + 1)
        return(cor(pval3, boot, method="spearman"))
    }
    spearmans <- sapply(1:length(samplePoints), function(i) {getCorrelation3(i)})
    spearmans3[,ci+1] <- spearmans
    write_lines(spearmans, paste(circ, "-spearman3.dat", sep=""))
    
    print("Done with circ")
    print(spearmans1)
    print(spearmans2)
    print(spearmans3)
}

write.table(spearmans1, file="spearman1.dat", row.names=FALSE, quote=FALSE)
write.table(spearmans2, file="spearman2.dat", row.names=FALSE, quote=FALSE)
write.table(spearmans3, file="spearman3.dat", row.names=FALSE, quote=FALSE)
