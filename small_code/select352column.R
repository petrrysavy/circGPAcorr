library('pracma')

res <- readRDS('hsa_circ_0000228.RDS')

evaluatePoint <- function(point) {
    print(paste("Evaluate point", point))
    samplePoints <- c(1:99, logspace(2, 6, 901))
    samplePoints <- as.integer(samplePoints)
    index352 <- which(samplePoints == point)

    pval <- unlist(sapply(colnames(res), function(key) as.numeric(unlist(res[,key]$pvalue))))
    nona <- !is.na(pval)

    pval3 <- pval[pval<0.05]
    nona3 <- !is.na(pval3) & !is.na(names(pval3))
    pval3 <- pval3[nona3]
    pval3 <- (pval3 * 1e6 + 1) / 1000001

    boot <- unlist(sapply(names(pval3), function(key) as.numeric(unlist(res[,key]$bootstrap[index352]))))
    boot <- (boot * point + 1) / (point + 1)

    diff <- abs(boot - pval3) / pval3

    print("All data")
    print(max(diff))
    print(mean(diff))
    print(sd(diff))
    print(median(diff))

    pval3 <- sort(pval3)
    boot <- boot[names(pval3)]

    pval3x <- pval3[21:length(pval3)]
    bootx <- boot[21:length(pval3)]

    diff <- abs(bootx - pval3x) / pval3x

    print("Without 20 first")
    print(max(diff))
    print(mean(diff))
    print(sd(diff))
    print(median(diff))
}

evaluatePoint(352)
evaluatePoint(1000000)
