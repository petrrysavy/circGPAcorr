library(xlsx)
library(readr)
library(stringr)

circs <- c("hsa_circ_0005232", "hsa_circ_0124919", "hsa_circ_0000994", "hsa_circ_0010023", "hsa_circ_0008101", "hsa_circ_0002059", "hsa_circ_0003889", "hsa_circ_0001523", "hsa_circ_0002212", "hsa_circ_0078150", "hsa_circ_0000551", "hsa_circ_0046760", "hsa_circ_0002711")
circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0004624", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")

sum <- 0
count <- 0
min <- 1e10
max <- 0
bootTsum <- 0
pvalTsum <- 0

for(circ in circs) {
    res <- read.xlsx(paste(circ, "-short.xlsx", sep=""), sheetIndex=1)
    pval <- as.numeric(res$pvalue)
    pvalT <- as.numeric(res$pvalueTime)
    bootT <- as.numeric(res$bootstrapTime)
    paste("outliers")
    print(sum(bootT > 300)) # outlier, RCI probably put the job away 
    bootTinlier <- bootT < 300
    bootT <- bootT[pval < 1.0 & bootTinlier]
    pvalT <- pvalT[pval < 1.0 & bootTinlier]
    bootTsum <- bootTsum + sum(bootT)
    pvalTsum <- pvalTsum + sum(pvalT)

    
    speedup <- bootT/pvalT
    sum <- sum + sum(speedup)
    min <- min(min, min(speedup))
    max <- max(max, max(speedup))
    count <- count + length(speedup)
    
    circ <- str_replace_all(circ, "_", "-")
    write_lines(speedup, paste(circ, "-pvaluespeedup.dat", sep=""))
}

print(sum/count)
# 2955.752
print(min)
print(max)
#15.88527
#6177351

print("Bootstrap overall time vs pvalue overall time")
print(bootTsum)
print(pvalTsum)
