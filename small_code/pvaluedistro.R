library(xlsx)
library(readr)
library(stringr)

circs <- c("hsa_circ_0005232", "hsa_circ_0124919", "hsa_circ_0000994", "hsa_circ_0010023", "hsa_circ_0008101", "hsa_circ_0002059", "hsa_circ_0003889", "hsa_circ_0001523", "hsa_circ_0002212", "hsa_circ_0078150", "hsa_circ_0000551", "hsa_circ_0046760", "hsa_circ_0002711")
circs <- c("hsa_circ_0000005", "hsa_circ_0002816", "hsa_circ_0000006", "hsa_circ_0001897", "hsa_circ_0004624", "hsa_circ_0024604", "hsa_circ_0000228", "hsa_circ_0001540", "hsa_circ_0044708", "hsa_circ_0003583")


for(circ in circs) {
    res <- read.xlsx(paste(circ, "-short.xlsx", sep=""), sheetIndex=1)
    pval <- res$pvalue
    circ <- str_replace_all(circ, "_", "-")
    write_lines(pval, paste(circ, "-pvaluedistro.dat", sep=""))
}
