# reads a xlsx file provided by the circGPA run and provides it in a format that can be processed by
# the enrichmentmap plugin of the cytoscape program
# Usage:
# 1) take the hsa_circ_xxxxxxx-emap.tsv file and truncate it so that it contains terms that should be
#    visualized (i.e., top 100, alternatively more, but it saves time)
# 2) Open Cytoscape
# 3) Apps > EnrichmentMap
# 4) click +
# 5) fill in Analysis Type : Generic/...
#            Enrichments   : the hsa_circ_xxxxxxx-emap.tsv
#            GMT           : c5.all.v7.4.symbols.gmt (from the MSIGDB)
# 7) set FDR cutoff
# 6) study the interactions


args = commandArgs(trailingOnly=TRUE)
circ <- args[1]

library(openxlsx)
data <- read.xlsx(paste(circ, ".xlsx", sep=""))
data <- data[order(data$pvalue),]
data <- cbind(data$Row.names, data$goTerms, data$pvalue, +1)
colnames(data) <- c("GO.ID", "Description", "p.Val", "Phenotype")
write.table(data, paste(circ, "-emap.tsv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
