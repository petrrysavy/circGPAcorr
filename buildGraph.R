 # BiocManager::install("miRBaseConverter")
# BiocManager::install("multiMiR")
 
library(miRBaseConverter)

#library(msigdbr)
library(GO.db)
library(org.Hs.eg.db)
library("biomaRt")
library(httr)
httr::set_config(config(ssl_verifypeer = FALSE))

library(multiMiR)

#library(Rfast) # for binary search - no way, only double

#source("spidermirdownload.R")
source("circinteractome.R")
source("mirnaGO.R")

#all_gene_sets = msigdbr(species = "Homo sapiens")
#allegs = get("GO:0070104", org.Hs.egGO2ALLEGS)
#genes = unlist(mget(allegs,org.Hs.egSYMBOL))

# see https://support.bioconductor.org/p/124462/
listGenes <- function() {
    tryCatch({
        mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        all_coding_genes <- getBM(attributes = c( "hgnc_symbol"), filters = c("biotype"), values = list(biotype="protein_coding"), mart = mart)
        return(sort(all_coding_genes$hgnc_symbol))
    }, error=function(cond) {
        print(cond)
        return(readRDS("all_coding_genes.RDS"))
    })
}

listMiRNA <- function() {
    miRNATab=getMiRNATable(version="v22",species="hsa")
    return(sort(unique(miRNATab$Mature1)))
}

multiMirGeneInteract <- function() {
    #map <- read.csv("mirna-mrna-interactions.csv")
    map <- read.csv("mirna-mrna-interactions.csv")
    map <- subset(map, select = -c(X) )
    return(map)
}

if (file.exists("graph.RData")) {
    load("graph.RData")
} else {
    all_gene_list <- listGenes()
    all_mirna_list <- listMiRNA()
    # save(all_gene_list, all_mirna_list, file = "graph.RData")
}

geneInteract <- multiMirGeneInteract()
circ_mirna_map <- loadMap()
mirna_go_map <- loadMapMiRNAGoTerms()

buildIm <- function() {
    Im <- matrix(0, nrow=length(all_gene_list), ncol=length(all_mirna_list))
    for(i in 1:nrow(geneInteract)) {
        miRNA <- geneInteract[i,1]
        gene <- geneInteract[i,2]
        
        miRNAI <- which(all_mirna_list == miRNA)
        geneI <- which(all_gene_list == gene)
        
        if(length(miRNAI) == 0 | length(geneI) == 0) { next }
        
        Im[geneI[1], miRNAI[1]] <- 1
    }
    return(Im)
}

buildMuCirc <- function(circRNA) {
    interactions <- circ_mirna_map[circ_mirna_map$circRNA == circRNA,]
    return(as.numeric(all_mirna_list %in% interactions$miRNAs))
}

buildGOMu <- function(GOterm) {
    goid <- names(Term(GOTERM)[which(longnames == GOterm)])
    mirnas <- mirna_go_map$Mature1[mirna_go_map$goTerm == GOterm | mirna_go_map$goTerm == goid]
    return(as.numeric(all_mirna_list %in% mirnas))
}

buildGOM <- function(GOterm) {
    if(GOterm %in% names(m_list)) {
        genes <- unlist(m_list[names(m_list) == GOterm])
    } else {
        allegs = get(GOterm, org.Hs.egGO2ALLEGS)
        genes = unlist(mget(allegs,org.Hs.egSYMBOL))
    }
    return(as.numeric(all_gene_list %in% genes))
}
