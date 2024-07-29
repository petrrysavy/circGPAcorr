source('buildGraph.R')

# rows of the expression matrix are named by the circRNAs
# columns of the expression matrix are samples, this order must be fixed!
# buildRhoMuCirc <- function(circRNA, circExpressionMatrix, miRNAExpressionMatrix) {
#     circExpression <- as.integer(circExpressionMatrix[which(rownames(circExpressionMatrix) == circrna),])
# 
#     rhoCirc <- sapply(all_mirna_list, function(mirna) {
#         index <- which(rownames(miRNAExpressionMatrix) == mirna)
#         if(length(index) == 0) { return(0.0) }
#         row <- as.integer(miRNAExpressionMatrix[index[1],])
#         return(abs(cor(circExpression, row, method = c("spearman"))))
#     })
#     
#     rhoCirc[is.na(rhoCirc)] <- 0.0
#     return(as.numeric(rhoCirc))
# }
# 
# buildRhoMMuImHadammar <- function(geneExpressionMatrix, miRNAExpressionMatrix, Im) {
#     rhoMuImHadammar <- Im
# 
#     for(geneIndex in 1:nrow(rhoMuImHadammar)) {
#         index <- which(rownames(geneExpressionMatrix) == all_gene_list[geneIndex])
#         if(length(index) == 0) { rhoMuImHadammar[geneIndex,] = 0 }
#         else {
#             geneExpression <- as.integer(geneExpressionMatrix[index[1],])
#             for(miRNAIndex in 1:ncol(rhoMuImHadammar)) {
#                 if(rhoMuImHadammar[geneIndex, miRNAIndex] != 0) {
#                     index <- which(rownames(miRNAExpressionMatrix) == all_mirna_list[miRNAIndex])
#                     rhoMuImHadammar[geneIndex, miRNAIndex] = rhoMuImHadammar[geneIndex, miRNAIndex] * (if(length(index) == 0) { 0.0 }
#                         else { abs(cor(geneExpression, as.integer(miRNAExpressionMatrix[index[1],]), method = c("spearman"))) })
#                 }
#             }
#         }
#     }
#     
#     rhoMuImHadammar[is.na(rhoMuImHadammar)] <- 0.0
#     return(rhoMuImHadammar)
# }

buildRhoMuCirc <- function(circRNA, circExpressionMatrix, miRNAExpressionMatrix) {
    circExpression <- as.integer(circExpressionMatrix[which(rownames(circExpressionMatrix) == circRNA),])

    rhoCirc <- sapply(all_mirna_list, function(mirna) {
        index <- which(rownames(miRNAExpressionMatrix) == mirna)
        if(length(index) == 0) { return(0.0) }
        row <- as.integer(miRNAExpressionMatrix[index[1],])
        return(max(0.0, -(cor(circExpression, row, method = c("spearman")))))
    })
    
    rhoCirc[is.na(rhoCirc)] <- 0.0
    return(as.numeric(rhoCirc))
}

buildRhoMMuImHadammar <- function(geneExpressionMatrix, miRNAExpressionMatrix, Im) {
    rhoMuImHadammar <- Im

    for(geneIndex in 1:nrow(rhoMuImHadammar)) {
        index <- which(rownames(geneExpressionMatrix) == all_gene_list[geneIndex])
        if(length(index) == 0) { rhoMuImHadammar[geneIndex,] = 0 }
        else {
            geneExpression <- as.integer(geneExpressionMatrix[index[1],])
            for(miRNAIndex in 1:ncol(rhoMuImHadammar)) {
                if(rhoMuImHadammar[geneIndex, miRNAIndex] != 0) {
                    index <- which(rownames(miRNAExpressionMatrix) == all_mirna_list[miRNAIndex])
                    rhoMuImHadammar[geneIndex, miRNAIndex] = rhoMuImHadammar[geneIndex, miRNAIndex] * (if(length(index) == 0) { 0.0 }
                        else { max(0.0, -(cor(geneExpression, as.integer(miRNAExpressionMatrix[index[1],]), method = c("spearman")))) })
                }
            }
        }
    }
    
    rhoMuImHadammar[is.na(rhoMuImHadammar)] <- 0.0
    return(rhoMuImHadammar)
}
