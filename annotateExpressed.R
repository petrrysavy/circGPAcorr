source('annotateCorrelated.R')
source('correlationMatrix.R')
Rcpp::sourceCpp("bootstrap2.cpp")

cleanExpressionVector <- function(expression, signToKeep) {
    expression[sign(expression) == -signToKeep] = 0.0
    expression[is.na(expression)] <- 0.0
    return(abs(expression))
}

# This script should calculate annotation of a circRNA with GO terms together with correlations
annotateCircExpressed <- function(circRNA, goTermList, circLogFold, miRNAResultTable, geneResultTable, goLower = -1, goUpper = Inf, ...) {
    mirnaExpressionVector <- miRNAResultTable[all_mirna_list, "log2FoldChange"]
    mirnaExpressionVector <- cleanExpressionVector(mirnaExpressionVector, -sign(circLogFold))
    print(mirnaExpressionVector)
    
    mrnaExpressionVector <- geneResultTable[all_gene_list, "log2FoldChange"]
    mrnaExpressionVector <- cleanExpressionVector(mrnaExpressionVector, sign(circLogFold))
    print(mrnaExpressionVector)

    muCirc <- buildMuCirc(circRNA)
    return(sapply(goTermList, function(goTerm) {
        tryCatch({
            print(paste("Probing GO term:", goTerm))
            goMu <- buildGOMu(goTerm)
            goM <- buildGOM(goTerm)
            goSize <- sum(goMu) + sum(goM)
            if(goSize > goUpper | goSize < goLower) { return(NULL) }
            return(annotateVectorizedExpressed(muCirc, goMu, goM, Im, mrnaExpressionVector, mirnaExpressionVector, ...))
        }, error=function(cond) {
            print(cond)
            return(NULL)
        })
    }))
}

# goMu - 0/1 vector of is mu in goTerm
# goM the same is mRNA in goTerm
# muCirc - 0/1/more vecotor of circrna neighborhood
# Im - matrix with interactions between muRNA nad miRNA
# test:
# gomu <- c(0,1,1); gom <- c(1,1,1,0,0);  mucirc <- c(1,1,1); im <- matrix(c(0,1,1,1,1,1,0,0,1,1,0,0,1,1,0), nrow=5, ncol=3, byrow=TRUE); annotateVectorized(mucirc, gomu, gom, im)
annotateVectorizedExpressed <- function(muCirc, goMu, goM, Im, mrnaExpressionVector, mirnaExpressionVector, ...) {
    time <- as.numeric(system.time({
        #ImmuCirc <- Im %*% muCirc
        weightsMu <- muCirc * mirnaExpressionVector
        wieghtedIm <- diag(mrnaExpressionVector) %*% Im
        #weightsM <- ImmuCirc * mrnaExpressionVecto
        # TEST:
        # wieghtedIm  <- matrix(c(1,2,0,0,2,1,3,0,0,4), nrow=2, byrow=TRUE)
        # weightsMu <- c(0,1,0,2,2)
        # [1]  7 10
        #TODO - tady je možná chyba!!!!!
        weightsM <- sapply(1:nrow(wieghtedIm), function(m) {sum(wieghtedIm[m,weightsMu > 0])+sum(weightsMu[wieghtedIm[m,] > 0])})
        score <- dot(weightsMu, goMu) + dot(weightsM, goM)
        print(score)
        pvalueTime <- as.numeric(system.time({
            pvalueVec <- pValueBounded(weightsMu, goMu, goM, weightsM, score, ...) #pValueFull(muCirc, goMu, goM, ImmuCirc, score)
        })['elapsed'])
        bootstrapTime <- as.numeric(system.time({
            bootstrap <- pValueBootstrapReal(weightsMu, goMu, goM, weightsM, score)
            #bootstrap <- NaN #pValueSimctestCorrelated(weightsMu, goMu, goM, weightsM, score) # NaN #
        })['elapsed'])
        expectedScore <- dot(weightsMu, uniformize(goMu)) + dot(weightsM, uniformize(goM))
        scorenorm <- score / expectedScore
        goSize <- sum(goM) + sum(goMu)
    })['elapsed'])
    return(new("annotationResultCorrelated",
        score=score, pvalueUB=pvalueVec$pvalueUB, pvalueLB=pvalueVec$pvalueLB, trials=pvalueVec$trials, scorenorm=scorenorm, expectedScore=expectedScore, time=time, pvalueTime=pvalueTime, bootstrap=bootstrap, bootstrapTime=bootstrapTime, goSize=goSize, accepted=pvalueVec$accepted
    ))
}
