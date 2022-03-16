# This script should calculate annotation of a circRNA with GO terms

library(polynom)
library(geometry) # only for the dot function
library(tictoc)
source('buildGraph.R')
source('pValue.R')

setClass("annotationResult",
    slots=list(
        score="numeric",
        pvalue="numeric",
        scorenorm="numeric",
        expectedScore="numeric",
        time="numeric",
        pvalueTime="numeric",
        bootstrap="numeric",
        bootstrapTime="numeric",
        goSize="numeric"
    )
)

Im <- buildIm()

annotateCirc <- function(circRNA, goTermList) {
    muCirc <- buildMuCirc(circRNA)
    return(sapply(goTermList, function(goTerm) {
        tryCatch({
            print(paste("Probing GO term:", goTerm))
            goMu <- buildGOMu(goTerm)
            goM <- buildGOM(goTerm)
            return(annotateVectorized(muCirc, goMu, goM, Im))
        }, error=function(cond) {
            print(cond)
            return(NULL)
        })
    }))
}

annotate <- function(circRNA, goTerm) {
    return(annotateVectorized(buildMuCirc(circRNA), buildGOMu(goTerm), buildGOM(goTerm), Im))
}

# dataVector = vector with counts how many times each muRNA/mRNA is connected to the circRNA
buildPolynomial <- function(dataVector) {
    tbl <- as.data.frame(table(factor(dataVector, levels = 0:max(dataVector))))
    return(polynomial(coef = tbl[,2]))
}

# uniform probability vector for goMu or goM - vec(1)*sum(arr) / len(arr)
uniformize <- function(goVec) {
    return(rep(1, length(goVec)) * sum(goVec) / length(goVec))
}

pValueSimple <- function(muCirc, goMu, goM, ImmuCirc, score) {
    # get the p-value, convert the muCirc vector to the polynomial
    muCircPoly <- buildPolynomial(muCirc)
    ImMuCircPoly <- buildPolynomial(ImmuCirc)
    # now calculate power of those to get the number of ways
    generatingPoly <- muCircPoly^sum(goMu) * ImMuCircPoly^sum(goM)
    generatingCoefs <- coefficients(generatingPoly)
    # the p-value is how many of the combinations gave higher score than we have
    return(sum(generatingCoefs[(score+1):length(generatingCoefs)]) / sum(generatingCoefs)) # score + 1 as coefficients are x^0 + x^1 indexed from 1
}

# goMu - 0/1 vector of is mu in goTerm
# goM the same is mRNA in goTerm
# muCirc - 0/1/more vecotor of circrna neighborhood
# Im - matrix with interactions between muRNA nad miRNA
# test:
# gomu <- c(0,1,1); gom <- c(1,1,1,0,0);  mucirc <- c(1,1,1); im <- matrix(c(0,1,1,1,1,1,0,0,1,1,0,0,1,1,0), nrow=5, ncol=3, byrow=TRUE); annotateVectorized(mucirc, gomu, gom, im)
annotateVectorized <- function(muCirc, goMu, goM, Im) {
    time <- as.numeric(system.time({
        ImmuCirc <- Im %*% muCirc
        score <- dot(muCirc, goMu) + dot(ImmuCirc, goM) # ignore the halfs, I love whole numbers
        pvalueTime <- as.numeric(system.time({
            pvalue <- pValueFull(muCirc, goMu, goM, ImmuCirc, score)
        })['elapsed'])
        bootstrapTime <- as.numeric(system.time({
            bootstrap <- NaN #pValueBootstrap(muCirc, goMu, goM, ImmuCirc, score) # NaN #
        })['elapsed'])
        expectedScore <- dot(muCirc, uniformize(goMu)) + dot(ImmuCirc, uniformize(goM))
        scorenorm <- score / expectedScore
        goSize <- sum(goM) + sum(goMu)
    })['elapsed'])
    # return
    return(new("annotationResult",
        score=score, pvalue=pvalue, scorenorm=scorenorm, expectedScore=expectedScore, time=time, pvalueTime=pvalueTime, bootstrap=bootstrap, bootstrapTime=bootstrapTime, goSize=goSize
    ))
}
