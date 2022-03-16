library(Rcpp)
sourceCpp("pvalue.cpp")
sourceCpp("bootstrap.cpp")
library(polynom)

pValueFull <- function(muCirc, goMu, goM, ImmuCirc, score) {
    if(score == 0) { return(1) }

    coefsMu <- getCounts(muCirc, sum(goMu))
    coefsM <- getCounts(ImmuCirc, sum(goM))
    
    muPoly <- polynomial(coef = coefsMu)
    mPoly <- polynomial(coef = coefsM)
    
    #print(muPoly)
    #print(mPoly)
    
    generatingCoefs <- coefficients(muPoly * mPoly)
    
    return(sum(generatingCoefs[(score+1):length(generatingCoefs)]) / sum(generatingCoefs))
}

pValueBootstrap <- function(muCirc, goMu, goM, ImmuCirc, score) {
    if(score == 0) { return(1) }
    
    return(pvalue_bootstrap(muCirc, sum(goMu), ImmuCirc, sum(goM), score, 1000000))
}

getCounts <- function(pathsCounts, y) {
    table <- as.data.frame(table(factor(pathsCounts)))
    weights <- as.integer(as.vector(table$Var1))
    counts <- as.integer(as.vector(table$Freq))
    #print("get counts")
    #print(weights)
    #print(counts)
    #print(y)
    return(generating_poly_general(weights, counts, y))
}

# library(Rcpp); sourceCpp('pvalue.cpp'); generating_poly_general(c(0,1), c(1750, 11), 10)
