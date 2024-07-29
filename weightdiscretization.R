source("weightclustering.R")

safeNLargestSum <- function(weights, n) {
    if(n == 0) return(0)
    weights <- sort(weights, decreasing=TRUE)
    return(sum(weights[1:min(length(weights), n)]))
}

weightDiscretizationUB <- function(weightsM, goMSize, weightsMu, goMuSize, score, limit=20000) {
#     mVals <- sortAndRemoveZeros(weightsM)
#     muVals <- sortAndRemoveZeros(weightsMu)
#     print(mVals)
#     print(muVals)
#     
#     maxim <- sum(mVals$weights[1:min(length(mVals$weights), goMSize)]) + sum(muVals$weights[1:min(length(muVals$weights), goMuSize)])
#     step <- maxim/limit * 0.8 # 0.8 is reserve for rounding errors - TODO:maybe I should use some estimate/expected value for the margin ... (vals round up)
#     
#     return(list(weightsM = ceiling(weightsM / step), 
#                 weightsMu = ceiling(weightsMu / step),
#                 score = floor(score / step)))

    maxim <- (safeNLargestSum(weightsM, goMSize) + safeNLargestSum(weightsMu, goMuSize))
    step <- maxim/limit * 0.8 # 0.8 is reserve for rounding errors - TODO:maybe I should use some estimate/expected value for the margin ... (vals round up)
    
    weightsM = ceiling(weightsM / step)
    weightsMu = ceiling(weightsMu / step)
    score = floor(score / step)
    valid = (safeNLargestSum(weightsM, goMSize) + safeNLargestSum(weightsMu, goMuSize)) >= score
    
    return(list(weightsM = weightsM, weightsMu = weightsMu, score = score, valid = valid))
}

weightDiscretizationLB <- function(weightsM, goMSize, weightsMu, goMuSize, score, limit=20000) {
    #mVals <- sortAndRemoveZeros(weightsM)
    #muVals <- sortAndRemoveZeros(weightsMu)
    #print(mVals)
    #print(muVals)
    
    maxim <- (safeNLargestSum(weightsM, goMSize) + safeNLargestSum(weightsMu, goMuSize))
    step <- maxim/limit * 0.8 # 0.8 is reserve for rounding errors - TODO:maybe I should use some estimate/expected value for the margin ... (vals round up)
    
    weightsM = floor(weightsM / step)
    weightsMu = floor(weightsMu / step)
    score = ceiling(score / step)
    valid = (safeNLargestSum(weightsM, goMSize) + safeNLargestSum(weightsMu, goMuSize)) >= score
    
    return(list(weightsM = weightsM, weightsMu = weightsMu, score = score, valid = valid))
}
