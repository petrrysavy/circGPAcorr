source('annotate.R')
source('correlationMatrix.R')
Rcpp::sourceCpp("bootstrap2.cpp")

library(simctest)

setClass("annotationResultCorrelated",
    slots=list(
        score="numeric",
        pvalueUB="numeric",
        pvalueLB="numeric",
        trials="numeric",
        scorenorm="numeric",
        expectedScore="numeric",
        time="numeric",
        pvalueTime="numeric",
        bootstrap="numeric",
        bootstrapTime="numeric",
        goSize="numeric",
        accepted="logical"
    )
)
print.annotationResultCorrelated <- function(x) { cat(paste("Normalized score:", x@scorenorm, "\n")); cat(paste("p-value in (", x@pvalueLB, ", ", x@pvalueUB, ")\n")); }

# This script should calculate annotation of a circRNA with GO terms together with correlations
annotateCircCorrelated <- function(circRNA, goTermList, goLower = -1, goUpper = Inf, circRNAExpressionMatrix, miRNAExpressionMatrix, geneExpressionMatrix, ...) {
    if(!(identical(colnames(circRNAExpressionMatrix), colnames(miRNAExpressionMatrix)) & identical(colnames(miRNAExpressionMatrix), colnames(geneExpressionMatrix)))) { stop("Sample names must be the same!") }
    if(! circRNA %in% rownames(circRNAExpressionMatrix)) { stop("Unknown circrna") }

    muCirc <- buildMuCirc(circRNA)
    rhoMuCirc <- buildRhoMuCirc(circRNA, circRNAExpressionMatrix, miRNAExpressionMatrix)
    # the following line takes way tooooooooo long
    rhoMMuImHadammar <- buildRhoMMuImHadammar(geneExpressionMatrix, miRNAExpressionMatrix, Im)
    return(sapply(goTermList, function(goTerm) {
        tryCatch({
            print(paste("Probing GO term:", goTerm))
            goMu <- buildGOMu(goTerm)
            goM <- buildGOM(goTerm)
            goSize <- sum(goMu) + sum(goM)
            if(goSize > goUpper | goSize < goLower) { return(NULL) }
            return(annotateVectorizedCorrelated(muCirc, goMu, goM, Im, rhoMuCirc, rhoMMuImHadammar, ...))
        }, error=function(cond) {
            print(cond)
            return(NULL)
        })
    }))
}

annotateCircCorrelated2 <- function(circRNA, goTermList, goLower = -1, goUpper = Inf, circRNAExpressionMatrix, miRNAExpressionMatrix, geneExpressionMatrix, fdrLevel = 1.0, circGPAdf = NULL, ...) {
    if(!(identical(colnames(circRNAExpressionMatrix), colnames(miRNAExpressionMatrix)) & identical(colnames(miRNAExpressionMatrix), colnames(geneExpressionMatrix)))) { stop("Sample names must be the same!") }
    if(! circRNA %in% rownames(circRNAExpressionMatrix)) { stop("Unknown circrna") }

    muCirc <- buildMuCirc(circRNA)
    rhoMuCirc <- buildRhoMuCirc(circRNA, circRNAExpressionMatrix, miRNAExpressionMatrix)
    rhoMMuImHadammar <- buildRhoMMuImHadammar(geneExpressionMatrix, miRNAExpressionMatrix, Im)
    print("graph created together with correlations")
    # <<- accesses the outer scope! - done for thresholding
    refutedNum <- 0
    return(sapply(goTermList, function(goTerm) {
        tryCatch({
            if(!is.null(circGPAdf) && (goTerm %in% rownames(circGPAdf)) && !is.na(circGPAdf[goTerm,"pvalue"]) && circGPAdf[goTerm,"pvalue"] == 1.0) {
                print(paste("GoTerm", goTerm,"skipped due to circGPA results"))
                refutedNum <<- refutedNum + 1
                return(new("annotationResultCorrelated",
                    score=0.0, pvalueUB=1.0, pvalueLB=1.0, trials=0, scorenorm=0.0, expectedScore=0.0, time=0, pvalueTime=0, bootstrap=NaN, bootstrapTime=NaN, goSize=as.numeric(circGPAdf[goTerm,"goSize"]), accepted=FALSE
                ))
            }
        
            print(paste("Probing GO term:", goTerm))
            goMu <- buildGOMu(goTerm)
            goM <- buildGOM(goTerm)
            goSize <- sum(goMu) + sum(goM)
            if(goSize > goUpper | goSize < goLower) { return(NULL) }
            pvalueThreshold <- if (fdrLevel >= 1.0) 1.0 else fdrLevel * (length(goTermList) - refutedNum) / length(goTermList)
            print(paste("pvalueThreshold", pvalueThreshold))
            res <- annotateVectorizedCorrelated(muCirc, goMu, goM, Im, rhoMuCirc, rhoMMuImHadammar, pvalueThreshold = pvalueThreshold, ...)
            if(!res@accepted){ refutedNum <<- refutedNum + 1 }
            return(res)
        }, error=function(cond) {
            print(cond)
            refutedNum <<- refutedNum + 1
            return(NULL)
        })
    }))
}

annotateCorrelated <- function(circRNA, goTerm, circRNAExpressionMatrix, miRNAExpressionMatrix, geneExpressionMatrix) {
    return(annotateVectorizedCorrelated(buildMuCirc(circRNA), buildGOMu(goTerm), buildGOM(goTerm), Im,
        buildRhoMuCirc(circRNA, circRNAExpressionMatrix, miRNAExpressionMatrix), buildRhoMMuImHadammar(geneExpressionMatrix, miRNAExpressionMatrix, Im)))
}

pValueBootstrapCorrelated <- function(weightsMu, goMu, goM, weightsM, score) {
    if(score == 0) { return(1) }
    
    return(pvalue_bootstrap(weightsMu, sum(goMu), weightsM, sum(goM), score, 1000000))
}

pValueSimctestCorrelated <- function(weightsMu, goMu, goM, weightsM, score) {
    if(score == 0) { return(1) }
    
    sampler <- new(BootstrapSampler, weightsMu, sum(goMu), weightsM, sum(goM), score)
    simctst <- simctest(function() sampler$bootstrapSample())
    return(simctst@p.value)
}

pValueBounded <- function(weightsMu, goMu, goM, weightsM, score, pvalueThreshold=1.0, trialsMin=100, trialsMax=100000, ...) {
#     print(score)
    if(score == 0) { return(list(pvalueLB = 1.0, pvalueUB = 1.0, trials = 0, accepted = FALSE)) }

    trials <- if(pvalueThreshold < 1.0) trialsMin else trialsMax
    
    while(TRUE) {
        discretization <- weightDiscretizationLB(weightsM, sum(goM), weightsMu, sum(goMu), score, limit=trials)
        # feed the discretized values to circGPA
        # this test is here to prevent underflow - we cannot calculate the pValue if there is no way to reach score
        pvalueLB <- if(!discretization$valid) 0.0
                    else pValueFull(discretization$weightsMu, goMu, goM, discretization$weightsM, discretization$score)
        
        if(is.na(pvalueLB)) {
            print(discretization$score)
            print(sum(goMu))
            print(sum(goM))
            #print(discretization$weightsM)
            #print(discretization$weightsMu)
            print(sum(discretization$weightsM))
            print(sum(discretization$weightsMu))
            print(sum(is.na(discretization$weightsM)))
            print(sum(is.na(discretization$weightsMu)))
        }
        
        print(paste("trials:", trials))
        print(paste("pvalueLB:", pvalueLB))
        
        if(pvalueLB > pvalueThreshold) return(list(pvalueLB = pvalueLB, pvalueUB = 1.0, trials = trials, accepted = FALSE))
        
        if(trials == trialsMax) {
            discretization <- weightDiscretizationUB(weightsM, sum(goM), weightsMu, sum(goMu), score, limit=trials)
            pvalueUB <- pValueFull(discretization$weightsMu, goMu, goM, discretization$weightsM, discretization$score)
            return(list(pvalueLB = pvalueLB, pvalueUB = pvalueUB, trials = trials, accepted = TRUE))
        }
        trials <- if(trials * 10 > trialsMax) trialsMax else trials * 10
    }
}

# goMu - 0/1 vector of is mu in goTerm
# goM the same is mRNA in goTerm
# muCirc - 0/1/more vecotor of circrna neighborhood
# Im - matrix with interactions between muRNA nad miRNA
# test:
annotateVectorizedCorrelated <- function(muCirc, goMu, goM, Im, rhoMuCirc, rhoMMuImHadammar, ...) {
    time <- as.numeric(system.time({
        rhoMuMuCircHadammar <- muCirc * rhoMuCirc
        weightsM <- rhoMMuImHadammar %*% rhoMuMuCircHadammar
        score <- dot(rhoMuMuCircHadammar, goMu) + dot(weightsM, goM) # ignore the halfs, I love whole numbers
#         print(score)
        pvalueTime <- as.numeric(system.time({
            pvalueVec <- pValueBounded(rhoMuMuCircHadammar, goMu, goM, weightsM, score, ...) #pValueFull(muCirc, goMu, goM, ImmuCirc, score)
        })['elapsed'])
        bootstrapTime <- as.numeric(system.time({
            bootstrap <- pValueBootstrapReal(weightsMu, goMu, goM, weightsM, score)
            #bootstrap <- NaN #pValueSimctestCorrelated(rhoMuMuCircHadammar, goMu, goM, weightsM, score) # NaN #
        })['elapsed'])
        expectedScore <- dot(rhoMuMuCircHadammar, uniformize(goMu)) + dot(weightsM, uniformize(goM))
        scorenorm <- score / expectedScore
        goSize <- sum(goM) + sum(goMu)
    })['elapsed'])
    return(new("annotationResultCorrelated",
        score=score, pvalueUB=pvalueVec$pvalueUB, pvalueLB=pvalueVec$pvalueLB, trials=pvalueVec$trials, scorenorm=scorenorm, expectedScore=expectedScore, time=time, pvalueTime=pvalueTime, bootstrap=bootstrap, bootstrapTime=bootstrapTime, goSize=goSize, accepted=pvalueVec$accepted
    ))
}
