#include "../bootstrap.h"

#include <Rcpp.h>

// this is a partial copy-paste of the bootstrap method to include logging of p-values at selected times
// should be used only for experimental plot in the paper
// [[Rcpp::export]]
NumericVector pvalue_bootstrap_logging(IntegerVector weightsMu, int goMuSize, IntegerVector weightsM, int goMSize, int score, IntegerVector testPoints) {
    const int muSize = weightsMu.size();
    const int mSize = weightsM.size();
    bool* goMu = new bool[muSize];
    bool* goM = new bool[mSize];
    
    const int trials = testPoints[testPoints.size() - 1];
    int trialsI = 0;
    NumericVector retval(testPoints.size());
    
    int positive = 0;
    for(int i = 1; i <= trials; ++i) {
        memset(goMu, 0, muSize * sizeof(bool));
        memset(goM, 0, mSize * sizeof(bool));
        
        sample(goMu, muSize, goMuSize);
        sample(goM, mSize, goMSize);
        
        const int sample = dot(weightsMu, goMu) + dot(weightsM, goM);
        if(sample >= score) ++positive;
        
        if(testPoints[trialsI] == i) {
            retval[trialsI] = ((double) positive) / i;
            ++trialsI;
        }
    }
    
     delete[] goMu;
     delete[] goM;
    
    return retval;
}
