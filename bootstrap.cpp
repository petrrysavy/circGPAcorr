#include "bootstrap.h"

#include <Rcpp.h>
#include <cstdlib>
#include <cstring>

void sample(bool* arr, int size, int n) {
    int pos;
    while(n != 0) {
        do { pos = (rand() % size); } while(arr[pos]);
        arr[pos] = true;
        --n;
    }
}

int dot(IntegerVector counts, bool* index) {
    int product = 0;
    for(IntegerVector::iterator i = counts.begin(); i != counts.end(); ++i, ++index)
        product += (*i) * ((int) (*index));
    return product;
}

// test: pvalue_bootstrap(c(1,1,1), 2, c(2,3,1,1,2), 3, 8, 1000)
// [[Rcpp::export]]
double pvalue_bootstrap(IntegerVector weightsMu, int goMuSize, IntegerVector weightsM, int goMSize, int score, int trials) {
    const int muSize = weightsMu.size();
    const int mSize = weightsM.size();
    bool* goMu = new bool[muSize];
    bool* goM = new bool[mSize];
    
    int positive = 0;
    for(int i = 0; i < trials; ++i) {
        memset(goMu, 0, muSize * sizeof(bool));
        memset(goM, 0, mSize * sizeof(bool));
        
        sample(goMu, muSize, goMuSize);
        sample(goM, mSize, goMSize);
        
        const int sample = dot(weightsMu, goMu) + dot(weightsM, goM);
        if(sample >= score) ++positive;
    }
    
     delete[] goMu;
     delete[] goM;
    
    return ((double) positive + 1) / trials + 1;
}
