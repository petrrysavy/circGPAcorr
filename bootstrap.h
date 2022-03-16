#ifndef GUARD_bootstrap_h
#define GUARD_bootstrap_h

#include<Rcpp.h>

using namespace std;
using namespace Rcpp;

void sample(bool* arr, int size, int n);

int dot(IntegerVector counts, bool* index);

double pvalue_bootstrap(IntegerVector weightsMu, int goMuSize, IntegerVector weightsM, int goMSize, int score, int trials);

#endif
