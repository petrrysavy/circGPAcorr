#include <Rcpp.h>
#include <bits/stdc++.h>
#include <cstring>
#include <float.h>
#include <limits.h>
using namespace std;
using namespace Rcpp;

//#include <boost/multiprecision/mpfr.hpp>
//using namespace boost::multiprecision;

#define MIN(a,b) (((a)<(b)) ? a : b)
#define DOUBLE_TYPE long double

// weights sorted! smalledst to biggest
int maxX(IntegerVector weights, IntegerVector counts, int maxY) {
    int maxX = 0;
    for(int i = weights.size() - 1; i >= 0 && maxY > 0; --i) {
        maxX += weights[i] * MIN(counts[i], maxY);
        maxY -= counts[i];
    }
    return maxX;
}

int max(IntegerVector vector) {
    int max = INT_MIN;
    for(IntegerVector::iterator i = vector.begin(); i != vector.end(); ++i)
        if(*i > max) max = *i;
    return max;
}

DOUBLE_TYPE max(DOUBLE_TYPE* arr, const int size) {
    const DOUBLE_TYPE* stopPtr = arr + size;
    DOUBLE_TYPE max = LDBL_MIN;
    for(; arr < stopPtr; ++arr)
        if((*arr) > max) max = (*arr);
    return max;
}

DOUBLE_TYPE* binomials(const int maxN, const int maxK) {
    const int ncol = maxK + 1;
    DOUBLE_TYPE* binoms = new DOUBLE_TYPE[(maxN + 1) * ncol]();
    for(int n = 0; n <= maxN; n++) {
        binoms[n * ncol] = 1;
        const int Kstop = MIN(n, maxK);
        for(int k = 1; k <= Kstop; k++) {
            binoms[n * ncol + k] = binoms[(n - 1) * ncol + k - 1] + binoms[(n-1) * ncol + k];
        }
    }
    //Rprintf("maximum : %f\n", max(binoms, (maxN + 1) * ncol));
    //Rprintf("maximum log : %f\n", (double) log10(max(binoms, (maxN + 1) * ncol)));
    //Rprintf("maximum log : %d %d\n", max(binoms, (maxN + 1) * ncol) > LDBL_MAX, max(binoms, (maxN + 1) * ncol) > DBL_MAX);

    return binoms;
}

DOUBLE_TYPE* binomialsForCounts(IntegerVector counts, int maxy) {
    return binomials(max(counts), maxy);
}

bool overflowPrevention(DOUBLE_TYPE* current, DOUBLE_TYPE* next, int size, DOUBLE_TYPE maxCount) {
    const DOUBLE_TYPE maxValue = max(next, size);
    // Rprintf("%f %f\n", maxCount, maxValue);
    if(maxValue * maxCount > DBL_MAX) { // worst case, ovely pesimistic
        const DOUBLE_TYPE* currentStopPtr = current + size;
        for(; current != currentStopPtr; ++current, ++next) {
            *current /= maxCount;
            *next /= maxCount;
        }
        return true;
    }
    
    return false;
    
    // scan the array to find whether multiplication of current will lead to overflow
    /*double* currentStopPtr = current + size;
    for(double* currentPtr = current; currentPtr != currentStopPtr; ++currentPtr)
        if((*currentPtr) * maxCount > DOUBLE_MAX) {
            overflowRisk = true;
            break;
        }*/
}

void multiply(DOUBLE_TYPE* current, DOUBLE_TYPE* next, int nrow, int ncol, int x, int y, DOUBLE_TYPE count) {
    // i goes from 0 to nrow - y (exclusive) in current and from y to nrow (exclusive) in next
    for(int i = y; i < nrow; i++) {
        // multiply row i-y in current and store to next i
        // j goes from 0 to ncol - x (exclusive) in current and from x to ncol (exclusive) in next
        DOUBLE_TYPE* currentPtr = current + ((i - y) * ncol);
        DOUBLE_TYPE* currentStopPtr = currentPtr + ncol - x;
        DOUBLE_TYPE* nextPtr = next + (i * ncol) + x;
        for(; currentPtr < currentStopPtr; ++currentPtr, ++nextPtr)
            *nextPtr += (*currentPtr) * count;
    }
}

/* weights and counts must be sorted from the smallest to the biggest! */

// [[Rcpp::export]]
NumericVector generating_poly_general(IntegerVector weights, IntegerVector counts, int maxy) {
    
    //Rprintf("long double size %d %d %f\n", sizeof(long double), LDBL_MAX_EXP, LDBL_MAX);
//    mpfr_float a = 2;
//    mpfr_float::default_precision(1000);
//    Rprintf("%s\n", mpfr_float::default_precision());
    
    //Rprintf("long double %f\n", std::numeric_limits<long double>::min());
    const int maxx = maxX(weights, counts, maxy);
    const int nrow = maxy + 1;
    const int ncol = maxx + 1;
    
    DOUBLE_TYPE* current = new DOUBLE_TYPE[nrow * ncol](); // initialized to zero
    DOUBLE_TYPE* next = new DOUBLE_TYPE[nrow * ncol](); // initialized to zero
    const DOUBLE_TYPE* binomials = binomialsForCounts(counts, maxy);
    
    current[0] = 1; // one possibility to use zero balls with zero weight
    next[0] = 1;
    
    int size = weights.size();
    for(int i = 0; i < size; ++i) { // for each weight and count ...
        int weight = weights[i];
        int count = counts[i];
        // ... multiply by respective powers of (1 + x^weight)^count
        int maxPower = MIN(count, maxy); // no point in too big powers - ytrim
        
        // this is overflow prevention, we calcualte a pesimistic estimate of the multiplicative factor
        // and if it is too big, we adjust the current and next arrays appropriately
        DOUBLE_TYPE maxCount = 1.0; // one because of the addition
        for(int pow = 1; pow <= maxPower; pow++)
            maxCount += binomials[count * nrow + pow];
        overflowPrevention(current, next, nrow * ncol, maxCount);
        
        // and finally do the multiplication as needed
        for(int pow = 1; pow <= maxPower; pow++)
            multiply(current, next, nrow, ncol, weight * pow, pow, binomials[count * nrow + pow]);
        
        // we are done with one weight, move values from next to current and start multiplying with a new polynomial
        // this will do the *1 multiplication for us easily
        std::memcpy(current, next, nrow * ncol * sizeof(DOUBLE_TYPE)); // dest, src, size
    }
    
    // and return the result
    DOUBLE_TYPE maxValue = max(next + maxy * ncol, ncol);
    if(maxValue < DBL_MAX) maxValue = 1.0;
    
    NumericVector retval(ncol);
    for(int i = 0; i < ncol; ++i)
        retval[i] = next[maxy * ncol + i] / maxValue;
    
    delete[] binomials;
    delete[] current;
    delete[] next;

    //for(int i = 0; i < ncol; i++)
    //    Rprintf("%f ", retval[i]);
    //Rprintf("\n");
    
    return retval;
}
