#include <stdio.h>

#include <mkl.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"

#ifndef INCLUDED_PRINTUTILS_H
#include "printUtils.h"
#endif

#define DIE(x) fprintf(stderr, x); exit(1)
#define MIN(a,b) a < b ? a : b
#define MAX(a,b) a > b ? a : b

int tid;
int numThreads;
int threadCounter;
#pragma omp threadprivate(threadCounter, tid)
int nActions1, nActions2;
float * payoffsA, * payoffsB;
int * allActions1, * allActions2;

int readGame(char * gameData);

int kSubsetsItHelper(int * arr, int n, int k);
void kSubsetsItOMP(int k, void (*f) (int *, int *, int)) {
    int arr1[k], arr2[k];
    for (int i = 0; i < k; ++i) arr1[i] = i, arr2[i] = i;
    do {
        do {
            if (threadCounter % numThreads == tid) {
                f(arr1, arr2, k);
                threadCounter++;
            } else {
                // printf("Thread %d skipping...\n", tid);
                threadCounter++;
            }
        } while (kSubsetsItHelper(arr2, nActions2, k));
    } while (kSubsetsItHelper(arr1, nActions1, k));
}

// Given payoff table and index, returns val at that index or errors out
// Column major ordering
float getPayoff(float * payoffs, int i, int j) {
    if (i >= nActions1 || j >= nActions2) {
        printf("%d, %d\n",i , j);
        DIE("Error: getPayoff: index out of bounds\n");
    } else {
        return payoffs[j * nActions1 + i];
    }
}

void buildMat(int * acc1, int * acc2, int suppSize, float * mat, bool isB) {
    // both payoffsA and payoffsB are col major
    //
    int k = 0;
    if (isB) {
        for (int i = 0; i < suppSize; ++i) {
            for (int j = 0; j < suppSize; ++j) 
                mat[k++] = getPayoff(payoffsB, acc1[i], acc2[j]);
            mat[k++] = 1.; 
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    } else {
        for (int j = 0; j < suppSize; ++j) {
            for (int i = 0; i < suppSize; ++i) 
                mat[k++] = getPayoff(payoffsA, acc1[i], acc2[j]);
            mat[k++] = 1.;
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    }
}

void buildFullStrat(int * ac, float * stratWeights, int sizeSubSet,
                    int sizeStrat, float * strategy) 
{
    for (int i = 0; i < sizeStrat; ++i) strategy[i] = 0.;
    for (int i = 0; i < sizeSubSet; ++i) strategy[ac[i]] = stratWeights[i];
}

void nashEq(int * acc1, int * acc2, int suppSize) {
    int matSize = suppSize+1;
    float matA[matSize * matSize];
    int nrhs = 1;
    int numEqs = matSize;
    int lda = numEqs;
    int ldb = lda;
    int ipiv[numEqs];

    float vecX[matSize], vecY[matSize];
    for (int i = 0; i < matSize; ++i) {
        vecX[i] = 0.; 
        vecY[i] = 0.;
    }
    vecX[matSize-1] = 1.; vecY[matSize-1] = 1.;

    // Get prob distribution for support in player 1's strategy
    // so that player 2 is indifferent towards 1's strategy
    buildMat(acc1, acc2, suppSize, matA, true);

    if (0!=LAPACKE_sgesv(LAPACK_COL_MAJOR, numEqs, nrhs, matA, lda, ipiv, vecX, lda)) {
        // printf("No linear system solution\n");
        return;
    }
    // Check that prob distribution is all non-negative weights
    for (int i = 0; i < suppSize; ++i) {
        if (vecX[i] < 0.) {
            // printf("found negative vecX\n");
            return;
        }
    }
    
    // Get prob distribution for support in player 2's strategy
    // so that player 1 is indifferent towards 2's strategy
    buildMat(acc1, acc2, suppSize, matA, false);
    if (0!=LAPACKE_sgesv(LAPACK_COL_MAJOR, numEqs, nrhs, matA, lda, ipiv, vecY, lda)){
        return;
    }
    // Check that prob distribution is all non-negative weights
    for (int i = 0; i < suppSize; ++i) {
        if (vecY[i] < 0.) {
            // printf("found negative vecY\n");
            return;
        }
    }
    float strat1[nActions1], strat2[nActions2];
    buildFullStrat(acc1, vecX, suppSize, nActions1, strat1); 
    buildFullStrat(acc2, vecY, suppSize, nActions2, strat2); 
    float xB[nActions2], Ay[nActions1];
    for (int i = 0; i < nActions1; ++i) Ay[i] = 0.;
    for (int i = 0; i < nActions2; ++i) xB[i] = 0.;
     
    int m = nActions1, n = 1, k = nActions2;
    lda = nActions1, ldb = nActions2;
    int ldc = nActions1;
    float alpha = 1., beta = 1.; 
   
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha,
                payoffsA, lda, strat2, ldb, beta, Ay, ldc);

    m = 1, n = nActions2, k = nActions1;
    lda = nActions1, ldb = nActions1, ldc = 1;

    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, 
                strat1, lda, payoffsB, ldb, beta, xB, ldc);

    float maxAy = Ay[0], maxxB = xB[0];
    for (int i = 0; i < nActions1; ++i) {
        maxAy = MAX(Ay[i], maxAy);
    }
    for (int i = 0; i < nActions2; ++i) {
        maxxB = MAX(xB[i], maxxB);
    }

    for (int i = 0; i < suppSize; ++i) {
        if (Ay[acc1[i]] < (maxAy - fabs(maxAy / 1000.))) {
            // printf("Failed maxAy test: %.5f < %.5f\n", Ay[acc1[i]], maxAy);
            return;
        }
    }

    for (int i = 0; i < suppSize; ++i) {
        if (xB[acc2[i]] < (maxxB - fabs(maxxB / 1000.))) {
            // printf("Failed maxxB test: %.5f < %.5f\n", xB[acc2[i]], maxxB);
            return;
        }
    }
    
    // printSolution(strat1, strat2, nActions1, nActions2); 
}

int main(int argc, char * argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Incorrect command line args\n");
        return 1;
    } else {
        readGame(argv[1]);
        numThreads = atoi(argv[2]); // Set Number of threads
        double wt1, wt2, ct1, ct2;
        timing(&wt1, &ct1);
        // Init set of all actions
        int allActs1[nActions1], allActs2[nActions2];
        for (int i = 0; i < nActions1; ++i) allActs1[i] = i;
        for (int i = 0; i < nActions2; ++i) allActs2[i] = i;
        allActions1 = allActs1, allActions2 = allActs2;

        int maxSupport = MIN(nActions1, nActions2);
        omp_set_num_threads(numThreads);
        // numThreads = 8; 
        #pragma omp parallel shared(maxSupport)
        {   
            threadCounter = 0, tid = omp_get_thread_num();
            for (int i = 1; i <= maxSupport; ++i) {
                kSubsetsItOMP(i, nashEq);
            }
        }
        
        free(payoffsA); free(payoffsB);
        timing(&wt2, &ct2);
        printf("%.5f\n", wt2 - wt1);
        return 0;
    }
}

