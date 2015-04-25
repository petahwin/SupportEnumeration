#include <mkl.h>
#include <stdio.h>

int main() {
    float x[3] = {1,1,1};
    float y[3] = {2,2,2};
    //float dotProd = cblas_sdot(3, x, 1, y, 1);
    //printf("%f\n", dotProd);
    lapack_int nrhs = 1; // Cols of B, always 1 in our case
    lapack_int numEqs = 3; // Size of side of matrix, or number of linear equations
    lapack_int lda = numEqs;
    lapack_int ldb = lda;
    float matA[] = {3.,7.,1., 4.,2.,1., 1.,1.,0.};
    float vecB[] = {0.,0.,1.};
    lapack_int ipiv[3] = {0,0,0};

    // Apparently the function can only work using LAPACK_COL_MAJOR order; 
    // investigate
    lapack_int info = LAPACKE_sgesv(LAPACK_COL_MAJOR, numEqs, nrhs, matA, 
                                    lda, ipiv, vecB, lda);

    // Can use status code to identify failures in computation; only 
    // status code=0 is good
    printf("status code:%d\n",info);
    for (int i=0; i < ldb; ++i) {
        printf("%f\n", vecB[i]);
    }
    // Test situation where solution is not possible

    // Now test inner product
    // 14 args
    // Matrix to be multiplied should be in col maj order
    float sqA[] = {3.,2.,4., 2.,2.,4.};
    float arrB[] = {1.,1.};
    float matC[] = {0.,0.};
    int rowsA = 3, colsA = 2;
    int m = rowsA, n = 1, k = colsA; 
    int aLd = rowsA, bLd = colsA, cLd = rowsA;
    float alpha = 1., beta = 1.;

    // First try A * b
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, sqA, aLd, 
                arrB, bLd, beta, matC, cLd);
    
    printf("A*b\n"); 
    for (int i=0; i < rowsA; ++i) {
        printf("%f\n", matC[i]);
    }
    // Then try b_T * A
    float arrB_T[] = {1.,1.,1.};
    m = 1, n = colsA, k = rowsA;
    aLd = rowsA;
    bLd = rowsA; // For Op(b), not b itself
    cLd = 1;
    matC[0] = 0.; matC[1] = 0.; matC[2] = 0.;

    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, arrB_T, aLd, 
                sqA, bLd, beta, matC, cLd);
    
    printf("b_T*A\n"); 
    for (int i=0; i < colsA; ++i) {
        printf("%f\n", matC[i]);
    }
    
    return 0;
}

