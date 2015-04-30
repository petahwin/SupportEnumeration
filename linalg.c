/* big just stores max of whatever
sum is just an accumulator
temp is just used for swapping
i, j and k are just counter values
I assume that imax has to do with the index at which the max is found, the argmax if you will
I don't know what dum is for
vv stores scaling of each row
d is a parity value; useless
*/

// In:  Column major square matrix array A
//      Number of equations in A
//      ipiv[numEqs] to store permutation of A
// Out: Filled ipiv
//      LU decomposition of A
// Return: 0 on success, 1 on failure

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "printUtils.h"

#define MAX(a,b) a > b ? a : b

int croutLU(float * A, int n, int * ipiv) {
    int argmaxi;
    float max, acc, temp;
    float vv[n];

    for (int i = 0; i < n; ++i) {     /* Get implicit pivot scaling info */ 
        max = 0.0;                    /* for each row                    */
        for (int j = 0; j < n; ++j) {
            max = MAX(max, fabs(A[j * n + i]));
        }
        if (max == 0.0) return 1;     // Singular matrix
        vv[i] = 1.0/max;              // Storing scaling info
    }
    
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) { // Computing Upper Triangular elmts
            acc = A[j * n + i];
            for (int k = 0; k < i; ++k) acc -= A[k * n + i] * A[j* n + k];
            A[j * n + i] = acc;
        }

        max = 0.0;

        for (int i = j; i < n; ++i) { // Computing Lower Triangular elmts
            acc = A[j * n + i];
            for (int k = 0; k < j; ++k) acc -= A[k * n + i] * A[j * n + k];
            A[j * n + i] = acc;
            
            if (vv[i]*fabs(acc) >= max) { /* Determining largest lower triangular */
                max = vv[i]*fabs(acc);    /* elmt of given column as pivot */
                argmaxi = i;
            }
        }

        // if (j != argmaxi) {         // Consider removing to avoid branching
            for (int k = 0; k < n; ++k) {
                temp = A[k * n + argmaxi];
                A[k * n + argmaxi] = A[k * n + j];
                A[k * n + j] = temp;
            }
            vv[argmaxi] = vv[j];
        //}
        ipiv[j] = argmaxi;

        if (A[j * n + j] == 0.0) return 1; // Singular matrix, should break NOW
        
        // if (j != n) {   // Consider removing to avoid branching
            acc = 1.0/(A[j * n + j]);
            for (int i = j + 1; i < n; ++i) A[j * n + i] *= acc;
        // }
    }
    return 0;
}

void backsub(float * A, int n, int * ipiv, float * b) {
    int indexPiv, ii = 0;
    float acc;
    for (int i = 0; i < n; ++i) {
        // printf("ii: %d\n", ii);
        indexPiv = ipiv[i];
        acc = b[indexPiv];
        b[indexPiv] = b[i];
        // if (ii >= 0) {  // Isn't this always true?
            for (int j = ii; j < i; ++j) acc -= A[j * n + i] * b[j];
        // } else if (acc) ii = i;

        b[i] = acc;
    }

    for (int i = n - 1; i >= 0; --i) {
        acc = b[i];
        for (int j = i+1; j < n; ++j) acc -= A[j * n + i] * b[j];
        b[i] = acc/A[i * n + i]; 
    }
}

int gesv(float * A, int n, float * b) {
    int info;
    int ipiv[n];

    if (! (info = croutLU(A, n, ipiv)) ) { // Successful factorization
        // printVectorI(ipiv, 3);
        backsub(A, n, ipiv, b); 
    }
    
    return info;
}

// M x k, and k x N
void MatVecGemm(int m, int vecLen, float * A, float * b, float * c) {
    float dotProd;
    for (int i = 0; i < m; ++i) {
        dotProd = 0.0;
        for (int j = 0; j < vecLen; ++j) {
            dotProd += A[j * m + i] * b[j];
        }
        c[i] = dotProd;
    }
}

void VecMatGemm(int vecLen, int n, float * b_t, float * A, float * c) {
    float dotProd;
    for (int j = 0; j < n; ++j) {
        dotProd = 0.0;
        for (int i = 0; i < vecLen; ++i) {
            dotProd += b_t[i] * A[j * vecLen + i];
        }
        c[j] = dotProd;
    }
}

float comb(int n, int k) {
    float prod = 1.;
    for (int i = 0; i < k; ++i) prod *= (n - i)/(float)(k - i);
}

int main() {
    // float matA[] = {3.,7.,1., 4.,2.,1., 1.,1.,0.};
    // float vecB[] = {0.,0.,1.};
    float matA[] = {0.00300,5.291, 59.14,-6.130};
    float vecB[] = {59.17, 46.78};
    int numEqs = 2;
    int info;

    printMatrix(matA, numEqs, numEqs);
    printVectorF(vecB, numEqs);

    if (! (info = gesv(matA, numEqs, vecB)) ) {
        // printMatrix(matA, 3, 3);
        printVectorF(vecB, numEqs);
    } else {
        printf("Singular Matrix\n");
    }

    float sqA[] = {3.,2.,4., 2.,2.,4.};
    float arrB[] = {1.,1.};
    float arrB_t[] = {1.,1.,1.};
    float c[3];

    // First try A * b
    MatVecGemm(3, 2, sqA, arrB, c);
    printVectorF(c, 3);
    // Then b_t * A
    VecMatGemm(3,2, arrB_t, sqA, c);
    printVectorF(c, 2);
    printf("%d choose %d = %f\n", 5, 3, comb(5,3));
    return 0;
}

