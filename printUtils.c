#include <stdio.h>

void printVectorI(int * v, int size) {
    printf("[");
    for (int i=0; i<size; ++i) printf("%d ", v[i]);
    printf("]\n");
}

void printVectorF(float * v, int size) {
    printf("[");
    for (int i=0; i<size; ++i) printf("%.5f ", v[i]);
    printf("]\n");
}

// Print matrix; matrix in column maj order
void printMatrix(float * mat, int rows, int cols) {
    printf("[\n");
    for (int i = 0; i < rows; ++i) {
        printf("\t[ ");
        for (int j = 0; j < cols; ++j) {
            printf("%f ", mat[j*rows + i]); 
        }
        printf("],\n");
    }
    printf(" ]\n");
}

// Debug function to determine correctness of subset enumeration
void printPair(int * acc1, int * acc2, int k) {
    printf("Acc1: ");
    for (int i = 0; i < k; ++i) {
        printf("%d ", acc1[i]); 
    }
    printf("Acc2: ");
    for (int i = 0; i < k; ++i) {
        printf("%d ", acc2[i]); 
    }
    printf("\n");
}

void printSolution(float * strat1, float * strat2, int n1, int n2) {
    printf("SOLUTION: ");
    for (int i = 0; i < n1; ++i) {
        printf("%f ", strat1[i]);
    }
    printf("| ");
    for (int i = 0; i < n2; ++i) {
        printf("%f ", strat2[i]);
    }
    printf("\n");
}

