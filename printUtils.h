#ifndef INCLUDED_PRINTUTILS_H
#define INCLUDED_PRINTUTILS_H

/* Misc utilities to print vectors and matrices for debugging */
void printVectorI (int * v, int size);
void printVectorF (int * v, int size);
void printMatrix(float * mat, int rows, int cols);
void printPair(int * actions1, int * actions2, int k);

#endif

