#define FP double

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define DIE(msg)    exit (fprintf (stderr, "%s\n", msg))
#define MAX(a,b)    (a > b)? a : b

// Calculates A * B = C matrix using the CPU, single threaded
// Using the "kij" variant of the matrix multiplication
/*
K is columns of A/Rows of B
I is rows of A/C
J is columns of B/C
Fix element of A
Fix row of B 
*/
void cpu_matrixmult(FP * A, FP * B, FP * C, int n, int p, int m) {
    FP r;
    for (int k = 0; k < p; k++) {
        for (int i = 0; i < n; i++) { // Rows of A, rows of C
            r = A[i*p + k];
            for (int j = 0; j < m; j++) C[i*m+j] += r * B[k*m+j];
        }
    }
}

void printMatrix(FP * C, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%f\t", C[i * m + j]);
        }
        printf("\n");
    }
}

// Holds all cmd line args for this programs
typedef struct CmdLineOpts_t {
    int n;  // rows of A, C
    int m;  // cols of B, C
    int p;  // cols of A, rows of B
    int bx; // block numRows of threads
    int by; // block numCols of threads
    int gx; // grid numRows of blocks
    int gy; // grid numCols of blocks
} CmdLineOpts_t;

CmdLineOpts_t parseCmdLineArgs(int argc, char * argv[]) {
    if (argc < 4) {
        DIE("Parse Error: usage: ./Task1CPU <n> <m> <p>");
    } 
    CmdLineOpts_t options;
    options.n = atoi(argv[1]);
    options.m = atoi(argv[2]); 
    options.p = atoi(argv[3]); 

    if (options.n <= 0 || options.m <= 0 || options.p <= 0) {
        DIE("Parse Error: Matrix dimensions must be positive");
    } else {
        return options;
    }
}

int main(int argc, char * argv[]) {
    CmdLineOpts_t options = parseCmdLineArgs(argc, argv);

    int i,j;

    FP * h_A, * h_B, * h_C; // host matrices

    int sizeA = options.n * options.p * sizeof(FP); 
    int sizeB = options.p * options.m * sizeof(FP); 
    int sizeC = options.n * options.m * sizeof(FP);

    if (! (h_A = (FP*) malloc(sizeA)) ) DIE("Failed malloc");
    if (! (h_B = (FP*) malloc(sizeB)) ) DIE("Failed malloc");
    if (! (h_C = (FP*) malloc(sizeC)) ) DIE("Failed malloc");
    
    srand(12345);

    /* INIT MATRICES */
    for (i = 0; i < options.n; ++i) {
        for (j = 0; j < options.p; ++j) {
            h_A[i*options.p+j] = (FP) rand() / (FP) RAND_MAX;
        }
    }

    for (i = 0; i < options.p; ++i) {
        for (j = 0; j < options.m; ++j) {
            h_B[i*options.m+j] = (FP) rand() / (FP) RAND_MAX;
        }
    }

    for (i = 0; i < options.n * options.m; ++i) {
        h_C[i] = (FP) 0;
    }

    cpu_matrixmult(h_A, h_B, h_C, 
        options.n, options.p, options.m); // do calculation on host

    printMatrix(h_C, options.n, options.m);
    
    free(h_A);
    free(h_B);
    free(h_C);

    return 0;
}

