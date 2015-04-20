#define FP double

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <stdbool.h>
#define DIE(msg)    exit (fprintf (stderr, "%s\n", msg))
#define MAX(a,b)    (a > b)? a : b

// CUDA Error macros attributed to
// http://choorucode.com/2011/03/02/how-to-do-error-checking-in-cuda/

#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckLastError()    __cudaCheckError( __FILE__, __LINE__ )

void __cudaSafeCall(cudaError err, const char *file, const int line ) {
    if ( cudaSuccess != err ) {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
    return;
}
 
void __cudaCheckError( const char *file, const int line ) {
    cudaThreadSynchronize();
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err ) {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
    return;
}

// In: matrix a, b, dimensions x and y
// Out: matrix c
// Calculates A * B = C matrix using the GPU
__global__ void gpu_matrixmult(FP * A, FP * B, FP * C, int n, int p, int m) {
    int col = threadIdx.x + blockDim.x * blockIdx.x;
    int row = threadIdx.y + blockDim.y * blockIdx.y;

    if(row < n && col < m) {
        int index = row * m + col;
        FP cval = 0.;
        for (int indexa = row*p, indexb = col, i = 0; i < p; 
        ++indexa, indexb += m, ++i) 
        {
            cval += A[indexa]*B[indexb];
        }
        C[index] = cval;
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
    if (argc < 8) {
        DIE("Parse Error: usage: ./Task1GPU <n> <m> <p> <bx> <by> <gx> <gy>");
    } 
    CmdLineOpts_t options;
    options.n = atoi(argv[1]);
    options.m = atoi(argv[2]); 
    options.p = atoi(argv[3]); 
    options.bx = atoi(argv[4]); 
    options.by = atoi(argv[5]); 
    options.gx = atoi(argv[6]); 
    options.gy = atoi(argv[7]);

    if (options.n <= 0 || options.m <= 0 || options.p <= 0) {
        DIE("Parse Error: Matrix dimensions must be positive");
    } else if (options.bx <= 0 || options.by <= 0) {
        DIE("Parse Error: Block dimensions must be positive");
    } else if (options.gx <= 0 || options.gy <= 0) {
        DIE("Parse Error: Grid dimensions must be positive");
    } else if (options.bx * options.by > 1024) {
        DIE("Parse Error: Too many threads in block");
    } else if (options.by * options.gy < options.n) {
        DIE("Parse Error: no. of threads in y dim < no. of rows in matrix");
    } else if (options.bx * options.gx < options.m) {
        DIE("Parse Error: no. of threads in x dim < no. of cols in matrix");
    } else {
        return options;
    }
}

int main(int argc, char * argv[]) {
    /*----------------READ IN COMMAND LINE ARGS-----------------*/
    CmdLineOpts_t options = parseCmdLineArgs(argc, argv);

    /*------------------SET PARAMS, MATRICES--------------------*/
     
    FP * h_A, * h_B, * h_C; // host matrices
    FP * d_A, * d_B, * d_CIn; // device matrices
    FP * d_COut;

    int sizeA = options.n * options.p * sizeof(FP); 
    int sizeB = options.p * options.m * sizeof(FP); 
    int sizeC = options.n * options.m * sizeof(FP);

    dim3 Block(options.bx, options.by);
//    dim3 Grid(options.gx,options.gy);
    int gridx = MAX(1, (options.m+options.bx-1)/options.bx);
    int gridy = MAX(1, (options.n+options.by-1)/options.by);
    dim3 Grid(gridx,gridy);
    
    if (! (h_A = (FP*) malloc(sizeA)) ) DIE("Failed malloc");
    if (! (h_B = (FP*) malloc(sizeB)) ) DIE("Failed malloc");
    if (! (h_C = (FP*) malloc(sizeC)) ) DIE("Failed malloc");
    if (! (d_COut = (FP*) malloc(sizeC)) ) DIE("Failed malloc");
    
    srand(12345);

    /* INIT MATRICES */
    for (int i = 0; i < options.n; ++i) {
        for (int j = 0; j < options.p; ++j) {
            h_A[i*options.p+j] = (FP) rand() / (FP) RAND_MAX;
        }
    }

    for (int i = 0; i < options.p; ++i) {
        for (int j = 0; j < options.m; ++j) {
            h_B[i*options.m+j] = (FP) rand() / (FP) RAND_MAX;
        }
    }

    for (int i = 0; i < options.n * options.m; ++i) {
        h_C[i] = (FP) 0;
    }

    /*------------------COMPUTATION ON DEVICE--------------------*/
    cudaSafeCall(cudaMalloc((void**)&d_A, sizeA)); 
    cudaSafeCall(cudaMalloc((void**)&d_B, sizeB)); 
    cudaSafeCall(cudaMalloc((void**)&d_CIn, sizeC));
    
     
    cudaSafeCall(cudaMemcpy(d_A, h_A, sizeA, cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy(d_B, h_B, sizeB, cudaMemcpyHostToDevice));
    
    gpu_matrixmult<<<Grid,Block>>>(d_A, d_B, d_CIn,
        options.n, options.p, options.m);

    cudaCheckLastError();
    
    cudaSafeCall(cudaMemcpy(d_COut, d_CIn, sizeC, cudaMemcpyDeviceToHost));  
    
    printMatrix(d_COut, options.n, options.m);
    
    free(h_A);
    free(h_B);
    free(h_C);
    free(d_COut);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_CIn);
    
    return 0;
}

