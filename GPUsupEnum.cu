#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include "linalg.h"

#define NTB 512 
#define DIE(msg) exit(fprintf(stderr, "Error: %s\n", msg))
#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckLastError()    __cudaCheckError( __FILE__, __LINE__ )

__constant__ float d_payoffsA[400];
__constant__ float d_payoffsB[400];
__constant__ int d_nActions1;
__constant__ int d_nActions2;

float * payoffsA, payoffsB;
int nActions1, nActions2;
char * supportArray, * d_supportArray;

// CUDA Error macros attributed to
// http://choorucode.com/2011/03/02/how-to-do-error-checking-in-cuda/

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

// Fill up all of the support combinations before sending to device
int kSubsetsChrItHelper(char * arr, int n, int k);


void kSubsetsPopulate(int k) {
    char arr[k];
    for (int i = 0; i < k; ++i) arr[i] = (char)i;
    int i = 0;
    do {
        memcpy(supportArray + (i++) * k, arr, k);
    } while (kSubsetsChrItHelper(arr, nActions1, k));
    
    cudaSafeCall(cudaMemcpy(d_supportArray, supportArray, 
        sizeAttheMoment, cudaMemcpyHostToDevice));
}

// Loads data from file
int readGame(char * gameData);

// Max elmts found at center
int maxSupportElmts() {
    int centerSupport = nActions1 % 2 == 0 ? 
                        nActions1 / 2 : (int)ceil(nActions1/2.);
    int maxNumSupports = (int)comb(nActions1, centerSupport);

    return maxNumSupports * centerSupport;
}

void initDeviceConstants() {
    // Load payoff matrices and constants to GPU
    cudaSafeCall(cudaMemcpyToSymbol(d_payoffsA, payoffsA, 
                 sizeof(float * nActions1 * nActions2)));
    cudaSafeCall(cudaMemcpyToSymbol(d_payoffsB, payoffsB,
                 sizeof(float * nActions1 * nActions2)));
    cudaSafeCall(cudaMemcpyToSymbol(d_nActions1, nActions1, sizeof(int)));
    cudaSafeCall(cudaMemcpyToSymbol(d_nActions2, nActions2, sizeof(int)));
    
    free(payoffsA); free(payoffsB);
}

__device__ void buildMat(char * acc1, char * acc2, int suppSize, float * mat, 
                         bool isB) 
{
    // both payoffsA and payoffsB are col major
    int k = 0;
    if (isB) {
        for (int i = 0; i < suppSize; ++i) {
            for (int j = 0; j < suppSize; ++j) 
                mat[k++] = d_payoffsB[acc2[j] * nActions1 + acc1[i]];
            mat[k++] = 1.; 
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    } else {
        for (int j = 0; j < suppSize; ++j) {
            for (int i = 0; i < suppSize; ++i) 
                mat[k++] = d_payoffsA[acc2[j] * nActions1 + acc1[i]];
            mat[k++] = 1.;
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    }
}

__device__ void buildFullStrat(char * acc, float * stratWeights, int sizeSubSet,
                               int sizeStrat, float * strategy)
{
    for (int i = 0; i < sizeStrat; ++i) strategy[i] = 0.;
    for (int i = 0; i < sizeSubSet; ++i) strategy[acc[i]] = stratWeights[i];
}

__global__ void gpu_nashEq(char * supArray, int k) {
    extern __shared__ char p1Support[]; 
    for (int i = 0; i < k; ++i) {
        p1Support[] = supArray[k * blockIdx.x + i];
    }

    int tid = threadIdx.x;
    int numSupports = comb();
    char * p2Support;

    int matSize = k+1;
    float matA[matSize * matSize];
    
    // Cycle as many times as it takes to perform one's work
    for (int i = 0; i * k + tid < numSupports; ++i) {
        p2Support = supArray + i*k+tid; 
    }
    
    // Assign associated segment of array with correct block
    // Allocate EXACTLY enough blocks, no need to turn off
    // Merely need to add additional cycles per block
    // No. cycles = number of supports / number of threads in block
}

int main() {
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("No. Devices: %d\n", devCount);
    if (argc < 2) {
        DIE("Incorrect command line args"); 
    } else {
        // Pull payoff info and action spaces of both players from file
        readGame(argv[1]);

        initDeviceConstants();

        // Allocate space to store all strategy combinations
        int maxSupportElmts = maxSupportElmts();
        char supArr[maxSupportElmts]; supportArray = supArr;
        cudaSafeCall(cudaMalloc((void**)&d_supportArray, maxSupportElmts)); 
        
        for (int i = 1; i <= nActions1; ++i) {
            kSubsetsPopulate(i);    // Send all subsets to device of size k
            /* Compute optimal number of blocks; I guess just the number
             * of combinations of a given size support? */
            int numBlocks = comb(nActions1, i);

            
            gpu_nashEq<<<numBlocks, NTB, i>>>(d_supportArray, i); 
        }

        cudaFree(d_supportArray);
    }

    return 0;
}

