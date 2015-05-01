#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <cuda.h>

// #ifndef INCLUDED_LINALG_H
// #include "linalg.h"
// #endif

// #define NTB 512
#define MAXSIZE 20
#define MAXMATSIZE 17
#define MAXSUPPORTSIZE 16
 
#define DIE(msg) exit(fprintf(stderr, "Error: %s\n", msg))
#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckLastError()    __cudaCheckError( __FILE__, __LINE__ )

__constant__ float d_payoffsA[400];
__constant__ float d_payoffsB[400];
__constant__ int d_nActions1;
__constant__ int d_nActions2;
__constant__ int NTB;

float * payoffsA, * payoffsB;
int nActions1, nActions2;
int * supportArray; 
int * d_supportArray;

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
//////////////////////////////////////////////////////////////////////

__device__ int croutLU(float * A, int n, int * ipiv) {
    int argmaxi;
    float max, acc, temp;
    float vv[MAXSIZE];

    for (int i = 0; i < n; ++i) {     /* Get implicit pivot scaling info */ 
        max = 0.0;                    /* for each row                    */
        for (int j = 0; j < n; ++j) {
            max = fmaxf(max, fabsf(A[j * n + i]));
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

__device__ void backsub(float * A, int n, int * ipiv, float * b) {
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

__device__ int gesv(float * A, int n, float * b) {
    int info;
    int ipiv[MAXSIZE];

    if (! (info = croutLU(A, n, ipiv)) ) { // Successful factorization
        // printVectorI(ipiv, 3);
        backsub(A, n, ipiv, b); 
    }
    
    return info;
}

// M x k, and k x N
__device__ void MatVecGemm(int m, int vecLen, float * A, float * b, float * c) {
    float dotProd;
    for (int i = 0; i < m; ++i) {
        dotProd = 0.0;
        for (int j = 0; j < vecLen; ++j) {
            dotProd += A[j * m + i] * b[j];
        }
        c[i] = dotProd;
    }
}

__device__ void VecMatGemm(int vecLen, int n, float * b_t, float * A, float * c) {
    float dotProd;
    for (int j = 0; j < n; ++j) {
        dotProd = 0.0;
        for (int i = 0; i < vecLen; ++i) {
            dotProd += b_t[i] * A[j * vecLen + i];
        }
        c[j] = dotProd;
    }
}

__host__ __device__ double comb(int n, int k) {
    double prod = 1.;
    for (int i = 0; i < k; ++i) prod *= (n - i)/(double)(k - i);
    return prod;
}
/////////////////////////////////////////////////////////////////////////////

// Fill up all of the support combinations before sending to device
int kSubsetsItHelper(int * arr, int n, int k) {
    int finished = 0;
    int changed = 0;
    int i, j;
    if (k > 0) {
        for (i = k - 1; !finished && !changed; --i) {
            if (arr[i] < (n-1) - (k-1) + i) {
                ++arr[i];
                if (i < k - 1) {
                    for (j = i + 1; j < k; ++j) {
                        arr[j] = arr[j-1] + 1;
                    }
                }
                changed = 1;
            }
            finished = i == 0;
        }
        if (!changed) {
            for (i = 0; i < k; ++i) {
                arr[i] = i;
            }
        }
    }
    return changed;
}

void kSubsetsPopulate(int k) {
    int arr[k];
    for (int i = 0; i < k; ++i) arr[i] = i;
    int i = 0;
    do {
        memcpy(supportArray + (i++) * k, arr, k*sizeof(int));
    } while (kSubsetsItHelper(arr, nActions1, k));
   
    int size = sizeof(int) * i * k; 
    /*for (int j = 0; j < i; ++j) {
        for (int l = 0; l < k; ++l) {
            printf("%d ", supportArray[j*k + l]);
        }
        printf("\n");
    }
    return;*/
    cudaSafeCall(cudaMemcpy(d_supportArray, supportArray, 
        size, cudaMemcpyHostToDevice));
}

// Loads data from file
int readGame(char * gameData) {
  FILE * fp;
  char line[60];
  if ( (fp = fopen(gameData, "r")) ) {
//    int nActions1 = 1, nActions2 = 1;
    float payoff1 = 0., payoff2 = 0.;

    /* Get no. actions for player 1 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player1: %d", &nActions1);
    
    /* Get no. actions for player 2 */
    fgets(line, sizeof(line), fp);
    sscanf(line, "player2: %d", &nActions2);

    payoffsA = (float *)malloc(nActions1 * nActions2 * sizeof(float));
    payoffsB = (float *)malloc(nActions1 * nActions2 * sizeof(float));

    int i = 0;
    /* Get all of the payoffs */
    while (fgets(line, sizeof(line), fp)) {
      sscanf(line, "%f %f", &payoff1, &payoff2);
      payoffsA[i] = payoff1, payoffsB[i] = payoff2;
      ++i;
    }
    fclose(fp);
    return 1;
  } else {
    return 0;  
  }
}

// Max elmts found at center
int maxSupportElmts() {
    int max = 0;
    for (int i = nActions1 / 2; i <= nActions1; ++i) {
        int numSupports = (int)comb(nActions1, i);
        if (numSupports * i > max) max = numSupports * i;
    }
    /*int centerSupport = nActions1 % 2 == 0 ? 
                        nActions1 / 2 : (int)ceil(nActions1/2.);
    int maxNumSupports = (int)comb(nActions1, centerSupport);

    return maxNumSupports * centerSupport;*/

    return max + 1;
}

void initDeviceConstants() {
    // Load payoff matrices and constants to GPU
    cudaSafeCall(cudaMemcpyToSymbol(d_payoffsA, payoffsA, 
                 sizeof(float) * nActions1 * nActions2));
    cudaSafeCall(cudaMemcpyToSymbol(d_payoffsB, payoffsB,
                 sizeof(float) * nActions1 * nActions2));
    cudaSafeCall(cudaMemcpyToSymbol(d_nActions1, &nActions1, sizeof(int)));
    cudaSafeCall(cudaMemcpyToSymbol(d_nActions2, &nActions2, sizeof(int)));
    
    free(payoffsA); free(payoffsB);
}

__device__ void buildMat(int * acc1, int * acc2, int suppSize, float * mat, 
                         bool isB) 
{
    // both payoffsA and payoffsB are col major
    int k = 0;
    if (isB) {
        for (int i = 0; i < suppSize; ++i) {
            for (int j = 0; j < suppSize; ++j) 
                mat[k++] = d_payoffsB[acc2[j] * d_nActions1 + acc1[i]];
            mat[k++] = 1.; 
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    } else {
        for (int j = 0; j < suppSize; ++j) {
            for (int i = 0; i < suppSize; ++i) 
                mat[k++] = d_payoffsA[acc2[j] * d_nActions1 + acc1[i]];
            mat[k++] = 1.;
        }
        for (int i = 0; i < suppSize; ++i) mat[k++] = 1.;
        mat[k] = 0.;
    }
}

__device__ void buildFullStrat(int * acc, float * stratWeights, int sizeSubSet,
                               int sizeStrat, float * strategy)
{
    for (int i = 0; i < sizeStrat; ++i) strategy[i] = 0.;
    for (int i = 0; i < sizeSubSet; ++i) strategy[acc[i]] = stratWeights[i];
}

__device__ void printSolution(float * strat1, float * strat2, int n1, int n2) {
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

__device__ void doNothing() {
    return;
}

__global__ void gpu_nashEq(int * supArray, int k) {
    int tid = threadIdx.x;
    int numSupports = comb(d_nActions1, k);
    extern __shared__ int p1Support[]; 
    // int * p1Support = supArray + (k * blockIdx.x);
    for (int i = 0; i < k; ++i) {
        p1Support[i] = supArray[k * blockIdx.x + i];
    }

    int * p2Support;

    int matSize = k + 1;
    float matA[MAXMATSIZE*MAXMATSIZE];
    float vecX[MAXMATSIZE], vecY[MAXMATSIZE];
    
    // Cycle as many times as it takes to perform one's work
    for (int i = 0; (i * NTB + tid) < numSupports; ++i) {
        p2Support = supArray + i*NTB*k+tid*k;

        for (int j = 0; j < matSize; ++j) vecX[j] = 0., vecY[j] = 0.;
        
        vecX[matSize-1] = 1., vecY[matSize-1] = 1.;
        
        buildMat(p1Support, p2Support, k, matA, true);
        if (! (gesv(matA, matSize, vecX)) ) {
            float min = vecX[0];
            for (int i = 0; i < k; ++i) min = fminf(min, vecX[i]);
            if (min >= -0.00001) {
                float strat1[MAXSUPPORTSIZE];
                float xB[MAXSUPPORTSIZE];
                buildFullStrat(p1Support, vecX, k, d_nActions1, strat1);
                VecMatGemm(d_nActions1, d_nActions2, strat1, d_payoffsB, xB);
                float maxxB = xB[0];
                for (int i = 0; i < d_nActions2; ++i) maxxB = fmaxf(xB[i], maxxB);
            
                if (xB[p2Support[0]] >= (maxxB - fabs(maxxB / 1000.))) {
                    
                    buildMat(p1Support, p2Support, k, matA, false);
                    if (! (gesv(matA, matSize, vecY)) ) {
                        min = vecY[0];
                        for (int i = 0; i < k; ++i) min = fminf(min, vecY[i]);
                        if (min >= 0.) {
                            float strat2[MAXSUPPORTSIZE];
                            float Ay[MAXSUPPORTSIZE];
                            buildFullStrat(p2Support, vecY, k, d_nActions2, strat2);
                            MatVecGemm(d_nActions1, d_nActions2, d_payoffsA, 
                                      strat2, Ay);
                            float maxAy = Ay[0];
                            for (int i = 0; i < d_nActions1; ++i) 
                                maxAy = fmaxf(Ay[i], maxAy);
                            if(Ay[p1Support[0]] >= (maxAy - fabsf(maxAy / 1000.)))
                                // printf
                                //for (int i = 0; i < NTB; ++i) {
                                //    if (i == tid) 
                                    doNothing();
                                        // printSolution(strat1, strat2, d_nActions1, d_nActions2);
                                //}
                        }
                    }
                }
            }
        }
    }
    // Assign associated segment of array with correct block
    // Allocate EXACTLY enough blocks, no need to turn off
    // Merely need to add additional cycles per block
    // No. cycles = number of supports / number of threads in block
}

int main(int argc, char * argv[]) {
    int devCount;
    cudaGetDeviceCount(&devCount);
    // printf("No. Devices: %d\n", devCount);
    if (argc < 2) {
        DIE("Incorrect command line args"); 
    } else {
        // Pull payoff info and action spaces of both players from file
        readGame(argv[1]);
        int NTBs = 512;
        if (argc == 3) {
            NTBs = atoi(argv[2]);
        }
        cudaSafeCall(cudaMemcpyToSymbol(NTB, &NTBs, sizeof(int)));
        // printf("Threads per block: %d\n", NTBs);

        initDeviceConstants();

        // Allocate space to store all strategy combinations
        int maxSupports = maxSupportElmts();
        int supArr[maxSupports]; supportArray = supArr;
        cudaSafeCall(cudaMalloc((void**)&d_supportArray, maxSupports*sizeof(int))); 
        
        for (int i = 1; i <= nActions1; ++i) {
            double numBlocks = comb(nActions1, i);
            // printf("nB : %f\n", numBlocks);
            kSubsetsPopulate(i);    // Send all subsets to device of size k
            /* Compute optimal number of blocks; I guess just the number
             * of combinations of a given size support? */
            // continue;
            
            gpu_nashEq<<<numBlocks, NTBs, i*(sizeof(int))>>>(d_supportArray, i); 
        }

        cudaFree(d_supportArray);
    }

    return 0;
}

