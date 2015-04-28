/* All credit attributed to Rampersaud and Wayne State University */

#include <Windows.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <limits.h>
#include <time.h>

#define TRUE  1
#define FALSE 0
#define swap(af,bf) (((af)!=(bf))? af^=bf^=af^=bf:af)
#define defRow 5
#define defCol 5
#define epsilon 1.0e-4
#define MAX_THREADS 1024

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

using namespace std;

////////////////////////////////////////////// DEVICE & HOST COMPUTE //////////////////////////////////////////////
__device__	int return_max(int *a, int size);
__host__ __device__	bool next(int n_, int *p_);
__device__ void ludcmp_array(float *a, int n, int *indx, float *d);
__device__ void bcksub_array(float *a, int n, int *indx, float *b);
__device__ bool validityCheck(float *, int);
////////////////////////////////////////////// DEVICE & HOST COMPUTE //////////////////////////////////////////////

////////////////////////////////////////////// PURE STRATEGY COMPUTE //////////////////////////////////////////////
__global__ void pureStrategyCompute	(	
										float* glo_array_payA, 
										float* glo_array_payB,
										int* glo_principalContainer,
										int *glo_returnEquilibrium
									);

cudaError_t pureCalcWrapper			(	
										float* __array_payA, 
										float* __array_payB,
										int glo_principalCount,
										int* glo_principalContainer,
										FILE* _f
									);
////////////////////////////////////////////// PURE STRATEGY COMPUTE //////////////////////////////////////////////


////////////////////////////////////////////// MIXED STRATEGY COMPUTE //////////////////////////////////////////////
__global__ void mixedStrategyCompute(	
									float* glo_array_payA, 
									float* glo_array_payB,
									int* glo_principalContainer,
									int* glo_support_iter,
									float* glo_returnEquilibrium
									);
cudaError_t mixedCalcWrapper		(
									float* __array_payA, 
									float* __array_payB, 
									int __principalCount, 
									int *__principalContainer, 
									int __support_iter,
									FILE* __f
									);

////////////////////////////////////////////// MIXED STRATEGY COMPUTE //////////////////////////////////////////////
 
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
 
struct timezone
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
// Definition of a gettimeofday function

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
// Define a structure to receive the current Windows filetime
  FILETIME ft;
 
// Initialize the present time to 0 and the timezone to UTC
  unsigned __int64 tmpres = 0;
  static int tzflag = 0;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
// Convert to microseconds by dividing by 10
    tmpres /= 10;
 
// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
// in seconds from Jan 1 1601.
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
 
// Finally change microseconds to seconds and place in the seconds value. 
// The modulus picks up the microseconds.
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 
  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
  
// Adjust for the timezone west of Greenwich
      tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
  return 0;
}

int main(int argc, char *argv[])
{
	double time_start, time_end;
	struct timeval tv; 
	struct timezone tz;
	gettimeofday(&tv, &tz); 
	time_start = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0; 

		FILE *file = fopen("parallelNashGPU.asc", "w");        
        int i, j;
        int support_iter;
		float **payoffA;
        float **payoffB;
        float *array_payA;
        float *array_payB;
   
        int _p1[defRow];                        // = {0,0,1};					//Support tracker. for _p1 support = x, ... 
		
		int principalCount;
		int *principalContainer;
		int *principalSupport;  
		int _pcCount[defRow];
		int pc_iterate;
		
// Get Matrix Entries from Directory PlayerA.txt & PlayerB.txt into arrays 
        {        
                FILE* input;
                payoffA                         = (float**) malloc(defRow * sizeof(float*));
                payoffB                         = (float**) malloc(defRow * sizeof(float*));
                
                for(i = 0; i < defRow; i++)
                {
                         payoffA[i]                 = (float*) malloc(defCol * sizeof(float));
                         payoffB[i]                 = (float*) malloc(defCol * sizeof(float));
                }
        
        for(i = 0; i < defRow; i++)
        
                        input = fopen("PlayerA.txt", "r");
                        for(j = 0; j < defRow; j++)
                                for(i = 0; i < defCol; i++)
                                        fscanf(input, "%f",&payoffA[j][i]);
                        fclose(input);

                        input = fopen("PlayerB.txt", "r");
                        for(j = 0; j < defRow; j++)
                                for(i = 0; i < defCol; i++)
                                        fscanf(input, "%f",&payoffB[j][i]);
                        fclose(input);
                        
                array_payA = (float *) malloc((float) (defCol*defRow) * sizeof(float));
                array_payB = (float *) malloc((float) (defCol*defRow) * sizeof(float));                        
                     
                for(i = 0; i < defRow; i++)
                {
                        for(j = 0; j < defCol; j++)
                        {        
                                array_payA[i*defCol + j] = (payoffA[i][j]);
                                array_payB[i*defCol + j] = (payoffB[i][j]);
                        }
                } 
		}  
// ........................................ CALCULATING PURE NASH .............................................//     
//Build Support Set        
        for(int i = 0; i < defRow-1; i++)
                {
					_p1[i] = 0;
					_pcCount[i] = 0;
				}
		
		_p1[defRow-1] = 1;
		_pcCount[defCol-1] = 1;

//Build Principal Component Support 
			principalCount = 0;
			do	{
					principalCount++;
				}while(next(defRow,_pcCount));  
			principalSupport = (int *)malloc((int) ((principalCount) * sizeof(int)));
			for(i = 0; i < defCol; i++)
				principalSupport[i] = 0;

			principalContainer = (int *)malloc((int) ((principalCount*defRow) * sizeof(int)));  

			for(i = 0; i < (principalCount*defRow); i++)
				principalContainer[i] = 0;
			for(i = 0; i < defCol; i++)
                principalSupport[i] = _p1[i];

			pc_iterate = 0;
			do	{
					for(int ps_iterate = 0; ps_iterate < defCol; ps_iterate++)
						{
							principalContainer[pc_iterate*defRow + ps_iterate] = principalSupport[ps_iterate];				
						}	
						pc_iterate++;
				}while(next(defRow,principalSupport));
	cudaError_t cudaStatus = pureCalcWrapper(array_payA, array_payB, principalCount, principalContainer, file);	
    if (cudaStatus != cudaSuccess)
	{
        fprintf(stderr, "pureCalcWrapper failed!");
		system("pause");
        exit(0);
    }
	// ........................................ CALCULATING PURE NASH .............................................//

	// ........................................ CALCULATING MIXED NASH ............................................//
for(support_iter = 2; support_iter <= defRow; support_iter++)         
{        
        fprintf(file,"\n------------------------------------------------");
		fprintf(file,"\nSupport \t%d\n", support_iter);
		fprintf(file,"------------------------------------------------\n");
        {
                // Build Support Set
                for(i = 0; i < defRow-support_iter; i++)
                        {
							_p1[i] = 0;
							_pcCount[i] = 0;
						}
                for(i = defRow - support_iter; i < defRow; i++)
                        {
							_p1[i] = 1;
							_pcCount[i] = 1;
						}

		}
		
//Build Principal Component Support 
			principalCount = 0;
			do	{
					principalCount++;
				}while(next(defRow,_pcCount)); 
			principalSupport = (int *)malloc(((principalCount) * sizeof(int)));

			for(int z1 = 0; z1 < principalCount; z1++)
				principalSupport[z1] = 0;
		
			principalContainer = (int *)malloc((principalCount*defRow) * sizeof(int));  

			for(int z2 = 0; z2 < (principalCount*defRow); z2++)
				principalContainer[z2] = 0;

			for(int z3 = 0; z3 < defCol; z3++)
                principalSupport[z3] = _p1[z3];

			pc_iterate = 0;
			do	{
					for(int ps_iterate = 0; ps_iterate < defRow; ps_iterate++)
						{
							principalContainer[pc_iterate*defRow + ps_iterate] = principalSupport[ps_iterate];
						}	
						pc_iterate++;
				}while(next(defRow,principalSupport)); 
	cudaStatus = mixedCalcWrapper(array_payA, array_payB, principalCount, principalContainer, support_iter, file);	
    if (cudaStatus != cudaSuccess)
	{
        printf("\n mixedCalcWrapper on HOST failed!\n");
		system("pause");
        exit(0);
    }
}
// ........................................ CALCULATING MIXED NASH ............................................//
	gettimeofday(&tv, &tz); 
	time_end = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0; 
	printf("%lf\n", time_end - time_start); 
	system("pause");
  return 0;
}

cudaError_t mixedCalcWrapper(float* __array_payA, float* __array_payB, int __principalCount, int *__principalContainer, int __support_iter, FILE* __f)
{
	int s1,s3,numOfThreads, numOfBlocks;
	cudaError_t cudaStatus = cudaSetDevice(0);

	// Transfer Variables inputted from CPU Wrapper Function 
	
	float *dev_array_payA = 0;
	float *dev_array_payB = 0;
	int* dev_principalContainer = 0;
	float*dev_returnEquilibrium = 0;	
		float *returnEquilibrium  = (float*) malloc(__principalCount*__principalCount*(defRow + defCol + 1) * sizeof(float));
	int * dev_support = 0;
		int *wrapperSupport = (int*) malloc(3 *sizeof(int));
		wrapperSupport[0] = __support_iter;
		wrapperSupport[1] = 0;
		wrapperSupport[2] = (int) ceil((float) __principalCount / MAX_THREADS);
////////////////////////////////////////////// Mixed Strategy CUDA Malloc  //////////////////////////////////////////////	

	cudaStatus = cudaMalloc((void**)&dev_array_payA, defCol*defRow*sizeof(float));
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payA malloc failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMalloc((void**)&dev_array_payB, defCol*defRow*sizeof(float));
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payB malloc failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMalloc((void**)&dev_principalContainer, __principalCount*defRow*sizeof(int));
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_principalContainer malloc failed!");
			system("pause");
			exit(0);
		}

	cudaStatus = cudaMalloc((void**)&dev_returnEquilibrium, __principalCount*__principalCount*(defRow + defCol + 1)*sizeof(float));
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_returnEquilibrium malloc failed!");
			system("pause");
			exit(0);
		}

		cudaStatus = cudaMalloc((void**)&dev_support, 3*sizeof(int));
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_support malloc failed!");
			system("pause");
			exit(0);
		}
		
////////////////////////////////////////////// Mixed Strategy CUDA Malloc ///////////////////////////////////////////////

////////////////////////////////////////////// Mixed Strategy CUDA Memcopy //////////////////////////////////////////////

	cudaStatus = cudaMemcpy(dev_array_payA,	__array_payA,	defCol*defRow*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payA memcpy failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMemcpy(dev_array_payB,	__array_payB,	defCol*defRow*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMemcpy(dev_principalContainer,	__principalContainer, __principalCount*defRow*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}
	printf("\nSupport : %d\t%d\n",__support_iter,__principalCount*__principalCount*(defRow + defCol + 1));
	cudaStatus = cudaMemcpy(dev_returnEquilibrium,	returnEquilibrium, __principalCount*__principalCount*(defRow + defCol + 1)*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}		

		cudaStatus = cudaMemcpy(dev_support, wrapperSupport, 3*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}		

////////////////////////////////////////////// Mixed Strategy CUDA Execution //////////////////////////////////////////////	

	if(__support_iter == defRow)
		{
			numOfThreads = defRow;
			numOfBlocks	 = 1;		
		}
	else if(__support_iter > 1 && __support_iter < 9)
		{
			numOfThreads = __principalCount;
			numOfBlocks = __principalCount;
		}
	else if(__support_iter > 9 && __support_iter < 13)
		{
			numOfThreads = __principalCount*2;
			numOfBlocks = __principalCount;
		}
	else
		{
			numOfThreads = MAX_THREADS;
			numOfBlocks = __principalCount;			
		}

	dim3 dimBlock(numOfThreads,1,1);
	dim3 dimGrid(numOfBlocks,1,1);
	mixedStrategyCompute<<<dimGrid,dimBlock>>>(dev_array_payA, dev_array_payB, dev_principalContainer, dev_support, dev_returnEquilibrium);
    cudaStatus = cudaDeviceSynchronize();
 
	if (cudaStatus != cudaSuccess) 
	{
        printf("mixedStrategyCompute returned error code %d after launching addKernel!\n", cudaStatus);
		system("pause");
		exit(0);
	}
	if(defRow < 13)
	{
		cudaStatus = cudaMemcpy(returnEquilibrium, dev_returnEquilibrium, (__principalCount*__principalCount)*(defRow + defCol + 1)*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("mixed dev_returnEquilibrium memcpy failed!");
			system("pause");
			exit(0);
		}	

			for(s1 = 0; s1 < (__principalCount*__principalCount); s1++)
					{                
						if(returnEquilibrium[s1*(defCol+defRow+1)] ==  1.0)
							{               			
									fprintf(__f,"\nP_1:\t(");
									for(s3 = 0; s3 < defRow; s3++)
										fprintf(__f," %10.10f ",returnEquilibrium[s1*(defCol+defRow+1) + s3 + 1]);
                                        fprintf(__f,")\nP_2:\t(");
                                    for(s3 = defRow; s3 < (defCol+defRow); s3++)
										fprintf(__f," %10.10f ",returnEquilibrium[s1*(defCol+defRow+1) + s3 + 1]);
                                        fprintf(__f,")\n");
                            }
					} 

	}
	return cudaStatus;
}

__global__ void mixedStrategyCompute(float* glo_array_payA, float *glo_array_payB, int *glo_principalContainer, int *glo_support, float *glo_returnEquilibrium)
{
///////////////////////////////////////// LOCAL VARIABLES /////////////////////////////////////////
	int i;
	int j;
	int k;
	int s1;
	int s2;
	
	float _payoffA[defRow*defRow*defRow]; 					
    float _payoffB[defRow*defRow*defRow]; 					
	float _x[defRow];
	float _y[defRow];
	int _idx_A[defRow];
	int _idx_B[defRow];
	int support_iter = glo_support[0];

	float placeHolderA;
	float placeHolderB;

	register int loop;
	register int iterate;
	float _Xmax; 
	float _Ymax;
	register bool proceed;
	register bool checkNoPlay;

	float loc_x[defRow];
	float loc_y[defRow];
	float loc_Ay[defRow];
	float loc_xB[defRow];

///////////////////////////////////////// LOCAL VARIABLES /////////////////////////////////////////
for(int masterloop = 1; masterloop <= glo_support[2]; masterloop++)
{
//initAll
for(int init = 0; init < defRow; init++)
{
	loc_x[init] = 0.0;
	loc_y[init] = 0.0;
	loc_Ay[init] = 0.0;
	loc_xB[init] = 0.0;
}
			for(i = 0; i < defRow*defRow; i++)
			{
				_payoffA[i] = 0.0;
				_payoffB[i] = 0.0;
			}

			for(i = 0; i < defRow; i++)
			{
				_x[i] = 0.0;
				_y[i] = 0.0;
				_idx_A[i] = 0;
				_idx_B[i] = 0;
			}
{	  
loop = 0;			
iterate = 0;       
	for(i = defRow-1; i > 0; i--)
        {
			if(glo_principalContainer[masterloop*(threadIdx.x*defRow + i)] == 1)
            {                               
				for(k = 0; k <= i - 1; k++)
					{
						if(glo_principalContainer[masterloop*(threadIdx.x*defRow + k)] == 1)
						{	
							for(j = 0; j < defCol; j++)
							{
								if(glo_principalContainer[masterloop*(blockIdx.x*defRow + j)] == 1)
								{	
									_payoffB[iterate] = glo_array_payB[j*defCol + k] - glo_array_payB[j*defCol + i];
									iterate++;
								}				
						}	
					}
				}
			}		
        }
}
			
	for(i = 0; i < support_iter; ++i)
		{
			_x[i] = 0.0;		
			_payoffB[(support_iter-1)*support_iter + i] = 1.0;
		}
		
		_x[(support_iter-1)] = 1.0;
	
			ludcmp_array(_payoffB, support_iter, _idx_B, &placeHolderB);					// device function
			bcksub_array(_payoffB, support_iter, _idx_B, _x);								// device function
			loop = 0;	
				{
					for(i = 0; i < defCol; i++)
						{
							if(glo_principalContainer[masterloop*(blockIdx.x*defRow + i)] == 0)
								{
									loc_x[i] = _x[i]*( glo_principalContainer[masterloop*(blockIdx.x*defRow + i)]);
								}
							else
								{
									loc_x[i] = _x[loop]*( glo_principalContainer[masterloop*(blockIdx.x*defRow + i)]);
									loop++;								
								}
						}
				
				for(s1= 0; s1 < defCol; s1++)
                    {
						for(s2= 0; s2 < defRow; s2++)
                        {                                        
                            loc_xB[s1] += __fmul_rn(loc_x[s2], glo_array_payB[s2*defCol + s1]);
                        }
                    }

					_Xmax = -10000000.00;
                    s2 = 0;
					for(s1 = 0; s1 < defCol; s1++)
					{
						if(_Xmax < loc_xB[s1])
                        {
							s2    = s1;
                            _Xmax = loc_xB[s1];
                        }
					}	
						proceed = FALSE;

						if(glo_principalContainer[masterloop*(threadIdx.x*defRow + s2)] == 1)
                        {
							proceed = TRUE;
						}					
						else
						{
							proceed = FALSE; 
						}
				}
__syncthreads();
if(proceed)
{
loop = 0;		
iterate = 0;       
	for(i = defCol-1; i > 0; i--)
        {
			if(glo_principalContainer[masterloop*(blockIdx.x*defRow + i)] == 1)
            {                               
				for(k = 0; k <= i-1; k++)
					{
						if(glo_principalContainer[masterloop*(blockIdx.x*defRow + k)] == 1)
						{	
							for(j = 0; j < defCol; j++)
							{
								if(glo_principalContainer[masterloop*(threadIdx.x*defRow + j)] == 1)
								{
									_payoffA[iterate] = glo_array_payA[k*defCol + j] - glo_array_payA[i*defCol + j];
									iterate++;									
								}				
						}	
					}
				}
			}		
        }
	__syncthreads();

	for(i = 0; i < support_iter; ++i)
		{
			_y[i] = 0.0;		
			_payoffA[(support_iter-1)*support_iter + i] = 1.0;
		}
					
		_y[(support_iter-1)] = 1.0;

			ludcmp_array(_payoffA, support_iter, _idx_A, &placeHolderA);		//device function
			bcksub_array(_payoffA, support_iter, _idx_A, _y);					//device function

				{	
					for(i = 0; i < defRow; i++)
						{
							if(glo_principalContainer[masterloop*(threadIdx.x*defRow + i)] == 0)
								{
									loc_y[i] = (float) _y[i] * glo_principalContainer[masterloop*(threadIdx.x*defRow + i)];
								}
							else
								{
									loc_y[i] = (float) _y[loop]*glo_principalContainer[masterloop*(threadIdx.x*defRow + i)];
									loop++;								
								}
						}

				for(s1 = 0; s1 < defRow; s1++)
                    {
                        for(s2 = 0; s2 < defCol; s2++)
                        {                                        
							loc_Ay[s1] += __fmul_rn(loc_y[s2], glo_array_payA[s1*defCol + s2]);
						}
                    }			

					_Ymax = 0;
                    s2 = 0;

                    for(s1 = 0; s1 < defRow; s1++)
                    {
						if(_Ymax < loc_Ay[s1])
                        {
							s2   = s1;
                            _Ymax = loc_Ay[s1];
                        }
					}	
						proceed = FALSE;
						if(glo_principalContainer[masterloop*(blockIdx.x*defRow + s2)] == 1)
                        {	
							proceed = TRUE;
						}					
						else
						{
							proceed = FALSE;
						}
			}
}

if(proceed)
{
// No play check
	checkNoPlay = TRUE;
	placeHolderA = 0.0;
	{

		for(s1 = 0; s1 < defCol; s1++)		// Checks Negative
			{
				if(loc_x[s1] < 0)
					checkNoPlay = FALSE;
			}			
		for(s1 = 0; s1 < defRow; s1++)		// Checks Negative
			{
				if(loc_y[s1] < 0)
					checkNoPlay = FALSE;
			}

		for(s1 = 0; s1 < defRow; s1++)		// Checks Equal to 1; no Pure plays
			{
				if(loc_x[s1] == 1.00)
					checkNoPlay = FALSE;
			}
		
		for(s1 = 0; s1 < defRow; s1++)		// Checks Equal to 1; no Pure plays
			{
				if(loc_y[s1] == 1.00)
					checkNoPlay = FALSE;
			}
		for(s1 = 0; s1 < defRow; s1++)		// Checks Small Probability Plays ( less than 1 in a million)
			{
				if(glo_principalContainer[masterloop*(blockIdx.x*defRow + s1)] == 1 && loc_x[s1] < epsilon)
					checkNoPlay = FALSE;
			}

		for(s1 = 0; s1 < defRow; s1++)		// Checks Small Probability Events ( less than 1 in a million)
			{
				if(glo_principalContainer[masterloop*(threadIdx.x*defRow + s1)] == 1 && loc_y[s1] < epsilon)
					checkNoPlay = FALSE;
			}

		__syncthreads();
}			
if(checkNoPlay)
	{		
		glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)] = (float) 1.00;
	}
else
	{
		glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)] = (float) 0.00;
	}

if(glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)] == (float) 1.00)
{
	for(s1= 0; s1 < defCol; s1++)
    {
        glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)+(s1+1)] = loc_x[s1];
    }	
	for(s1= defRow; s1 < (defCol+defRow); s1++)
	{
		glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)+(s1+1)] = loc_y[(s1 - defRow)];
    }
}
__syncthreads();
}	   
}	//masterloop
// ........................................ CALCULATING MIXED NASH ............................................/                        
}

cudaError_t pureCalcWrapper(float* __array_payA, float* __array_payB, int __principalCount, int* __principalContainer, FILE* _f)
{
	int s1,s2,s3;
	cudaError_t cudaStatus;

	// Transfer Variables inputted from CPU Wrapper Function 
	float *dev_array_payA = 0;
	float *dev_array_payB = 0;
	int* dev_principalContainer = 0;
	int *dev_returnEquilibrium = 0;	
	int *returnEquilibrium  = (int *) malloc(__principalCount*__principalCount*(defRow + defCol + 1) * sizeof(int));
{
		fprintf(_f,"\n------------------------------------------------");
        fprintf(_f,"\nPayoff A\n");
        fprintf(_f,"------------------------------------------------\n");								
                for(s1 = 0; s1 < defRow; ++s1)
                {
                        for(s2 = 0; s2 < defCol; ++s2)
                        {
                                fprintf(_f,"%f ", __array_payA[s1*defCol + s2]);
                        }
                        fprintf(_f,"\n");
                }
        fprintf(_f,"\n------------------------------------------------");
        fprintf(_f,"\nPayoff B\n");
        fprintf(_f,"------------------------------------------------\n");				

                for(s1 = 0; s1 < defRow; ++s1)
                {
                        for(s2 = 0; s2 < defCol; ++s2)
                        {
                                fprintf(_f,"%f ", __array_payB[s1*defRow + s2]);
                        }
                        fprintf(_f,"\n");
                }
        fprintf(_f,"\n------------------------------------------------");
        fprintf(_f,"\nSupport 1\n");
        fprintf(_f,"------------------------------------------------\n");	
}

////////////////////////////////////////////// Pure Strategy CUDA Malloc  //////////////////////////////////////////////	

	cudaStatus = cudaMalloc((void**)&dev_array_payA, defCol*defRow*sizeof(float));
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payA malloc failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMalloc((void**)&dev_array_payB, defCol*defRow*sizeof(float));
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payB malloc failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMalloc((void**)&dev_principalContainer, __principalCount*defRow*sizeof(int));
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payB malloc failed!");
			system("pause");
			exit(0);
		}
		printf("\n%d\n", __principalCount*__principalCount*(defRow + defCol + 1));
	cudaStatus = cudaMalloc((void**)&dev_returnEquilibrium, __principalCount*__principalCount*(defRow + defCol + 1)*sizeof(int));
	if (cudaStatus != cudaSuccess)
		{
			cudaError_t error = cudaGetLastError();
			printf("pure dev_returnEquilibrium malloc failed!");
			printf("%s\n", cudaGetErrorString(cudaStatus));
			system("pause");
			exit(0);
		}
////////////////////////////////////////////// Pure Strategy CUDA Malloc ///////////////////////////////////////////////

////////////////////////////////////////////// Pure Strategy CUDA Memcopy //////////////////////////////////////////////

	cudaStatus = cudaMemcpy(dev_array_payA,	__array_payA,	defCol*defRow*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payA memcpy failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMemcpy(dev_array_payB,	__array_payB,	defCol*defRow*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMemcpy(dev_principalContainer,	__principalContainer, __principalCount*defRow*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}
	cudaStatus = cudaMemcpy(dev_returnEquilibrium,	returnEquilibrium, __principalCount*__principalCount*(defRow + defCol + 1)*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_array_payB memcpy failed!");
			system("pause");
			exit(0);
		}		
//		exit(0);
////////////////////////////////////////////// Pure Strategy CUDA Execution //////////////////////////////////////////////	
	dim3 dimBlock(__principalCount,1,1);
	dim3 dimGrid(__principalCount,1,1);
	//cudaEventRecord(start, 0);
	pureStrategyCompute<<<dimGrid,dimBlock>>>(dev_array_payA, dev_array_payB, dev_principalContainer, dev_returnEquilibrium);
	//cudaEventRecord(end, 0);
//	cudaEventSynchronize(end);
//	cudaEventElapsedTime(&time, start, end);
	if (cudaStatus != cudaSuccess) 
	{
		printf("pureStrategyCompute failed!");
			system("pause");
			exit(0);
	}

	cudaStatus = cudaMemcpy(returnEquilibrium, dev_returnEquilibrium, (__principalCount*__principalCount)*(defRow + defCol + 1)*sizeof(int), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			printf("pure dev_returnEquilibrium memcpy failed!");
			system("pause");
			exit(0);
		}	
			for(s1 = 0; s1 < (__principalCount*__principalCount); s1++)
					{                
						if(returnEquilibrium[s1*(defCol+defRow+1)] == 1)
							{                                                                                                                                        
								fprintf(_f,"\nP_1:\t(");
									for(s3 = 0; s3 < defRow; s3++)
										fprintf(_f," %d ",returnEquilibrium[s1*(defCol+defRow+1) + s3 + 1 ]);
                                        fprintf(_f,")\nP_2:\t(");
                                    for(s3 = defRow; s3 < (defCol+defRow); s3++)
										fprintf(_f," %d ",returnEquilibrium[s1*(defCol+defRow+1) + s3 + 1 ]);
                                        fprintf(_f,")\n");
                            }
					}
//	fprintf (_f,"\nPure Strategy Execution Time :\t%f ms\n",time);	
	
//	cudaFree(dev_array_payA);
//    cudaFree(dev_array_payB);
//    cudaFree(dev_principalContainer);
  //  cudaFree(dev_returnEquilibrium);

	return cudaStatus;
}

__device__ int return_max(int *a, int size)
{
        int iter;
        int max_hold = INT_MIN;
        int max_index = INT_MIN;
        
        for(iter = 0; iter < size; iter++)
                {        
                        if(a[iter] >= max_hold)
                                {
                                        max_hold = a[iter];
                                        max_index = iter;
                                }
                }
        return max_index;                
}

__global__ void	pureStrategyCompute(float* glo_array_payA, float *glo_array_payB, int *glo_principalContainer, int *glo_returnEquilibrium)
{               
///////////////////////////////////////// LOCAL VARIABLES /////////////////////////////////////////
	int s1,s2,s3;
	int ident_x = 0;
	int ident_y = 0;
	int loc_Ay[defRow];
	int loc_xB[defCol];
///////////////////////////////////////// LOCAL VARIABLES /////////////////////////////////////////
//initAll
for(int init = 0; init < defCol; init++)
{
	loc_Ay[init] = 0.0;
	loc_xB[init] = 0.0;
}
                        {        // Calculate Ay                
                            for(s1 = 0; s1 < defRow; s1++)
								loc_Ay[s1] = 0;
                            
							for(s1 = 0; s1 < defRow; s1++)
                                        {
                                        for(s2 = 0; s2 < defCol; s2++)
                                                {                                        
                                                        loc_Ay[s1] += glo_principalContainer[threadIdx.x*defRow + s2] * glo_array_payA[s1*defCol + s2];			
                                                }
										}

							for(s3 = 0; s3 < defCol; s3++)
                                        {
                                                if(glo_principalContainer[threadIdx.x*defRow + s3] == 1)
                                                {        
                                                        ident_y = s3;
                                                }
                                        }
                        }

                        {        // Calculate xB        
                                for(s1 = 0; s1 < defCol; (s1)++)
                                        loc_xB[s1] = 0;
                                for(s1= 0; s1 < defCol; s1++)
                                        {
                                        for(s2= 0; s2 < defRow; s2++)
                                                {                                        
                                                        loc_xB[s1] += glo_principalContainer[blockIdx.x*defRow + s2] * glo_array_payB[s2*defCol + s1];
                                                }
                                        }                                
                                for(s3 = 0; s3 < defRow; s3++)
                                        {
                                                if(glo_principalContainer[blockIdx.x*defRow + s3] == 1)
                                                {
                                                        ident_x = s3;
                                                }
                                        }
                        }          

				if((return_max(loc_Ay,defRow)*defCol + ident_y) == (ident_x*defCol + return_max(loc_xB,defCol)))
						{
							glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)] = 1;
						}
                else
						{
							glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)] = 0;
                        }	

				for(s1= 0; s1 < defRow; s1++)
                        {
                                glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)+(s1+1)] = glo_principalContainer[threadIdx.x*defRow + s1];
                        }
                for(s1= defRow; s1 < (defCol+defRow); s1++)
                        {
                                glo_returnEquilibrium[(blockIdx.x*defRow + threadIdx.x)*(defCol+defRow+1)+(s1+1)] = glo_principalContainer[(blockIdx.x*defRow + s1 - defRow)];
                        }                                                
}

__host__ __device__ bool next(int n_, int *p_)
{
        const int n1 = n_ - 1;
        int i = n1;
        int s;
        int r; 
        int j;
        do { --i; } while ( p_[i] >= p_[i+1] );
        if ( (long)i<0 ) 
                return FALSE;
        j = n1;
        while ( p_[i] >= p_[j] ) { --j; }
        swap(p_[i], p_[j]);
        r = n1;
        s = i + 1;
        while ( r > s ) { swap(p_[r], p_[s]); --r; ++s; }
                return TRUE;
}

__device__ bool validityCheck(float *a,int size)
{
	int s4;
	for(s4 = 0; s4 < size; s4++)
		{
			if(a[s4] < float (0) || a[s4] > float (1))
				return FALSE;
		}		
	return TRUE;
}

__device__ void bcksub_array(float *a, int n, int *indx, float *b)
{
	int i,j,ip,ii=0;
	float sum;
 
	for(i=0;i<n;i++) 
	{
		ip = indx[i];
		sum	= b[ip];
		b[ip] = b[i];
		if(ii>=0) 
			for(j = ii; j < i; j++) 
				sum -= a[i*n + j]*b[j];
		else if(sum) 
				ii=i;	/* a nonzero element encounted */
		b[i]=sum;
	}
 
	for(i=n-1;i>=0;i--) 
	{
		sum = b[i];
		for(j = i + 1; j < n; j++) 
			sum-=a[i*n + j]*b[j];
		b[i]=sum/a[i*n + i];
	}
}

__device__ void ludcmp_array(float *a, int n, int *indx, float *d)
{
	int i, imax, j, k;
	float big, dum, sum, temp;
	float vv[defCol];

	*d = 1.0;
	for(i = 0; i < n; i++)
	{
		big = 0.0;
		for(j = 0; j < n; j++)
		{
			if((temp = fabs(a[i*n + j])) > big) 
				big = temp;
		}
		if(big == 0.0) 
		{
			printf("Singular matrix in routine ludcmp.");
			exit(1);
		}
		vv[i] = 1.0 / big;	//vv correct for both
	}

	for(j = 0; j < n; j++)
	{
		for(i = 0; i < j; i++)
		{
			sum = a[i*n + j];
			for(k = 0; k < i; k++)
			{
				sum -= a[i*n + k] * a[k*n + j];
			}
				a[i*n + j] = sum;								
		}
		
		big = 0.0;
		
		for(i = j; i < n; i++)
		{
			sum = a[i*n + j];
			for(k = 0; k < j; k++) 
				sum -= a[i*n + k] * a[k*n + j];
			a[i*n + j] = sum;


			if((dum = vv[i] * fabs(sum)) >= big)
			{	
				big = dum;
				imax = i;
			}
		}		
		if(j != imax)
		{	
			for(k = 0; k < n; k++)
			{
				dum = a[imax*n + k];
				a[imax*n + k] = a[j*n + k];
				a[j*n + k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}

		indx[j] = imax;
		if(a[j*n + j] == 0.0) 
			a[j*n + j] = epsilon;

		if(j != n)
		{
			dum = 1.0 / (a[j*n + j]);

			for(i = j + 1; i < n; i++) 
				a[i*n + j] *= dum;
		}
	}
}
