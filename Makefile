CUDAPATH = /home/apps/fas/GPU/cuda_6.0.37

CC = icc
NVCC = $(CUDAPATH)/bin/nvcc
CFLAGS = -mkl=sequential -g -O3 -xHost -fno-alias -Wall -Werror -pedantic -std=c99 -openmp #-profile-functions 

# NVCC Compiler-specific flags (phitest supports only 2.x capability)
NVCCFLAGS = -I$(CUDAPATH)/include -O3 -rdc=true 

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

all: SerialsupEnum OMPsupEnum GPUsupEnum

SerialsupEnum: SerialsupEnum.o readData.o subsets.o printUtils.o \
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

OMPsupEnum: OMPsupEnum.o readData.o subsets.o printUtils.o \
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

MagmasupEnum: MagmasupEnum.o
	$(CC) -o $@ $(MCFLAGS) $^

GPUsupEnum: GPUsupEnum.o 
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

%.o : %.cu
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) -c $< -lm

clean:
	rm -f *.o SerialsupEnum OMPsupEnum GPUsupEnum 

