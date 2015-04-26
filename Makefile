CUDAPATH = /home/apps/fas/GPU/cuda_6.5.14
MAGMAPATH = /home/apps/fas/GPU/MAGMA/magma-1.6.1/

CC = icc
NVCC = $(CUDAPATH)/bin/nvcc
CFLAGS = -mkl -g -O3 -xHost -fno-alias -std=c99 -openmp #-profile-functions 

MCFLAGS = -fPIC -openmp -L$(MAGMAPATH)/lib -lmagma -lcublas -lcudart -I$(MAGMAPATH)/include

# NVCC Compiler-specific flags (phitest supports only 2.x capability)
# override NVCCFLAGS += -I$(CUDAPATH)/include -O3
NVCCFLAGS = -I$(CUDAPATH)/include -O3

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

SerialsupEnum: SerialsupEnum.o readData.o\
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

OMPsupEnum: OMPsupEnum.o readData.o\
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

MagmasupEnum: MagmasupEnum.o
	$(CC) -o $@ $(MCFLAGS) $^

GPUpairs: GPUpairs.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

MagmasupEnum.o : MagmasupEnum.c
	$(CC) -o $@ $(MCFLAGS) -c $<

%.o : %.cu
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o serial GPUpairs SerialsupEnum OMPsupEnum

