CUDAPATH = /home/apps/fas/GPU/cuda_6.0.37
CC = icc
NVCC = $(CUDAPATH)/bin/nvcc
CFLAGS = -mkl -debug -g -O3 -xHost -fno-alias -std=c99 #-profile-functions 

# NVCC Compiler-specific flags (phitest supports only 2.x capability)
override NVCCFLAGS += -I$(CUDAPATH)/include -O3
LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

serial: serial.o readData.o\
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

GPUpairs: GPUpairs.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

%.o : %.cu
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o serial GPUpairs

