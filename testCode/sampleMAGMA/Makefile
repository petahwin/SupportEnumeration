MAGMAPATH = /home/apps/fas/GPU/MAGMA/magma-1.6.1/

CC = icc
NVCC = $(CUDAPATH)/bin/nvcc

MCFLAGS = -DADD_ -DMAGMA_SETAFFINITY -DHAVE_CUBLAS -DMIN_CUDA_ARCH=200\
 -fPIC -openmp -L$(MAGMAPATH)/lib -lmagma -lcublas -lcudart -I$(MAGMAPATH)/include -std=c99

magmatest: magmatest.o
	$(CC) -o $@ $(MCFLAGS) $^

%.o : %.c
	$(CC) $(MCFLAGS) -c $<

clean:
	rm -f *.o magmatest

