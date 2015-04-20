CC = icc
CFLAGS = -g -O3 -openmp -xHost -fno-alias -std=c99

serial:	serial.o vertexNode.o readData.o fastqueue.o\
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

parallel: parallel.o vertexNode.o readData.o fastqueue.o\
/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o serial parallel
