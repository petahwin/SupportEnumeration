#ifndef INCLUDED_VERTEXNODE_H
#include "vertexNode.h"
#endif

#ifndef INCLUDED_QUEUE_H
#include "queue.h"
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"

#define INF -1

double wtime0, wtime1, cptime;

// Returns 0 on error, 1 on success
// In param: graph/source data files
// Out param: adjList (populated adjacency list), NAdjList (no. vertices)
// Out param: sourceArr (vertex ids of sources), NsourceArr(no. sources)
int readGraphSource(char * graphData, char * sourceData, 
    vertexNode *** adjList, int * NAdjList, int ** sourceArr, int * NSourceArr);

// implement Moore's alg here?
double mooreShortestPaths(vertexNode ** adjList, int N, int src) {
    double qtime0, qtime1, qtime;

    timing(&wtime0, &cptime);
    int vi, vj, tid, dist[N+1], numThreads;
    char isDoneMask[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    
    omp_lock_t distLock[N+1];
    
    // Need the locks as well

    Queue * q = makeQueue();
    vertexNode * vjPtr;
    // Locks for ALL of the vertices

    timing(&wtime1, &cptime);

    printf("time to init vars %9.4f\n", wtime1 - wtime0);

    timing(&wtime0, &cptime);
    #pragma omp parallel private(vi, vj, vjPtr, tid, qtime1, qtime0, qtime) \
    shared(distLock, q, dist, isDoneMask, adjList, N, src, numThreads) 
    {
        Queue * locInQ = makeQueue(), * locOutQ = makeQueue();

        qtime = 0.;
        tid = omp_get_thread_num();
        omp_init_lock(&(distLock[0]));
        #pragma omp for
        for (int i = 1; i <= N; ++i) {
            dist[i] = INF;
            omp_init_lock(&(distLock[i]));
        }

        #pragma omp master
        {
            numThreads = omp_get_num_threads();
            for (int i = 0; i < numThreads; ++i) {
                isDoneMask[i] = 1;
            }
            dist[src] = 0;
            timing(&qtime0, &cptime);
            enqueue(q, src);
            timing(&qtime1, &cptime);
            timing(&wtime1, &cptime);
            printf("time to assign vars %9.4f\n", wtime1 - wtime0);
        
            timing(&wtime0, &cptime);
            printf("Starting with %d threads, src=%d\n", numThreads, src);
        }

        #pragma omp barrier
        
        int numNodes = 0;
        while (true) {
            bool emptyQueue = false, kill = false;
            
            if (size(q) == 0) {
                isDoneMask[tid] = 1;
                bool allDone = true;
                for(int i = 0; i < numThreads; ++i) {
                    if (isDoneMask[i] == 0) allDone = false;
                }
                if (allDone) {
                    kill = true;
                } else {
                    emptyQueue = true;
                }

            } else {
                #pragma omp critical (queueOp)
                {
                    if (size(q) == 0) {
                        emptyQueue = true;
                    } else {
                        isDoneMask[tid] = 0;
                        if (size(q) <= numThreads) {
                            numNodes++;
                            enqueue(locInQ, dequeue(q));
                        } else {
                            int share = size(q) / numThreads;
                            numNodes += share;
                            for(int i = 0; i < share; ++i) {
                                enqueue(locInQ, dequeue(q));
                            }
                        }
                    }
                }
            }

            if (kill) break;
            if (emptyQueue) continue;
            
            while (vi = dequeue(locInQ)) { 
            
            vjPtr = adjList[vi]; // get first edge out of vi
            
            // Loop over edges out of vi
            while (vjPtr) {
                vj = vjPtr->vertex; // get vertex no.

                // Relax distance and add to queue if pass
                if (dist[vj] == INF || dist[vi] + vjPtr->weight < dist[vj]) {
                    // do the check again in critical section for given vj
                    // Use lock for dist[vj]
                    omp_set_lock(&(distLock[vj]));
                    {
                        if (dist[vj] == INF || dist[vi] + vjPtr->weight < dist[vj]) {
                            dist[vj] = dist[vi] + vjPtr->weight;
                            enqueue(locOutQ, vj);
                        }
                    }
                    omp_unset_lock(&(distLock[vj]));
                }
                vjPtr = vjPtr->next;

            } 
            }
            
            #pragma omp critical (queueOp)
            {
                int sizeQ = size(locOutQ);
                for(int i = 0; i < sizeQ; ++i) {
                    enqueue(q, dequeue(locOutQ));
                }
            }
        }
          
        #pragma omp for
        for (int i = 0; i <=N; ++i) {
            omp_destroy_lock(&(distLock[i]));
        }

        freeQueue(locInQ); freeQueue(locOutQ);

        #pragma omp master
        {
            timing(&wtime1, &cptime);
            printf("\nFor source %d:\n", src);
            printf("\tElapsed time: %9.4f seconds\n", wtime1 - wtime0);   
            printf("\tDest Vertex\tDistance\n");
            printf("\t___________\t________\n");
            for (int i = 1, j = (N < 10001)? N : 10001; i <= j; i += 1000) {
                printf("\t%8d\t%8d\n", i, dist[i]);
            }
            printf("\t%8d\t%8d\n", N, dist[N]);
        }
        
        #pragma omp barrier
        printf("tid %d: processed %d nodes\n", tid, numNodes);
    }

    freeQueue(q);

    return wtime1 - wtime0;
}

int main(int argc, char ** argv) {
    vertexNode ** adjList;
    int * sources;
    int NadjList, Nsources;
    double totalTime = 0.;
    // Read in all file data, store in adjacency list
    //
    timing(&wtime0, &cptime);
    if (argc != 3) {
        fprintf(stderr, "Usage: ./serial graphFile sourceFile\n");
        return 1;
    } else if (!readGraphSource(argv[1], argv[2], &adjList, &NadjList,
        &sources, &Nsources)) 
    {
        fprintf(stderr, "Error reading files\n");
        return 1;
    }

    timing(&wtime1, &cptime);
    totalTime += wtime1 - wtime0;

    printf("Elapsed time for init and I/O:\t%9.4f seconds\n", totalTime);

    for (int i = 0; i < Nsources; ++i) {
        totalTime += mooreShortestPaths(adjList, NadjList, sources[i]);
    }

    printf("\nTotal Elapsed time: %9.4f seconds\n", totalTime);

    freeAdjList(adjList, NadjList);
    free(sources);
    
    return 0;
}
