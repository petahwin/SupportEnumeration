#!/usr/bin/python

import sys
import argparse
from scipy import misc



sizeFloat = 4

def memThread(n, k):
    supports = 2*k
    strats = 2*n
    matrix = (k+1)**2
    rhs = k + 1
    # Also some unaccounted for device calls and other counters
    return sizeFloat * (supports + strats + matrix + rhs)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("nActions", type=int)
    args = parser.parse_args()
    n = args.nActions

    maxMem = 0.0
    maxSunkSupportsMem = 0.0
    maxSunkThreadsMem = 0.0

    memPayoffs = 2 * (n**2)
    argmax = 0
    maxThreads = 0.0

    for i in range (1, n + 1):
        numSupports = misc.comb(n, i)
        numThreads = numSupports*512 # assuming that each thread handles one pair
        
        memSupports = numSupports * i * sizeFloat
        memThreads = memThread(n, i) * numThreads
        totalMem = memThreads + memSupports + memPayoffs
        if (totalMem > maxMem): 
            maxMem = totalMem
            maxSunkSupportsMem = memSupports
            maxSunkThreadsMem = memThreads
            argmax = i
            maxThreads = numThreads

    print "For {} choose {}:".format(n, argmax)
    print "maxMem = {}".format(maxMem / 1024.0**3)
    print "maxThreads = {}".format(maxThreads)
    print "sunk Support mem = {}".format(maxSunkSupportsMem)
    print "sunk Threads mem = {}".format(maxSunkThreadsMem / 1024.**3)

main()

