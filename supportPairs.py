#!/usr/bin/python

import sys
import argparse
import numpy
from scipy import misc

parser = argparse.ArgumentParser()
parser.add_argument("nActions1", type=int)
parser.add_argument("nActions2", type=int)
args = parser.parse_args()
nActions1, nActions2 = args.nActions1, args.nActions2

maxSupport = min(nActions1, nActions2)

sum = 0
totalMem = 0
maxMem = 0
totalElem = 0
maxElem = 0
sizeInt = 4

for i in range(1, maxSupport+1):
    ac1, ac2 = misc.comb(nActions1, i), misc.comb(nActions2, i) # no. of subsets
    elem1, elem2 = ac1 * i, ac2 * i
    maxElem = max(elem1, elem2, maxElem)
    totalElem += elem1 + elem2
    # mem1, mem2 = elem1 * sizeInt, elem2 * sizeInt
    # maxMem = max(mem1,mem2,maxMem)
    # totalMem += mem1 + mem2
    print "{} c {} = {} |\t{} c {} = {}".format(nActions1, i, ac1, nActions2, i, ac2)
    sum += ac1 * ac2

print "total pairs: {}".format(sum)
print "max mem = {}".format(maxElem * sizeInt)
print "total mem = {}".format(totalElem * sizeInt)
print "max elem = {}".format(maxElem)
print "total elem = {}".format(totalElem)

