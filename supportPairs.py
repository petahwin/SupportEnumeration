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
for i in range(1, maxSupport+1):
    sum += misc.comb(nActions1,i) * misc.comb(nActions2,i)

print "total pairs: {}".format(sum)

