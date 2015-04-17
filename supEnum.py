#!/usr/bin/python

import sys
import argparse

try:
    import numpy
except:
    print "numpy is not installed"

print "NumPy is installed: verion " + numpy.__version__


# def kSubsetsHelper(k, items, acc, index):
#     # (n-1 c k) + (n - 1 c k - 1) 
#     if (len(items) - index) < k:
#         return
#     elif k == 0: 
#         print acc
#     else:
#         kSubsetsHelper(k, items, acc, index + 1) # do not include item at index
#         acc.append(items[index])
#         kSubsetsHelper(k - 1, items, acc, index + 1) # include item at index
#         acc.pop()
# 
# def kSubsets(items, k):
#     kSubsetsHelper(k, items, [], 0) 
#

numPairs = 0 

def kSubsetsHelper(k, kCur, items1, items2, acc1, acc2, index, proc, f):
    if proc:
        (acc, items) = acc2, items2
    else:
        acc, items = acc1, items1

    if (len(items) - index) < kCur:
        return
    elif kCur == 0:
        if proc:
             f (acc1, acc2)
        else:
            kSubsetsHelper(k, k, items1, items2, acc1, [], 0, True, f)
    else:
        kSubsetsHelper(k, kCur, items1, items2, acc1, acc2, index + 1, proc, f)
        acc.append(items[index])
        kSubsetsHelper(k, kCur - 1, items1, items2, acc1, acc2, index + 1, proc, f)
        acc.pop()

def kSubsets(items1, items2, k, f):
    kSubsetsHelper(k, k, items1, items2, [], [], 0, False, f) 

def printPair (set1, set2):
    global numPairs
    numPairs += 1
    print((set1,set2))

def parseArgs():
  parser = argparse.ArgumentParser()
  parser.add_argument("gameFile", type=file)
  args = parser.parse_args()
  parseFile(args.gameFile)

def parseFile (fObj):
  for line in fObj:
    line = line.strip()
    if not line:
      continue
    else:
      line = line.split()
      if line[0][0].isdigit():
	continue
      elif line[0] == "Players:":
	print "players = " + line[1]
      elif line[0] == "Actions:":
	print "actions = " + line[1] + ", " + line[2]
      else:
	continue

def main ():
    parseArgs()
    items1 = ['a', 'b', 'c', 'd']
    items2 = ['e', 'f', 'g', 'h']
    global numPairs
    # for i in range(1,5):
    #     numPairs = 0
    #     kSubsets(items1, items2, i, printPair)
    #     print (numPairs)

main()

# I would have to maintain 2 

