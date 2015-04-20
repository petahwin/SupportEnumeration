#!/usr/bin/python

import sys
import argparse
import numpy

payOffTableA, payOffTableB = [], []
nAction1, nAction2 = 1, 1

def buildMat(ac1, ac2, payoffT, rowOrCol):
  row, col = True, False
  mat = [] 
  # Build Ax = b equation to solve
  matSize = 1 + len (ac1)
  # Build equations one by one
  if rowOrCol == col: # for eq 3
    for j in ac2:
      equationAcc = []
      for i in ac1:
        equationAcc.append(payoffT[i][j])
      equationAcc.append(1)
      mat.append(equationAcc)
  else:
    for i in ac1: # for eq 5
      equationAcc = []
      for j in ac2:
        equationAcc.append(payoffT[i][j])
      equationAcc.append(1)
      mat.append(equationAcc)

  # Row to ensure probabilities sum to 1
  lastRow = [1 for x in range(matSize)]
  lastRow[matSize-1] = 0
  mat.append(lastRow)
  vectorB = [0 for x in range(matSize)]
  vectorB[matSize-1] = 1
  return numpy.array(mat),numpy.array(vectorB)

def buildFullStrat(ac, strat, nActions):
  stratF = [0 for x in range(nActions)] 
  i = 0
  for x in ac:
    stratF[x] = strat[i]
    i += 1
  return stratF

def nashEq(ac1, ac2):
  matA, vecB = buildMat(ac1,ac2,payOffTableB, False)
  vecX = numpy.delete(numpy.linalg.lstsq(matA, vecB)[0], -1)
  for x in vecX:
    if x < 0.0: 
      # print "negative elmt in X; next iteration"
      return

  matA, vecB = buildMat(ac1,ac2,payOffTableA, True)
  vecY = numpy.delete(numpy.linalg.lstsq(matA, vecB)[0], -1)
  for y in vecY:
    if y < 0.0: 
      # print "negative elmt in Y; next iteration"
      return

  strat1 = buildFullStrat(ac1, vecX, nAction1)
  strat2 = buildFullStrat(ac2, vecY, nAction2)
  Ay = numpy.dot(payOffTableA, strat2)
  xB = numpy.dot(strat1, payOffTableB)
  maxAy, maxxB = reduce(max, Ay), reduce(max, xB)
  
  # check for best response
  for x in ac1:
    if Ay[x] < (maxAy - abs(maxAy / 1000.)): 
      # print "MaxAy = {}".format(maxAy)
      # print "Ay[x] = {}".format(Ay[x])
      # print "not best response in Ay; next iteration"
      return
  for y in ac2:
    if xB[y] < (maxxB - abs(maxxB / 1000.)): 
      # print "MaxxB = {}".format(maxxB)
      # print "xB[y] = {}".format(xB[y])
      # print "not best response in xB; next iteration"
      return
  print "SOLUTION: {}".format((strat1, strat2))

def parseArgs():
  parser = argparse.ArgumentParser()
  parser.add_argument("gameFile", type=file)
  args = parser.parse_args()
  parseFile(args.gameFile)

def parseFile (fObj):
  global payOffTableA, payOffTableB
  global nAction1, nAction2

  for line in fObj:
    line = line.strip()
    if not line:
      continue
    else:
      line = line.split()
      if line[0] == "Actions:":
	nAction1, nAction2 = int(line[1]), int(line[2])
	payOffTableA = [[] for x in range(nAction1)]
	payOffTableB = [[] for x in range(nAction1)]
      elif line[0][0].isdigit() or line[0][0] == '-':
	i = 0
	for item in range(0, len(line), 2):
	  index = i % nAction1
	  payOffTableA[index].append(float(line[item]))
	  payOffTableB[index].append(float(line[item+1]))
	  i+=1
      else:
	continue

def kSubsetsHelper(k, kCur, items1, items2, acc1, acc2, index, proc, f):
  if proc:
    acc, items = acc2, items2
  else:
    acc, items = acc1, items1
  if (len (items) - index) < kCur:
    return
  elif kCur == 0:
    if proc:
      # print "testing {}".format((acc1,acc2))
      f (acc1, acc2)
    else:
      kSubsetsHelper(k,k, items1, items2, acc1, [], 0, True, f)
  else:
    kSubsetsHelper(k,kCur, items1, items2, acc1, acc2, index+1, proc,f)
    acc.append(items[index])
    kSubsetsHelper(k,kCur-1,items1,items2,acc1,acc2,index+1,proc,f)
    acc.pop()

def kSubsets(items1, items2, k, f):
  kSubsetsHelper(k,k,items1,items2,[],[],0,False,f)

def main():
  parseArgs()
  items1 = range(nAction1)
  items2 = range(nAction2)
  for i in range(1, min(nAction1,nAction2)+1):
    kSubsets(items1,items2,i,nashEq)

main()

