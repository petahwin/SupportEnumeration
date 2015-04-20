#!/usr/bin/python

import sys
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("gameFile", type=file)
  args = parser.parse_args()
  fOut = open(args.gameFile.name + ".easy", 'w')
  for line in args.gameFile:
    line = line.strip()
    if not line:
      continue
    else:
      line = line.split()
      if line[0] == "Actions:":
	fOut.write("player1: " + line[1]+"\n")
	fOut.write("player2: " + line[2]+"\n")
      elif line[0][0].isdigit() or line[0][0] == '-':
	for item in range(0, len(line), 2):
	  fOut.write(line[item] + " " + line[item+1] + "\n")
      else:
	continue

main()

