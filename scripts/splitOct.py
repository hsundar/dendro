#! /usr/bin/env python

import string
from sys import argv, exit

## Usage is intended as ...
##
##  split octFile numProcs splitPrefix

if len(argv) < 2:
  print "Usage: splitOct octFile npes splitPrefix"
  exit(0);

fin = open(argv[1], "r")
dim_dep = fin.readline()
n = int(fin.readline())
p = int(argv[2]);
nlocal = n/p;

for i in range(p):
  ## filename
  fname = argv[3] + repr(i) + '_' +  repr(p) + ".ot"
  ## open file
  fout = open(fname, "w")
  ## write out the header info ...
  fout.write(dim_dep)
  if i == p-1:
    nlocal = n - nlocal*i
  fout.write(str(nlocal)+"\n")
  for j in range(nlocal):
    tstr = fin.readline()
    fout.write(tstr)
  fout.close()

fin.close()
