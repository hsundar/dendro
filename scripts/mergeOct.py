#! /usr/bin/env python

import string
from sys import argv, exit

## Usage is intended as ...
##
##  merge prefix npes saveFile

if len(argv) < 4:
  print "Usage: mergeOct prefix npes saveFile [suffix] "
  exit(0);

p = int(argv[2]);

ntotal = 0;
## first compute the total number of nodes.

suff = "";

if len(argv) > 4:
   suff = argv[4]

for i in range(p):
  fname = argv[1] + repr(i)+ '_' +  repr(p) +  suff + ".ot"
  fin = open (fname, "r")
  dim_dep = fin.readline()
  if dim_dep == "":
     n = 0
  else:
     n = int(fin.readline())
  ntotal += n
  fin.close();

## write out the header ...
fout = open(argv[3], 'w')
fout.write(dim_dep)
fout.write(str(ntotal)+'\n')
## write out the nodes ...
for i in range(p):
  fname = argv[1] + repr(i)+ '_' +  repr(p) +  suff + ".ot"
  fin = open (fname, "r")
  dim_dep = fin.readline()
  if dim_dep == "":
     n = 0
  else:
     n = int(fin.readline())
  for j in range(n):
    tstr = fin.readline()
    fout.write(tstr);
  fin.close()
fout.close()

