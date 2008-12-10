#! /usr/bin/env python

import string
from sys import argv, exit

## Usage is intended as ...
##
##  merge prefix npes saveFile

if len(argv) < 5:
  print "Usage: mergeOct prefix npes timeSteps saveFilePrefix"
  exit(0);

p = int(argv[2]);


## first compute the total number of nodes.

timeSteps = int(argv[3]);

   
for k in range(timeSteps):
    ntotal = 0
    for i in range(p):
        fname = argv[1] + '_' + repr(i) + '_' +   repr(k) + ".ot"
        fin = open (fname, "r")
        dim_dep = fin.readline()
        n = int(fin.readline())
        ntotal += n
        fin.close();

## write out the header ...
    fname2 = argv[4] + '_' + repr(k) + ".ot"
    fout = open(fname2, 'w')
    fout.write(dim_dep)
    fout.write(str(ntotal)+'\n')
## write out the nodes ...
    for i in range(p):
        fname = argv[1] + '_' + repr(i)+ '_' + repr(k) + ".ot"
        fin = open (fname, "r")
        dim_dep = fin.readline()
        n = int(fin.readline())
        for j in range(n):
            tstr = fin.readline()
            fout.write(tstr);
        fin.close()
    fout.close()

