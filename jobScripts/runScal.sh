#!/bin/bash

export inp=gt 
export complx=0
export incInt=0
export runOpt=2
export readPts=1
export readOct=0
export maxD=30
export writeP=1
export writeB=1
export dim=3
export maxNumPts=1
export incCorner=1

mpirun -np $1 ./runScal ${inp} ${complx} ${incInt} ${runOpt} ${readPts} ${readOct} ${maxD} ${writeP} ${writeB} ${dim} ${maxNumPts} ${incCorner} >& con.${inp}.txt


