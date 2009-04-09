#!/bin/sh
#PBS -l nodes=4:ppn=2
#PBS -l walltime=20:00
#PBS -j oe
#PBS

numprocs=8
regLev=7
maxD=29
dim=3
maxNum=1
incCor=1
compressLut=0
mgLoadFac=1.5

prefixPath=/opt/openmpi/1.2.4

cd ${DENDRO_DIR}

mpirun -prefix $prefixPath -nooversubscribe -np ${numprocs} ./tstMgReg $regLev $maxD ${dim} ${incCor} ${compressLut} ${mgLoadFac} >& tstMgRg.${regLev}.${numprocs}.txt

status=$?
exit ${status}

