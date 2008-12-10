#!/bin/sh
#PBS -l nodes=4:ppn=2
#PBS -l walltime=50:00
#PBS -j oe
#PBS

numprocs=8
inp=p0.25M
maxD=29
solveU=1
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=0
mgLoadFac=1.5

prefixPath=/opt/openmpi/1.2.4

cd ${HOME}/otk

mpirun -prefix $prefixPath -nooversubscribe -np ${numprocs} ./elasticitySolver ${inp} ${maxD} ${solveU} ${writeB} ${dim} ${maxNum} ${incCor} ${compressLut} ${mgLoadFac} >& ${HOME}/otk/nvlm.${inp}.${numprocs}.txt

status=$?
exit ${status}

