#!/bin/sh
#PBS -l nodes=2:ppn=2
#PBS -l walltime=50:00
#PBS -j oe
#PBS

numprocs=3
inp=p1M
maxD=30
solveU=1
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=0
mgLoadFac=1.5

prefixPath=/opt/openmpi/1.2.4

cd ${DENDRO_DIR}

mpirun -prefix $prefixPath -nooversubscribe -np ${numprocs} ./tstMg ${inp} ${maxD} ${solveU} ${writeB} ${dim} ${maxNum} ${incCor} ${compressLut} ${mgLoadFac}  >& tstMg.${inp}.${numprocs}.txt

status=$?
exit ${status}

