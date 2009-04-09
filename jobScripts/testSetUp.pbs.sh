#!/bin/sh
#PBS -l nodes=2:ppn=2
#PBS -l walltime=20:00
#PBS -j oe
#PBS

numprocs=0
tmpfile=/tmp/${PBS_JOBCOOKIE}.dat
echo ${tmpfile}
echo ${PBS_NODEFILE}

for s in `sort < ${PBS_NODEFILE} | uniq`
do
echo $s
echo $s slots=2 >> ${tmpfile}
numprocs=$((numprocs + 2))
done

inp=rg4M
maxD=29
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=1
mgLoadFac=1.5

prefixPath=/ronaldo/opt/openmpi/1.1.1

cd ${DENDRO_DIR}

mpirun -prefix $prefixPath -hostfile ${tmpfile} -nooversubscribe -np ${numprocs} ./tstMgSetUp ${inp} ${maxD} ${writeB} ${dim} ${maxNum} ${incCor} ${compressLut} ${mgLoadFac} >& setUpTest.${inp}.${numprocs}.txt

status=$?
rm -f ${tmpfile}
exit ${status}

