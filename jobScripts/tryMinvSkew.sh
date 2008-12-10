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

consistentRhs=1
maxD=29
dim=3
maxNum=1
incCor=1
compressLut=1
mgLoadFac=1.2
regLev=6
skew=0
noRandom=0


prefixPath=/ronaldo/opt/openmpi/1.1.1

cd ${HOME}/otk

mpirun -prefix $prefixPath -hostfile ${tmpfile} -nooversubscribe -np ${numprocs} ./tstMinvSkew $regLev ${skew} $consistentRhs  ${noRandom} $maxD ${dim} ${incCor} ${compressLut} ${mgLoadFac} >& MinvSkew.${regLev}.${skew}.${noRandom}.txt

status=$?
rm -f ${tmpfile}
exit ${status}

