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

inp=p1M
solveU=1
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=1
mgLoadFac=1.2
consistentRhs=1
maxD=29
noRandom=1

prefixPath=/ronaldo/opt/openmpi/1.1.1

cd ${DENDRO_DIR}

mpirun -prefix $prefixPath -hostfile ${tmpfile} -nooversubscribe -np ${numprocs} ./tstMinv ${inp} 1 ${noRandom} ${maxD} 1 0 3 1 1 1 1 1.2 20 >& Minv.${inp}.${maxD}.${noRandom}.txt

status=$?
rm -f ${tmpfile}
exit ${status}

