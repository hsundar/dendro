#!/bin/csh
#PBS -l size=1
#PBS -l walltime=15:00
#PBS -j eo
#PBS -r n
#PBS -m bea
#PBS -M hsundar@gmail.com,rahul.sampath@gmail.com

unalias cp
unalias rm

set echo

set inp = p1M
set dim = 3
set maxD = 29
set solveU = 0
set writeB = 0
set maxNum = 1
set incCor = 1
set numLoops = 5
set compressLut = 1

cd ${SCRATCH}

pbsyod -size ${PBS_NPROCS} ./tstMatVec $inp $maxD $solveU $writeB $dim $maxNum $incCor $numLoops ${compressLut} >& ./${inp}.${PBS_NPROCS}.${compressLut}.txt
 
mv ./${inp}.${PBS_NPROCS}.${compressLut}.txt ${HOME}/otk/Results/${inp}.${PBS_NPROCS}.${compressLut}.txt

