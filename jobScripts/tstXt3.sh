#!/bin/csh
#PBS -l size=1:1
#PBS -l walltime=5:00
#PBS -j eo
#PBS -r n
#PBS -m bea
#PBS -M rahul.sampath@gmail.com

unalias cp
unalias rm

set echo

set inp = p1M
set cmplx = 1
set incInt = 1
set runOpt = 2
set readP = 1
set readO = 0
set dim = 3
set maxD = 30
set writeP = 0
set writeO = 0
set maxNum = 1
set incCor = 1

cd ${SCRATCH}

pbsyod -size ${PBS_NPROCS} ./runScal $inp $cmplx $incInt $runOpt $readP $readO $maxD $writeP $writeO $dim $maxNum $incCor  >& ./${inp}.${PBS_NPROCS}.${maxD}.txt

mv ./${inp}.${PBS_NPROCS}.${maxD}.txt ${HOME}/${inp}.${PBS_NPROCS}.${maxD}.txt

