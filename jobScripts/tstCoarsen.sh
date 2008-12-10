inp=$2
cmplx=1
incInt=1
runOpt=2
readP=1
readO=0
dim=3
maxD=29
writeP=0
writeO=1
maxNum=1
incCor=1

mpirun -np $1 ./chkCoarsen $inp $cmplx $incInt $runOpt $readP $readO $maxD $writeP $writeO $dim $maxNum $incCor  

