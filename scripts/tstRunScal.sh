inp=p2M
cmplx=1
incInt=1
runOpt=2
readP=1
readO=0
dim=3
maxD=30
writeP=0
writeO=0
maxNum=1
incCor=1

mpirun -np 2 ./runScal $inp $cmplx $incInt $runOpt $readP $readO $maxD $writeP $writeO $dim $maxNum $incCor  

