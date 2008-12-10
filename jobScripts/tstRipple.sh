
inp=$2
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
chkBail=1
rePart=1
rippleType=3

mpirun -np $1 ./tstRipple $inp $cmplx $incInt $runOpt $readP $readO $maxD $writeP $writeO $dim $maxNum $incCor $chkBail $rePart $rippleType

