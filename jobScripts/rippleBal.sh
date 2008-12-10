
inp=$2
writeO=0
incCor=1
chkBail=0
rePart=0
rippleType=2

mpirun -np $1 ./rippleBal $inp $writeO $incCor $chkBail $rePart $rippleType

