
inp=$2
writeO=0
incCor=1

mpirun -np $1 ./justBal $inp $writeO $incCor  

