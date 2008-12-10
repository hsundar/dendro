
inp=$2
writeO=1

mpirun -np $1 ./justCoarsen $inp $writeO   

