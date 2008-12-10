inp=$2
maxD=$3
dim=3
maxNum=1

mpirun -np $1 ./downSample ${inp} ${maxD} ${dim} ${maxNum} 


