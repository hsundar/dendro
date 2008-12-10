inp=$2
maxD=$3
dim=3
incCor=1
compressLut=1
mgLoadFac=$4

mpirun -np $1 ./tstSetUp2 ${inp} ${maxD} ${dim} ${incCor} ${compressLut} ${mgLoadFac} 


