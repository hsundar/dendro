regLev=$2
dim=3
maxD=$3
compressLut=1
mgLoadFac=$4

mpirun -np $1 ./tstRPreg ${regLev} ${dim} ${maxD} ${compressLut} ${mgLoadFac} 

