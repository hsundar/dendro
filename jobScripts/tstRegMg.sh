regLev=$2
maxD=$3
dim=3
maxNum=1
incCor=1
compressLut=0
mgLoadFac=1.2

mpirun -np $1 ./tstMgReg $regLev $maxD ${dim} ${incCor} ${compressLut} ${mgLoadFac} 


