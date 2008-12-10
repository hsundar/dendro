inp=$2
maxD=$3
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=1
mgLoadFac=$4

mpirun -np $1 ./tstMgSetUp ${inp} ${maxD} ${writeB} ${dim} ${maxNum} ${incCor} ${compressLut} ${mgLoadFac} 


