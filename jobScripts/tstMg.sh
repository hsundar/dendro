inp=$2
maxD=$3
solveU=1
writeB=0
dim=3
maxNum=1
incCor=1
compressLut=0
mgLoadFac=1.2

mpirun -np $1 ./tstMg ${inp} ${maxD} ${solveU} ${writeB} ${dim} ${maxNum} ${incCor} ${compressLut} ${mgLoadFac}


