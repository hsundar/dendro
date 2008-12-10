#!/bin/bash
for i in $(find . -type f -name '*\.inp')
do
    src=$i
    tgtPrefix=$(echo $i | sed -e "s/.inp//")
    for ((j=0, k=0; j < $1; j = j+1000, k++))
    do	
        tgt=${tgtPrefix}_$k.inp
        cp $src $tgt
    done
done 
