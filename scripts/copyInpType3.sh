#!/bin/bash
for i in $(find . -type f -name '*\.inp')
do
    src=$i
    tgtPrefix=$(echo $i | sed -e "s/.inp//")
    for ((j=0; j < $1; j++))
    do	
        tgt=${tgtPrefix}_$j.inp
        cp $src $tgt
    done
done 
