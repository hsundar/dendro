#!/bin/bash
for i in $(find . -type f -name '*\.cpp')
do
    src=$i
    tgt=$(echo $i | sed -e "s/cpp/C/")
    mv $src $tgt
done 
