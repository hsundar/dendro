#!/bin/bash
for i in $(find . -type f -name '*_*\.inp')
do
    rm -f $i    
done 
