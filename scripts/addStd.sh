#!/bin/bash
for file in $(find . -type f -name '*\.C')
do
 cat $file |sed "s/cout/std::cout/g" > $file.1
 cat $file.1 |sed "s/cin/std::cin/g" > $file.2
 cat $file.2 |sed "s/endl/std::endl/g" > $file.3
 cat $file.3 |sed "s/ifstream/std::ifstream/g" > $file.4
 cat $file.4 |sed "s/ofstream/std::ofstream/g" > $file.5
 cat $file.5 |sed "s/cerr/std::cerr/g" > $file.6
 mv $file.6 $file
 rm $file.?
done 
