#!/bin/bash
for file in $(find . -type f -name '*\.C')
do
 cat $file |sed "s/std::cout/cout/g" > $file.1
 cat $file.1 |sed "s/std::cin/cin/g" > $file.2
 cat $file.2 |sed "s/std::endl/endl/g" > $file.3
 cat $file.3 |sed "s/std::ifstream/ifstream/g" > $file.4
 cat $file.4 |sed "s/std::ofstream/ofstream/g" > $file.5
 cat $file.5 |sed "s/std::cerr/cerr/g" > $file.6
 mv $file.6 $file
 rm $file.?
done 
