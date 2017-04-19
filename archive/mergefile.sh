#!/bin/bash

# for file in `ls -d */`; do cd $file; rm *.ms ; mv ./bak/*.ms ./; cd .. ; done
# for file in `ls -d */`; do cd $file; rm *_???.ms ; cd .. ; done

for i in `for file in \`ls *.ms\`; do echo ${file%_*.*.ms}; done | uniq`; do mkdir $i ; done
for file in `ls *.ms`; do mv $file ./${file%_*.*.ms}; done

# first call the two command above, then call this file
# and can get ms and mt ft vt file

for directory in `ls -d */`;
do
    cd $directory
    mkdir bak
    for file in `ls *.ms`;
    do
        cat $file >> ${file%.*.ms}.ms
        mv $file ./bak/
    done
    wc *.ms | awk '$4 != "total" && $1 != 292 {print "ERROR " $1 " " $4}'
    for msfile in `ls *.ms`;
    do
        ~/GeneFinding/bioinformatic/archive/parse.py $msfile
    done
    wc *.?t | awk '$4 != "total" && $1 != 92  {print "ERROR " $0 }'
    cd ..
done

