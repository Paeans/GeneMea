#!/bin/bash
for l in 11 ;
do
    for k in 0.60 ;
    do
        echo `date +%H:%M:%S` "start" $l "on cutoff" $k
        for j in `seq 0 24`; do ./gammadist.py $(( j * 12 )) 12 $l $k | grep -v 'calc' & done;
        #for i in `seq 0 9`; do ./gammadist.py $(( i * 30 )) 30 11 $k | grep -v 'calc' & done; 
        wait `ps -e | grep python | awk '{print $1}'`
        cp -r ../GeneData/MS ../GeneData/MS.${l}.20161222.C${k/./}
        rm -r ../GeneData/MS
        mkdir ../GeneData/MS
        echo `date +%H:%M:%S` "end" $l "on cutoff" $k
    done
done
