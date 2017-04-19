for k in 0.81 0.82 0.83 0.84 0.86 0.87 0.88 0.89 ; 
do
    for i in 1 2 3 7 10 12; do for j in `seq 0 4`; do ./gammadist.py $(( j * 59 )) 59 $i $k | grep -v 'calc' & done; done; 
    for i in `seq 0 9`; do ./gammadist.py $(( i * 30 )) 30 11 $k | grep -v 'calc' & done; 
    wait `ps -e | grep python | awk '{print $1}'`
    cp -r ../GeneData/MS ../GeneData/MS.20160730.C${k/./}
done
