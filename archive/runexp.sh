for k in 0.00 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 ;
do
    for j in `seq 0 24`; do ./gammadist.py $(( j * 12 )) 12 7 $k | grep -v 'calc' & done;
    #for i in `seq 0 9`; do ./gammadist.py $(( i * 30 )) 30 11 $k | grep -v 'calc' & done; 
    wait `ps -e | grep python | awk '{print $1}'`
    cp -r ../GeneData/MS ../GeneData/MS.20160801.C${k/./}
    rm -r ../GeneData/MS
    mkdir ../GeneData/MS
done
