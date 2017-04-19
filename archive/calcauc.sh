#!/bin/bash
echo "fileID = fopen('result.txt', 'w');"
for directory in `ls -d */`;
do
    cd $directory;
    for file in `ls *.ms *.mt`;
    do
        
        filename=${file%.*}
        echo "test = ["
        #cat $file
        [[ ${file: -1: 1} == s ]] && awk '{print "["$1", "$5"]; "}' $file
        [[ ${file: -1: 1} == t ]] && awk '{print "["$3", "$2"]; "}' $file
        echo "];"
        
        echo "[X, Y, T, AUC] = perfcurve(test(:, 2), test(:, 1), 1);"
        echo "fprintf(fileID, '%s %s %f\n', '${filename: -3: 3}', '$file', AUC);"
    done;
    cd .. ;
done;
echo "fclose(fileID);"
