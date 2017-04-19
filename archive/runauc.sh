#!/bin/bash

biodir=~/GeneFinding/bioinformatic/archive

function getauc {
    . $biodir/mergefile.sh
    . $biodir/calcauc.sh > result.m
    matlab -nodisplay -nodesktop -nojvm < result.m
}


for directory in `ls -d */`;
do
    cd $directory
    #. $biodir/mergefile.sh
    #. $biodir/calcauc.sh > result.m
    #matlab -nodisplay -nodesktop -nojvm < result.m &
    getauc &
    cd ..
done

wait `ps -e | grep getauc | awk '{print $1}'`
