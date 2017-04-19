#!/bin/bash

python main.py e
# generate sub sequence of gene data
# for every peak in chr?

step=20
# define how many motifs will each process execute
for j in `seq 0 19`
do
    python main.py c $(( j * step )) $step &
done
# calc the score list of motif from $i to ($i + $step)
# on the gene sub sequences of each peak

step=1
# define how many peaks will each process execute
for j in `seq 0 14`
do
    python main.py g $(( j * step )) $step &
done
# calc freq score list of peak $i to ($i + $step)

for j in `seq 0 14`
do
    python main.py m $(( j * step )) $step &
done
# calc the centdist score of each motif

for j in `seq 0 14`
do
    python main.py g $(( j * step )) $step  && python main.py m $(( i * step )) $step &
done
######
# the three steps need to be execute sequencely
# the 2nd part need the result of 1st part
# and the 3rd part need the result of 2nd part
# execute each part SEPRATELY and SEQUENCELY
######

wc -l file.name

sort -k1 -nr file.name

#while read p; do echo -n ${p}":"; cat TTFF.txt | grep $p | awk '{printf $1}'; echo -n "; "; done < res.txt
#while read p; do echo -n \"${p}\"":"; cat TTFF.txt | grep $p | awk '{printf "\""$1"\""}'; echo -n "; "; done < res.txt
#cat GSM.ms | awk '{print $4}' | while read p; do echo -n ${p}": "; cat tfgroup.txt | grep $p | awk '{print $1}'; done
#cat GSM.ms | awk '{print $4}' | while read p; do echo -n ${p}": "; cat tfgroup.txt | grep $p | awk '{print $1}'; done | awk '{print $2}'
#awk '!x[$0]++' #unique without sorting. SO AMAZING!

#echo -n "[" && while read -r -a p; do echo -n "["${p[0]}", "; awk '$2 =~ (^|\|)"${p[3]}"($|\|) && "NANOG SOX" =~ (^|[[:space:]])"$1"($|[[:space:]]) {printf "1], "} $2 =~ (^|\|)"${p[3]}"($|\|) && ! "NANOG SOX" =~ (^|[[:space:]])"$1"($|[[:space:]]) {printf "0], "}' tfgroup.txt; done < msfile && echo "]" > result.txtecho -n "["


#while read -r -a p;
#do
#	echo -n "["${p[0]}", "
##	awk '$2 =~ (^|\|)"${p[3]}"($|\|), "NANOG SOX" =~ (^|[[:space:]])"$1"($|[[:space:]]) {printf "1], "} \
##	     $2 =~ (^|\|)"${p[3]}"($|\|), ! "NANOG SOX" =~ (^|[[:space:]])"$1"($|[[:space:]]) {printf "0], "}' ../tfgroup.txt
##	awk '/"${p[3]}"/ ~ "$2" {printf $1}' ../tfgroup.txt
#	grep "${p[3]}" ../tfgroup.txt | (
#		read -r -a q
#		[[ "NANOG SOX ERE OCT STAT" =~ (^|[[:space:]])"${q[0]}"($|[[:space:]]) ]] && echo -n "1], " || echo -n "0], " )
#done < 1.1.nok.GSM288345_ES_Nanog.ms
#echo "]"
