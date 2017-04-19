i=1
c=1
resultdir=""
echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
j=1
for file in `ls *.ms`;
do
#    [ i%3 == 0 ] && i=$(( i + 1 ))
#    [ i%3 == 2 ] && i=$(( i + 2 ))

    resultdir=${file::-7}
    for t in `seq 1 3`;
    do
    [[ $i -gt 12 ]] && { echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf');"
        i=1
        c=$(( c + 1 ))
        echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
    }

    filename=${file%.*}
    echo "X1 = 1:"$( awk '$5 == 1' $file | wc -l )";"
    echo "test1 = ["
    #cat $file
    [[ $t -eq 1 ]] && sort -k1,1n $file | awk '$5 == 1 {print "["$1", "$5"]; "}'
    [[ $t -eq 2 ]] && sort -k7,7n $file | awk '$5 == 1 {print "["$7", "$5"]; "}'
    [[ $t -eq 3 ]] && sort -k8,8n $file | awk '$5 == 1 {print "["$8", "$5"]; "}'
    echo "];"
    echo "X2 = 1:"$( awk '$5 == 0' $file | wc -l )";"
    echo "test2 = ["
    [[ $t -eq 1 ]] && sort -k1,1n $file | awk '$5 == 0 {print "["$1", "$5"]; "}'
    [[ $t -eq 2 ]] && sort -k7,7n $file | awk '$5 == 0 {print "["$7", "$5"]; "}'
    [[ $t -eq 3 ]] && sort -k8,8n $file | awk '$5 == 0 {print "["$8", "$5"]; "}'
    echo "];"

    echo "ax = subplot(4, 3, "$(( i ))");"
    echo "hold on"
    echo "scatter(X2, test2(:, 1), 3, 'g', 'filled', 'o');"
    echo "scatter(X1, test1(:, 1), 5, 'r', 'filled', 'o');"
    echo "hold off"
    #echo "axis equal"
    #echo "legend('Y', 'N');"
    [[ $t -eq 1 ]] && echo "xlabel(''); ylabel('mscore');"
    [[ $t -eq 2 ]] && echo "xlabel(''); ylabel('fscore');"
    [[ $t -eq 3 ]] && echo "xlabel(''); ylabel('vscore');"
    echo "title('${file//_/\\_}', 'Fontsize', 5);"

    i=$(( i + 1 ))
    done
done


echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf')"
