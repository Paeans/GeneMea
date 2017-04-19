i=1
c=1
resultdir=""
echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
for file in `ls *.ms *.mt`;
do
    resultdir=${file::-7}
    [[ $i -gt $(( 4 * 2 )) ]] && { echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf');"
        i=1
        c=$(( c + 1 ))
        echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
    }
    
    filename=${file%.*}
    echo "test = ["
    #cat $file
    [[ ${file: -1: 1} == s ]] && awk '{print "["$1", "$5"]; "}' $file
    [[ ${file: -1: 1} == t ]] && awk '{print "["$3", "$2"]; "}' $file
    echo "];"
    
    echo "[X, Y, T, AUC] = perfcurve(test(:, 2), test(:, 1), 1);"
    echo "AUC"
    
    echo "subplot(4, 2, "${i}");"
    echo "plot(X, Y, 'g');"
    echo "title(['${file//_/\\_} ', strcat('\\fontsize{10}\color{red}', num2str(AUC))], 'Fontsize', 5);"
    
    i=$(( i + 1 ))
done

echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf');"
echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
c=$(( c + 1 ))
i=1
j=1
for file in `ls *.mt *.ft *.vt`;
do
#    [ i%3 == 0 ] && i=$(( i + 1 ))
#    [ i%3 == 2 ] && i=$(( i + 2 ))

    resultdir=${file::-7}
    for t in `seq 1 5`;
    do
    [[ $j -gt 3 ]] && { echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf');"
        i=1
        j=1
        c=$(( c + 1 ))
        echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);"
    }

    filename=${file%.*}
    echo "test1 = ["
    #cat $file
    [[ $t -eq 1 ]] && awk '$2 == 1 {print "["$5", "$6"]; "}' $file
    [[ $t -eq 2 ]] && awk '$2 == 1 {print "["$3", "$4"]; "}' $file
    [[ $t -eq 3 ]] && awk '$2 == 1 {print "["$3", "$6"]; "}' $file
    [[ $t -eq 4 ]] && awk '$2 == 1 {print "["$4", "$6"]; "}' $file
    [[ $t -eq 5 ]] && awk '$2 == 1 {print "["$7", "$6"]; "}' $file
    echo "];"

    echo "test2 = ["
    [[ $t -eq 1 ]] && awk '$2 == 0 {print "["$5", "$6"]; "}' $file
    [[ $t -eq 2 ]] && awk '$2 == 0 {print "["$3", "$4"]; "}' $file
    [[ $t -eq 3 ]] && awk '$2 == 0 {print "["$3", "$6"]; "}' $file
    [[ $t -eq 4 ]] && awk '$2 == 0 {print "["$4", "$6"]; "}' $file
    [[ $t -eq 5 ]] && awk '$2 == 0 {print "["$7", "$6"]; "}' $file
    echo "];"

    echo "ax = subplot(5, 3, "$(( i ))");"
    echo "hold on"
    echo "scatter(test2(:, 1), test2(:, 2), 3, 'g', 'filled', 'o');"
    echo "scatter(test1(:, 1), test1(:, 2), 5, 'r', 'filled', 'o');"
    echo "hold off"
    #echo "axis equal"
    #echo "legend('Y', 'N');"
    [[ $t -eq 1 ]] && echo "xlabel('mean'); ylabel('std');"
    [[ $t -eq 2 ]] && echo "xlabel('max'); ylabel('min');"
    [[ $t -eq 3 ]] && echo "xlabel('max'); ylabel('std');"
    [[ $t -eq 4 ]] && echo "xlabel('min'); ylabel('std');"
    [[ $t -eq 5 ]] && echo "xlabel('median'); ylabel('std');"
    echo "title('${file//_/\\_}', 'Fontsize', 5);"

    i=$(( i + 3 ))
    [[ $i -gt 15 ]] && { j=$(( j + 1 ))
       i=$(( j ))
    }
    done
done


echo "print('"${resultdir}"_"$(printf "%02d" $c)".pdf', '-dpdf')"
