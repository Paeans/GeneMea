i=1
c=1
echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1])"
for file in `ls *.ms`;
do
[[ $i -gt $(( 6 * 3 )) ]] && { echo "print('result_"${c}".pdf', '-dpdf')"
i=1
c=$(( c + 1 ))
echo "figure('PaperType', 'a4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1])" 
}
filename=${file%.*}
echo "test = ["
#cat $file
awk '{print "["$1", "$5"]; "}' $file
echo "];"
echo "[X, Y, T, AUC] = perfcurve(test(:, 2), test(:, 1), 1);"
echo "AUC"
#echo "figure"
echo "subplot(6, 3, "${i}")"
echo "plot(X, Y, 'g')"
echo "title(strcat('${filename:6:3}', {'  '}, num2str(AUC)))"
#echo "print('$filename.pdf', '-dpdf')"
i=$(( i + 1 ))
done

echo "print('result_"${c}".pdf', '-dpdf')"
