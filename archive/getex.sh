#!/bin/bash

for file in `ls -d */`;
do
    cd $file
    ../create.sh | matlab -nodisplay &> /dev/null
    mv *.pdf ../
    cd ..
done > /dev/null

echo "\\documentclass[a4paper]{article}"
echo "\\usepackage[multidot]{grffile}"
echo "\\usepackage{pdfpages}"
echo "\\begin{document}"
for pdffile in `ls *.pdf`;
do
    echo "\\includepdf[pages=-]{${pdffile}}"
done
echo "\\end{document}"
