rm *.aux
rm -rf `biber --cache`
if pdflatex -shell-escape Plots.tex
then
  if biber Plots.bcf
    then
      pdflatex -shell-escape Plots.tex
      rm *.bbl *.bcf *.blg *.log *.out *.run.xml *.aux *.auxlock 
    else
      rm *.aux
  fi
else
  rm *.aux
fi 
