TEXISRC=./docsrc/mockgals.texinfo
JAVASCRIPTS=./docsrc/javascripts.txt
HTMLINDEX=./mockgals/index.html

#Make the mockgals folder and the HTML:
rm -rf ./mockgals
texi2any --html --css-ref=./manualstyle.css  $TEXISRC
mkdir ./mockgals/mockgals-figures/
./addjavascript $JAVASCRIPTS ./mockgals/*.html
./correctindextop $HTMLINDEX
cp ./docsrc/manualstyle.css ./mockgals/
cp ./mockgals-figures/*.jpg ./mockgals/mockgals-figures/
cp ./mockgals-figures/*.png ./mockgals/mockgals-figures/
cp ./mockgals-figures/*.gif ./mockgals/mockgals-figures/

#Make the PDF:
texi2pdf $TEXISRC
rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr *.cps *.fns
mv mockgals.pdf ./mockgals/

#Make the three other verions:
texi2any $TEXISRC
texi2any --docbook $TEXISRC
texi2any --plaintext -o mockgals.txt $TEXISRC
mv mockgals.info mockgals.txt mockgals.xml ./mockgals/

#Also copy the source file:
cp ./docsrc/mockgals.texinfo ./mockgals/


