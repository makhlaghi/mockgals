texi2pdf mockgals.texinfo
rm *.aux *.cp *.fn *.ky *.log *.pg *.toc *.tp *.vr
texi2any mockgals.texinfo
texi2any --docbook mockgals.texinfo
texi2any --plaintext -o mockgals.txt mockgals.texinfo

rm -rf ./mockgals 
texi2any --html --split=chapter mockgals.texinfo
