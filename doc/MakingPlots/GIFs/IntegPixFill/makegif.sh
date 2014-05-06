fitstojpg -al -icentinteg.fits -f0 -g0 -q500
convert -delay 50 centinteg_{2..33}.jpg a.gif
convert a.gif -resize 200x200 s7_integpixfill.gif
rm a.gif *.jpg
