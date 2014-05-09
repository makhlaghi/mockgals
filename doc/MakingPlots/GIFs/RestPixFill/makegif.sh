fitstojpg -al -irestpixfill.fits -f0 -g0 -q500
convert -delay 1 restpixfill_{1..500}.jpg a.gif
convert a.gif -resize 200x200 s7_restpixfill.gif
rm a.gif *.jpg
