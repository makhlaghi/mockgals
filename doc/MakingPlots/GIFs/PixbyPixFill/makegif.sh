fitstojpg -alt -itmp.fits -f0
convert -delay 1 tmp_{2..400}.jpg a.gif
convert a.gif -resize 300x300 s7_pixfill.gif
rm a.gif *.jpg
