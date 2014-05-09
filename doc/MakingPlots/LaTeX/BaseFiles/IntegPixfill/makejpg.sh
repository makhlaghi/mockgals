for i in {0..2}
  do
    for j in {0..7}
       do
        fitstojpg -l -icentinteg.fits -e$(($i*8+$j)) -o $i\_$j.jpg -f0 -g0 -q200
       done
  done
