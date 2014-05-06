for i in {0..2}
  do
    for j in {0..10}
       do
        fitstojpg -l -icentinteg.fits -e$(($i*11+$j)) -o $i\_$j.jpg -f0 -g0 -q500
       done
  done
