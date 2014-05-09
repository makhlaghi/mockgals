for i in {0..2}
  do
    for j in {0..10}
       do
        fitstojpg -l -irestpixfill.fits -e$(( ($i*11+$j)*13 )) -o $i\_$j.jpg -f0 -g0 -q500
       done
  done
fitstojpg -l -irestpixfill.fits -e444 -o 2\_10.jpg -f0 -g0 -q500
