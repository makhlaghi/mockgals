for i in {0..7}
  do
    for j in {0..19}
       do
        fitstojpg -l -itmp.fits -e$(($i*20+$j+1)) -o $i\_$j.jpg -f0 -g0	
       done
  done
