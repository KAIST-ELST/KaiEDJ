set xra [  0.0000:   8.2759]
set yra [  0.0000: 657.7469]
set ylabel "Energy (meV)" 
set xtics ("G"    0.0000,"N"    1.5409,"P"    2.6427,"G"    4.5315,"H"    6.7268,"N"    8.2759)
plot "magnon.dat" u 1:2 w l 
replot "vline.dat" u 1:2 w l lc 0 
pause -1
