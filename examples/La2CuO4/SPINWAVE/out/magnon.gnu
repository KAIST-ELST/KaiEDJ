set xra [  0.0000:   3.3949]
set yra [  0.0000: 660.0945]
set ylabel "Energy (meV)" 
set xtics ("P"    0.0000,"M"    0.5811,"X"    1.4035,"P"    1.9880,"G"    2.5691,"X"    3.3949)
plot "magnon.dat" u 1:2 w l 
replot "magnon.dat" u 1:3 w l 
replot "vline.dat" u 1:2 w l lc 0 
pause -1
