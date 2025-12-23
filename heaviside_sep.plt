set terminal pngcairo size 600,800
set output "heaviside_sep.png"
set multiplot layout 3,1 title "Heaviside and relatives"
set xlabel ""
set ylabel "h(t)"
plot "heaviside_data.dat" using 1:2 with lines notitle
set ylabel "int_h(t)"
plot "heaviside_data.dat" using 1:3 with lines notitle
set xlabel "t"
set ylabel "dh(t)"
plot "heaviside_data.dat" using 1:4 with lines notitle
unset multiplot
