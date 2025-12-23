set terminal pngcairo size 800,400
set output "heaviside_all.png"
set key left top
set xlabel "t"
set ylabel "value"
plot "heaviside_data.dat" using 1:2 with lines title "Heaviside", "heaviside_data.dat" using 1:3 with lines title "Integral", "heaviside_data.dat" using 1:4 with lines title "Derivative"
