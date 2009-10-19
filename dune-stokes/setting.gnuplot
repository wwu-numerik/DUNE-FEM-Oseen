set xlabel "row"
set ylabel "col"
set yrange [] reverse
set terminal postscript enhanced color solid lw 2 "Times-Roman" 20
set output "mat_W.ps"
plot "mat_W.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_E.ps"
plot "mat_E.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_R.ps"
plot "mat_R.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_Z.ps"
plot "mat_Z.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_Y.ps"
plot "mat_Y.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_X.ps"
plot "mat_X.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1
set output "mat_M-inv.ps"
plot "mat_M-inv.gnuplot" using 2:1 with points lt 1 lw 1 pt 1 ps 0.1 
