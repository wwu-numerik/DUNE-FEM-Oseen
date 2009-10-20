xtics_num = 10
ytics_num = 10

set xlabel "row"
set ylabel "col"
set terminal postscript enhanced color solid lw 2 "Times-Roman" 20

m_name = "mat_W"
input_f = sprintf("%s.gnuplot",m_name) 
output_f = sprintf("%s.ps",m_name) 
plot input_f using 2:1 with dots lt 1 lw 1 
xmax = GPVAL_DATA_X_MAX + 1
ymax = GPVAL_DATA_Y_MAX + 1
set ytics out autofreq 0,floor(ymax/ytics_num)
set xtics out autofreq 0,floor(xmax/xtics_num)
set yrange [-1:ymax] reverse
set xrange [-1:xmax]

#set title sprintf("%d",ax)
set output output_f

replot
cur_ind = cur_ind + 1
if(cur_ind<mat_num) reread
q
set output "mat_E.ps"
plot "mat_E.gnuplot" using 2:1 with dots lt 1 lw 1
set output "mat_R.ps"
plot "mat_R.gnuplot" using 2:1 with dots lt 1 lw 1
set output "mat_Z.ps"
plot "mat_Z.gnuplot" using 2:1 with dots lt 1 lw 1
set output "mat_Y.ps"
plot "mat_Y.gnuplot" using 2:1 with dots lt 1 lw 1
set output "mat_X.ps"
plot "mat_X.gnuplot" using 2:1 with dots lt 1 lw 1
set output "mat_M-inv.ps"
plot "mat_M-inv.gnuplot" using 2:1 with dots lt 1 lw 1 
