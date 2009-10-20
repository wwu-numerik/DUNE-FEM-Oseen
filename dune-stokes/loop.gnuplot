m_name = "mat_".mat_string[cur_ind:cur_ind]
input_f = sprintf("%s.gnuplot",m_name) 

set output #don't plot in prev ps after first run
plot input_f using 2:1 with dots lt 1 lw 1 

xmax = GPVAL_DATA_X_MAX + 1
ymax = GPVAL_DATA_Y_MAX + 1

set ytics out autofreq 0,floor(ymax/ytics_num)
set xtics out autofreq 0,floor(xmax/xtics_num)
set yrange [-1:ymax] reverse
set xrange [-1:xmax]

#set title sprintf("%d",ax)

output_f = sprintf("%s.ps",m_name) 
set output output_f

replot
cur_ind = cur_ind + 1
if(cur_ind<mat_num+1) reread

