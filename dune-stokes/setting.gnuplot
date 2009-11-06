mat_string = "WXERZYM"
mat_num = 7
cur_ind = 1
xtics_num = 10
ytics_num = 10

set xlabel "row"
set ylabel "col"
set terminal postscript enhanced color solid lw 2 "Times-Roman" 20
call "loop.gnuplot"
