set multiplot layout 2, 1 title 'TIME' font "Helvetica,20"
set size ratio 0.333 
set xr [0:3]
set yr [0:1]
set pm3d 
set palette model RGB defined ( 0 'blue', 1 'orange', 2 'yellow')
unset surface
set view map
set size ratio 0.333
set tics font "Helvetica,20"
set xlabel 'x' font "Helvetica,20"
set ylabel 'y' font "Helvetica,20" rotate by 0
set title 'S' font "Helvetica,20"
sp [0:3][0:1] 'S_Ra50alpha0_20000.txt' t ''
#
set palette model RGB defined ( 0 'blue', 1 'red', 2 'orange', 3 'yellow', 4 'white' )
set tics font "Helvetica,20"
set xlabel 'x' font "Helvetica,20"
set ylabel 'y' font "Helvetica,20" rotate by 0
set title 'T' font "Helvetica,20"
sp [0:3][0:1] 'T_Ra50alpha0_20000.txt' t ' '
unset multiplot
#    EOF






