#terminal options (last number = fontsize)

set terminal postscript eps enhanced 20
set output "2ds.eps"


#view options
unset key
set view map
unset surface

#tics
set xtics 50
set ytics 50
set cbtics 0.2

#labels
set xlabel "{/Symbol w}_1/2{/Symbol p}c [cm^{-1}]"
set ylabel "{/Symbol w}_3/2{/Symbol p}c [cm^{-1}]"

#color palette
set palette defined (0 "midnight-blue", 1 "dark-cyan", 2 "yellow", 3 "red")

#contours
set contour
unset clabel
set cntrparam levels incremental -1,0.1,1
set cntrparam cubicspline

#fake filled contours
#set palette maxcolors 10

#border and background
set border linewidth 2
set grid noxtics noytics front
set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back
set object 1 rect fc rgb "midnight-blue" fillstyle solid 1.0

#plot ranges
set cbrange [-1:1]
set xrange [1400:1600]
set yrange [1400:1600]
#interpolated color surface
set pm3d
set pm3d interpolate 4,4

#the plot
splot '2D.par.dat' u 1:2:($4/1000) w l lt 3 lw 3

#convert eps to pdf
#! epstopdf 2ds.eps
