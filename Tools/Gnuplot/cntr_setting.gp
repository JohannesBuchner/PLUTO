#
# FILE: cntr_setting.gp
#
# Set a default settings for making 2D contour plots with gnuplot.
#
# Last modified: Nov 30, 2014 by A. Mignone (mignone@ph.unito.it) 
#

print "> Load default contour setting"

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

set contour base
set view map                         
set cntrparam levels 10                   
unset surface
unset clabel  #  disable coloring
                       
set title "  "
set xlabel "x"
set ylabel "y"
unset key
set style data lines # use lines when doing contour 

# set margins

set lmargin at screen 0.15
set rmargin at screen 0.9
set bmargin at screen 0.1
set tmargin at screen 0.925

