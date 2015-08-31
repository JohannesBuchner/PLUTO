#
# FILE: pm3d_setting.gp
#
# Set a default settings for viewing PLUTO binary data files
# with Gnuplot using the pm3d style. 
#
# Last modified: Nov 30, 2014 by A. Mignone (mignone@ph.unito.it) 
#

print "> Load default pm3d setting"

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

unset contour
set title "  "
set xlabel "x"
set ylabel "y"
#set key top left
unset key
set surface
set pm3d map

# -- Set margins --

set lmargin at screen 0.15
set rmargin at screen 0.8
set bmargin at screen 0.1
set tmargin at screen 0.93

set palette defined
