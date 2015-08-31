#
# File: grid.gp
#
# Gnuplot script to read "grid.out" and store the information into 
# gnuplot variables without calling shell or system commands.
#
# At the end of the script, the following variables will be available
#
# - xbeg, ybeg:  initial coordinates of the computational domain;
# - xend, yend:  final   coordinates of the computational domain;
# - dx, dy:      mesh spacings
# - nx, ny:      number of points
#
# Last Modified:  23 Apr, 2014 by A. Mignone (mignone@ph.unito.it)
# 

set term push    # remembers the current terminal setting
set term unknown

# ------------------------------------------------------
# Read grid information from 'grid.out' by reading nx, 
# xbeg = the minimum of 2nd column
# xend = max of 3rd column (xbeg) and so on.       
# Trick reference:
# http://stackoverflow.com/questions/11211339/gnuplot-store-one-number-from-data-file-into-variable
# ------------------------------------------------------

unset pm3d
set xrange[0:1]
set yrange[0:1]
plot 'grid.out' index 0 every 1:1:0:0:0:0
nx = int(GPVAL_DATA_Y_MIN)
plot 'grid.out' index 0 every 1:1:1:0:1:0 using 2
xbeg = GPVAL_DATA_Y_MIN
plot 'grid.out' index 0 every 1:1:nx:0:nx:0 using 3
xend = GPVAL_DATA_Y_MAX

i = nx + 1 
plot 'grid.out' index 0 every 1:1:i:0:i:0
ny = int(GPVAL_DATA_Y_MIN)
plot 'grid.out' index 0 every 1:1:(i+1):0:(i+1):0 using 2
ybeg = GPVAL_DATA_Y_MIN
plot 'grid.out' index 0 every 1:1:(i+ny):0:(i+ny):0 using 3
yend = GPVAL_DATA_Y_MAX

set xrange [ * : * ] noreverse nowriteback # default shown by 'show xrange'
set yrange [ * : * ] noreverse nowriteback 

set term pop # restore previous terminal

# ------------------------------------------------------
# Compute grid spacing (uniform) and domain size
# ------------------------------------------------------

dx = (xend - xbeg)/nx
dy = (yend - ybeg)/ny
Lx = (xend - xbeg)
Ly = (yend - ybeg)

# ------------------------------------------------------
# Do some printing 
# ------------------------------------------------------

print "> Reading grid.out:"
print "  xbeg = ",xbeg,"; xend = ",xend,"; nx = ",nx, "; Lx = ",Lx
print "  ybeg = ",ybeg,"; yend = ",yend,"; ny = ",ny, "; Ly = ",Ly

#set size ratio (yend-ybeg)/(xend-xbeg)

