#
# File: macros.gp
#
# Gnuplot script to define some useful macros for:
#
# - binary array format specifications
# - 1D slices
# - colormaps
# 
# Note: it requires the data size (dsize=4 or 8) to be set.
#
# Last Modified: 30 Nov, 2014 by A. Mignone (mignone@ph.unito.it)
# 

icut = 0;
jcut = 0;
dformat = (dsize == 4 ? "%f":"%lf")

if (dsize == 4) print "> Data size is single precision"
if (dsize == 8) print "> Data size is double precision"

# ----------------------------------------------------------------------
# Define the BINARR macro for plotting 2D single or double 
# precision data files.
# These macros can be used together with splot to 
# provide the necessary size and grid information, e.g.,
#
# gnuplot> splot "data.0002.dbl"  @BINARR
# ----------------------------------------------------------------------

print "> Setting macros @BINARR"

set macro  # activate macro expansion.
str1 = sprintf("bin array=%dx%d format='%s' ",nx,ny, dformat)
str2 = sprintf("dx=dx dy=dy origin=(xbeg, ybeg, 0.0) ");
str3 = sprintf("skip=(nx*ny*dsize*nvar)");
BINARR = str1.str2.str3    # concatenate strings

# ----------------------------------------------------------------------
# Define the XSLICE and YSLICE macros for plotting 1D profiles along,
# respectively, rows or columns of 2D binary data files.
#
# gnuplot> jcut = 2; plot "data.0002.dbl"  @XSLICE
# gnuplot> icut = 5; plot "data.0002.dbl"  @YSLICE
# ----------------------------------------------------------------------

print "> Setting macros @XSLICE, @YSLICE"

str1 = sprintf("bin array=%d format='%s' ",nx*ny,dformat)
str2 = sprintf("dx=dx origin=(-jcut*Lx+0.5*dx, 0.0) ");
str3 = sprintf("skip=(nx*ny*dsize*nvar) every 1:1:(nx*jcut):0:(nx*jcut+nx-1):0");
XSLICE = str1.str2.str3

str1 = sprintf("bin array=%d format='%s' ",nx*ny,dformat)
str2 = sprintf("dx=dy/nx "); # two consecutive points in y are spaced by nx zones in x
str3 = sprintf("skip=(nx*ny*dsize*nvar) every nx::icut:0:(icut+(ny-1)*nx):0");
YSLICE = str1.str2.str3      # concatenate strings

# ----------------------------------------------------------------------
# Define a few colormap macros
# ----------------------------------------------------------------------

HOT = "rgbformulae 22,13,-31"
RED = "rgbformulae 21,22,23"
RYG = 'model RGB defined ( 0 "red", 0.5 "yellow", 1 "green" )'
RAINBOW = "rgbformulae 33,13,10 " # rainbow (blue-green-yellow-red)
JET = "defined (0  0.0 0.0 0.5, \
                1  0.0 0.0 1.0, \
                2  0.0 0.5 1.0, \
                3  0.0 1.0 1.0, \
                4  0.5 1.0 0.5, \
                5  1.0 1.0 0.0, \
                6  1.0 0.5 0.0, \
                7  1.0 0.0 0.0, \
                8  0.5 0.0 0.0 )"

GREEN_CT = 'model RGB defined (0 "black", 1 "slateblue1", 2 "white")'
print "  - available colormap macros: HOT, RED, RYG, RAINBOW, JET, GREEN_CT"

