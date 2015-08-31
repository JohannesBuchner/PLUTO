import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/HD/Sod/'
nlinf = pp.nlast_info(w_dir=wdir)

D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D.

f1 = figure()
ax1 = f1.add_subplot(111)
plot(D.x1,D.rho,'r',D.x1,D.prs,'k',D.x1,D.vx1,'g')
xlabel(r'x')
ylabel(r'$\rho$ [red], P [black], $V_{\rm x}$ [green]')
title(r'Sod shock Tube test')
axis([0.0,1.0,-0.2,1.2])
savefig('sod_1.pdf')
show()
