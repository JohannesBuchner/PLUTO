import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/HD/Rayleigh_Taylor/'

D0 = pp.pload(0,w_dir=wdir)
D1 = pp.pload(1,w_dir=wdir) # Loading the data into a pload object D.
D2 = pp.pload(2,w_dir=wdir)

## SMART WAY##
I = pp.Image()
I.multi_disp(D0.rho,D1.rho,D2.rho,x1=D0.x1,x2=D0.x2,Ncols=3,label1=3*['x'],
             label2=3*['y'],title=[r'$\tau=0$',r'$\tau=1$',r'$\tau=2$'],
             cbar=(True,'vertical','each'),figsize=[12,7])

##BRUTE FORCE WAY##
## f1 = figure(figsize=[12,7])
## ax1 = f1.add_subplot(131)
## pcolormesh(D0.x1,D0.x2,D0.rho.T)
## colorbar()
## ax1.set_xlabel(r'x')
## ax1.set_ylabel(r'y')
## ax1.axis([-0.5,0.5,0.0,4.0])
## ax1.set_aspect('equal')
## ax1.set_title(r'$\tau$ = 0')


## ax2 = f1.add_subplot(132)
## pcolormesh(D1.x1,D1.x2,D1.rho.T)
## colorbar()
## ax2.set_xlabel(r'x')
## ax2.set_ylabel(r'y')
## ax2.axis([-0.5,0.5,0.0,4.0])
## ax2.set_aspect('equal')
## ax2.set_title(r'$\tau$ = 1')

## ax3 = f1.add_subplot(133)
## pcolormesh(D2.x1,D2.x2,D2.rho.T)
## colorbar()
## ax3.set_xlabel(r'x')
## ax3.set_ylabel(r'y')
## ax3.axis([-0.5,0.5,0.0,4.0])
## ax3.set_aspect('equal')
## ax3.set_title(r'$\tau$ = 2')

savefig('RayleighTaylor_multi.png')
show()

