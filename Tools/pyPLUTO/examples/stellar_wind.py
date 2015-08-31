import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

#To run this example it is suggested to get data in 2D using pluto_01.ini and set the data in flt datatype instead of dbl.h5

plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/HD/Stellar_Wind/'
nlinf = pp.nlast_info(w_dir=wdir,datatype='float')

D = pp.pload(nlinf['nlast'],w_dir=wdir,datatype='float') # Loading the data into a pload object D.

I = pp.Image()
I.pldisplay(D, log10(D.rho),x1=D.x1,x2=D.x2,label1='x',label2='y',
            title=r'Log Density $\rho$ [Stellar Wind]',cbar=(True,'vertical'),figsize=[8,12])

# Code to plot arrows. --> Spacing between the arrow can be adjusted by modifying the newdims tuple of conrid function.
T = pp.Tools()
newdims = 2*(20,)
Xmesh, Ymesh = meshgrid(D.x1.T,D.x2.T)
xcong = T.congrid(Xmesh,newdims,method='linear')
ycong = T.congrid(Ymesh,newdims,method='linear')
velxcong = T.congrid(D.vx1.T,newdims,method='linear')
velycong = T.congrid(D.vx2.T,newdims,method='linear')
gca().quiver(xcong, ycong, velxcong, velycong,color='w')

savefig('stellar_wind_1.png')
show()
