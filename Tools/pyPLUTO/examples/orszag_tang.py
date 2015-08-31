import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/MHD/Orszag_Tang/'
nlinf = pp.nlast_info(w_dir=wdir)

D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D.
I = pp.Image()
I.pldisplay(D, D.rho,x1=D.x1,x2=D.x2,label1='x',label2='y',title=r'Density $\rho$ [Orszag Tang test]',cbar=(True,'vertical'))
savefig('orszag_tang_1.png') # Only to be saved as either .png or .jpg
show()
