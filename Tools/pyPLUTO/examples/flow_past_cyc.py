import os
import matplotlib.pyplot as plt
import pyPLUTO as pp
import numpy as np

I = pp.Image()
plutodir = os.environ['PLUTO_DIR']
wdir = plutodir+'/Test_Problems/HD/Viscosity/Flow_Past_Cylinder/'

try:
    D = pp.pload(10,datatype='hdf5',level=3, w_dir=wdir)

except IndexError:
    print 'data.0010.hdf5 not found .. Loading data.0000.hdf5 at level 0'
    D = pp.pload(0,datatype='hdf5',w_dir=wdir) 

#Plot 2D AMR R-Phi Polar Data on Cartesian X-Y Plane. 
I.pldisplay(D, D.rho, x1=D.x1, x2=D.x2, 
            polar=[True, True], fignum=1, 
            cbar=[True, 'horizontal'],figsize=[8,8],
            label1='X', label2='Y',title=r'Density [$\rho$]')
ax1 = plt.gca()
ax1.axis([-10, 50,-20.0, 20.0])

#Using the OplotBox routine to overplot sublevels.
I.oplotbox(D.AMRLevel,lrange=[0,3],geom=D.geometry) 

plt.savefig('amr_flowcyc.png')
