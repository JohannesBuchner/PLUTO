# -*- coding: utf-8 -*-
import os
import sys
import array
import numpy as np
import scipy.ndimage
import scipy.interpolate
from scipy.interpolate import UnivariateSpline
from matplotlib.pyplot import *
from matplotlib.mlab import *

####### Check for h5py to Read AMR data ######
#try:
#    import h5py as h5
#    hasH5 = True
#except ImportError:
#    hasH5 = False
hasH5 = False

def curdir():
	""" Get the current working directory.
	"""
        curdir = os.getcwd()+'/'
        return curdir

def get_nstepstr(ns):
	""" Convert the float input *ns* into a string that would match the data file name.

	**Inputs**:

	 ns -- Integer number that represents the time step number. E.g., The ns for data.0001.dbl is 1.\n

	**Outputs**:

	 Returns the string that would be used to complete the data file name. E.g., for data.0001.dbl, ns = 1 and pyPLUTO.get_nstepstr(1) returns '0001'
	 """
	nstepstr = str(ns)
	while len(nstepstr) < 4:
           nstepstr= '0'+nstepstr
	return nstepstr

def nlast_info(w_dir=None,datatype=None):
	""" Prints the information of the last step of the simulation as obtained from out files

	**Inputs**:
	
	  w_dir -- path to the directory which has the dbl.out(or flt.out) and the data\n
	  datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

        **Outputs**:
	
	  This function returns a dictionary with following keywords - \n

	  nlast -- The ns for the last file saved.\n
	  time -- The simulation time for the last file saved.\n
	  dt -- The time step dt for the last file. \n
	  Nstep -- The Nstep value for the last file saved.


	**Usage**:
	
	  In case the data is 'float'.

	  ``wdir = /path/to/data/directory``\n
	  ``import pyPLUTO as pp``\n
	  ``A = pp.nlast_info(w_dir=wdir,datatype='float')``	
	"""
	if w_dir is None: w_dir=curdir()
	if datatype == 'float':
		fname_v = w_dir+"flt.out"
	elif datatype == 'vtk':
		fname_v = w_dir+"vtk.out"
	else:
		fname_v = w_dir+"dbl.out"
	last_line = file(fname_v,"r").readlines()[-1].split()
	nlast = int(last_line[0])
	SimTime =  float(last_line[1])
	Dt = float(last_line[2])
	Nstep = int(last_line[3])
	    
	print "------------TIME INFORMATION--------------"
	print 'nlast =',nlast
	print 'time  =',SimTime
	print 'dt    =', Dt
	print 'Nstep =',Nstep
	print "-------------------------------------------"
	    
	return {'nlast':nlast,'time':SimTime,'dt':Dt,'Nstep':Nstep}
    

class pload(object):
    def __init__(self, ns, w_dir=None, datatype=None, level = 0, x1range=None, x2range=None, x3range=None):
        """Loads the data.
	
        **Inputs**:
	  
          ns -- Step Number of the data file\n
	  w_dir -- path to the directory which has the data files\n
          datatype -- Datatype (default = 'double')
	  
        **Outputs**:
          
          pyPLUTO pload object whose keys are arrays of data values.

	"""
        self.NStep = ns
        self.Dt = 0.0

        self.n1 = 0
        self.n2 = 0
        self.n3 = 0

        self.x1 = []
        self.x2 = []
        self.x3 = []
        self.dx1 = []
        self.dx2 = []
        self.dx3 = []
        
        self.x1range = x1range
        self.x2range = x2range
        self.x3range = x3range

        self.NStepStr = str(self.NStep)
        while len(self.NStepStr) < 4:
            self.NStepStr = '0'+self.NStepStr

        if datatype is None:
            datatype = "double"
        self.datatype = datatype

        if ((not hasH5) and (datatype == 'hdf5')):
            print 'To read AMR hdf5 files with python'
            print 'Please install h5py (Python HDF5 Reader)'
            return

        self.level = level

        if w_dir is None:
            w_dir = os.getcwd() + '/'
        self.wdir = w_dir
        
        Data_dictionary = self.ReadDataFile(self.NStepStr)
        for keys in Data_dictionary:
            object.__setattr__(self, keys, Data_dictionary.get(keys))

    def ReadTimeInfo(self, timefile):
        """ Read time info from the outfiles.

	**Inputs**:
	  
	  timefile -- name of the out file which has timing information. 

	"""

        if (self.datatype == 'hdf5'):
            fh5 = h5.File(timefile,'r') 
            self.SimTime = fh5.attrs.get('time')
            #self.Dt = 1.e-2 # Should be erased later given the level in AMR
            fh5.close()
        else:
            ns = self.NStep
            f_var = open(timefile, "r")
            tlist = []
            for line in f_var.readlines():
                tlist.append(line.split())
            self.SimTime = float(tlist[ns][1])
            self.Dt = float(tlist[ns][2])

    def ReadVarFile(self, varfile):
        """ Read variable names from the outfiles.

	**Inputs**:
	  
	  varfile -- name of the out file which has variable information. 

	"""
        if (self.datatype == 'hdf5'):
            fh5 = h5.File(varfile,'r') 
            self.filetype = 'single_file' 
            self.endianess = '>' # not used with AMR, kept for consistency
            self.vars = []
            for iv in range(fh5.attrs.get('num_components')):
                self.vars.append(fh5.attrs.get('component_'+str(iv)))
            fh5.close()
        else:
            vfp = open(varfile, "r")
            varinfo = vfp.readline().split()
            self.filetype = varinfo[4]
            self.endianess = varinfo[5]
            self.vars = varinfo[6:]
            vfp.close()

    def ReadGridFile(self, gridfile):
        """ Read grid values from the grid.out file.

	**Inputs**:
	  
	  gridfile -- name of the grid.out file which has information about the grid. 

	"""
        xL = []
        xR = []
        nmax = []
        gfp = open(gridfile, "r")
        for i in gfp.readlines():
            if len(i.split()) == 1:
                try:
                    int(i.split()[0])
                    nmax.append(int(i.split()[0]))
                except:
                    pass
                
            if len(i.split()) == 3:
                try:
                    int(i.split()[0])
                    xL.append(float(i.split()[1]))
                    xR.append(float(i.split()[2]))
                except:
                    if (i.split()[1] == 'GEOMETRY:'):
                        self.geometry=i.split()[2]
                    pass
                
        self.n1, self.n2, self.n3 = nmax
        n1 = self.n1
        n1p2 = self.n1 + self.n2
        n1p2p3 = self.n1 + self.n2 + self.n3
        self.x1 = np.asarray([0.5*(xL[i]+xR[i]) for i in range(n1)])
        self.dx1 = np.asarray([(xR[i]-xL[i]) for i in range(n1)])
        self.x2 = np.asarray([0.5*(xL[i]+xR[i]) for i in range(n1, n1p2)])
        self.dx2 = np.asarray([(xR[i]-xL[i]) for i in range(n1, n1p2)])
        self.x3 = np.asarray([0.5*(xL[i]+xR[i]) for i in range(n1p2, n1p2p3)])
        self.dx3 = np.asarray([(xR[i]-xL[i]) for i in range(n1p2, n1p2p3)])


        # Stores the total number of points in '_tot' variable in case only
        # a portion of the domain is loaded. Redefine the x and dx arrays
        # to match the requested ranges
        self.n1_tot = self.n1 ; self.n2_tot = self.n2 ; self.n3_tot = self.n3
        if (self.x1range != None):
            self.n1_tot = self.n1
            self.irange = range(abs(self.x1-self.x1range[0]).argmin(),abs(self.x1-self.x1range[1]).argmin()+1)
            self.n1  = len(self.irange)
            self.x1  = self.x1[self.irange]
            self.dx1 = self.dx1[self.irange]
        else:
            self.irange = range(self.n1)
        if (self.x2range != None):
            self.n2_tot = self.n2
            self.jrange = range(abs(self.x2-self.x2range[0]).argmin(),abs(self.x2-self.x2range[1]).argmin()+1)
            self.n2  = len(self.jrange)
            self.x2  = self.x2[self.jrange]
            self.dx2 = self.dx2[self.jrange]
        else:
            self.jrange = range(self.n2)
        if (self.x3range != None):
            self.n3_tot = self.n3
            self.krange = range(abs(self.x3-self.x3range[0]).argmin(),abs(self.x3-self.x3range[1]).argmin()+1)
            self.n3  = len(self.krange)
            self.x3  = self.x3[self.krange]
            self.dx3 = self.dx3[self.krange]
        else:
            self.krange = range(self.n3)
        self.Slice=(self.x1range != None) or (self.x2range != None) or (self.x3range != None)


        # Create the xr arrays containing the edges positions
        # Useful for pcolormesh which should use those
        self.x1r = np.zeros(len(self.x1)+1) ; self.x1r[1:] = self.x1 + self.dx1/2.0 ; self.x1r[0] = self.x1r[1]-self.dx1[0]
        self.x2r = np.zeros(len(self.x2)+1) ; self.x2r[1:] = self.x2 + self.dx2/2.0 ; self.x2r[0] = self.x2r[1]-self.dx2[0]
        self.x3r = np.zeros(len(self.x3)+1) ; self.x3r[1:] = self.x3 + self.dx3/2.0 ; self.x3r[0] = self.x3r[1]-self.dx3[0]


        prodn = self.n1*self.n2*self.n3
        if prodn == self.n1:
            self.nshp = (self.n1)
        elif prodn == self.n1*self.n2:
            self.nshp = (self.n2, self.n1)
        else:
            self.nshp = (self.n3, self.n2, self.n1)


    def DataScanVTK(self, fp, n1, n2, n3, endian, dtype):
        """ Scans the VTK data files. 
	
        **Inputs**:
	
	 fp -- Data file pointer\n
	 n1 -- No. of points in X1 direction\n
	 n2 -- No. of points in X2 direction\n
	 n3 -- No. of points in X3 direction\n
	 endian -- Endianess of the data\n
	 dtype -- datatype 
	 
        **Output**:
          
          Dictionary consisting of variable names as keys and its values. 

	"""
        ks = []
        vtkvar = []
        while True:
            l = fp.readline()
            try:
                l.split()[0]
            except IndexError:
                pass
            else:
                if l.split()[0] == 'SCALARS':
                    ks.append(l.split()[1])
                elif l.split()[0] == 'LOOKUP_TABLE':
                    A = array.array(dtype)
                    fmt = endian+str(n1*n2*n3)+dtype
                    nb = np.dtype(fmt).itemsize 
                    A.fromstring(fp.read(nb))
                    if (self.Slice):
                        darr = np.zeros((n1*n2*n3))
                        indxx = np.sort([self.n3_tot*self.n2_tot*k + j*self.n2_tot + i for i in self.irange for j in self.jrange for k in self.krange])
                        if (sys.byteorder != self.endianess):
                            A.byteswap()
                        for ii,iii in enumerate(indxx):
                            darr[ii] = A[iii]
                        vtkvar_buf = [darr]
                    else:
                        vtkvar_buf = np.frombuffer(A,dtype=np.dtype(fmt))
                    vtkvar.append(np.reshape(vtkvar_buf,self.nshp).transpose())
                else:
                    pass
            if l == '':
                break
        
        vtkvardict = dict(zip(ks,vtkvar))
        return vtkvardict
            
    def DataScanHDF5(self, fp, myvars, ilev):
        """ Scans the Chombo HDF5 data files for AMR in PLUTO. 
        
        **Inputs**:
		
          fp     -- Data file pointer\n
          myvars -- Names of the variables to read\n
          ilev   -- required AMR level
		
        **Output**:
        
          Dictionary consisting of variable names as keys and its values. 

        **Note**:

          Due to the particularity of AMR, the grid arrays loaded in ReadGridFile are overwritten here.
        
        """
        # Read the grid information
        dim = fp['Chombo_global'].attrs.get('SpaceDim')
        nlev = fp.attrs.get('num_levels')
        il = min(nlev-1,ilev)
        lev  = []
        for i in range(nlev):
            lev.append('level_'+str(i))	    
        freb = np.zeros(nlev,dtype='int')
        for i in range(il+1)[::-1]:
            fl = fp[lev[i]]
            if (i == il):
                pdom = fl.attrs.get('prob_domain')
                dx = fl.attrs.get('dx')
                dt = fl.attrs.get('dt')
                ystr = 1. ; zstr = 1. ; logr = 0
                try:
                    geom = fl.attrs.get('geometry')
                    logr = fl.attrs.get('logr')
                    if (dim == 2):
                        ystr = fl.attrs.get('g_x2stretch')
                    elif (dim == 3):
                        zstr = fl.attrs.get('g_x3stretch')
                except:
                    print 'Old HDF5 file, not reading stretch and logr factors'
                freb[i] = 1
                x1b = fl.attrs.get('domBeg1')
                if (dim == 1):
                    x2b = 0
                else:
                    x2b = fl.attrs.get('domBeg2')
                if (dim == 1 or dim == 2):
                    x3b = 0
                else:
                    x3b = fl.attrs.get('domBeg3')
                jbeg = 0 ; jend = 0 ; ny = 1
                kbeg = 0 ; kend = 0 ; nz = 1
                if (dim == 1):
                    ibeg = pdom[0] ; iend = pdom[1] ; nx = iend-ibeg+1
                elif (dim == 2):
                    ibeg = pdom[0] ; iend = pdom[2] ; nx = iend-ibeg+1
                    jbeg = pdom[1] ; jend = pdom[3] ; ny = jend-jbeg+1
                elif (dim == 3):
                    ibeg = pdom[0] ; iend = pdom[3] ; nx = iend-ibeg+1
                    jbeg = pdom[1] ; jend = pdom[4] ; ny = jend-jbeg+1
                    kbeg = pdom[2] ; kend = pdom[5] ; nz = kend-kbeg+1
            else:
                rat = fl.attrs.get('ref_ratio')
                freb[i] = rat*freb[i+1]
        
        dx0 = dx*freb[0]

        ## Allow to load only a portion of the domain
        if (self.x1range != None):
            if logr == 0:
                self.x1range = self.x1range-x1b
            else:
                self.x1range = [log(self.x1range[0]/x1b),log(self.x1range[1]/x1b)]
            ibeg0 = min(self.x1range)/dx0 ; iend0 = max(self.x1range)/dx0
            ibeg  = max([ibeg, int(ibeg0*freb[0])]) ; iend = min([iend,int(iend0*freb[0]-1)])
            nx = iend-ibeg+1
        if (self.x2range != None):
            self.x2range = (self.x2range-x2b)/ystr
            jbeg0 = min(self.x2range)/dx0 ; jend0 = max(self.x2range)/dx0
            jbeg  = max([jbeg, int(jbeg0*freb[0])]) ; jend = min([jend,int(jend0*freb[0]-1)])
            ny = jend-jbeg+1
        if (self.x3range != None):
            self.x3range = (self.x3range-x3b)/zstr
            kbeg0 = min(self.x3range)/dx0 ; kend0 = max(self.x3range)/dx0
            kbeg  = max([kbeg, int(kbeg0*freb[0])]) ; kend = min([kend,int(kend0*freb[0]-1)])
            nz = kend-kbeg+1
	    
        ## Create uniform grids at the required level
        if logr == 0:
            x1 = x1b + (ibeg+np.array(range(nx))+0.5)*dx
        else:
            x1 = x1b*(exp((ibeg+np.array(range(nx))+1)*dx)+exp((ibeg+np.array(range(nx)))*dx))*0.5
        
        x2 = x2b + (jbeg+np.array(range(ny))+0.5)*dx*ystr
        x3 = x3b + (kbeg+np.array(range(nz))+0.5)*dx*zstr
        if logr == 0:
            dx1 = np.ones(nx)*dx
        else:
            dx1 = x1b*(exp((ibeg+np.array(range(nx))+1)*dx)-exp((ibeg+np.array(range(nx)))*dx))
        dx2 = np.ones(ny)*dx*ystr
        dx3 = np.ones(nz)*dx*zstr

        # Create the xr arrays containing the edges positions
        # Useful for pcolormesh which should use those
        x1r = np.zeros(len(x1)+1) ; x1r[1:] = x1 + dx1/2.0 ; x1r[0] = x1r[1]-dx1[0]
        x2r = np.zeros(len(x2)+1) ; x2r[1:] = x2 + dx2/2.0 ; x2r[0] = x2r[1]-dx2[0]
        x3r = np.zeros(len(x3)+1) ; x3r[1:] = x3 + dx3/2.0 ; x3r[0] = x3r[1]-dx3[0]
        NewGridDict = dict([('n1',nx),('n2',ny),('n3',nz),\
                            ('x1',x1),('x2',x2),('x3',x3),\
                            ('x1r',x1r),('x2r',x2r),('x3r',x3r),\
                            ('dx1',dx1),('dx2',dx2),('dx3',dx3),\
                            ('Dt',dt)])
	    
        # Variables table
        nvar = len(myvars)
        vars = np.zeros((nx,ny,nz,nvar))
        
        LevelDic = {'nbox':0,'ibeg':ibeg,'iend':iend,'jbeg':jbeg,'jend':jend,'kbeg':kbeg,'kend':kend}
        AMRLevel = [] 
        AMRBoxes = np.zeros((nx,ny,nz))
        for i in range(il+1):
            AMRLevel.append(LevelDic.copy())
            fl = fp[lev[i]]
            data = fl['data:datatype=0']
            boxes = fl['boxes']
            nbox = len(boxes['lo_i'])
            AMRLevel[i]['nbox'] = nbox
            ncount = 0L
            AMRLevel[i]['box']=[]
            for j in range(nbox): # loop on all boxes of a given level
                AMRLevel[i]['box'].append({'x0':0.,'x1':0.,'ib':0L,'ie':0L,\
                                           'y0':0.,'y1':0.,'jb':0L,'je':0L,\
                                           'z0':0.,'z1':0.,'kb':0L,'ke':0L})
                # Box indexes
                ib = boxes[j]['lo_i'] ; ie = boxes[j]['hi_i'] ; nbx = ie-ib+1
                jb = 0 ; je = 0 ; nby = 1
                kb = 0 ; ke = 0 ; nbz = 1
                if (dim > 1):
                    jb = boxes[j]['lo_j'] ; je = boxes[j]['hi_j'] ; nby = je-jb+1
                if (dim > 2):
                    kb = boxes[j]['lo_k'] ; ke = boxes[j]['hi_k'] ; nbz = ke-kb+1
                szb = nbx*nby*nbz*nvar
                # Rescale to current level
                kb = kb*freb[i] ; ke = (ke+1)*freb[i] - 1
                jb = jb*freb[i] ; je = (je+1)*freb[i] - 1
                ib = ib*freb[i] ; ie = (ie+1)*freb[i] - 1
	    		
                # Skip boxes lying outside ranges
                if ((ib > iend) or (ie < ibeg) or \
                    (jb > jend) or (je < jbeg) or \
                    (kb > kend) or (ke < kbeg)):
                    ncount = ncount + szb
                else:

                    ### Read data
                    q = data[ncount:ncount+szb].reshape((nvar,nbz,nby,nbx)).T
                    
                    ### Find boxes intersections with current domain ranges
                    ib0 = max([ibeg,ib]) ; ie0 = min([iend,ie])
                    jb0 = max([jbeg,jb]) ; je0 = min([jend,je])
                    kb0 = max([kbeg,kb]) ; ke0 = min([kend,ke])
	            
                    ### Store box corners in the AMRLevel structure
                    if logr == 0:
                        AMRLevel[i]['box'][j]['x0'] = x1b + dx*(ib0)
                        AMRLevel[i]['box'][j]['x1'] = x1b + dx*(ie0+1)
                    else:
                        AMRLevel[i]['box'][j]['x0'] = x1b*exp(dx*(ib0))
                        AMRLevel[i]['box'][j]['x1'] = x1b*exp(dx*(ie0+1))
                    AMRLevel[i]['box'][j]['y0'] = x2b + dx*(jb0)*ystr
                    AMRLevel[i]['box'][j]['y1'] = x2b + dx*(je0+1)*ystr
                    AMRLevel[i]['box'][j]['z0'] = x3b + dx*(kb0)*zstr
                    AMRLevel[i]['box'][j]['z1'] = x3b + dx*(ke0+1)*zstr
                    AMRLevel[i]['box'][j]['ib'] = ib0 ; AMRLevel[i]['box'][j]['ie'] = ie0 
                    AMRLevel[i]['box'][j]['jb'] = jb0 ; AMRLevel[i]['box'][j]['je'] = je0 
                    AMRLevel[i]['box'][j]['kb'] = kb0 ; AMRLevel[i]['box'][j]['ke'] = ke0 
                    AMRBoxes[ib0-ibeg:ie0-ibeg+1, jb0-jbeg:je0-jbeg+1, kb0-kbeg:ke0-kbeg+1] = il
                    
                    ### Extract the box intersection from data stored in q
                    cib0 = (ib0-ib)/freb[i] ; cie0 = (ie0-ib)/freb[i]
                    cjb0 = (jb0-jb)/freb[i] ; cje0 = (je0-jb)/freb[i]
                    ckb0 = (kb0-kb)/freb[i] ; cke0 = (ke0-kb)/freb[i]
                    q1 = np.zeros((cie0-cib0+1, cje0-cjb0+1, cke0-ckb0+1,nvar))
                    q1 = q[cib0:cie0+1,cjb0:cje0+1,ckb0:cke0+1,:]
                    
                    # Remap the extracted portion
                    if (dim == 1):
                        new_shape = (ie0-ib0+1,1)
                    elif (dim == 2):
                        new_shape = (ie0-ib0+1,je0-jb0+1)
                    else:
                        new_shape = (ie0-ib0+1,je0-jb0+1,ke0-kb0+1)
                        
                    stmp = list(new_shape)
                    while stmp.count(1) > 0:
                        stmp.remove(1)
                    new_shape = tuple(stmp)
                    
                    myT = Tools()
                    for iv in range(nvar):
                        vars[ib0-ibeg:ie0-ibeg+1,jb0-jbeg:je0-jbeg+1,kb0-kbeg:ke0-kbeg+1,iv] = \
                            myT.congrid(q1[:,:,:,iv].squeeze(),new_shape,method='linear',minusone=True).reshape((ie0-ib0+1,je0-jb0+1,ke0-kb0+1))
                    ncount = ncount+szb
	    
        h5vardict={}
        for iv in range(nvar):
            h5vardict[myvars[iv]] = vars[:,:,:,iv].squeeze()
        AMRdict = dict([('AMRBoxes',AMRBoxes),('AMRLevel',AMRLevel)])
        OutDict = dict(NewGridDict)
        OutDict.update(AMRdict)
        OutDict.update(h5vardict)
        return OutDict

    
    def DataScan(self, fp, n1, n2, n3, endian, dtype, off=None):
	""" Scans the data files in all formats. 
        
        **Inputs**:
          
          fp -- Data file pointer\n
          n1 -- No. of points in X1 direction\n
          n2 -- No. of points in X2 direction\n
          n3 -- No. of points in X3 direction\n
          endian -- Endianess of the data\n
          dtype -- datatype, eg : double, float, vtk, hdf5\n
          off -- offset (for avoiding staggered B fields) 
	 
        **Output**:
         
          Dictionary consisting of variable names as keys and its values. 

	"""
        if off is not None:
            off_fmt = endian+str(off)+dtype
            nboff = np.dtype(off_fmt).itemsize
            fp.read(nboff)

        n1_tot = self.n1_tot ; n2_tot = self.n2_tot; n3_tot = self.n3_tot

        A = array.array(dtype)
        fmt = endian+str(n1_tot*n2_tot*n3_tot)+dtype
        nb = np.dtype(fmt).itemsize 
        A.fromstring(fp.read(nb))
        
        if (self.Slice):
            darr = np.zeros((n1*n2*n3))
            indxx = np.sort([n3_tot*n2_tot*k + j*n2_tot + i for i in self.irange for j in self.jrange for k in self.krange])
            if (sys.byteorder != self.endianess):
                A.byteswap()
            for ii,iii in enumerate(indxx):
                darr[ii] = A[iii]
            darr = [darr]
        else:
            darr = np.frombuffer(A,dtype=np.dtype(fmt))
        
        return np.reshape(darr[0],self.nshp).transpose()

    def ReadSingleFile(self, datafilename, myvars, n1, n2, n3, endian,
                       dtype, ddict):
        """Reads a single data file, data.****.dtype.
	
        **Inputs**:	

          datafilename -- Data file name\n
	  myvars -- List of variable names to be read\n
          n1 -- No. of points in X1 direction\n
          n2 -- No. of points in X2 direction\n
          n3 -- No. of points in X3 direction\n
          endian -- Endianess of the data\n
          dtype -- datatype\n
          ddict -- Dictionary containing Grid and Time Information
          which is updated
	 
        **Output**:

          Updated Dictionary consisting of variable names as keys and its values.
	"""
        if self.datatype == 'hdf5':
            fp = h5.File(datafilename,'r')
        else:
            fp = open(datafilename, "rb")
        
        print "Reading Data file : %s"%datafilename
        
        if self.datatype == 'vtk':
            vtkd = self.DataScanVTK(fp, n1, n2, n3, endian, dtype)
            ddict.update(vtkd)
        elif self.datatype == 'hdf5':
            h5d = self.DataScanHDF5(fp,myvars,self.level)
            ddict.update(h5d)      
        else:
            for i in range(len(myvars)):
                if myvars[i] == 'bx1s':
                    ddict.update({myvars[i]: self.DataScan(fp, n1, n2, n3, endian,
                                                           dtype, off=n2*n3)})
                elif myvars[i] == 'bx2s':
                    ddict.update({myvars[i]: self.DataScan(fp, n1, n2, n3, endian,
                                                           dtype, off=n3*n1)})
                elif myvars[i] == 'bx3s':
                    ddict.update({myvars[i]: self.DataScan(fp, n1, n2, n3, endian,
                                                           dtype, off=n1*n2)})
                else:
                    ddict.update({myvars[i]: self.DataScan(fp, n1, n2, n3, endian,
                                                           dtype)})
        
        
        fp.close()

    def ReadMultipleFiles(self, nstr, dataext, myvars, n1, n2, n3, endian,
                          dtype, ddict):
        """Reads a  multiple data files, varname.****.dataext.
	
        **Inputs**:
	  
          nstr -- File number in form of a string\n
	  dataext -- Data type of the file, e.g., 'dbl', 'flt' or 'vtk' \n
          myvars -- List of variable names to be read\n
          n1 -- No. of points in X1 direction\n
          n2 -- No. of points in X2 direction\n
          n3 -- No. of points in X3 direction\n
          endian -- Endianess of the data\n
          dtype -- datatype\n
          ddict -- Dictionary containing Grid and Time Information
          which is updated.
	 
        **Output**:
          
          Updated Dictionary consisting of variable names as keys and its values.
	
	"""
        for i in range(len(myvars)):
            datafilename = self.wdir+myvars[i]+"."+nstr+dataext
            fp = open(datafilename, "rb")
            if self.datatype == 'vtk':
                ddict.update(self.DataScanVTK(fp, n1, n2, n3, endian, dtype))
            else:
                ddict.update({myvars[i]: self.DataScan(fp, n1, n2, n3, endian,
                                                       dtype)})
            fp.close()

    def ReadDataFile(self, num):
        """Reads the data file generated from PLUTO code.

	**Inputs**:
	
	  num -- Data file number in form of an Integer.

        **Outputs**:
	
	  Dictionary that contains all information about Grid, Time and 
	  variables.

	"""
        gridfile = self.wdir+"grid.out"
        if self.datatype == "float":
            dtype = "f"
            varfile = self.wdir+"flt.out"
            dataext = ".flt"
        elif self.datatype == "vtk":
            dtype = "f"
            varfile = self.wdir+"vtk.out"
            dataext=".vtk"
        elif self.datatype == 'hdf5':
            dtype = 'd'
            dataext = '.hdf5'
            nstr = num
            varfile = self.wdir+"data."+nstr+dataext
        else:
            dtype = "d"
            varfile = self.wdir+"dbl.out"
            dataext = ".dbl"
        
        self.ReadVarFile(varfile)
        self.ReadGridFile(gridfile)
        self.ReadTimeInfo(varfile)
        nstr = num
        if self.endianess == 'big':
            endian = ">"
        elif self.datatype == 'vtk':
            endian = ">"
        else:
            endian = "<"

        D = [('NStep', self.NStep), ('SimTime', self.SimTime), ('Dt', self.Dt),
             ('n1', self.n1), ('n2', self.n2), ('n3', self.n3),
             ('x1', self.x1), ('x2', self.x2), ('x3', self.x3),
             ('dx1', self.dx1), ('dx2', self.dx2), ('dx3', self.dx3),
             ('endianess', self.endianess), ('datatype', self.datatype),
             ('filetype', self.filetype)]
        ddict = dict(D)

        if self.filetype == "single_file":
            datafilename = self.wdir+"data."+nstr+dataext
            self.ReadSingleFile(datafilename, self.vars, self.n1, self.n2,
                                self.n3, endian, dtype, ddict)
        elif self.filetype == "multiple_files":
            self.ReadMultipleFiles(nstr, dataext, self.vars, self.n1, self.n2,
                                   self.n3, endian, dtype, ddict)
        else:
            print "Wrong file type : CHECK pluto.ini for file type."
            print "Only supported are .dbl, .flt, .vtk, .hdf5"
            sys.exit()

        return ddict
      

class Tools(object):
	"""
	
	This Class has all the functions doing basic mathematical
	operations to the vector or scalar fields.
	It is called after pyPLUTO.pload object is defined.
	
	"""

	def deriv(self,Y,X=None):
		"""
		Calculates the derivative of Y with respect to X.

		**Inputs:**

		  Y : 1-D array to be differentiated.\n
		  X : 1-D array with len(X) = len(Y).\n

		  If X is not specified then by default X is chosen to be an equally spaced array having same number of elements
		  as Y.

                **Outputs:**

                  This returns an 1-D array having the same no. of elements as Y (or X) and contains the values of dY/dX.
		
		"""
		n = len(Y)
		n2 = n-2
		if X==None : X = np.arange(n)
		Xarr = np.asarray(X,dtype='float')
		Yarr = np.asarray(Y,dtype='float')
		x12 = Xarr - np.roll(Xarr,-1)   #x1 - x2
		x01 = np.roll(Xarr,1) - Xarr    #x0 - x1
		x02 = np.roll(Xarr,1) - np.roll(Xarr,-1) #x0 - x2
		DfDx = np.roll(Yarr,1) * (x12 / (x01*x02)) + Yarr * (1./x12 - 1./x01) - np.roll(Yarr,-1) * (x01 / (x02 * x12))
		# Formulae for the first and last points:

		DfDx[0] = Yarr[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - Yarr[1] * x02[1]/(x01[1]*x12[1]) + Yarr[2] * x01[1]/(x02[1]*x12[1])
		DfDx[n-1] = -Yarr[n-3] * x12[n2]/(x01[n2]*x02[n2]) + Yarr[n-2]*x02[n2]/(x01[n2]*x12[n2]) - Yarr[n-1]*(x02[n2]+x12[n2])/(x02[n2]*x12[n2])

		return DfDx
	
	def Grad(self,phi,x1,x2,dx1,dx2,polar=False):
		""" This method calculates the gradient of the 2D scalar phi.

                **Inputs:**

                  phi -- 2D scalar whose gradient is to be determined.\n
		  x1 -- The 'x' array\n
		  x2 -- The 'y' array\n
                  dx1 -- The grid spacing in 'x' direction.\n
                  dx2 -- The grid spacing in 'y' direction.\n
                  polar -- The keyword should be set to True inorder to estimate the Gradient in polar co-ordinates. By default it is set to False.
                
		**Outputs:**

                  This routine outputs a 3D array with shape = (len(x1),len(x2),2), such that [:,:,0] element corresponds to the gradient values of phi wrt to x1 and [:,:,1] are the gradient values of phi wrt to x2.
 
		"""
		(n1, n2) = phi.shape 
		grad_phi = np.zeros(shape=(n1,n2,2))
		h2 = np.ones(shape=(n1,n2))
		if polar == True:
			for j in range(n2):
				h2[:,j] = x1
		
		for i in range(n1):
			scrh1 = phi[i,:]
			grad_phi[i,:,1] = self.deriv(scrh1,x2)/h2[i,:]
		for j in range(n2):
			scrh2 = phi[:,j]
			grad_phi[:,j,0] = self.deriv(scrh2,x1)

		return grad_phi

	def Div(self,u1,u2,x1,x2,dx1,dx2,geometry=None):
		""" This method calculates the divergence of the 2D vector fields u1 and u2.

                **Inputs:**

                  u1 -- 2D vector along x1 whose divergence is to be determined.\n
                  u2 -- 2D vector along x2 whose divergence is to be determined.\n
		  x1 -- The 'x' array\n
		  x2 -- The 'y' array\n
                  dx1 -- The grid spacing in 'x' direction.\n
                  dx2 -- The grid spacing in 'y' direction.\n
                  geometry -- The keyword *geometry* is by default set to 'cartesian'. It can be set to either one of the following : *cartesian*, *cylindrical*, *spherical* or *polar*. To calculate the divergence of the vector fields, respective geometric corrections are taken into account based on the value of this keyword.

                **Outputs:**

		  A 2D array with same shape as u1(or u2) having the values of divergence.

		"""
		(n1, n2) = u1.shape
		Divergence = np.zeros(shape=(n1,n2))
		du1 = np.zeros(shape=(n1,n2))
		du2 = np.zeros(shape=(n1,n2))

		A1 = np.zeros(shape=n1)
		A2 = np.zeros(shape=n2)

		dV1 = np.zeros(shape=(n1,n2))
		dV2 = np.zeros(shape=(n1,n2))

		if geometry == None : geometry = 'cartesian'
		
		#------------------------------------------------
		#  define area and volume elements for the
		#  different coordinate systems
		#------------------------------------------------

		if geometry == 'cartesian' :
			A1[:] = 1.0
			A2[:] = 1.0
			dV1   = np.outer(dx1,A2)
			dV2   = np.outer(A1,dx2)

		if geometry == 'cylindrical' :
			A1 = x1
			A2[:] = 1.0
			dV1 = np.meshgrid(x1*dx1,A2)[0].T*np.meshgrid(x1*dx1,A2)[1].T
			for i in range(n1) : dV2[i,:] = dx2[:]
		
		if geometry == 'polar' :
			A1    = x1
			A2[:] = 1.0
			dV1   = np.meshgrid(x1,A2)[0].T*np.meshgrid(x1,A2)[1].T
			dV2   = np.meshgrid(x1,dx2)[0].T*np.meshgrid(x1,dx2)[1].T

		if geometry == 'spherical' :
			A1 = x1*x1
			A2 = np.sin(x2)
			for j in range(n2): dV1[:,j] = A1*dx1
			dV2   = np.meshgrid(x1,np.sin(x2)*dx2)[0].T*np.meshgrid(x1,np.sin(x2)*dx2)[1].T

		# ------------------------------------------------
		#              Make divergence
		# ------------------------------------------------
		
		
		for i in range(1,n1-1):
			du1[i,:] = 0.5*(A1[i+1]*u1[i+1,:] - A1[i-1]*u1[i-1,:])/dV1[i,:]
		for j in range(1,n2-1):
			du2[:,j] = 0.5*(A2[j+1]*u2[:,j+1] - A2[j-1]*u2[:,j-1])/dV2[:,j]

		Divergence = du1 + du2
		return Divergence



	def RTh2Cyl(self,R,Th,X1,X2):
		""" This method does the transformation from spherical coordinates to cylindrical ones.
		
                **Inputs:**

                  R - 2D array of spherical radius coordinates.\n
	          Th - 2D array of spherical theta-angle coordinates.\n
		  X1 - 2D array of radial component of given vector\n
		  X2 - 2D array of thetoidal component of given vector\n
                
		**Outputs:**

                  This routine outputs two 2D arrays after transformation.
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
                  ``TH,R=np.meshgrid(D.x2,D.x1)``\n
		  ``Br,Bz=ppt.RTh2Cyl(R,TH,D.bx1,D.bx2)``

                D.bx1 and D.bx2 should be vectors in spherical coordinates. After transformation (Br,Bz) corresponds to vector in cilindrical coordinates.

 
		"""
		Y1=X1*np.sin(Th)+X2*np.cos(Th)
		Y2=X1*np.cos(Th)-X2*np.sin(Th)
		return Y1,Y2




	def myInterpol(self,RR,N):
		""" This method interpolates (linear interpolation) vector 1D vector RR to 1D N-length vector. Useful for stretched grid calculations. 
		
                **Inputs:**

                  RR - 1D array to interpolate.\n
		  N  - Number of grids to interpolate to.\n
				  
		**Outputs:**

                  This routine outputs interpolated 1D array to the new grid (len=N).
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
                  ``x=linspace(0,1,10) #len(x)=10``\n
		  ``y=x*x``\n
		  ``Ri,Ni=ppt.myInterpol(y,100) #len(Ri)=100``

                  Ri - interpolated numbers;
		  Ni - grid for Ri
 
		"""	
		
		NN=np.linspace(0,len(RR)-1,len(RR))
		spline_fit=UnivariateSpline(RR,NN,k=3,s=0)
		
		RRi=np.linspace(RR[0],RR[-1],N)
		NNi=spline_fit(RRi)
		NNi[0]=NN[0]+0.00001
		NNi[-1]=NN[-1]-0.00001
		return RRi,NNi
		
	def getUniformGrid(self,r,th,rho,Nr,Nth):
		""" This method transforms data with non-uniform grid (stretched) to uniform. Useful for stretched grid calculations. 
		
                **Inputs:**

		  r  - 1D vector of X1 coordinate (could be any, e.g D.x1).\n
		  th - 1D vector of X2 coordinate (could be any, e.g D.x2).\n
		  rho- 2D array of data.\n
		  Nr - new size of X1 vector.\n
		  Nth- new size of X2 vector.\n
				  
		**Outputs:**

                  This routine outputs 2D uniform array Nr x Nth dimension
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
		  ``X1new, X2new, res = ppt.getUniformGrid(D.x1,D.x2,D.rho,20,30)``

                  X1new - X1 interpolated grid len(X1new)=20
		  X2new - X2 interpolated grid len(X2new)=30
		  res   - 2D array of interpolated variable
 
		"""	

		Ri,NRi=self.myInterpol(r,Nr)
		Ra=np.int32(NRi);Wr=NRi-Ra

		YY=np.ones([Nr,len(th)])
		for i in range(len(th)):
		      YY[:,i]=(1-Wr)*rho[Ra,i] + Wr*rho[Ra+1,i]

		THi,NTHi=self.myInterpol(th,Nth)
		THa=np.int32(NTHi);Wth=NTHi-THa

		ZZ=np.ones([Nr,Nth])
		for i in range(Nr):
		      ZZ[i,:]=(1-Wth)*YY[i,THa] + Wth*YY[i,THa+1]

		return Ri,THi,ZZ
	
        def sph2cyl(self,D,Dx,rphi=None,theta0=None):
		""" This method transforms spherical data into cylindrical applying interpolation. Works for stretched grid as well, transforms poloidal (R-Theta) data by default. Fix theta and set rphi=True to get (R-Phi) transformation.
				
                **Inputs:**

                  D  - structure  from 'pload' method.\n
		  Dx - variable to be transformed (D.rho for example).\n
				  
		**Outputs:**

                  This routine outputs transformed (sph->cyl) variable and grid.
				  
                **Usage:**

                  ``import pyPLUTO as pp``\n
		  ``import numpy as np``\n
                  ``D = pp.pload(0)``\n
		  ``ppt=pp.Tools()``\n
		  ``R,Z,res = ppt.sph2cyl(D,D.rho.transpose())``

                  R - 2D array with cylindrical radius values
		  Z - 2D array with cylindrical Z values
		  res - 2D array of transformed variable
 
		"""	
		
		if rphi is None or rphi == False:
		    rx=D.x1
		    th=D.x2		    
		else:
                    rx=D.x1*np.sin(theta0)
		    th=D.x3
		    
		rx,th,Dx=self.getUniformGrid(rx,th,Dx.T,200,200)
		Dx=Dx.T
		
		if rphi is None or rphi == False:
                    
                    r0=np.min(np.sin(th)*rx[0])
                    rN=rx[-1]
                    dr=rN-r0
                    z0=np.min(np.cos(th)*rN)
                    zN=np.max(np.cos(th)*rN)
                    dz=zN-z0
                    dth=th[-1]-th[0]
                    rl=np.int32(len(rx)*dr/(rx[-1]-rx[0]))  
                    zl=np.int32(rl* dz/dr)
                    thl=len(th)
                    r=np.linspace(r0, rN, rl)
                    z=np.linspace(z0, zN, zl)
		else:
                    r0=np.min([np.sin(th)*rx[0] , np.sin(th)*rx[-1]])
                    rN=np.max([np.sin(th)*rx[0] , np.sin(th)*rx[-1]])
                    dr=rN-r0
                    z0=np.min(np.cos(th)*rN)
                    zN=np.max(np.cos(th)*rN)
                    dz=zN-z0
                    dth=th[-1]-th[0]
                    rl=np.int32(len(rx)*dr/(rx[-1]-rx[0]))  
                    zl=np.int32(rl* dz/dr)
                    thl=len(th)
                    r=np.linspace(r0, rN, rl)
                    z=np.linspace(z0, zN, zl)
                
                R,Z = np.meshgrid(r, z)
		Rs = np.sqrt(R*R + Z*Z)
		
		
		Th = np.arccos(Z/Rs)
		kv_34=find(R<0)
		Th.flat[kv_34]=2*np.pi - Th.flat[kv_34]
		
		
		ddr=rx[1]-rx[0]
		ddth=th[1]-th[0]
		
		Rs_copy=Rs.copy()
		Th_copy=Th.copy()
				
		nR1=find(Rs<rx[0])  
		Rs.flat[nR1]=rx[0] 
		nR2=find(Rs>rN)
		Rs.flat[nR2]=rN
		
		nTh1=find(Th>th[-1])
		Th.flat[nTh1]=th[-1]
		nTh2=find(Th<th[0])
		Th.flat[nTh2]=th[0]
		
		
		ra = ((len(rx)-1.001)/(np.max(Rs.flat)-np.min(Rs.flat)) *(Rs-np.min(Rs.flat)))  
		tha = ((thl-1.001)/dth *(Th-th[0]))  

		rn = np.int32(ra)
		thn = np.int32(tha)
		dra=ra-rn
		dtha=tha-thn
		w1=1-dra
		w2=dra
		w3=1-dtha
		w4=dtha
		lrx=len(rx)
		NN1=np.int32(rn+thn*lrx)
		NN2=np.int32((rn+1)+thn*lrx)
		NN3=np.int32(rn+(thn+1)*lrx)
		NN4=np.int32((rn+1)+(thn+1)*lrx)
		n=np.transpose(np.arange(0,np.product(np.shape(R))))
		DD=Dx.copy()
		F=R.copy()
		F.flat[n]=w1.flat[n]*(w3.flat[n]*Dx.flat[NN1.flat[n]] + w4.flat[n]*Dx.flat[NN3.flat[n]] )+\
		    w2.flat[n]*(w3.flat[n]*Dx.flat[NN2.flat[n]] + w4.flat[n]*Dx.flat[NN4.flat[n]] )
		    
		nR1=find(Rs_copy<rx[0]-ddr/1.5)
		nR2=find(Rs_copy>rN+ddr/1.5)
		nTh1=find(Th_copy>th[-1]+ddth/1.5)
		nTh2=find(Th_copy<th[0]-ddth/1.5)

		nmask=np.concatenate((nR1,nR2,nTh1,nTh2))
		F.flat[nmask]=np.nan;
		return R,Z,F
        

	def congrid(self, a, newdims, method='linear', centre=False, minusone=False):
	    """
	    Arbitrary resampling of source array to new dimension sizes.
	    Currently only supports maintaining the same number of dimensions.
	    To use 1-D arrays, first promote them to shape (x,1).

	    Uses the same parameters and creates the same co-ordinate lookup points
	    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
	    routine of the same name.

	    **Inputs:**

	      a -- The 2D array that needs resampling into new dimensions.\n
	      newdims -- A tuple which represents the shape of the resampled data\n
	      method -- This keyword decides the method used for interpolation.\n
	        neighbour - closest value from original data\n
	        nearest and linear - uses n x 1-D interpolations using scipy.interpolate.interp1d
	        (see Numerical Recipes for validity of use of n 1-D interpolations)\n
	        spline - uses ndimage.map_coordinates\n
	      centre -- This keyword decides the positions of interpolation points.\n
	        True - interpolation points are at the centres of the bins\n
	        False - points are at the front edge of the bin\n
	      minusone -- This prevents extrapolation one element beyond bounds of input array\n
	        For example- inarray.shape = (i,j) & new dimensions = (x,y)\n
	        False - inarray is resampled by factors of (i/x) * (j/y)\n
	        True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)\n

            **Outputs:**

              A 2D array with resampled data having a shape corresponding to newdims. 
	    
	    """
	    if not a.dtype in [np.float64, np.float32]:
		a = np.cast[float](a)

	    m1 = np.cast[int](minusone)
	    ofs = np.cast[int](centre) * 0.5
	    old = np.array( a.shape )
	    ndims = len( a.shape )
	    if len( newdims ) != ndims:
		print "[congrid] dimensions error. " \
		      "This routine currently only support " \
		      "rebinning to the same number of dimensions."
		return None
	    newdims = np.asarray( newdims, dtype=float )
	    dimlist = []

	    if method == 'neighbour':
		for i in range( ndims ):
		    base = np.indices(newdims)[i]
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		cd = np.array( dimlist ).round().astype(int)
		newa = a[list( cd )]
		return newa

	    elif method in ['nearest','linear']:
		# calculate new dims
		for i in range( ndims ):
		    base = np.arange( newdims[i] )
		    dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
				    * (base + ofs) - ofs )
		# specify old dims
		olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

		# first interpolation - for ndims = any
		mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
		newa = mint( dimlist[-1] )

		trorder = [ndims - 1] + range( ndims - 1 )
		for i in range( ndims - 2, -1, -1 ):
		    newa = newa.transpose( trorder )

		    mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
		    newa = mint( dimlist[i] )

		if ndims > 1:
		    # need one more transpose to return to original dimensions
		    newa = newa.transpose( trorder )

		return newa
	    elif method in ['spline']:
		oslices = [ slice(0,j) for j in old ]
		oldcoords = np.ogrid[oslices]
		nslices = [ slice(0,j) for j in list(newdims) ]
		newcoords = np.mgrid[nslices]

		newcoords_dims = range(np.rank(newcoords))
		#make first index last
		newcoords_dims.append(newcoords_dims.pop(0))
		newcoords_tr = newcoords.transpose(newcoords_dims)
		# makes a view that affects newcoords

		newcoords_tr += ofs

		deltas = (np.asarray(old) - m1) / (newdims - m1)
		newcoords_tr *= deltas

		newcoords_tr -= ofs

		newa = scipy.ndimage.map_coordinates(a, newcoords)
		return newa
	    else:
		print "Congrid error: Unrecognized interpolation type.\n", \
		      "Currently only \'neighbour\', \'nearest\',\'linear\',", \
		      "and \'spline\' are supported."
		return None


        








class Image(object):
    ''' This Class has all the routines for the imaging the data
	and plotting various contours and fieldlines on these images.
	CALLED AFTER pyPLUTO.pload object is defined
	'''
    def pldisplay(self, D, var,**kwargs):
        """ This method allows the user to display a 2D data using the 
        matplotlib's pcolormesh.

        **Inputs:**

          D   -- pyPLUTO pload object.\n
          var -- 2D array that needs to be displayed.
        
        *Required Keywords:*

          x1 -- The 'x' array\n
          x2 -- The 'y' array
        
        *Optional Keywords:*

          vmin -- The minimum value of the 2D array (Default : min(var))\n
          vmax -- The maximum value of the 2D array (Default : max(var))\n
          title -- Sets the title of the image.\n
          label1 -- Sets the X Label (Default: 'XLabel')\n
          label2 -- Sets the Y Label (Default: 'YLabel')\n
          polar -- A list to project Polar data on Cartesian Grid.\n
            polar = [True, True] -- Projects r-phi plane.\n
            polar = [True, False] -- Project r-theta plane.\n
            polar = [False, False] -- No polar plot (Default)\n
          cbar -- Its a tuple to set the colorbar on or off. \n
            cbar = (True,'vertical') -- Displays a vertical colorbar\n
            cbar = (True,'horizontal') -- Displays a horizontal colorbar\n
            cbar = (False,'') -- Displays no colorbar.
         
        **Usage:**
          
          ``import pyPLUTO as pp``\n
          ``wdir = '/path/to/the data files/'``\n
          ``D = pp.pload(1,w_dir=wdir)``\n
          ``I = pp.Image()``\n
          ``I.pldisplay(D, D.v2, x1=D.x1, x2=D.x2, cbar=(True,'vertical'),\
          title='Velocity',label1='Radius',label2='Height')``
        """
        x1 = kwargs.get('x1')
        x2 = kwargs.get('x2')
        var = var.T

        f1 = figure(kwargs.get('fignum',1), figsize=kwargs.get('figsize',[10,10]),
                    dpi=80, facecolor='w', edgecolor='k')
        ax1 = f1.add_subplot(111)
        ax1.set_aspect('equal')        

        if kwargs.get('polar',[False,False])[0]:
            xx, yy = self.getPolarData(D,kwargs.get('x2'),rphi=kwargs.get('polar')[1])
            pcolormesh(xx,yy,var,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        else:
            ax1.axis([np.min(x1),np.max(x1),np.min(x2),np.max(x2)])
            pcolormesh(x1,x2,var,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        
        title(kwargs.get('title',"Title"),size=kwargs.get('size'))
        xlabel(kwargs.get('label1',"Xlabel"),size=kwargs.get('size'))
        ylabel(kwargs.get('label2',"Ylabel"),size=kwargs.get('size'))
        if kwargs.get('cbar',(False,''))[0] == True:
            colorbar(orientation=kwargs.get('cbar')[1])
    
    
    def multi_disp(self,*args,**kwargs):
        mvar = []
        for arg in args:
            mvar.append(arg.T)    		
        
        xmin = np.min(kwargs.get('x1'))
        xmax = np.max(kwargs.get('x1'))		
        ymin = np.min(kwargs.get('x2'))
        ymax = np.max(kwargs.get('x2'))
        mfig = figure(kwargs.get('fignum',1),figsize=kwargs.get('figsize',[10,10]))
        Ncols = kwargs.get('Ncols')
        Nrows = len(args)/Ncols
        mprod = Nrows*Ncols
        dictcbar=kwargs.get('cbar',(False,'','each'))

        for j in range(mprod):
            mfig.add_subplot(Nrows,Ncols,j+1)
            pcolormesh(kwargs.get('x1'),kwargs.get('x2'), mvar[j])
            axis([xmin,xmax,ymin,ymax])
            gca().set_aspect('equal')
                
            xlabel(kwargs.get('label1',mprod*['Xlabel'])[j])
            ylabel(kwargs.get('label2',mprod*['Ylabel'])[j])
            title(kwargs.get('title',mprod*['Title'])[j])
            if (dictcbar[0] == True) and (dictcbar[2] =='each'):
                colorbar(orientation=kwargs.get('cbar')[1])
            if dictcbar[0] == True and dictcbar[2]=='last':
                if (j == np.max(range(mprod))):colorbar(orientation=kwargs.get('cbar')[1])
         
    def oplotbox(self, AMRLevel, lrange=[0,0], cval=['b','r','g','m','w','k'],\
                     islice=-1, jslice=-1, kslice=-1,geom='CARTESIAN'):
        """ 
        This method overplots the AMR boxes up to the specified level. 

        **Input:**

          AMRLevel -- AMR object loaded during the reading and stored in the pload object
        
        *Optional Keywords:*

          lrange     -- [level_min,level_max] to be overplotted. By default it shows all the loaded levels\n
          cval       -- list of colors for the levels to be overplotted.\n
          [ijk]slice -- Index of the 2D slice to look for so that the adequate box limits are plotted. 
                        By default oplotbox considers you are plotting a 2D slice of the z=min(x3) plane.\n
          geom       -- Specified the geometry. Currently, CARTESIAN (default) and POLAR geometries are handled.
        """

        nlev = len(AMRLevel)
        lrange[1] = min(lrange[1],nlev-1)
        npl  = lrange[1]-lrange[0]+1
        lpls = [lrange[0]+v for v in range(npl)]  
        cols = cval[0:nlev]
        # Get the offset and the type of slice
        Slice = 0 ; inds = 'k'
        xx = 'x' ; yy ='y'
        if (islice >= 0):
            Slice = islice + AMRLevel[0]['ibeg'] ; inds = 'i'
            xx = 'y' ; yy ='z'
        if (jslice >= 0):
            Slice = jslice + AMRLevel[0]['jbeg'] ; inds = 'j'
            xx = 'x' ; yy ='z'
        if (kslice >= 0):
            Slice = kslice + AMRLevel[0]['kbeg'] ; inds = 'k'
            xx = 'x' ; yy ='y'
            
        # Overplot the boxes
        hold(True)
        for il in lpls:
            level = AMRLevel[il]
            for ib in range(level['nbox']):
                box = level['box'][ib]
                if ((Slice-box[inds+'b'])*(box[inds+'e']-Slice) >= 0):
                    if (geom == 'CARTESIAN'):
                        x0 = box[xx+'0'] ; x1 = box[xx+'1'] 
                        y0 = box[yy+'0'] ; y1 = box[yy+'1'] 
                        plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],color=cols[il])
                    elif (geom == 'POLAR') or (geom == 'SPHERICAL'):
                        dn = np.pi/50.
                        x0 = box[xx+'0'] ; x1 = box[xx+'1'] 
                        y0 = box[yy+'0'] ; y1 = box[yy+'1']
                        if y0 == y1:
                            y1 = 2*np.pi+y0 - 1.e-3
                        xb = np.concatenate([
                                [x0*np.cos(y0),x1*np.cos(y0)],\
                                x1*np.cos(np.linspace(y0,y1,num=int(abs(y0-y1)/dn) )),\
                                [x1*np.cos(y1),x0*np.cos(y1)],\
                                x0*np.cos(np.linspace(y1,y0,num=int(abs(y0-y1)/dn)))])
                        yb = np.concatenate([
                                [x0*np.sin(y0),x1*np.sin(y0)],\
                                x1*np.sin(np.linspace(y0,y1,num=int(abs(y0-y1)/dn))),\
                                [x1*np.sin(y1),x0*np.sin(y1)],\
                                x0*np.sin(np.linspace(y1,y0,num=int(abs(y0-y1)/dn)))])
                        plot(xb,yb,color=cols[il])

        hold(False)
         
    def field_interp(self,var1,var2,x,y,dx,dy,xp,yp):
        """ This method interpolates value of vector fields (var1 and var2) on field points (xp and yp).
        The field points are obtained from the method field_line.

        **Inputs:**
		
          var1 -- 2D Vector field in X direction\n
          var2 -- 2D Vector field in Y direction\n
          x -- 1D X array\n
          y -- 1D Y array\n
          dx -- 1D grid spacing array in X direction\n
          dy -- 1D grid spacing array in Y direction\n
          xp -- field point in X direction\n
          yp -- field point in Y direction\n
          
        **Outputs:**
        
          A list with 2 elements where the first element corresponds to the interpolate field 
          point in 'x' direction and the second element is the field point in 'y' direction.  

        """
        q=[]
        U = var1
        V = var2
        i0 = np.abs(xp-x).argmin()
        j0 = np.abs(yp-y).argmin()
        scrhUx = np.interp(xp,x,U[:,j0])
        scrhUy = np.interp(yp,y,U[i0,:])
        q.append(scrhUx + scrhUy - U[i0,j0])
        scrhVx = np.interp(xp,x,V[:,j0])
        scrhVy = np.interp(yp,y,V[i0,:])
        q.append(scrhVx + scrhVy - V[i0,j0])
        return q
    
    def field_line(self,var1,var2,x,y,dx,dy,x0,y0):
        """ This method is used to obtain field lines (same as fieldline.pro in PLUTO IDL tools).
    
        **Inputs:**
    
          var1 -- 2D Vector field in X direction\n
          var2 -- 2D Vector field in Y direction\n
          x -- 1D X array\n
          y -- 1D Y array\n
          dx -- 1D grid spacing array in X direction\n
          dy -- 1D grid spacing array in Y direction\n
          x0 -- foot point of the field line in X direction\n
          y0 -- foot point of the field line in Y direction\n
          
        **Outputs:**
    
          This routine returns a dictionary with keys - \n
          qx -- list of the field points along the 'x' direction.
          qy -- list of the field points along the 'y' direction.
          
        **Usage:**
          
          See the myfieldlines routine for the same.          
        """
        xbeg = x[0] - 0.5*dx[0]
        xend = x[-1] + 0.5*dx[-1]
        
        ybeg = y[0]  - 0.5*dy[0]
        yend = y[-1] + 0.5*dy[-1]
    
        inside_domain = x0 > xbeg and x0 < xend and y0 > ybeg and y0 < yend
    
        MAX_STEPS = 5000
        xln_fwd = [x0]
        yln_fwd = [y0]
        xln_bck = [x0]
        yln_bck = [y0]
        rhs = []
        k = 0
    
        while inside_domain == True:
            R1 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])
            dl = 0.5*np.max(np.concatenate((dx,dy)))/(np.sqrt(R1[0]*R1[0] + R1[1]*R1[1] + 1.e-14))
            xscrh = xln_fwd[k] + 0.5*dl*R1[0]
            yscrh = yln_fwd[k] + 0.5*dl*R1[1]
            
            R2 = self.field_interp(var1,var2,x,y,dx,dy,xln_fwd[k],yln_fwd[k])
            
            xln_one = xln_fwd[k] + dl*R2[0]
            yln_one = yln_fwd[k] + dl*R2[1]
            
            xln_fwd.append(xln_one)
            yln_fwd.append(yln_one)
            inside_domain = xln_one > xbeg and xln_one < xend and yln_one > ybeg and yln_one < yend
            inside_domain = inside_domain and (k < MAX_STEPS-3)
            k = k + 1
    
    
        k_fwd = k
        qx = np.asarray(xln_fwd[0:k_fwd])
        qy = np.asarray(yln_fwd[0:k_fwd])
        flines={'qx':qx,'qy':qy}
        
        return flines
      
      
    def myfieldlines(self,Data,x0arr,y0arr,stream=False,**kwargs):
        """ This method overplots the magnetic field lines at the footpoints given by (x0arr[i],y0arr[i]).
        
        **Inputs:**
        
          Data -- pyPLUTO.pload object\n
          x0arr -- array of x co-ordinates of the footpoints\n
          y0arr -- array of y co-ordinates of the footpoints\n
          stream -- keyword for two different ways of calculating the field lines.\n
          True -- plots contours of rAphi (needs to store vector potential)\n
          False -- plots the fieldlines obtained from the field_line routine. (Default option)
          
        *Optional Keywords:*
		
          colors -- A list of matplotlib colors to represent the lines. The length of this list should be same as that of x0arr.\n
          lw -- Integer value that determines the linewidth of each line.\n
          ls -- Determines the linestyle of each line.

        **Usage:**
          
          Assume that the magnetic field is a given as **B** = B0$\hat{y}$. 
          Then to show this field lines we have to define the x and y arrays of field foot points.\n 
          
          ``x0arr = linspace(0.0,10.0,20)``\n
          ``y0arr = linspace(0.0,0.0,20)``\n
          ``import pyPLUTO as pp``\n
          ``D = pp.pload(45)``\n
          ``I = pp.Image()``\n
          ``I.myfieldlines(D,x0arr,y0arr,colors='k',ls='--',lw=1.0)``
        """
	       
        if len(x0arr) != len(y0arr) : print "Input Arrays should have same size"
        QxList=[]
        QyList=[]
        StreamFunction = []
        levels =[]
        if stream == True:
            X, Y = np.meshgrid(Data.x1,Data.x2.T)
            StreamFunction = X*(Data.Ax3.T)
            for i in range(len(x0arr)):
                nx = np.abs(X[0,:]-x0arr[i]).argmin()
                ny = np.abs(X[:,0]-y0arr[i]).argmin()
                levels.append(X[ny,nx]*Data.Ax3.T[ny,nx])
			
            contour(X,Y,StreamFunction,levels,colors=kwargs.get('colors'),linewidths=kwargs.get('lw',1),linestyles=kwargs.get('ls','solid'))
        else:
            for i in range(len(x0arr)):
                QxList.append(self.field_line(Data.bx1,Data.bx2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qx'))
                QyList.append(self.field_line(Data.bx1,Data.bx2,Data.x1,Data.x2,Data.dx1,Data.dx1,x0arr[i],y0arr[i]).get('qy'))
                plot(QxList[i],QyList[i],color=kwargs.get('colors'))
                axis([min(Data.x1),max(Data.x1),min(Data.x2),max(Data.x2)])

    def getSphData(self,Data,w_dir=None,datatype=None,**kwargs):
        """This method transforms the vector and scalar  fields from Spherical co-ordinates to Cylindrical.

        **Inputs**:
        
          Data -- pyPLUTO.pload object\n
          w_dir -- /path/to/the/working/directory/\n
          datatype -- If the data is of 'float' type then datatype = 'float' else by default the datatype is set to 'double'.

        *Optional Keywords*:
	    
	  rphi -- [Default] is set to False implies that the r-theta plane is transformed. If set True then the r-phi plane is transformed.\n
          x2cut -- Applicable for 3D data and it determines the co-ordinate of the x2 plane while r-phi is set to True.\n
          x3cut -- Applicable for 3D data and it determines the co-ordinate of the x3 plane while r-phi is set to False.
        
        """
 
        Tool = Tools()
        key_value_pairs = []
        allvars = []
        if w_dir is None: w_dir = curdir()		    
        for v in Data.vars:
            allvars.append(v)
            
        if kwargs.get('rphi',False)==True:
            R,TH = np.meshgrid(Data.x1,Data.x3)
            if Data.n3 != 1:
                for variable in allvars:
                    key_value_pairs.append([variable,getattr(Data,variable)[:,kwargs.get('x2cut',0),:].T])
			
                SphData = dict(key_value_pairs)
                if ('bx1' in allvars) or ('bx2' in allvars):
                    (SphData['b1c'],SphData['b3c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx3'))
                    allvars.append('b1c')
                    allvars.append('b3c')
                if ('vx1' in allvars) or ('vx2' in allvars):
                    (SphData['v1c'],SphData['v3c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx3'))
                    allvars.append('v1c')
                    allvars.append('v3c')
            else:
                print "No x3 plane for 2D data"
        else:
            R,TH = np.meshgrid(Data.x1,Data.x2)
            if Data.n3 != 1:
                for variable in allvars:
                    key_value_pairs.append([variable,getattr(Data,variable)[:,:,kwargs.get('x3cut',0)].T])
                SphData = dict(key_value_pairs)
                if ('bx1' in allvars) or ('bx2' in allvars):
                    (SphData['b1c'],SphData['b2c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx2'))
                    allvars.append('b1c')
                    allvars.append('b2c')
                if ('vx1' in allvars) or ('vx2' in allvars):
                    (SphData['v1c'],SphData['v2c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx2'))
                    allvars.append('v1c')
                    allvars.append('v2c')
            else:
                for variable in allvars:
                    key_value_pairs.append([variable,getattr(Data,variable)[:,:].T])
                SphData = dict(key_value_pairs)
                if ('bx1' in allvars) or ('bx2' in allvars):
                    (SphData['b1c'],SphData['b2c']) = Tool.RTh2Cyl(R,TH,SphData.get('bx1'),SphData.get('bx2'))
                    allvars.append('b1c')
                    allvars.append('b2c')
                if ('vx1' in allvars) or ('vx2' in allvars):
                    (SphData['v1c'],SphData['v2c']) = Tool.RTh2Cyl(R,TH,SphData.get('vx1'),SphData.get('vx2'))
                    allvars.append('v1c')
                    allvars.append('v2c')
            
        for variable in allvars:
            if kwargs.get('rphi',False)==True:
                R,Z,SphData[variable]= Tool.sph2cyl(Data,SphData.get(variable),rphi=True,theta0=Data.x2[kwargs.get('x2cut',0)])
            else:
                if Data.n3 != 1:
                    R,Z,SphData[variable] = Tool.sph2cyl(Data,SphData.get(variable),rphi=False)
                else:
                    R,Z,SphData[variable] = Tool.sph2cyl(Data,SphData.get(variable),rphi=False)

        return R,Z,SphData
    
    
    def getPolarData(self, Data, ang_coord, rphi=False):
        """To get the Cartesian Co-ordinates from Polar.
        
        **Inputs:**
        
          Data -- pyPLUTO pload Object\n
          ang_coord -- The Angular co-ordinate (theta or Phi)
         
        *Optional Keywords:*
        
          rphi -- Default value FALSE is for R-THETA data, 
          Set TRUE for R-PHI data.\n

        **Outputs**:
        
          2D Arrays of X, Y from the Radius and Angular co-ordinates.\n
          They are used in pcolormesh in the Image.pldisplay functions.
        """
        D = Data
        if ang_coord is D.x2:
            x2r = D.x2r
        elif ang_coord is D.x3:
            x2r = D.x3r
        else:
            print "Angular co-ordinate must be given"
            
        rcos = np.outer(np.cos(x2r), D.x1r)
        rsin = np.outer(np.sin(x2r), D.x1r)        
        
        if rphi:
            xx = rcos
            yy = rsin
        else:
            xx = rsin
            yy = rcos
            
        return xx, yy

    def pltSphData(self,Data,w_dir=None,datatype=None,**kwargs):
        """This method plots the transformed data obtained from getSphData using the matplotlib's imshow
        
        **Inputs:**
        
          Data -- pyPLUTO.pload object\n
          w_dir -- /path/to/the/working/directory/\n
          datatype -- Datatype.

        *Required Keywords*:
          
          plvar -- A string which represents the plot variable.\n

	*Optional Keywords*:

          logvar -- [Default = False] Set it True for plotting the log of a variable.\n
          rphi -- [Default = False - for plotting in r-theta plane] Set it True for plotting the variable in r-phi plane. 

        """
              
        if w_dir is None: w_dir=curdir()
        R,Z,SphData = self.getSphData(Data,w_dir=w_dir,datatype=datatype,**kwargs)
        
        extent=(np.min(R.flat),max(R.flat),np.min(Z.flat),max(Z.flat))
        dRR=max(R.flat)-np.min(R.flat)
        dZZ=max(Z.flat)-np.min(Z.flat)


        isnotnan=-np.isnan(SphData[kwargs.get('plvar')])
        maxPl=max(SphData[kwargs.get('plvar')][isnotnan].flat)
        minPl=np.min(SphData[kwargs.get('plvar')][isnotnan].flat)
        normrange=False
        if minPl<0:
            normrange=True
        if maxPl>-minPl:
            minPl=-maxPl
        else:
            maxPl=-minPl	  
        if (normrange and kwargs.get('plvar')!='rho' and kwargs.get('plvar')!='prs'):
            SphData[kwargs.get('plvar')][-1][-1]=maxPl
            SphData[kwargs.get('plvar')][-1][-2]=minPl
    
        if (kwargs.get('logvar') == True):
            SphData[kwargs.get('plvar')] = np.log10(SphData[kwargs.get('plvar')])
	
        imshow(SphData[kwargs.get('plvar')], aspect='equal', origin='lower', cmap=cm.jet,extent=extent, interpolation='nearest')



	

	
		
		



	


      




		    

	    
    
	    
	    
	    
	    
	    
	    



	
        
	  






    
         
    







    
   
