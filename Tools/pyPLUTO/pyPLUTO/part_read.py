import sys
import os
import numpy as np
import array

class ploadparticles(object):
    def __init__(self, ns, wdir=None, datatype=None):
        self.ns = ns
        if wdir == None:
            self.wdir = os.getcwd()+'/'
        else:
            self.wdir = wdir
        
        if datatype == None:
            datatype = 'dbl'
            self.extsize = 8

        self.ext = datatype
        if self.ext == 'flt': self.extsize = 4
        if self.ext == 'vtk': self.extsize = 4

        self.Header = {}
        try:
            self.GetParticleDataFileName()
        except (IOError):
            print "!ERROR! - Particle Data file not found : %s"%self.fname
        else:
            if self.ext == 'vtk':
                Part_Data_dict = self.GetParticleDataVTK()
            else:
                self.ReadParticleFile()
                Part_Data_dict = self.GetParticleData()
            
            self.partdata_keys = []
            for keys in Part_Data_dict:
                self.partdata_keys.append(keys)
                object.__setattr__(self, keys, Part_Data_dict.get(keys))

    def GetParticleDataFileName(self):
        nsstr_ = str(self.ns)
        while len(nsstr_) < 4:
            nsstr_ = '0'+ nsstr_
        self.fname = self.wdir+"particles."+nsstr_+"."+self.ext


    def GetParticleFileHeader(self, line_):
        if line_.split()[1] != 'PLUTO':
            hlist = line_.split()[1:]
            if hlist[0] == 'field_names':
                vars_ = hlist[1:]
                t = tuple([hlist[0], vars_])
                self.Header.update(dict([t]))
            elif hlist[0] == 'field_dim':
                varsdim_ = hlist[1:]
                t = tuple([hlist[0], varsdim_])
                self.Header.update(dict([t]))
            else:
                t = tuple([''.join(hlist[:-1]), hlist[-1]])
                self.Header.update(dict([t]))
            
    def ReadParticleFile(self):
        print "Reading Particle Data file : %s"%self.fname
        fp_ = open(self.fname,'rb')
        for l in fp_.readlines():
            if l.split()[0] == '#':
            	self.GetParticleFileHeader(l)
            else:
            	break
        fp_.close()
        
        self.tot_fdim = np.sum(np.array(self.Header['field_dim'], dtype=np.int))
        HeadBytes_ = int(self.Header['nparticles'])*self.tot_fdim*self.extsize
        fp_ = open(self.fname,'rb')
        scrh_ = fp_.read()
        self.datastr = scrh_[len(scrh_)-HeadBytes_:]
        
        fp_.close()

    def GetParticleDataVTK(self):
        print "Reading Particle Data file : %s"%self.fname
        fp_ = open(self.fname,'rb')
        ks_ = []
        vtk_scalars_ = []
        while True:
            l = fp_.readline()
            try:
                l.split()[0]
            except IndexError:
                pass
            else:
                if l.split()[0] == "POINTS":
                    np_ = int(l.split()[1])
                    endianess = '>'
                    dtyp_ = np.dtype(endianess+'%df'%(3*np_))
                    nb_ = np.dtype(dtyp_).itemsize
                    posbuf_ = array.array('f')
                    posbuf_.fromstring(fp_.read(nb_))
                    pos_ = np.frombuffer(posbuf_, dtype=np.dtype(dtyp_))
                    
                    pos1_ = np.zeros(np_)
                    pos2_ = np.zeros(np_)
                    pos3_ = np.zeros(np_)
                    for i in range(np_):
                        pos1_[i] = pos_[0][3*i]
                        pos2_[i] = pos_[0][3*i+1]
                        pos3_[i] = pos_[0][3*i+2]

                elif l.split()[0] == "VECTORS":
                    endianess = '>'
                    dtyp_ = np.dtype(endianess+'%df'%(3*np_))
                    nb_ = np.dtype(dtyp_).itemsize
                    velbuf_ = array.array('f')
                    velbuf_.fromstring(fp_.read(nb_))
                    vel_ = np.frombuffer(velbuf_, dtype=np.dtype(dtyp_))

                    vel1_ = np.zeros(np_)
                    vel2_ = np.zeros(np_)
                    vel3_ = np.zeros(np_)
                    for i in range(np_):
                        vel1_[i] = vel_[0][3*i]
                        vel2_[i] = vel_[0][3*i+1]
                        vel3_[i] = vel_[0][3*i+2]

                elif l.split()[0] == "SCALARS":
                    ks_.append(l.split()[1])
                elif l.split()[0] == "LOOKUP_TABLE":
                    endianess = '>'
                    if ks_[-1] == "Identity":
                        dtyp_ = np.dtype(endianess+'%di'%np_)
                    else:
                        dtyp_ = np.dtype(endianess+'%df'%np_) 
                        
                    nb_ = np.dtype(dtyp_).itemsize
                    scalbuf_ = array.array('f')
                    scalbuf_.fromstring(fp_.read(nb_))
                    scal_ = np.frombuffer(scalbuf_, dtype=np.dtype(dtyp_))
                    
                    vtk_scalars_.append(scal_[0])
                else:
                    pass

            if l == '':
                break
        ks_ = ['id', 'color', 'tinj']
        vtkvardict = dict(zip(ks_,vtk_scalars_))
        tup_ = [('x1', pos1_),('x2', pos2_),('x3', pos3_),\
                ('vx1', vel1_),('vx2', vel2_),('vx3', vel3_)]
        vtkvardict.update(dict(tup_))
        return vtkvardict


    def GetParticleData(self):
        vars_ = self.Header['field_names']
        endianess = '<'
        if self.Header['endianity'] == 'big':endianess = '>'
        dtyp_ = np.dtype(endianess+self.ext[0])
        DataDict_ = self.Header
        np_ = int(DataDict_['nparticles'])
        data_ = np.fromstring(self.datastr,dtype=dtyp_)

        fdims = np.array(self.Header['field_dim'], dtype=np.int)
        indx = np.where(fdims == 1)[0]
        spl_cnt = len(indx)
        counter = 0
        
        A = np.array_split(data_, [self.tot_fdim*i for i in range(1,np_)])
        
        tup_ = []
        arr_ = np.zeros([np_, len(vars_)])
        
        for s in range(spl_cnt):
            for i in range(len(A)):
                arr_[i,s] = A[i][s]
            tup_.append((vars_[s], arr_[:,s]))
    
        try:
            scrh = vars_.index('vL')
        except ValueError:
            nflx_ = 0
            pass
        else:
            nflx_ = fdims[scrh]
            shkvL_ = np.zeros((np_, nflx_))
            shkvR_ = np.zeros((np_, nflx_))
    
            for n in range(np_):
                shkvL_[n,:] = A[n][spl_cnt:spl_cnt+nflx_]
                shkvR_[n,:] = A[n][spl_cnt+nflx_:spl_cnt+2*nflx_]

            tup_.append((vars_[spl_cnt+counter], shkvL_))
            counter += 1
            tup_.append((vars_[spl_cnt+counter], shkvR_))

        try:
            scrh = vars_.index('eng')
        except ValueError:
            pass
        else:
            nebins_ = fdims[scrh]
            eng_ = np.zeros((np_,nebins_))
            chi_ = np.zeros((np_,nebins_))

            for n in range(np_):
                eng_[n,:] = A[n][spl_cnt+2*nflx_:spl_cnt+2*nflx_+nebins_]
                chi_[n,:] = A[n][spl_cnt+2*nflx_ + nebins_:spl_cnt+2*nflx_ + 2*nebins_]

            if (nflx_ != 0): counter += 1
            tup_.append((vars_[spl_cnt+counter], eng_))
            counter += 1
            tup_.append((vars_[spl_cnt+counter], chi_))

        DataDict_.update(dict(tup_))
        return DataDict_





