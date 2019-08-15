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
                vars_ = hlist[1:]#[v.strip(' ') for v in hlist[1:]]
                t = tuple([hlist[0], vars_])
                self.Header.update(dict([t]))
            elif hlist[0] == 'field_dim':
                varsdim_ =hlist[1:]
                t = tuple([hlist[0], varsdim_])
                self.Header.update(dict([t]))
                self.Header['nflux'] = int(varsdim_[-3])
                self.Header['ebins'] = int(varsdim_[-1])
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
        
        cnt_ = 9
        HeadBytes_ = int(self.Header['nparticles'])*cnt_*self.extsize
        if int(self.Header['spectrum']):
            cnt_ += 6+2*int(self.Header['nflux'])+ 2*int(self.Header['ebins'])
            HeadBytes_ = int(self.Header['nparticles'])*cnt_*self.extsize

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

        spectrum_ = int(DataDict_['spectrum'])
        ids_   = np.zeros(np_)
        pos1_  = np.zeros(np_)
        pos2_  = np.zeros(np_)
        pos3_  = np.zeros(np_)
        vel1_  = np.zeros(np_)
        vel2_  = np.zeros(np_)
        vel3_  = np.zeros(np_)
        tinj_  = np.zeros(np_)
        col_   = np.zeros(np_)
        
        A = np.array_split(data_, [9*i for i in range(1,np_)])
        
        
        if spectrum_:
            nmic_      = np.zeros(np_)
            cmpratio_  = np.zeros(np_)
            shkflg_    = np.zeros(np_)
            shkgradp_  = np.zeros(np_)
            a_         = np.zeros(np_)
            c_         = np.zeros(np_)
            A = np.array_split(data_, [10*i for i in range(1,np_)])
            A = np.array_split(data_, [11*i for i in range(1,np_)])
            A = np.array_split(data_, [12*i for i in range(1,np_)])
            A = np.array_split(data_, [13*i for i in range(1,np_)])
            A = np.array_split(data_, [14*i for i in range(1,np_)])

            ebins_ = int(DataDict_['ebins'])
            nflx_  = int(DataDict_['nflux'])
                     
            shkvL_ = np.zeros((np_, nflx_))
            shkvR_ = np.zeros((np_, nflx_))

            specE_ = np.zeros((np_,ebins_))
            specN_ = np.zeros((np_,ebins_))
            A = np.array_split(data_, [(15+2*nflx_)*i for i in range(1,np_)])
            A = np.array_split(data_, [(15+2*nflx_+2*ebins_)*i for i in range(1,np_)])
            
    
        for n in range(np_):
            ids_[n]   = A[n][0]
            pos1_[n] = A[n][1]
            pos2_[n] = A[n][2]
            pos3_[n] = A[n][3]
            vel1_[n] = A[n][4]
            vel2_[n] = A[n][5]
            vel3_[n] = A[n][6]
            tinj_[n]  = A[n][7]
            col_[n]   = A[n][8]
            
            if spectrum_:
                nmic_[n]  = A[n][9]
                cmpratio_ = A[n][10]
                shkflag_  = A[n][11]
                shkgradp_ = A[n][12]
                a_        = A[n][13]
                b_        = A[n][14]
                shkvL_[n,:] = A[n][15:15+nflx_]
                shkvR_[n,:] = A[n][15+nflx_:15+2*nflx_]
                
                specE_[n,:] = A[n][15+2*nflx_:15+2*nflx_+ebins_]
                specN_[n,:] = A[n][15+2*nflx_ + ebins_:15+2*nflx_ + 2*ebins_]

        tup_ = [(vars_[0], ids_),(vars_[1], pos1_),(vars_[2], pos2_),(vars_[3], pos3_),\
                (vars_[4], vel1_),(vars_[5], vel2_),(vars_[6], vel3_),(vars_[7], tinj_),\
                (vars_[8], col_)]

        if spectrum_:
            tup_.append((vars_[9],nmic_))
            tup_.append((vars_[10], cmpratio_))
            tup_.append((vars_[11], shkflag_))
            tup_.append((vars_[12], shkgradp_))
            tup_.append((vars_[13], a_))
            tup_.append((vars_[14], c_))
            tup_.append((vars_[15], shkvL_))
            tup_.append((vars_[16], shkvR_))
            tup_.append((vars_[17], specE_))
            tup_.append((vars_[18], specN_))
        
        DataDict_.update(dict(tup_))
        return DataDict_
        




