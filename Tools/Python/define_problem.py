import os
import sys
import menu as menu
import pluto_files_IO as pfIO

class DefineProblem(object):
  def __init__(self, work_dir, pluto_dir, auto_update):
    """Defines the problem for the PLUTO code.

    This class essentially creates the definitions.h file
    based on the inputs given by the user from the PlutoSetup menu.
    In case the header file is already present, this class will
    read those default values to re-create the file.

    **Inputs**:
      1. work_dir = Path to PLUTO code working directory
      2. pluto_dir = Path to PLUTO source directory
      3. auto-update = Boolean that indicates auto-update of definitions.h file.

    **Output**:
      It generates a definitions.h file.

    """
    #Some class attributes that will be used in various class methods.
    self.work_dir         = work_dir
    self.ini_fname        = self.work_dir + '/init.c'
    self.pluto_dir        = pluto_dir
    self.krome_dir        = pluto_dir + '/Src/Cooling/KROME/krome/'
    self.auto_update      = auto_update
    self.def_fname        = self.work_dir + '/definitions.h'
    self.additional_files = []
    self.header_files     = []
    self.additional_flags = []
    self.pluto_path       = []
    self.def_file_list    = []

    self.ChkCompatiblity() #To Replace Old Keywords with New Ones

    # defining the PLUTO entries and its default values in lists.  
    self.entries = ['PHYSICS', 'DIMENSIONS', 'COMPONENTS', 'GEOMETRY',
                    'BODY_FORCE', 'FORCED_TURB','COOLING', 'RECONSTRUCTION', 
		    'TIME_STEPPING','DIMENSIONAL_SPLITTING', 'NTRACER', 
		    'USER_DEF_PARAMETERS']
    self.default = ['HD', '1', '1', 'CARTESIAN','NO', 'NO',
                    'NO','LINEAR','RK2',
                    'NO', '0', '0']

    # Creating a dictionary of flags that are invoked by giving arguments.
    flag_keys = ['WITH-CHOMBO', 'FULL', 'WITH-FD', 'WITH-SB', 'WITH-FARGO',
                 'WITH-PARTICLES','WITH-CR_TRANSPORT']
    #self.flag_dict = {key: False for key in flag_keys} DOESNT WORK WITH PYTHON 2.6
    self.flag_dict = {'WITH-CHOMBO':False, 'FULL':False, 'WITH-FD':False,
                      'WITH-SB':False, 'WITH-FARGO':False,'WITH-PARTICLES':False,
                      'WITH-CR_TRANSPORT':False}
    
    for arg in sys.argv:
      if arg[2:].upper() in flag_keys:
        self.flag_dict[arg[2:].upper()] = True
      elif arg[2:] == 'with-chombo:':
        self.flag_dict['WITH-CHOMBO'] = True
      else:
        pass 


    # Generates Full Option List.
    self.GenerateOptionsList()
    
    # !!!!!!!!!!! CR_TRANSPORT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # if (self.flag_dict['WITH-CR_TRANSPORT'] == True):
    #   self.options[self.entries.index('PHYSICS')] = ['CR_TRANSPORT']
    # else:
    #   i = self.entries.index('PHYSICS')
    #   self.options[i].remove('CR_TRANSPORT')
      

    # print self.options[self.entries.index('PHYSICS')]
    # print self.entries
    # print self.options
    #   
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #Updates Options, default based on FLAGS.
    if True in self.flag_dict.values():
      self.AfterFlagLists() 

    #Read the exsisting definition.h file or Browse the menu for Setting up problem.
    self.ReadOrBrowse(Ents=self.entries, Defs=self.default, Opts=self.options, MenuTitle="Setup problem")
    
    
    #Process the PHYSICS Modules.
    if self.default[self.entries.index('PHYSICS')] == 'HD':
      self.ProcessHDModule()
      self.ReadOrBrowse(Ents = self.entries_HD, Defs = self.default_HD, Opts = self.options_HD, MenuTitle = "HD Menu")
      self.eos = self.default_HD[self.entries_HD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'RHD':
      self.ProcessRHDModule()
      self.ReadOrBrowse(Ents = self.entries_RHD, Defs = self.default_RHD, Opts = self.options_RHD, MenuTitle = "RHD Menu")
      self.eos = self.default_RHD[self.entries_RHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'MHD':
      self.ProcessMHDModule()
      self.ReadOrBrowse(Ents = self.entries_MHD, Defs = self.default_MHD, Opts = self.options_MHD, MenuTitle = "MHD Menu")
      self.eos = self.default_MHD[self.entries_MHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'RMHD':
      self.ProcessRMHDModule()
      self.ReadOrBrowse(Ents = self.entries_RMHD, Defs = self.default_RMHD, Opts = self.options_RMHD, MenuTitle = "RMHD Menu")
      self.eos = self.default_RMHD[self.entries_RMHD.index('EOS')]

    if self.default[self.entries.index('PHYSICS')] == 'CR_TRANSPORT':
      self.ProcessCR_TransportModule()
      self.ReadOrBrowse(Ents = self.entries_CR_TRANSPORT,
                        Defs = self.default_CR_TRANSPORT,
                        Opts = self.options_CR_TRANSPORT, MenuTitle = "CR_Transport Menu")
      self.eos = self.default_CR_TRANSPORT[self.entries_CR_TRANSPORT.index('EOS')]
    
    # Process the KROME cooling module when required. 
    if self.default[self.entries.index('COOLING')] == 'KROME':
      self.ProcessKROMEModule()
      self.ReadOrBrowse(Ents = self.entries_KROME, Defs = self.default_KROME, Opts = self.options_KROME, MenuTitle = "KROME Menu")
      self.WriteKromeOpts()

    #UserDef Para and Constants
    self.udef_params = []
    self.udef_const = []
    self.udef_const_vals = []
    self.nparam = int(self.default[self.entries.index('USER_DEF_PARAMETERS')])

    #self.nconst = int(self.default[self.entries.index('USER_DEF_CONSTANTS')])
    if self.nparam > 0:
      self.ProcessUserDefPara()
    
    self.ProcessUserDefConst()

    # Write a List def_file_list which will be written as the header file.    
    self.WriteDefFileList()
    pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
    pf.List2File(self.def_file_list)
    
    self.AppendInitFile()


  def ChkCompatiblity(self):
    oldKeys_     = ['INTERPOLATION', 'MHD_FORMULATION','RESISTIVE_MHD']
    replaceKeys_ = ['RECONSTRUCTION','DIVB_CONTROL',   'RESISTIVITY']
    if (os.path.exists(self.work_dir+'/definitions.h')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')

      for i in range(len(oldKeys_)):
        pf.ReplaceWord(oldKeys_[i], replaceKeys_[i])
                

  def GenerateOptionsList(self):
    """Creates a default option list.

    This method of class DefineProblem will create a
    default valued options list for each entry in the
    entries list. These are essentially the options
    that will be browsed in the Pluto Setup Menu.
    """
    phylist = ['ADVECTION','HD','RHD','MHD','RMHD','CR_TRANSPORT']
    dimlist = ['1','2','3']
    comlist = ['1','2','3']
    geolist = ['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL']
    bfolist = ['NO','VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)']
    ftblist = ['NO','YES']
    coolist = ['NO','POWER_LAW','TABULATED','SNEq','MINEq','H2_COOL', 'KROME']
    intlist = ['FLAT','LINEAR','LimO3','WENO3','PARABOLIC']
    tmslist = ['EULER','RK2','RK3','HANCOCK','CHARACTERISTIC_TRACING']
    dislist = ['YES','NO']
    ntrlist = ['%d'%n for n in range(9)]
    udplist = ['%d'%n for n in range(32)]
    udclist = ['%d'%n for n in range(32)]

    self.options = [phylist, dimlist, comlist, geolist, bfolist,ftblist,
                    coolist, intlist, tmslist,
                    dislist, ntrlist, udplist,
                    udclist]

  def AfterFlagLists(self):
    """Modify options and default list based on command-line flags.

    This method is called after generation of default option list.
    It modifies the members of the options list and if required
    also the default list based on the conditions required by
    the flags set using system arguments.
    """
    if self.flag_dict['FULL']:
      self.options[self.entries.index('RECONSTRUCTION')] = ['FLAT','LINEAR','LimO3', 'WENO3','PARABOLIC', 'MP5']
    
    if self.flag_dict['WITH-CHOMBO']:
      self.options[self.entries.index('GEOMETRY')] = ['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL']
      self.options[self.entries.index('RECONSTRUCTION')] = ['FLAT','LINEAR','WENO3','PARABOLIC']
      self.options[self.entries.index('TIME_STEPPING')] = ['EULER','HANCOCK','CHARACTERISTIC_TRACING','RK2']
      self.default[self.entries.index('TIME_STEPPING')] = 'HANCOCK'
      self.options[self.entries.index('DIMENSIONAL_SPLITTING')] = ['NO']
      self.default[self.entries.index('DIMENSIONAL_SPLITTING')] = 'NO'
    
    if self.flag_dict['WITH-FARGO']:
      self.options[self.entries.index('PHYSICS')] = ['HD', 'MHD']
      self.options[self.entries.index('DIMENSIONS')] = ['2','3']
      self.default[self.entries.index('DIMENSIONS')] = '2'
      self.options[self.entries.index('DIMENSIONAL_SPLITTING')] = ['NO']
      self.default[self.entries.index('DIMENSIONAL_SPLITTING')] = 'NO'
        
    if self.flag_dict['WITH-FD']:
      self.options[self.entries.index('PHYSICS')] = ['HD', 'MHD']
      self.options[self.entries.index('GEOMETRY')] = ['CARTESIAN']
      self.options[self.entries.index('RECONSTRUCTION')] = ['WENO3_FD', 'WENOZ_FD', 'MP5_FD','LIMO3_FD']
      self.default[self.entries.index('RECONSTRUCTION')] = 'WENOZ_FD'
      self.options[self.entries.index('TIME_STEPPING')] = ['RK3','SSP_RK4']
      self.default[self.entries.index('TIME_STEPPING')] = 'RK3'

    if self.flag_dict['WITH-SB']:
      self.options[self.entries.index('PHYSICS')] = ['HD', 'MHD']
      self.default[self.entries.index('PHYSICS')] = 'MHD'
      self.options[self.entries.index('DIMENSIONS')] = ['2', '3']
      self.default[self.entries.index('DIMENSIONS')] = '2'
      self.options[self.entries.index('COMPONENTS')] = ['2', '3']
      self.default[self.entries.index('COMPONENTS')] = '2'
      self.options[self.entries.index('BODY_FORCE')] = ['VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)']
      self.default[self.entries.index('BODY_FORCE')] = 'VECTOR'


  def ReadOrBrowse(self, Ents=None, Opts=None, Defs=None, MenuTitle=None):
    """Reads or Browses the entries, options and defaults to create header file.

    This method either reads the already exsisting definitions.h file or browses
    throught the three lists which are provided as inputs.

    **Inputs**:
      1. Ents = List of entries. [None]
      2. Opts = List of options corresponding to each member in Ents [None]
      3. Defs = List of default value from Opts corresponding to each member in Ents [None]
      4. MenuTitle = Title for the Menu to be Browsed [None]
    """
    if (os.path.exists(self.work_dir+'/definitions.h')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
      pf.UpdateListFromFile(Ents, Defs)
      for i in range(len(Ents)):
        if Defs[i] not in Opts[i]:
          Defs[i] = Opts[i][0]
        else:
          pass

    # Provides Browsing options using the menu file in case of no automatic update flag.
    if self.auto_update == 0:
      selection = ''
      menu.SetTitle(MenuTitle)
      selection = menu.Browse(Ents, default=Defs, options=Opts)
      
  def ProcessHDModule(self):
    """
    Provides entries, options and defaults specific to Hydro Module.
    Also updates them accordingly if required by flags.
    """
    self.entries_HD = ['EOS', 'ENTROPY_SWITCH',
                       'THERMAL_CONDUCTION', 'VISCOSITY',
                       'ROTATING_FRAME']
    self.default_HD = ['IDEAL', 'NO',
                       'NO',    'NO',
                       'NO']
    self.options_HD = [['IDEAL','PVTE_LAW','ISOTHERMAL'],
#                           ['NO','YES'],
                        ['NO','SELECTIVE','ALWAYS','CHOMBO_REGRID'],
                       ['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_LEGENDRE'],
                       ['NO','EXPLICIT','SUPER_TIME_STEPPING','RK_LEGENDRE'],
                       ['NO','YES']]
  
    if self.flag_dict['WITH-CHOMBO']: # Chombo does not support STS at the
                                      # moment. Only explicit allowed with Chombo
      i =  self.entries_HD.index('THERMAL_CONDUCTION')
      self.options_HD[i] = ['NO','EXPLICIT']
      i =  self.entries_HD.index('VISCOSITY')
      self.options_HD[i] = ['NO','EXPLICIT']

  def ProcessRHDModule(self):
    """
    Provides entries, options and defaults specific to Relativistic
    Hydro Module. Also updates them accordingly if required by flags.
    """
    self.entries_RHD = ['EOS', 'ENTROPY_SWITCH']
    self.default_RHD = ['IDEAL', 'NO', 'NO']
    self.options_RHD = [['IDEAL','TAUB'],
#                            ['NO','YES']]
                          ['NO','SELECTIVE','ALWAYS','CHOMBO_REGRID']]

  def ProcessMHDModule(self):
    """
    Provides entries, options and defaults specific to Magneto-
    Hydro Module.Also updates them accordingly if required by flags.
    """
    self.entries_MHD = ['EOS', 'ENTROPY_SWITCH', 'DIVB_CONTROL',
                        'BACKGROUND_FIELD',
                        'AMBIPOLAR_DIFFUSION', 'RESISTIVITY', 'HALL_MHD',
                        'THERMAL_CONDUCTION', 'VISCOSITY',
                        'ROTATING_FRAME']
    self.default_MHD = ['IDEAL','NO','EIGHT_WAVES',
                        'NO',
                        'NO','NO','NO',
                        'NO','NO',
                        'NO']
    self.options_MHD = [['IDEAL','PVTE_LAW','ISOTHERMAL'],
#                            ['NO','YES'],
                        ['NO','SELECTIVE','ALWAYS','CHOMBO_REGRID'],
                        ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'],
                        ['NO','YES'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','EXPLICIT', 'SUPER_TIME_STEPPING','RK_LEGENDRE'],
                        ['NO','YES']]
    
    if self.flag_dict['WITH-CHOMBO']:
      indx_ = self.entries_MHD.index('DIVB_CONTROL')
      self.options_MHD[indx_] = ['NO','EIGHT_WAVES','DIV_CLEANING']
      indx_ = self.entries_MHD.index('RESISTIVITY')
      self.options_MHD[indx_] = ['NO','EXPLICIT']
      indx_ = self.entries_MHD.index('THERMAL_CONDUCTION')
      self.options_MHD[indx_] = ['NO','EXPLICIT']
      indx_ = self.entries_MHD.index('VISCOSITY')
      self.options_MHD[indx_] = ['NO','EXPLICIT']

    if self.flag_dict['WITH-FD']:
      indx_ = self.entries_MHD.index('DIVB_CONTROL')
      self.options_MHD[indx_] = ['NO','EIGHT_WAVES','DIV_CLEANING']

    if self.flag_dict['WITH-SB'] or self.flag_dict['WITH-FARGO']:
      indx_ = self.entries_MHD.index('DIVB_CONTROL')
      self.options_MHD[indx_] = ['CONSTRAINED_TRANSPORT']
      self.default_MHD[indx_] = 'CONSTRAINED_TRANSPORT'
  
  def ProcessRMHDModule(self):
    """
    Provides entries, options and defaults specific to Relativisitc
    Magneto-Hydro Module.Also updates them accordingly if required by flags.
    """
    self.entries_RMHD = ['EOS', 'ENTROPY_SWITCH','DIVB_CONTROL','RESISTIVITY']
    self.default_RMHD = ['IDEAL', 'NO', 'NO','NO']
    self.options_RMHD = [['IDEAL', 'TAUB'],
#                            ['NO','YES'],
                         ['NO','SELECTIVE','ALWAYS','CHOMBO_REGRID'],
                         ['NO','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'],
                         ['NO','EXPLICIT','IMEX','OP_SPLIT','NEW_SPLIT']]

    if self.flag_dict['WITH-CHOMBO']:
      indx_ = self.entries_RMHD.index('DIVB_CONTROL')
      self.options_RMHD[indx_] = ['NO','EIGHT_WAVES','DIV_CLEANING','NO']


  def ProcessCR_TransportModule(self):
    """
    Provides entries, options and defaults specific to Hydro Module.
    Also updates them accordingly if required by flags.
    """
    self.entries_CR_TRANSPORT = ['EOS','ANISOTROPY']
    self.default_CR_TRANSPORT = ['ISOTHERMAL','NO'] 
    self.options_CR_TRANSPORT = [['ISOTHERMAL'],['NO','YES']]

  def ProcessKROMEModule(self):
    """
    Provides entries, options, and defaults for the KROME cooling module.
    """
    netwrkfiles = os.listdir(self.krome_dir+'networks/')
    self.entries_KROME = ['NETWORK_FILE', 'USE_N','COOLING_TYPE', 'HEATING_TYPE', 'GAMMA_TYPE']
    self.default_KROME = [netwrkfiles[0], 'NO', 'NONE', 'NONE', 'DEFAULT']
    self.options_KROME = [netwrkfiles, ['NO', 'YES'], [ 'NONE', 'H2', 'ATOMIC', 'Z', 'H2,ATOMIC'], ['NONE', 'COMPRESS', 'PHOTO', 'XRAY', 'CHEM'], ['DEFAULT','FULL', 'ROT', 'VIB', 'EXACT', 'REDUCED']]

    if self.eos == 'IDEAL':
      self.options_KROME[4] = ['DEFAULT']



  def ProcessUserDefPara(self):
    """
    Sets the Userdefined parameters
    """
    self.udef_params = ['USER_PAR_%.2d'%i for i in range(self.nparam)]
    if (os.path.exists(self.work_dir+'/definitions.h')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
      scrh = pf.LocateString('parameters')
      k0   = scrh[0][0] + 2
      par_lines = pf.ReadLines(k0, k0 + self.nparam)
      for n in range(self.nparam):
        try:
          x = par_lines[n].split()
          x[0] == '#define'
        except IndexError:
          pass
        else:
          if (x[0] == "#define"):
            self.udef_params[n] = x[1]
          else:          
            break;
                    
    if self.auto_update == 0:
      menu.SetTitle ("User-defined Parameters")
      par_entries = ['%d'%i for i in range(self.nparam)]
      menu.Insert(par_entries,self.udef_params)

          
  def ProcessUserDefConst(self):
    """
    Sets the Userdefined Constants.
    """
    if (os.path.exists(self.work_dir+'/definitions.h')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
      old_beg_scrh = pf.LocateString('symbolic')
      try:
        old_beg_scrh[0][0]
      except IndexError:
        pass
      else:
        del_indx = pf.LocateString('USER_DEF_CONSTANTS')
        pf.DeleteLines(del_indx[0][0], del_indx[0][0])
        old_beg_scrh = pf.LocateString('symbolic')
        pf.ReplaceLine('/* [Beg] user-defined constants (do not change this line) */', old_beg_scrh[0][0])
        old_end_scrh = pf.LocateString('supplementary')
        pf.InsertLine('/* [End] user-defined constants (do not change this line) */', old_end_scrh[0][0] - 1)
          
      scrh_beg = pf.LocateString('[Beg]')
      k_beg   = scrh_beg[0][0]+1
      scrh_end = pf.LocateString('[End]')
      k_end   = scrh_end[0][0]-1
      const_lines = pf.ReadLines(k_beg, k_end)
      #print const_lines
      for n in range(len(const_lines)):
        x = const_lines[n].split()
        try:                
          x[0] == '#define'
        except IndexError: 
          pass
        else:
          if (x[0] == '#define'):
            self.udef_const.append(x[1])
            self.udef_const_vals.append(x[2])
          else:
            continue
          

  def NonUserFriendlyConst(self):
    """
        Sets the non-user friendly constants. : OBSOLETE
        
        This function only stores the Non-User Friendly constants
        whose values are different than their DEFAULT values and
        incorporates them in User Defined Constants.
    """
    self.non_usfr = []
    tmplist1 = ['INITIAL_SMOOTHING', 'WARNING_MESSAGES',
                'INTERNAL_BOUNDARY', 'SHOCK_FLATTENING']
    tmplist2     = ['NO','YES','NO','NO']
    tmplist2_def = ['NO','YES','NO','NO']
    
    
#       if self.flag_dict['WITH-CHOMBO']:
#           tmplist1 += ['CHOMBO_REF_VAR','CHOMBO_LOGR']
#           tmplist2 += ['ENG',       'NO']  

    if not self.flag_dict['WITH-FD']:
      tmplist1 = tmplist1 + ['CHAR_LIMITING', 'LIMITER']
      tmplist2 = tmplist2 + ['NO', 'DEFAULT']
      tmplist2_def = tmplist2_def + ['NO', 'DEFAULT']
        
    if 'DIVB_CONTROL' in self.mod_entries:
      divb_mode = self.mod_default[self.mod_entries.index('DIVB_CONTROL')]
      if divb_mode == 'CONSTRAINED_TRANSPORT':
        tmplist1 = tmplist1 + ['CT_EMF_AVERAGE', 'CT_EN_CORRECTION', 'ASSIGN_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['UCT_HLL', 'NO', 'NO']
        tmplist2_def = tmplist2_def + ['UCT_HLL', 'NO', 'NO']
      else:
        tmplist1 = tmplist1 + ['ASSIGN_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['NO']
        tmplist2_def = tmplist2_def + ['NO']
      
      if not self.flag_dict['WITH-CHOMBO']:
        tmplist1 = tmplist1 + ['UPDATE_VECTOR_POTENTIAL']
        tmplist2 = tmplist2 + ['NO']
        tmplist2_def = tmplist2_def + ['NO']

#        if 'HANCOCK' in self.default:
#            if (self.phymodule == 'RMHD'):
#                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
#                tmplist2 = tmplist2 + ['NO']
#           else:
#                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
#                tmplist2 = tmplist2 + ['YES']
    
    longword = max(len(w) for w in tmplist1)
    
    if (os.path.exists(self.work_dir+'/definitions.h')):
      pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
      scrh = pf.LocateString('supplementary')
      try:
        scrh[0]
      except IndexError:
        pass
      else:
        pf.UpdateListFromFile(tmplist1, tmplist2)
    
    for i in range(len(tmplist1)):
      if tmplist2[i] != tmplist2_def[i]:
        self.non_usfr = self.non_usfr + ['#define  '+tmplist1[i].ljust(28)+'   '+tmplist2[i]+'\n']

  def AppendAdditionalFiles(self):
    """
    Adds additional object files based on
    modular defintions and requirements. 
    """
    interp_mode = self.default[self.entries.index('RECONSTRUCTION')]

    if interp_mode == 'LINEAR':
      self.additional_files.append('plm_states.o')
    elif interp_mode == 'PARABOLIC':
      self.additional_files.append('ppm_states.o')
      self.additional_files.append('ppm_coeffs.o')
      self.header_files.append('ppm_coeffs.h')
    elif interp_mode in ['FLAT', 'LimO3', 'WENO3']:
      self.additional_files.append(interp_mode.lower()+'_states.o')
    else:
      pass
    
    if self.flag_dict['WITH-FD']:
      self.additional_files += ['fd_states.o', 'fd_reconstruct.o', 'fd_flux.o']

    if self.default[self.entries.index('COOLING')] not in ['NO', 'POWER_LAW', 'KROME']:
      self.additional_files += ['cooling_source.o','cooling_ode_solver.o']

    if self.phymodule == 'MHD' or self.phymodule == 'RMHD':
      self.additional_files.append('vec_pot_diff.o')
      if not self.flag_dict['WITH-CHOMBO']:
        self.additional_files.append('vec_pot_update.o')

    if self.flag_dict['WITH-CHOMBO']:
      if self.default[self.entries.index('TIME_STEPPING')] in ['EULER', 'RK2']:
        self.additional_files.append('PatchEuler.o')
        self.additional_files.append('update_stage.o')
      else:
        self.additional_files.append('PatchCTU.o')
    else:
      cmset = set(['CHARACTERISTIC_TRACING', 'HANCOCK']) & set(self.default)
      if len(cmset) != 0 and self.default[self.entries.index('DIMENSIONAL_SPLITTING')] == 'NO':
        self.additional_files.append('ctu_step.o')
      elif self.default[self.entries.index('TIME_STEPPING')] == 'SSP_RK4':
        self.additional_files.append('unsplit.ssprk.o')                
      else:
        self.additional_files.append('rk_step.o')
        self.additional_files.append('update_stage.o')
    
    if 'HANCOCK' in self.default:
      self.additional_files.append('hancock.o')

    if 'CHARACTERISTIC_TRACING' in self.default:
      self.additional_files.append('char_tracing.o')

    # -- Diffusion time stepping -- 
    parabolic_update = 0
    if 'SUPER_TIME_STEPPING' in self.mod_default:
      parabolic_update = 1
      self.additional_files += ['sts.o']
      
    if 'RK_LEGENDRE' in self.mod_default:
      parabolic_update = 1
      self.additional_files += ['rkl.o']

    if 'EXPLICIT' in self.mod_default:
      parabolic_update = 1

    if (parabolic_update):
      self.additional_files.append('parabolic_update.o')

      
  def AppendPlutoPathAndFlags(self):
    """
    Adds additional C flags and path to 'makefile' based on
    modular defintions and requirements. 
    """
    self.pluto_path.append(self.phymodule+'/')
    
    dis_eff = ['Dust','Thermal_Conduction', 'Viscosity']
    for de in dis_eff:
        if (     (de.upper() in self.mod_entries)
             and (self.mod_default[self.mod_entries.index(de.upper())] != 'NO')):
            self.pluto_path.append(de+'/')
    
    if self.phymodule == 'MHD' or self.phymodule == 'RMHD':
      divb_mode   = self.mod_default[self.mod_entries.index('DIVB_CONTROL')]
      resistivity = self.mod_default[self.mod_entries.index('RESISTIVITY')]

      if divb_mode == 'CONSTRAINED_TRANSPORT':
        self.pluto_path.append('MHD/CT/')
      elif divb_mode == 'DIV_CLEANING':
        self.pluto_path.append('MHD/GLM/')
      else:
          pass

      if self.phymodule == 'MHD' and \
        self.mod_default[self.mod_entries.index('AMBIPOLAR_DIFFUSION')] != 'NO':
        self.pluto_path.append('MHD/Ambipolar_Diffusion/')

      if self.phymodule == 'MHD' and \
        self.mod_default[self.mod_entries.index('HALL_MHD')] != 'NO':
        self.pluto_path.append('MHD/Hall_MHD/')

      if self.phymodule == 'MHD' and \
        self.mod_default[self.mod_entries.index('RESISTIVITY')] != 'NO':
        self.pluto_path.append('MHD/Resistivity/')


      if self.phymodule == 'RMHD' and \
        resistivity   != 'NO':
        self.pluto_path.append('RMHD/Resistivity/')
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Hot fix to reverse path list for debugging purposed.
        # This will have to be removed at some point later on...
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.pluto_path=list(reversed(self.pluto_path)) 

        if ('rk_step.o' in self.additional_files):
          if (resistivity != "EXPLICIT"):
            self.additional_files.remove('rk_step.o')
            self.additional_files.append('rk_step_imex.o')

           
    if self.flag_dict['WITH-FARGO']:
      self.pluto_path.append('Fargo/')
      self.additional_flags.append(' -DFARGO')

    if self.flag_dict['WITH-FD']:
      self.additional_flags.append(' -DFINITE_DIFFERENCE')

    if self.flag_dict['WITH-PARTICLES']:
      self.pluto_path.append('Particles/')
      self.additional_flags.append(' -DPARTICLES')

    if self.flag_dict['WITH-SB']:
      self.pluto_path.append('MHD/ShearingBox/')
      self.additional_flags.append(' -DSHEARINGBOX')
    
    bd_force = self.default[self.entries.index('BODY_FORCE')]
    for_turb = self.default[self.entries.index('FORCED_TURB')]
    if for_turb == 'YES':
      self.pluto_path.append('Forced_Turb/') 
    
    cool_mode = self.default[self.entries.index('COOLING')]
    if cool_mode != 'NO':
      if cool_mode == 'TABULATED':
        self.pluto_path.append('Cooling/TABULATED/')
      elif cool_mode == 'POWER_LAW':
        self.pluto_path.append('Cooling/Power_Law/')
      else:
        self.pluto_path.append('Cooling/'+ cool_mode +'/')

    if 'EOS' in self.mod_entries:
      if 'PVTE_LAW' in self.mod_default:
        tmp1 = 'PVTE'
      else:
        tmp1 = self.eos[0]+self.eos[1:].lower()
      self.pluto_path.append('EOS/'+tmp1+'/')

  def UpdatePlutoIni(self):
    """
    Updates pluto.ini file with values of UserDefined Parameters
    """
    pf   = pfIO.PlutoFiles(self.work_dir+'/pluto.ini')
    scrh = pf.LocateString('[Parameters]')
    try:
      scrh[0]
    except IndexError:
      print("Parameters keyword not found in pluto.ini")
      sys.exit()
    else:
      pass
    
    ipos = scrh[0][0] + 2
    tmplist1 = pf.ReadLines(ipos,ipos+512)
    paradict = {}
    cmms = []
    for x in tmplist1:
      if (len(x.split()) == 0): continue  # skip blank lines 
      paradict.update({x.split()[0]:x.split()[1]})
      try:
        cmmval = x.split()[2]
      except IndexError:
        cmms.append('')
        continue
      else:
        if cmmval == '#' or cmmval.startswith('#'):
          cmms.append(' '.join(x.split()[2:]))
        else:
          cmms.append('')

    for x in self.udef_params:
      if x in paradict.keys():
        pf.InsertLine(x.ljust(28) + paradict[x] +'  '+cmms[self.udef_params.index(x)]+'\n', ipos)
      else:
        try:
          cmms[self.udef_params.index(x)]
        except IndexError:
          pf.InsertLine(x.ljust(28) + '0.0' + '\n', ipos)
        else:
          pf.InsertLine(x.ljust(28) + '0.0'+'  '+cmms[self.udef_params.index(x)]+ '\n', ipos)
      ipos = ipos + 1
    pf.DeleteLines(ipos,ipos+100)
  

  def WriteKromeOpts(self):
    self.kromeoptstr = 'python '+self.krome_dir+'krome -interfaceC -n='+self.krome_dir+'networks/'+ self.default_KROME[self.entries_KROME.index('NETWORK_FILE')]
    if self.default_KROME[self.entries_KROME.index('USE_N')] == 'YES':
      self.kromeoptstr += ' -useN'
    else:
      self.kromeoptstr += ' -useX'

    if self.default_KROME[self.entries_KROME.index('COOLING_TYPE')] != 'NONE':
      self.kromeoptstr += ' -cooling='+ self.default_KROME[self.entries_KROME.index('COOLING_TYPE')]

    if self.default_KROME[self.entries_KROME.index('HEATING_TYPE')] != 'NONE':
      self.kromeoptstr += ' -heating='+ self.default_KROME[self.entries_KROME.index('HEATING_TYPE')]

    if self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')] != 'DEFAULT':
      self.kromeoptstr += ' -gamma='+ self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')]

    #remove all checks : 
    self.kromeoptstr += ' -noCheck -noSinkCheck'

  def WriteDefFileList(self):
    """
    Writes all modular entries, options, defaults into a list.
    """
    for x in self.entries:
        self.def_file_list.append('#define  '+x.ljust(28)+'   '+self.default[self.entries.index(x)]+'\n')

    self.def_file_list.append('\n/* -- physics dependent declarations -- */\n\n')
    self.phymodule = self.default[self.entries.index('PHYSICS')]
    
    tmp1 = ['entries_%s'%self.phymodule, 'default_%s'%self.phymodule]

    # print self.phymodule
    # sys.exit(0)       
 
    self.mod_entries = self.__getattribute__(tmp1[0])
    self.mod_default = self.__getattribute__(tmp1[1])

    for x in self.mod_entries:
      self.def_file_list.append('#define  '+x.ljust(28)+'   '+ self.mod_default[self.mod_entries.index(x)]+'\n')

    self.AppendAdditionalFiles()
    self.AppendPlutoPathAndFlags()

    try:
     self.entries_KROME
    except AttributeError:
     pass
    else:
     self.def_file_list.append('\n/* -- KROME Cooling Options -- */\n\n')
     for x in self.entries_KROME:
        self.def_file_list.append('#define  '+x.ljust(28)+'   '+self.default_KROME[self.entries_KROME.index(x)]+'\n')

    # always insert user-defined paramters and constants sections

    self.def_file_list.append('\n/* -- user-defined parameters (labels) -- */\n\n')
    for x in self.udef_params:
      self.def_file_list.append('#define  '+x.ljust(28)+'   '+'%d'%self.udef_params.index(x)+'\n')

    self.UpdatePlutoIni()

            
    self.def_file_list.append('\n/* [Beg] user-defined constants (do not change this line) */\n\n')
    self.NonUserFriendlyConst()
    if len(self.non_usfr) > 0:
      for x in self.non_usfr:
        self.def_file_list.append(x)
    
    for i in range(len(self.udef_const)):
      self.def_file_list.append('#define  '+self.udef_const[i].ljust(28)+'   '+self.udef_const_vals[i]+'\n')

    self.def_file_list.append('\n/* [End] user-defined constants (do not change this line) */\n')
    
    #self.def_file_list.append('\n/* -- supplementary constants (user editable) -- */ \n\n')
    


  def AppendInitFile(self):
    #Read init_domain.c from the Src/Templates folder.
    inidomainfl_name = self.pluto_dir + '/Src/Templates/init_domain.c'
    pf1 = pfIO.PlutoFiles(inidomainfl_name)

    #Convert the contents of file init_domain.c into a list
    inidomain = pf1.File2List()

    #Check if the init.c file exists in the work directory
    inifl_exits = os.path.exists(self.ini_fname)
    
    #If it does read the file and convert it contents into a list
    # Then break the list into part above the Init Function and part below the Init Function.
    if inifl_exits:
      pf = pfIO.PlutoFiles(self.ini_fname)
      nlist = pf.File2List()
      scrh1 = pf.LocateString('void Init (double *v, double x1, double x2, double x3)\n')
      scrh3 = pf.LocateString('void InitDomain (Data *d, Grid *grid)\n')
      
      if len(scrh1) == 0:
        for item in nlist:
          line = item.split()
          if any(x in line for x in ["void", "Init", "Init(double"]):
            lindx1 = nlist.index(item)
            break
      else:
        lindx1 = scrh1[0][0]
          
      # If the InitDomain function is not present in the init.c file then
      # read line by line from the Init Function till you get number
      # of open curly brackets equal to close curly brackets.

      if len(scrh3) == 0: 
        l1 = nlist[0:lindx1]
        l2 = nlist[lindx1:]
        opBracket = 0
        clBracket = 0
	addIndx   = 0
        for item in l2:
          line = item.split()
	  addIndx += 1
          if ((len(line) > 0) and (line[0] !='*')):
            for w in line:
              if "{" in list(w): opBracket += 1
              if "}" in list(w): clBracket += 1
          if (opBracket > 0):
            if(opBracket == clBracket):
               addline = addIndx
               break
        # Finally create a new list by adding the content of init_domain.c file
        NewInit = l1 + l2[:addline] + inidomain + l2[addline:] 
        # Create the new init.c file. 
        pf.List2File(NewInit) 
    
    
  
