import os
import sys
import menu
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
        self.work_dir = work_dir
        self.pluto_dir = pluto_dir
        self.krome_dir = pluto_dir + 'Lib/krome/'
        self.auto_update = auto_update
        self.def_fname = self.work_dir + '/definitions.h'
        self.additional_files = []
        self.header_files = []
        self.additional_flags = []
        self.pluto_path = []
        self.def_file_list = []

        # defining the PLUTO entries and its default values in lists.  
        self.entries = ['PHYSICS', 'DIMENSIONS', 'COMPONENTS', 'GEOMETRY',
                        'BODY_FORCE', 'COOLING', 'INTERPOLATION', 'TIME_STEPPING',
                        'DIMENSIONAL_SPLITTING', 'NTRACER', 'USER_DEF_PARAMETERS'
                        ]
        self.default = ['HD', '1', '1', 'CARTESIAN','NO',
                        'NO','LINEAR','RK2',
                        'YES', '0', '0',
                        ]

        # Creating a dictionary of flags that are invoked by giving arguments.
        flag_keys = ['WITH-CHOMBO', 'FULL', 'WITH-FD', 'WITH-SB', 'WITH-FARGO']
        #self.flag_dict = {key: False for key in flag_keys} DOESNT WORK WITH PYTHON 2.6
	self.flag_dict = {'WITH-CHOMBO':False, 'FULL':False, 'WITH-FD':False, 'WITH-SB':False, 'WITH-FARGO':False}
        
        for arg in sys.argv:
            if arg[2:].upper() in flag_keys:
                self.flag_dict[arg[2:].upper()] = True
            elif arg[2:] == 'with-chombo:':
                self.flag_dict['WITH-CHOMBO'] = True
            else:
                pass 

        # Generates Full Option List.
        self.GenerateOptionsList()
        
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

        if self.default[self.entries.index('COOLING')] == 'KROME':
            self.ProcessKROMEModule()
            self.ReadOrBrowse(Ents = self.entries_KROME, Defs = self.default_KROME, Opts = self.options_KROME, MenuTitle = "KROME Menu")
            self.WriteKromeOpts()
            #self.eos = self.default_KROME[self.entries_RMHD.index('EOS')]


        #UserDef Para and Constants
        self.udef_params = []
        self.udef_const = []
        self.udef_const_vals = []
        self.nparam = int(self.default[self.entries.index('USER_DEF_PARAMETERS')])
        #self.nconst = int(self.default[self.entries.index('USER_DEF_CONSTANTS')])
        if self.nparam > 0:
            self.ProcessUserDefPara()
        
        #if self.nconst > 0:
        self.ProcessUserDefConst()

        # Write a List def_file_list which will be written as the header file.    
        self.WriteDefFileList()
        pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
        pf.List2File(self.def_file_list)
    
    def GenerateOptionsList(self):
        """Creates a default option list.

        This method of class DefineProblem will create a
        default valued options list for each entry in the
        entries list. These are essentially the options
        that will be browsed in the Pluto Setup Menu.
        """
        phylist = ['ADVECTION','HD','RHD','MHD','RMHD']
        dimlist = ['1','2','3']
        comlist = ['1','2','3']
        geolist = ['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL']
        bfolist = ['NO','VECTOR', 'POTENTIAL', '(VECTOR+POTENTIAL)']
        coolist = ['NO','POWER_LAW','TABULATED','SNEq','MINEq','H2_COOL', 'KROME']
        #parlist = ['NO','YES']
        intlist = ['FLAT','LINEAR','LimO3','WENO3','PARABOLIC']
        tmslist = ['EULER','RK2','RK3','HANCOCK','CHARACTERISTIC_TRACING']
        dislist = ['YES','NO']
        ntrlist = ['%d'%n for n in range(5)]
        udplist = ['%d'%n for n in range(32)]
        udclist = ['%d'%n for n in range(32)]

        self.options = [phylist, dimlist, comlist, geolist, bfolist,
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
             self.options[self.entries.index('INTERPOLATION')] = ['FLAT','LINEAR','LimO3', 'WENO3','PARABOLIC', 'MP5']
        
        if self.flag_dict['WITH-CHOMBO']:
            self.options[self.entries.index('GEOMETRY')] = ['CARTESIAN','CYLINDRICAL','POLAR','SPHERICAL']
            self.options[self.entries.index('INTERPOLATION')] = ['FLAT','LINEAR','WENO3','PARABOLIC']
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
            self.options[self.entries.index('INTERPOLATION')] = ['WENO3_FD', 'WENOZ_FD', 'MP5_FD','LIMO3_FD']
            self.default[self.entries.index('INTERPOLATION')] = 'WENOZ_FD'
            self.options[self.entries.index('TIME_STEPPING')] = ['RK3']
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
        self.entries_HD = ['EOS', 'ENTROPY_SWITCH', 'THERMAL_CONDUCTION', 'VISCOSITY', 'ROTATING_FRAME']
        self.default_HD = ['IDEAL', 'NO', 'NO', 'NO', 'NO']
        self.options_HD = [['IDEAL','PVTE_LAW','ISOTHERMAL'], ['NO','YES'], ['NO','EXPLICIT','SUPER_TIME_STEPPING'],
                            ['NO','EXPLICIT','SUPER_TIME_STEPPING'], ['NO','YES']]
    
        if self.flag_dict['WITH-CHOMBO']:
            self.options_HD[2] = ['NO','EXPLICIT']
            self.options_HD[3] = ['NO','EXPLICIT']

    def ProcessRHDModule(self):
        """
        Provides entries, options and defaults specific to Relativistic
        Hydro Module. Also updates them accordingly if required by flags.
        """
        self.entries_RHD = ['EOS', 'ENTROPY_SWITCH', 'USE_FOUR_VELOCITY']
        self.default_RHD = ['IDEAL', 'NO', 'NO']
        self.options_RHD = [['IDEAL','TAUB'], ['NO','YES'],  ['NO','YES']]

    def ProcessMHDModule(self):
        """
        Provides entries, options and defaults specific to Magneto-
        Hydro Module.Also updates them accordingly if required by flags.
        """
        self.entries_MHD = ['EOS', 'ENTROPY_SWITCH', 'MHD_FORMULATION', 'BACKGROUND_FIELD',
                            'RESISTIVE_MHD', 'THERMAL_CONDUCTION', 'VISCOSITY', 'ROTATING_FRAME']
        self.default_MHD = ['IDEAL','NO','EIGHT_WAVES','NO','NO','NO','NO','NO']
        self.options_MHD = [['IDEAL','PVTE_LAW','ISOTHERMAL'], ['NO','YES'],
                            ['NONE','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT'],
                            ['NO','YES'],['NO','EXPLICIT', 'SUPER_TIME_STEPPING'],
                            ['NO','EXPLICIT', 'SUPER_TIME_STEPPING'],['NO','EXPLICIT', 'SUPER_TIME_STEPPING'],['NO','YES']]
        
        if self.flag_dict['WITH-CHOMBO']:
            self.options_MHD[2] = ['NONE','EIGHT_WAVES','DIV_CLEANING']
            self.options_MHD[4] = ['NO','EXPLICIT']
            self.options_MHD[5] = ['NO','EXPLICIT']
            self.options_MHD[6] = ['NO','EXPLICIT']

        if self.flag_dict['WITH-FD']:
            self.options_MHD[2] = ['NONE','EIGHT_WAVES','DIV_CLEANING']

        if self.flag_dict['WITH-SB'] or self.flag_dict['WITH-FARGO']:
            self.options_MHD[2] = ['CONSTRAINED_TRANSPORT']
            self.default_MHD[2] = 'CONSTRAINED_TRANSPORT'
    
    def ProcessRMHDModule(self):
        """
        Provides entries, options and defaults specific to Relativisitc
        Magneto-Hydro Module.Also updates them accordingly if required by flags.
        """
        self.entries_RMHD = ['EOS', 'ENTROPY_SWITCH','MHD_FORMULATION' ]
        self.default_RMHD = ['IDEAL', 'NO', 'NONE']
        self.options_RMHD = [['IDEAL', 'TAUB'], ['NO','YES'],
                           ['NONE','EIGHT_WAVES','DIV_CLEANING','CONSTRAINED_TRANSPORT']]

        if self.flag_dict['WITH-CHOMBO']:
            self.options_RMHD[2] = ['NONE','EIGHT_WAVES','DIV_CLEANING']

    def ProcessKROMEModule(self):
        """
        Provides entries, options, and defaults for the KROME cooling module.
        """
        netwrkfiles = os.listdir(self.krome_dir+'/networks/')
        self.entries_KROME = ['NETWORK_FILE', 'USE_N','COOLING_TYPE', 'HEATING_TYPE', 'GAMMA_TYPE']
        self.default_KROME = [netwrkfiles[0], 'NO', 'NONE', 'NONE', 'DEFAULT']
        self.options_KROME = [netwrkfiles, ['NO', 'YES'], [ 'NONE', 'H2', 'ATOMIC', 'Z'], ['NONE', 'COMPRESS', 'PHOTO', 'XRAY'], ['DEFAULT','FULL', 'ROT', 'VIB', 'EXACT', 'REDUCED']]
        
        
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
            scrh_beg = pf.LocateString('User_Constants_Beg')
            k_beg   = scrh_beg[0][0]+1
            scrh_end = pf.LocateString('User_Constants_End')
            k_end   = scrh_end[0][0]-1
            const_lines = pf.ReadLines(k_beg, k_end)
            print const_lines
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
        Sets the non-user friendly constants.
        """
        tmplist1 = ['INITIAL_SMOOTHING', 'WARNING_MESSAGES', 'PRINT_TO_FILE', 'INTERNAL_BOUNDARY', 'SHOCK_FLATTENING']
        tmplist2 = len(tmplist1)*['NO']
        
        if self.flag_dict['WITH-CHOMBO']:
            tmplist1 += ['CHOMBO_EN_SWITCH','CHOMBO_REF_VAR','CHOMBO_LOGR']
            tmplist2 += [              'NO',           'ENG',       'NO']  

        if not self.flag_dict['WITH-FD']:
            tmplist1 = tmplist1 + ['ARTIFICIAL_VISCOSITY', 'CHAR_LIMITING', 'LIMITER']
            tmplist2 = tmplist2 + ['NO', 'NO', 'DEFAULT']
            
        if 'MHD_FORMULATION' in self.mod_entries:
            divb_mode = self.mod_default[self.mod_entries.index('MHD_FORMULATION')]
            if divb_mode == 'CONSTRAINED_TRANSPORT':
                tmplist1 = tmplist1 + ['CT_EMF_AVERAGE', 'CT_EN_CORRECTION', 'ASSIGN_VECTOR_POTENTIAL']
                tmplist2 = tmplist2 + ['UCT_HLL', 'NO', 'YES']
            else:
                tmplist1 = tmplist1 + ['ASSIGN_VECTOR_POTENTIAL']
                tmplist2 = tmplist2 + ['NO']
            
            if not self.flag_dict['WITH-CHOMBO']:
                tmplist1 = tmplist1 + ['UPDATE_VECTOR_POTENTIAL']
                tmplist2 = tmplist2 + ['NO']

        if 'HANCOCK' in self.default:
            if (self.phymodule == 'RMHD'):
                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
                tmplist2 = tmplist2 + ['NO']
            else:
                tmplist1 = tmplist1 + ['PRIMITIVE_HANCOCK']
                tmplist2 = tmplist2 + ['YES']
        
        if 'SUPER_TIME_STEPPING' in self.mod_default:
            tmplist1 = tmplist1 + ['STS_nu']
            tmplist2 = tmplist2 + ['0.01']

        longword = max(len(w) for w in tmplist1)
        
        if (os.path.exists(self.work_dir+'/definitions.h')):
            pf = pfIO.PlutoFiles(self.work_dir+'/definitions.h')
            pf.UpdateListFromFile(tmplist1, tmplist2)
            
        self.non_usfr = ['#define  '+tmplist1[i].ljust(longword+3)+tmplist2[i]+'\n' for i in range(len(tmplist1))]

    def AppendAdditionalFiles(self):
        """
        Adds additional object files based on
        modular defintions and requirements. 
        """
        interp_mode = self.default[self.entries.index('INTERPOLATION')]

        if interp_mode == 'LINEAR':
            self.additional_files.append('plm_states.o')
        elif interp_mode == 'PARABOLIC':
            self.additional_files.append('ppm_states.o')
            self.additional_files.append('ppm_coeffs.o')
            self.header_files.append('ppm_coeffs.h')
        elif interp_mode in ['FLAT', 'LimO3', 'WENO3']:
            self.additional_files.append('states_'+interp_mode.lower()+'.o')
	else:
	    pass
        
        if self.flag_dict['WITH-FD']:
            self.additional_files += ['states_fd.o', 'fd_reconstruct.o', 'fd_flux.o']

        if self.default[self.entries.index('COOLING')] not in ['NO', 'POWER_LAW', 'KROME']:
            self.additional_files += ['cooling_source.o','cooling_ode_solver.o']

        if self.phymodule == 'MHD' or self.phymodule == 'RMHD':
            self.additional_files.append('vec_pot_diff.o')
            if not self.flag_dict['WITH-CHOMBO']:
                self.additional_files.append('vec_pot_update.o')

        if self.flag_dict['WITH-CHOMBO']:
            if self.default[self.entries.index('DIMENSIONS')] == '1':
                self.additional_files.append('PatchCTU.1D.o')
            elif self.default[self.entries.index('TIME_STEPPING')] in ['EULER', 'RK2']:
                self.additional_files.append('PatchEuler.3D.o')
            else:
                self.additional_files.append('PatchCTU.3D.o')
        else:
            cmset = set(['CHARACTERISTIC_TRACING', 'HANCOCK']) & set(self.default)
            if len(cmset) != 0 and self.default[self.entries.index('DIMENSIONAL_SPLITTING')] == 'NO':
                self.additional_files.append('ctu_update.o')
            else:
                self.additional_files.append('rk_update.o')
        
        if 'HANCOCK' in self.default:
            self.additional_files.append('hancock.o')

        if 'CHARACTERISTIC_TRACING' in self.default:
            self.additional_files.append('char_tracing.o')

        if 'SUPER_TIME_STEPPING' in self.mod_default:
            self.additional_files += ['sts.o', 'parabolic_rhs.o']

        if 'EXPLICIT' in self.mod_default:
            self.additional_files.append('parabolic_flux.o')
        
        
    def AppendPlutoPathAndFlags(self):
        """
        Adds additional C flags and path to 'makefile' based on
        modular defintions and requirements. 
        """
        self.pluto_path.append(self.phymodule+'/')

        dis_eff = ['Thermal_Conduction', 'Viscosity']
        for de in dis_eff:
            if de.upper() in self.mod_entries and self.mod_default[self.mod_entries.index(de.upper())] != 'NO':
                self.pluto_path.append(de+'/')
        
        if self.phymodule == 'MHD' or self.phymodule == 'RMHD':
            divb_mode = self.mod_default[self.mod_entries.index('MHD_FORMULATION')]
            if divb_mode == 'CONSTRAINED_TRANSPORT':
                self.pluto_path.append('MHD/CT/')
            elif divb_mode == 'DIV_CLEANING':
                self.pluto_path.append('MHD/GLM/')
            else:
                pass

            if self.phymodule == 'MHD' and self.mod_default[self.mod_entries.index('RESISTIVE_MHD')] != 'NO':
                self.pluto_path.append('MHD/Resistive/')

        if self.flag_dict['WITH-SB']:
            self.pluto_path.append('MHD/ShearingBox/')
            self.additional_flags.append(' -DSHEARINGBOX')
            
        if self.flag_dict['WITH-FARGO']:
            self.pluto_path.append('Fargo/')
            self.additional_flags.append(' -DFARGO')

        if self.flag_dict['WITH-FD']:
            self.additional_flags.append(' -DFINITE_DIFFERENCE')

        cool_mode = self.default[self.entries.index('COOLING')]
        if cool_mode != 'NO':
            if cool_mode == 'TABULATED':
                self.pluto_path.append('Cooling/Tab/')
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
        pf = pfIO.PlutoFiles(self.work_dir+'/pluto.ini')
        scrh = pf.LocateString('[Parameters]')
        try:
            scrh[0]
        except IndexError:
            print "Parameters keyword not found in pluto.ini"
            sys.exit()
        else:
            pass
        
        ipos = scrh[0][0] + 2
        tmplist1 = pf.ReadLines(ipos,100)
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
                pf.InsertLine(x.ljust(21) + paradict[x] +'  '+cmms[self.udef_params.index(x)]+'\n', ipos)
            else:
                try:
                    cmms[self.udef_params.index(x)]
                except IndexError:
                    pf.InsertLine(x.ljust(21) + '0.0' + '\n', ipos)
                else:
                    pf.InsertLine(x.ljust(21) + '0.0'+'  '+cmms[self.udef_params.index(x)]+ '\n', ipos)
            ipos = ipos + 1
        pf.DeleteLines(ipos,ipos+100)

    def WriteKromeOpts(self):
        self.kromeoptstr = 'python '+self.krome_dir+'krome -n='+self.krome_dir+'networks/'+ self.default_KROME[self.entries_KROME.index('NETWORK_FILE')]
        if self.default_KROME[self.entries_KROME.index('USE_N')]:
            self.kromeoptstr += ' -useN'

        if self.default_KROME[self.entries_KROME.index('COOLING_TYPE')] != 'NONE':
            self.kromeoptstr += ' -cooling='+ self.default_KROME[self.entries_KROME.index('COOLING_TYPE')]

        if self.default_KROME[self.entries_KROME.index('HEATING_TYPE')] != 'NONE':
            self.kromeoptstr += ' -heating='+ self.default_KROME[self.entries_KROME.index('HEATING_TYPE')]

        if self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')] != 'DEFAULT':
            self.kromeoptstr += ' -gamma='+ self.default_KROME[self.entries_KROME.index('GAMMA_TYPE')]

    def WriteDefFileList(self):
        """
        Writes all modular entries, options, defaults into a list.
        """
        for x in self.entries:
            self.def_file_list.append('#define  '+x.ljust(21)+'   '+self.default[self.entries.index(x)]+'\n')

        self.def_file_list.append('\n/* -- physics dependent declarations -- */\n\n')
        self.phymodule = self.default[self.entries.index('PHYSICS')]
        tmp1 = ['entries_%s'%self.phymodule, 'default_%s'%self.phymodule]
        self.mod_entries = self.__getattribute__(tmp1[0])
        self.mod_default = self.__getattribute__(tmp1[1])

        for x in self.mod_entries:
            self.def_file_list.append('#define  '+x.ljust(21)+'   '+ self.mod_default[self.mod_entries.index(x)]+'\n')

        self.AppendAdditionalFiles()
        self.AppendPlutoPathAndFlags()

        # always insert user-defined paramters and constants sections

        self.def_file_list.append('\n/* -- user-defined parameters (labels) -- */\n\n')
        for x in self.udef_params:
            self.def_file_list.append('#define  '+x.ljust(21)+'   '+'%d'%self.udef_params.index(x)+'\n')
    
        self.UpdatePlutoIni()
    
                
        self.def_file_list.append('\n/* -- User_Constants_Beg --  DO NOT REMOVE THIS LINE*/\n\n')
        for i in range(len(self.udef_const)):
            self.def_file_list.append('#define  '+self.udef_const[i].ljust(21)+'   '+self.udef_const_vals[i]+'\n')
        self.def_file_list.append('\n/* -- User_Constants_End -- DO NOT REMOVE THIS LINE*/\n\n')
        

        self.def_file_list.append('\n/* -- supplementary constants (user editable) -- */ \n\n')
        self.NonUserFriendlyConst()
        for x in self.non_usfr:
            self.def_file_list.append(x)
        
        
        
        
        
        
            
        

        
        
