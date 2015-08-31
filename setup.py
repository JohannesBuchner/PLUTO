import os
import sys
import shutil

try:
  os.environ['PLUTO_DIR']
except KeyError:
  print 'PLUTO_DIR not defined. Setting it to the Current Directory'
  pluto_directory = os.getcwd()
  pass
else:
  pluto_directory = os.environ['PLUTO_DIR']

sys.path.append(pluto_directory + '/Tools/Python/')
import menu
import configure
from make_problem import MakeProblem


def PlutoInterFace(pluto_dir, do_auto_update = False):
  work_dir = os.getcwd()
  interface_optval = ''
  interface_opts = ['Setup problem', 'Change makefile', 
                    'Auto-update','Save Setup','Quit']
    
  if do_auto_update:
    MakeProblem(work_dir, pluto_dir, 1, 1)
    interface_optval = interface_opts[-1]   # set to "Quit" so it'll skip next loop

  while interface_optval != interface_opts[-1]:
    menu.SetTitle("Python setup (Aug 2015)","Working dir: "+work_dir+"\nPLUTO dir  : "+pluto_dir)

    interface_optval = menu.Browse(interface_opts)

    if interface_optval == interface_opts[0]:
      if not os.path.exists(work_dir+'/init.c'):
        shutil.copy(pluto_dir+'/Src/Templates/init.c',work_dir+'/init.c')
    
      if not os.path.exists(work_dir+'/pluto.ini'):
        shutil.copy(pluto_dir+'/Src/Templates/pluto.ini',work_dir+'/pluto.ini')

      MakeProblem(work_dir, pluto_dir, 0, 1)
    
    if interface_optval == interface_opts[1]:
      MakeProblem(work_dir, pluto_dir, 1, 0)

    if interface_optval == interface_opts[2]:
      menu.Prompt('Press Enter to Update '+work_dir)
      MakeProblem(work_dir, pluto_dir, 1, 1)
      menu.Print ("Configuration up to date",sleep=0.75)
      break

    if interface_optval == interface_opts[3]: #Save Setup
      sys.exit()

  if (menu.CursesIsActive()): menu.RestoreScreen()
  print "\n> Done."
  sys.exit()

if __name__ == "__main__":   # starts from here
  auto_update = 0
  print "\n> Checking system architecture\n"
  configure.check(pluto_directory, 1)
  for x in sys.argv[1:]:     # check argument list
    if (x == "--get-arch"):
      sys.exit(1)
      break

    if (x == "--with-chombo" or x == "--with-chombo:"): 
      print "Enabling Chombo support for AMR"
      cmset = set(['--with-fd','--with-sb','--with-fargo']) & set(sys.argv)
      if len(cmset) != 0:
        print '! Incompatible modules, ',x,' + '.join(y for y in cmset) 
        sys.exit(1)
      break

    elif (x == "--with-sb"): 
      print "Enabling support for shearing box module"
      if '--with-fd' in sys.argv:
        print '! Incompatible modules, ',x,' +  --with-fd'
        sys.exit(1)
  
    elif (x == "--with-fd"): 
      print "Enabling support for finite difference module"

    elif (x == "--with-fargo"): 
      print "Enabling support for FARGO scheme"

    elif (x == "--no-curses"):
      print ""

    elif (x == "--auto-update"):
      auto_update = 1	

    elif (x == "--help" or x == "-help"):
      print "Usage: python $PLUTO_DIR/setup.py [options]\n" 
      print "Here [options] can be:\n"
      print " --with-sb       enable the shearing box module."
      print " --with-fd       enable the finite difference module."
      print " --with-fargo    enable the FARGO-MHD module"
      print " --with-chombo   enable support for adaptive mesh refinement."
      print "                 (AMR) module using the Chombo library."
      print " --no-curses     disable ncurses library and use a"
      print "                 simpler text-based menu."
      sys.exit(1)

    else:
      print "! Unrecognized option '",x,"'"
      sys.exit(1)

  print '\n> Loading PLUTO Interface...'
  
  if auto_update == 1:
    PlutoInterFace(pluto_directory,do_auto_update=True)
  else:
    PlutoInterFace(pluto_directory)
  
  
