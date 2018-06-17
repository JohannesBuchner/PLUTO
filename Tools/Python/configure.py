import os
import sys
import pluto_files_IO as pfIO


def check(pluto_dir, get_arch):

  work_dir = os.getcwd() 
  print ("\n>> Inspecting system configuration << \n")

  log_file     = work_dir+'/sysconf.out'
  python_dir   = pluto_dir+'/Tools/Python'
  param_file   = pluto_dir+'/Src/pluto.h'

  os.chdir (python_dir)
  log = []

 #  ---------------------
 #    get user name
 #  ---------------------

  try:
    user = os.getlogin()
  except OSError:
    user = 'unknown'

 #  ---------------------
 #    get system spec.
 #  ---------------------

  PLATFORM = os.uname()
 
 #  ---------------------
 #   find current 
 #   version of PLUTO
 #  ---------------------

  pf = pfIO.PlutoFiles(param_file)
  scrh = pf.LocateString("PLUTO_VERSION")
  ipos = scrh[0][0]
  vers = scrh[0][1].split()[-1]

 # ------------------------------
 #  print a short system summary 
 # ------------------------------
 
  print ("User.......................",user)
  print ("System name................",PLATFORM[0])
  print ("Node name..................",PLATFORM[1])
  print ("Release....................",PLATFORM[2])
  print ("Arch.......................",PLATFORM[4])
  print ("Byteorder..................",sys.byteorder)
  print ("Version....................",PLATFORM[3])
  print ("Working_dir................",work_dir)
  print ("PLUTO main dir.............",pluto_dir)
  print ("PLUTO version..............",vers)

  
# --------------------------------------
#  Check for gcc or another c compiler
# --------------------------------------

  compiler_list = ['gcc','cc','gcc2']
  COMPILER_NAME = ''
 
  for x in compiler_list:
    if (CHECK_FOR(x) == 'YES'):
      print ("C Compiler................. ",x)
      COMPILER_NAME = x
      break

  if (COMPILER_NAME == ''):
    print ("! Can not find a C compiler       !")

# -----------------
#  check for mpi 
# -----------------

  mpi_compiler_list = ['mpicc','mpiCC','mpcc_r','hcc','mpcc']
  MPI_COMPILER_NAME = ''

  for x in mpi_compiler_list:
    if (CHECK_FOR(x) == 'YES'):
      print ("MPI Compiler .............. ",x)
      MPI_COMPILER_NAME = x
      break
   
  if (MPI_COMPILER_NAME == ''):
    print ("MPI Compiler............... NOT FOUND")


  if (get_arch):
    print ("\n")
    print ("Proposed makefile names: \n")
    print ("> "+PLATFORM[0]+"."+PLATFORM[4]+"."+COMPILER_NAME+".defs")
    print ("> "+PLATFORM[0]+"."+PLATFORM[4]+"."+MPI_COMPILER_NAME+".defs")

# ---------------------------------------------------
#  Build log list, that will be compared to the 
#  sysconf.out file.
# ---------------------------------------------------

  log.append("USER           = "+user+'\n')
  log.append("WORKING_DIR    = "+work_dir+"\n")
  log.append("SYSTEM_NAME    = "+PLATFORM[0]+"\n")
  log.append("NODE_NAME      = "+PLATFORM[1]+"\n")
  log.append("RELEASE        = "+PLATFORM[2]+"\n")
  log.append("ARCH           = "+PLATFORM[4]+"\n")
  log.append("BYTE_ORDER     = "+sys.byteorder+"\n")
  log.append("VERSION        = "+PLATFORM[3]+"\n")
  log.append("PLUTO_DIR      = "+pluto_dir+'\n')
  log.append("PLUTO_VERSION  = "+vers+'\n')
  log.append("C_COMPILER     = "+COMPILER_NAME+'\n')
  log.append("MPI_C_COMPILER = "+MPI_COMPILER_NAME+'\n')

 #  ----------------------------
 #    check for online updates
 #  ----------------------------
 
  print ("\n> Checking for updates (canceled)...\n")
# try:
#   urllib.urlretrieve("http://plutocode.ph.unito.it/updates.txt","updates.txt")
#   scrh = file.word_find ("updates.txt","release")
#   ipos = scrh[0]
#   scrh = file.read_lines ("updates.txt", ipos, ipos + 1)
#   rels = string.split(scrh[0])
#   if (rels[1] != vers):
#     print "   ******************************************************* "
#     print "    A new version of PLUTO ("+rels[1]+") is available at"
#     print "    http://plutocode.oato.inaf.it"
#     print "   *******************************************************\n"
#     scrh = raw_input("> Press enter to continue.")
#     os.chdir(work_dir)   
# except:
#   print "! Connection not available\n"

 
# ------------------------------------------------
#  Compare the list 'log' with the file log_file;
#  
#   - if they match, no update is necessary, 
#               ==> return to main menu
#
#   - if they do not match or if log_file does not
#     exists, create a new one
#
# ------------------------------------------------

  if (os.path.exists(log_file)):
    pf = pfIO.PlutoFiles(log_file)
    scrh = pf.ReadLines(0,128)
    if (scrh[0:] == log[0:]):
      os.chdir(work_dir)
      return
    else:
      print ("\n> System configuration file is not up to date. Updating...")
  else:
    print ("\n> System configuration file not found, creating one...")
    pf = pfIO.PlutoFiles(log_file)
    pf.List2File(log)

# ------------
#  Make Tools
# ------------

# print (" > Making binary tools in "+bintools_dir+"...") 
# os.chdir(bintools_dir)
# os.system('make -s clean')
# os.system('make -s dbl2flt')
# os.system('make -s bin2ascii')
# if (HAVE_LIBPNG == 'YES'):
#   os.system('make -s bin2png')


# ---------------------------
#   Add important info here
# ---------------------------

# scrh = raw_input(" > Press enter to continue.")
  os.chdir(work_dir)
  return

######################################################
def CHECK_FOR (file):

#
#  find whether file can be found in 
#  the user's path
#
#
#######################################################

# 
# determine user's path
#
 
  scrh = os.getenv('PATH')
  path = scrh.split(':')

#
# search file 
#
 
  have_file = "NO"
  for x in path:
    if (os.path.exists(x+"/"+file)):
      have_file = "YES"
     
  return have_file
  
