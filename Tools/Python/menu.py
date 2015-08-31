################################################################
#
# menu.py provides terminal control capabilities to browse 
# through a user-defined menu.
# It employs the Python curses library.
# An alternative (text-based) menu is also possible if the
# --no-curses is given in the command line options.
#
# Last modified:
#
#   Tuesday June 22, 2010 by A. Mignone (mignone@ph.unito.it)
#
################################################################

import os, sys, traceback, time, string
have_curses = 1
for x in sys.argv:  # avoid curses library with the --no-curses option.
  if (x == "--no-curses"):
    have_curses = 0

if (have_curses == 1):
  import curses, curses.textpad

#####################################################
#
# The class gb contains global variables and 
# pointers to ease up accessibility in coding the 
# functions of this module 
#
#####################################################
class gb:
  scrn = None  # a pointer to a window object
  rbeg = 6     # starting row
  cbeg = 1     # starting column
  row  = 0     # the current row
  ncol = 1     # number of columns in the menu
  csep = 30    # separation between multi-columns
  title    = 'No Title'
  subtitle = ''
  init = 0     # a value of 0 means that curses functionality has not
               # been activated (i.e. curses.initscr() has never been called)

 
#####################################################
#
# Sets the title (but does not print it)
# that will be used by ShowMenu
#
#####################################################
def SetTitle (title, subtitle = ''):
  
  gb.title    = title
  gb.subtitle = subtitle

#####################################################
#
# Print a message on the screen 
#
#####################################################
def Print (message, sleep=0.7,row=1):

  for x in sys.argv:  # avoid curses library with the --no-curses option.
    if (x == "--no-curses" or gb.init == 0):
      Print_no_curses(message,sleep,row)
      return
  
  if (gb.scrn == None): return   # need this when using test_pluto.py script
  if (row == 1): gb.scrn.erase()
  gb.scrn.addstr(row,1,message, curses.A_BOLD)
  gb.scrn.refresh()
  time.sleep(sleep)

#####################################################
#
# Prompt a message, wait for any key to be pressed
#
#####################################################
def Prompt (message):

  for x in sys.argv:  # avoid curses library with the --no-curses option.
    if (x == "--no-curses"):
      Prompt_no_curses(message)
      return

  gb.scrn.erase()
  gb.scrn.addstr(1,1,message, curses.A_BOLD)
  gb.scrn.refresh()
  c = gb.scrn.getch()

#####################################################
#
#  Show menu using entries (1st column) and
#  default (2nd column).
#  row is an optional argument giving the line
#  to be highlighted (default is first line)
#   
#####################################################
def ShowMenu(entries, default, row=0):
  # display title

  gb.scrn.clear()
  gb.scrn.addstr(0,0,">> "+gb.title+" <<", curses.A_BOLD)
  gb.rbeg = 3
  if (len(gb.subtitle) > 1): 
    gb.scrn.addstr(2,0, gb.subtitle, curses.A_UNDERLINE)
    gb.rbeg = 6

  lastrow = gb.rbeg
  for ln in entries[0:len(entries)]:
    indx = entries.index(ln)
    gb.scrn.addstr(lastrow, gb.cbeg, ln)
    if (gb.ncol == 2): 
      gb.scrn.addstr(lastrow, gb.cbeg + gb.csep, default[indx])
    lastrow += 1

  if (row == 0 or gb.row < gb.rbeg or gb.row > lastrow): 
    gb.row = gb.rbeg  # initial position

  n = gb.row - gb.rbeg
  gb.scrn.addstr(gb.row, gb.cbeg        , entries[n], curses.A_REVERSE)
  if (gb.ncol == 2): gb.scrn.addstr(gb.row, gb.cbeg+gb.csep, default[n], curses.A_UNDERLINE)
  gb.scrn.refresh()


#####################################################
#
# Allow the cursor to move up and down in the list
#
#####################################################
def UpDown(entries, default, inc):
  tmp = gb.row + inc

  # ignore attempts to go off the edge of menu

  if tmp >= gb.rbeg and tmp < (gb.rbeg + len(entries)):
    # unhighlight the current line by rewriting it in default attributes
    gb.scrn.addstr(gb.row, gb.cbeg          , entries[gb.row-gb.rbeg])
    if (gb.ncol == 2): gb.scrn.addstr(gb.row, gb.cbeg + gb.csep, default[gb.row-gb.rbeg])
    # highlight the previous/next line
    gb.row = tmp
    c1 = entries[gb.row-gb.rbeg]
    if (gb.ncol == 2): c2 = default[gb.row-gb.rbeg]
     
    gb.scrn.addstr(gb.row, gb.cbeg          , c1, curses.A_REVERSE)
    if (gb.ncol == 2): gb.scrn.addstr(gb.row, gb.cbeg + gb.csep, c2, curses.A_UNDERLINE)
    gb.scrn.refresh()


#####################################################
# 
# Allow left/right keys to switch options in the 
# second column and change default values
#
#####################################################
def LeftRight(entries, default, options, inc):

  i    = gb.row - gb.rbeg 
  idef = options[i].index(default[i])
  nopt = len(options[i])
  if (inc > 0): idef = idef + 1
  if (inc < 0): idef = idef - 1
  if (idef < 0):     idef = nopt-1
  if (idef == nopt): idef = 0

  default[i] = options[i][idef]
  gb.scrn.addstr(gb.row, gb.cbeg          , entries[i], curses.A_REVERSE)
  gb.scrn.addstr(gb.row, gb.cbeg + gb.csep, default[i], curses.A_UNDERLINE)
  gb.scrn.clrtoeol()
  gb.scrn.refresh()

#####################################################
#
#  Browse a menu with entries (1st column) and
#  default (2nd column, optional)
#  Note: with Python > 2.5 we had some troubles
#        initializing curses more than once. 
#        For this reason we prefer to initialize
#        curses only at the beginning.
#        
#####################################################
def Browse(entries, default=[], options=[]):

  gb.ncol = 1
  if (len(default) > 0): gb.ncol = 2

  for x in sys.argv:  # avoid curses library with the --no-curses option.
    if (x == "--no-curses"):
      return Browse_no_curses(entries, default, options)

#
# window setup will be done just once.
# 
  if (gb.init == 0):
    gb.scrn = curses.initscr()
    curses.noecho()
    curses.cbreak()
    gb.scrn.keypad(1)
    gb.init = 1
  ShowMenu(entries, default)
  while True:
    # get user command
    c  = gb.scrn.getch()
    try:    cc = chr(c)
    except: cc = 0

    if (c == 10): 
#      RestoreScreen()
      return entries[gb.row-gb.rbeg]
    elif (cc == 'q'): 
      RestoreScreen()
#      curses.reset_shell_mode()
      sys.exit()
    elif (cc == 'u' or c == curses.KEY_UP):   UpDown(entries, default, -1)
    elif (cc == 'd' or c == curses.KEY_DOWN): UpDown(entries, default, 1)
    elif (gb.ncol > 1): 
      if (cc == 'r' or c == curses.KEY_RIGHT):
        LeftRight(entries, default, options, 1)
      elif (cc == 'l' or c == curses.KEY_LEFT):  
        LeftRight(entries, default, options, -1)

#####################################################
#
# Similar to Browse, but allow the user to directly
# input the default values by a reading a string
#
#####################################################
def Insert(entries, default):

  gb.ncol = 2

  for x in sys.argv:  # avoid curses library with the --no-curses option.
    if (x == "--no-curses"):
      return Insert_no_curses(entries, default)

  #
  # window setup will be done just once.
  # 

  if (gb.init == 0):
    gb.scrn = curses.initscr()
    curses.noecho()
    curses.cbreak()
    gb.scrn.keypad(1)
    gb.init = 1

#  entries = []
#  for n in range(num): entries.append(repr(n))

#  RestoreScreen()
#  print default
#  sys.exit()
   
  ShowMenu(entries, default)
  while True:
    c  = gb.scrn.getch()   # get user command
    try:    cc = chr(c)
    except: cc = 0

    if (c == 10): 
      return 
    elif (cc == 'q'): 
      RestoreScreen()
      sys.exit()
    elif (cc == 'u' or c == curses.KEY_UP):   UpDown(entries, default, -1)
    elif (cc == 'd' or c == curses.KEY_DOWN): UpDown(entries, default, 1)
    elif (cc == 'r' or c == curses.KEY_RIGHT):
      curses.echo()
      gb.scrn.addstr(gb.row,gb.cbeg+gb.csep,'                                 ')
      gb.scrn.addstr(gb.row,gb.cbeg+gb.csep,'NAME or VALUE > ',curses.A_UNDERLINE)
      new_name = gb.scrn.getstr()
      i = gb.row-gb.rbeg
      default.pop(i)
      default.insert(i,new_name)
      curses.noecho()
      gb.scrn.clrtoeol()
      ShowMenu(entries,default, gb.row)

   
#####################################################
#
# Restore screen back to shell functionality.
# Note that RestoreScreen should be followed by 
# sys.exit() in order to avoid troubleshooting 
# observed with Python > 2.5
#
#####################################################
def RestoreScreen():

  for x in sys.argv:  # avoid curses library with the --no-curses option.
    if (x == "--no-curses"):
      return

  curses.reset_shell_mode()
  curses.nocbreak()
  gb.scrn.keypad(0)
  curses.echo()
  curses.endwin()

if __name__ == '__browse__':
  try:
    browse()
  except:
    RestoreScreen()
    # print error message re exception
    traceback.print_exc()

#####################################################
#
# Return 1 if curses have been activated
#
#####################################################
def CursesIsActive():

  return gb.init
  
#####################################################
#
#  The next set of functions replicate the previous 
#  ones without using curses library.
#  They are intended to provide a simpler way select
#  options through a terminal-based replacement.
#
#####################################################
def Print_no_curses(message, sleep, row):

  global xglb
#  if (row == 1): os.system("clear")
  print message
  time.sleep(sleep)

######################################################
def Prompt_no_curses (message):
#
#
######################################################

  os.system("clear")
  print message
  q = raw_input()

######################################################
def Browse_no_curses(entries,  default, options):
#
#
######################################################

  q = "c"
  while (q != ''):
    os.system("clear")
    print ">> ",gb.title+"\n"
    for x in entries:
      i = entries.index(x)
      if (len(default) > 0): 
        print str(i).rjust(2),') ',x.ljust(28), default[i]
      else:
        print str(i).rjust(2),') ',x.ljust(28)

    print " "
    q = raw_input(">> choice ? ")
    if (q == ''):
      print "Enter"
    else:
      try:
        q = int(q)
        if (len(default) == 0): return entries[q]
      except:
        continue

      opt_list = ''
      for x in options[q]:
        i = options[q].index(x)
        opt_list += repr(i)+") "+x+"   "
    
      print "\n"+entries[q]+": ",opt_list
      c = raw_input(">> choice ["+default[q]+"] ? ")
      try: 
        c = int(c)
      except:
        continue

      default[q] = options[q][c]

  return
######################################################
def Insert_no_curses(entries, names):
#
#
######################################################

  q = "c"
  while (q != ''):
    os.system("clear")
    print ">> ",gb.title+"\n"
    for x in entries:
      i = entries.index(x)
      print str(i).rjust(2),') ',names[i].ljust(28)

    print " "
    q = raw_input(">> choice ? ")
    if (q == ''):
      print "Enter"
    else:
      try:
        q = int(q)
      except:
        continue

      newname = raw_input(">> new name ? ")
      names[q] = newname

  return
