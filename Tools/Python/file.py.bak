"""               file.py   
 
 Useful utilities for file editing/manipulation

              by Andrea Mignone 

 Last modified Feb 2nd, 2008
""" 

import string
import os
import sys

###########################################################
def openr (fname):
# 
#  Open file fname for reading
#

  try:
    fl = open (fname,'r')
  except IOError:
    print " ! file.openr: Cannot open file '"+fname+"'"
    sys.exit()

  return (fl)

###########################################################
def read (fname):
# 
#  Read a file fname and return a list 
#  containing all the lines.
#

  fl    = openr(fname)
  lines = fl.readlines()
  fl.close()

  return (lines)

###########################################################
def create_file(f_name, lines):

# 
#  Creat a file f_name, 
#  where lines is a list containing the lines
#

   fl = open(f_name,'w')
   for tmp in lines: fl.write(tmp)
   fl.close()


##########################################################
def delete_lines(fname, lbeg, lend):

#
# Delete lines from lbeg to lend in 
# file fname; lbeg = 0is considered the
# first line in the file
#
   
  scrh = read (fname)
  del scrh[lbeg:lend + 1]
  create_file (fname, scrh)

##########################################################
def read_lines(fname, lbeg, lend):

#
# Read lines from lbeg to lend in 
# file fname; lbeg = 0 is considered the
# first line in the file
#

  scrh  = read(fname)
  lscrh = len(scrh)
  if (lend > lscrh): lend = lscrh 
  return (scrh[lbeg:lend])

###########################################################
def insert (fname, string, lbeg):

# 
#  Insert strings "string" at line "lbeg" in file
#  fname. 
#

  scrh = read(fname)
  scrh[lbeg:lbeg] = [string]
  create_file (fname, scrh);

###########################################################
def replace_line(fname, newstring, nline):

# 
#  Replaces line number nline with newline
#

  scrh = read(fname)
  del scrh[nline]
  scrh[nline:1] = [newstring]

  create_file (fname, scrh);

###########################################################
def string_list(f_name, str):

#
# Return a list of all the lines in file f_name
# where str occurs at least once as a SUBSTRING.
#

  return find(f_name,str,action='string',want='list')

###########################################################
def string_find(f_name, str):

#
# Return a list of all the line numbers in file f_name
# where str occurs at least once as a SUBSTRING.
#

  return find(f_name,str,action='string',want='enum')

###########################################################
def word_list(f_name, str):

#
# Return a list of all the lines in file f_name
# where str occurs at least once as a SEPARATE WORD.
#


  return find(f_name, str,action='word',want='list')

###########################################################
def word_find(f_name, str):

#
# Return a list of all the line numbers in file f_name
# where str occurs at least once as a SEPARATE WORD.
#


  return find(f_name, str,action='word',want='enum')


###########################################################
def find(fname, str, action, want):

  scrh = read(fname)
     
  line_str = []
  line_num = []
  ipos     = -1
  for tmp in scrh:
    ipos = ipos + 1
    x    = string.split(tmp)
    if (action == 'string'):
 
      if (string.find(tmp,str) >= 0): 
        line_str.append(tmp) 
        line_num.append(ipos) 

    elif (action == 'word'):

      if (x.count(str) > 0): 
        line_str.append(tmp) 
        line_num.append(ipos)

      
  if (want == 'list' ): 
    return line_str
  elif (want == 'enum'):
    return line_num

###########################################################
def count_lines(fname):

#
# Return the number of lines in file f_name
#

  scrh = read(fname)
  nlines = 0
  for x in scrh:
    if (len(string.split(x)) != 0): 
      nlines = nlines + 1

  return nlines

###########################################################
def count_words(f_name, i_line, f_line = -1):

#
# Return a list with all the words 
# from i_line (included) to f_line;
# if f_line is not given, end of 
# file is assumed.
#

   try:
     fl = open(f_name,'r')
   except IOError:
     print "Error in count_words [file.py]:"
     print "Cannot open file '"+f_name+"'"
     sys.exit()
   
# go to line i_line
   
   tmp = []
   for i in range(i_line): fl.readline()
   if (f_line < 0):
      tmp = fl.readlines()
   else:
      for i in range(f_line-i_line+1):
         tmp.append(fl.readline())
   

   scrh = []
   for yy in tmp:
      this_line = string.split(yy)
      if (len(this_line) != 0): 
         for x in this_line:    
          scrh.append(x)

   return scrh

   
   
    
##########################################################
def swap_lines(fname, l1, l2):

#
# swap lines l1 and l2 in file fname
#

  scrh = read (fname)
  
  line_1 = scrh[l1]
  line_2 = scrh[l2]
  
  replace_line (fname, l1, line_2)
  replace_line (fname, l2, line_1)
