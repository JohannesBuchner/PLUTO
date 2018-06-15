#
# Go through all the .cpp and .H files in a directory.  Make a copy of each file
# (adding a suffix ".new"), in which after that last #include, insert
# '#include "NamespaceHeader.H"', and before the closing include guard, insert
# '#include "NamespaceFooter.H"'.

# For most files (over 90% of .H and .cpp files other than BoxTools)
# this is all it takes.  For the rest, you'll need to do some
# hand-editing to get the NamespaceHeader and NamespaceFooter in the
# right place.  You need to do such hand editing whenever the file doesn't
# conform to the following rules:
# 1. All forward class declarations and "using namespace std" (and the like)
#    have to go after NamespaceHeader.H.
# 2. Conditional compilation directives (e.g. #ifdef CH_USE_MEMORY_TRACKING,
#    #ifdef CH_MPI) can't hide the namespace header and footer.

# I've been running this script from the lib/src/* directories, using the
# following command:
#  for f in *.cpp *.H; do
#    echo "*** $f ***"
#    python ../../util/multidim/namespace_inserter.py $f
#  done
#
# Scan the output for problems.  If there are none, do a 
#   for f in *.cpp *.H; do mv ${f}.new $f; done
# and then try to build the library.  If it doesn't build usually it's
# because of a violation of rules (1) or (2) above.

import sys
import os
import re

def findLastInclude( filename_or_list, chombo_header_directory=False ):
    """
    Find the last #include in a file, or in a list of strings.  Return its
    line number (counting from zero).

    If optional arg chombo_header_directory is set to something, then in
    looking for that "last #include" we ignore #include's of files found in
    chombo_header_directory.  That's because we can count on those files having
    the proper header and footer that prevents namespace nesting (and which we
    of course can't count on from any other header files).

    If no #include's are found at all, returns -1.
    """
    headers = []
    if chombo_header_directory:
        headers = os.listdir( chombo_header_directory )
        #print "headers=", headers

    if type(filename_or_list) == type([]):
        lines = filename_or_list
    else:
        lines = open(filename).readlines()

    line_num = -1

    compiled_regex = re.compile('^# *include ')

    for i in range(0,len(lines)):
        line = lines[i]
        if compiled_regex.search( line ):
            if len(headers) > 0:
                tokenized_line = line[line.find('include'):].split()
                if len(tokenized_line) == 1:
                    sys.stderr.write("Error: #include is last thing on line!\n")
                    sys.exit(1)
                the_file = tokenized_line[1].strip('"<> ')
                if not the_file in headers:
                    line_num = i
            else:
                line_num = i
    return line_num


def blankoutComments( filename ):
    """
    Replace a file with blanks everywhere the file had C/C++ comments.
    Return the resulting file as a list of strings.
    """
    tmp_file = '/tmp/inserter.' + str(os.getpid())
    os.system(   'stripcomments.py infile=' + filename + '  language=cpp '
               + 'replace=True > ' + tmp_file )
    lstr = open(tmp_file).readlines()
    os.unlink( tmp_file )
    return lstr


def findIncludeGuards( filename_or_list ):
    """
    Return, as a tuple, the positions of the include guard (the second line of
    it -- the one that says #define INCLUDED_FOO_H) and the closing #endif.

    If include guards not found, return None.
    """
    if type(filename_or_list) == type([]):
        lines = filename_or_list
    else:
        lines = open(filename).readlines()

    compiled_regex = re.compile('^ *\t*#ifndef _[A-Z]|^ *\t*#ifndef BL_[A-Z]')

    foundit = False
    i = 0
    while not foundit and  i < len(lines):
        if compiled_regex.search( lines[i] ):
            foundit = True
            break
        else:
            i += 1

    if not foundit:
        return None

    guard_position = i+1

    #
    # Check that the line right after the presumed include guard looks like
    # what we expect.  If not, 
    guard_name = lines[guard_position-1].split()[1]
    if(  (len(lines) < guard_position+1)
      or (len(lines[guard_position].split()) < 2)
      or (lines[guard_position].split()[1] != guard_name) ):
        sys.stderr.write( "Include guard lines not consecutive, or I have\n"
            "mistaken something other than an include guard for one.\n"
            "Check file.  Exiting now...\n" )
        sys.exit(1)

    # Find last #endif in the file -- the #endif that matches the include guard.
    i = len(lines) - 1
    while lines[i].find('#endif') == -1:
        if i == guard_position:
            sys.stderr.write( "No matching #endif for include guard!?\n" )
            sys.exit(1)
        i -= 1
    closing_endif = i
    
    return (guard_position, closing_endif)


def printStats( filename ):
    lines = blankoutComments( filename )
    if len(lines) == 0:
        print "Empty file."
        return

    compiled_regex = re.compile('^ *\t*\n$')
    code_lines = filter(lambda x: not compiled_regex.search(x), lines)

    """
    print "len(lines)=", len(lines)
    print "len(code_lines)=", len(code_lines)
    """
    """
    print "******* code_lines:"
    for l in code_lines:
        print l[:-1]
    print "******* lines:"
    for l in lines:
        print l[:-1]
    """

    own_dir = os.path.dirname( filename )
    if own_dir == '':
        own_dir = '.'

    last_n = findLastInclude( lines, chombo_header_directory=own_dir )
    print last_n, str("(" + str(100.0*last_n/len(lines)) + "%)")

    include_guard_locations = findIncludeGuards( lines )
    print "Include guards found at lines", include_guard_locations
    if (include_guard_locations != None) and (last_n != -1):
        assert( last_n > include_guard_locations[0] )

    last_n = findLastInclude( code_lines, chombo_header_directory=own_dir )
    try:
        print last_n, str("(" + str(100.0*last_n/len(code_lines)) + "%)")
    except:
        print "Zero code lines"


def editFile( filename ):
    """
    Insert '#include "NamespaceHeader.H"' and '#include "NamespaceFooter.H"' at what we hope are
    the correct locations.

    We don't actually change the file; we create a new file for the edited
    version -- same filename + '.new'.
    """
    lines = blankoutComments( filename )
    if len(lines) == 0:
        print "Empty file."
        return

    edited_filename = filename + '.new'

    # If file already has '#include "NamespaceHeader.H"' in it, make the ".new" version
    # of it just a copy of the original.
    if len( filter( lambda i: i.find('NamespaceHeader.H') != -1, lines ) ) > 0:
        print "File already contains '#include NamespaceHeader.H'."
        os.system( 'cp ' + filename + ' ' + edited_filename )
        return

    own_dir = os.path.dirname( filename )
    if own_dir == '':
        own_dir = '.'
    last_n = findLastInclude( lines, chombo_header_directory=None ) #=own_dir
    print last_n, str("(" + str(100.0*last_n/len(lines)) + "%)")

    include_guard_locations = findIncludeGuards( lines )
    print "Include guards found at lines", findIncludeGuards( lines )
    if (include_guard_locations != None) and (last_n != -1):
        assert( last_n > include_guard_locations[0] )

    # Now do the editing.
    rawlines = open(filename).readlines()
    outfile = open( edited_filename, 'w' )

    if include_guard_locations:
        first_stop = max( last_n+1, include_guard_locations[0]+1 )
    else:
        first_stop = last_n+1
    for i in range(0, first_stop):
        outfile.write( rawlines[i] )

    outfile.write( '#include "NamespaceHeader.H"\n' )

    if include_guard_locations:
        second_stop = include_guard_locations[1]
    else:
        second_stop = len(rawlines)

    for i in range( first_stop, second_stop ):
        outfile.write( rawlines[i] )
    outfile.write( '#include "NamespaceFooter.H"\n' )
    for i in range( second_stop, len(rawlines) ):
        outfile.write( rawlines[i] )
 

if __name__ == '__main__':
    filename = sys.argv[1]
    printStats( filename )
    editFile( filename )
