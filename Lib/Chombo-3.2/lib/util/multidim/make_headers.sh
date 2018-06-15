#!/bin/sh

#
# Generates special header files that go into the install
# include directory.  These are:
# 1. Copies of all the .H files but with .2D and .3D suffixes
#    and dimension-specific include guards.
# 2. Generates, for each .H file, a .multidim analog which
#    #include's the .2D and .3D versions after setting CH_SPACEDIM
#    appropriately for each one.
# 3. Generates a gutted version of each .H file, which just
#    delegates to the .2D or .3D version of that header.
#
# This script is run from Chombo/lib
#

if test $# -ne 4; then
    echo "*********************************************************************"
    echo "Usage: $0 srcdir pkgincludedir mindim maxdim                         "
    echo "  where srcdir is, like, lib/src/BoxTools,                           "
    echo "  pkgincludedir is, like, prefix/include/Chombo,                     "
    echo "  mindim and maxdim are the dimensions we want to support (e.g. 1 3  "
    echo "*********************************************************************"
    exit 1
fi

#
# Assign cmd-line args to variables.
#
curdir=`pwd`
srcdir=$1
cd $srcdir
# (dfm 4/28/08) -- the original version of this script had this cd
# however, for the current usage pattern, it's no longer necessary
# (or even correct)
#cd ..
abs_srcdir=`pwd`
cd $curdir
pkgincludedir=$2
mindim=$3
maxdim=$4

utildir=`dirname $0`
mkdir -p $pkgincludedir/multidim

#
# Generate patterns for .H.multidim and .H.XD headers.
#
echo '#ifdef CH_SPACEDIM' >  MultidimHeader.pat
echo '#undef CH_SPACEDIM' >> MultidimHeader.pat
echo '#endif'             >> MultidimHeader.pat
echo "#if  CH_SPACEDIM == 0" >  HeaderRedirector.pat
echo "Ayeeee               " >> HeaderRedirector.pat
d=$mindim
while test $d -le $maxdim ; do
  echo "#define CH_SPACEDIM $d" >> MultidimHeader.pat
  echo "#ifdef CH_USE${d}D " >> MultidimHeader.pat
  echo '#include "HEADERFILE"'  >> MultidimHeader.pat
  echo '#endif'                  >> MultidimHeader.pat
  echo '#undef CH_SPACEDIM'     >> MultidimHeader.pat

  echo "#elif CH_SPACEDIM == $d"    >> HeaderRedirector.pat
  echo "#include \"HEADERFILE.${d}D"\" >> HeaderRedirector.pat

  d=`expr $d + 1`
done
echo '#endif' >> HeaderRedirector.pat


# 
# Install headers that don't need dim-specific versions
#
while read headerfile; do
    cp -p $abs_srcdir/$headerfile $pkgincludedir/multidim
done < $abs_srcdir/multidim/dim-independent-headers.txt

for f in $abs_srcdir/*.H; do
    bn=`basename $f`
    newFile=true
    if test -f $pkgincludedir/multidim/$bn; then
        if test installed_multidim.flag -nt $f; then
            newFile=false
            continue
        fi
    fi

    if $newFile; then 
        echo "*** installing $bn ***"
        grep `basename $f` $abs_srcdir/multidim/dim-independent-headers.txt > /dev/null
        if test $? -eq 0; then
            continue
        fi
        regex1='#ifndef  *_[A-Z][A-Z0-9_][A-Z0-9_]*_H_'
        regex2='#define  *_[A-Z][A-Z0-9_][A-Z0-9_]*_H_'
        if test `grep -c "$regex1" $f` -gt 1; then
            echo "Error: found more than one include-guard-like thing in $f."
            echo "Was looking for things matching regex '${regex1}'.  Please fix $f."
            exit 1
        fi
        line1=`grep -n "$regex1" $f | awk -F":" '{print $1}'`
        line2=`grep -n "$regex2" $f | awk -F":" '{print $1}'`
        if test `expr $line2 - $line1` -ne 1; then
            echo "Error: include guard pair are not on adjacent lines, in $f."
            echo "Please fix $f."
            exit 2
        fi
        
        bf=`basename $f`
        dim=$mindim
        while test $dim -le $maxdim ; do
            sed -e "s/\(${regex1}\)/\1${dim}D/" -e "s/\(${regex2}\)/\1${dim}D/" $f > $pkgincludedir/multidim/${bf}.${dim}D
            dim=`expr $dim + 1`
        done
        
        sed "s/HEADERFILE/$bf/" MultidimHeader.pat > $pkgincludedir/multidim/${bf}.multidim
        sed "s/HEADERFILE/$bf/" HeaderRedirector.pat > $pkgincludedir/multidim/$bf
    fi 
done

touch installed_multidim.flag
