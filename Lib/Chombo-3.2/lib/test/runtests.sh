#!/bin/sh
#
# Runs every executable that's been installed in $(prefix)/libexec/Chombo/tests,
# and reports its exit status.
#

for d in `find . -type d | egrep -v '^\.$|CVS'`; do
    (cd $d

    # Skip ccse tests, except facade (either that, or put the ccse inputs
    # files back in where the executables expect them).
    echo $d | grep ccse > /dev/null
    if test $? -eq 0 ; then
        echo $d | grep facade > /dev/null
        if test $! -ne 0 ; then
            continue
        fi
    fi
    
    echo "Entered directory $d"
    for f in `find . -type f -printf "%f "`; do
        perm=`find . -name $f -printf "%m"`
        if [ $perm = "755" ]; then
            ./$f
            echo "$f finished with status $?"
            rm -f `find . -name "*.hdf5"`
        fi
    done)
done
