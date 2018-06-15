#!/bin/sh

#
# Tests all header files
# Meant to be run from subdirectories of lib/src (e.g. BoxTools).
# Usage: "./test_all_headers.sh"
#
# Reads headers from a complete list in lib/src/*/multidim/tested_headers.txt,
# skipping any headers preceded by a '#'.
#

curdir=`dirname $0`
while read header; do
    echo $header | egrep '^#|^ *$' > /dev/null
    if test $? -eq 0; then
        continue
    fi
    echo "Testing $header ..."
    ./$curdir/test_header.sh $header
    if test $? -ne 0; then
        exit 1
    fi
done < ./multidim/tested_headers.txt
