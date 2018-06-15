#!/bin/sh

./sliceTest.*.ex >& /tmp/MultiDim_Slicing.out
cmp canonical.out /tmp/MultiDim_Slicing.out
if test $? -eq 0; then
    echo "Success"
else
    echo "Failure: diff canonical.out /tmp/MultiDim_Slicing.out"
fi
