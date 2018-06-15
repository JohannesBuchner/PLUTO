#!/bin/csh
#note the backquote
set filelist = `ls ./*.MPI.ex` 
##module unload hdf5
##module load hdf5/parallel
set pout0 = "pout.0"
set pout1 = "pout.1"
foreach file1 ($filelist)
 set file2 = `echo $file1 | sed 's/2d\.Linux\.64\.mpiCC\.gfortran\.DEBUG\.MPI\.ex/.out/'`
## echo $file1 $file2
if(-e $file2) then
  echo "rm $file2"
  rm $file2
endif
if(-e $pout0) then
  echo "rm $pout0"
  rm $pout0
endif
if(-e $pout1) then
  echo "rm $pout1"
  rm $pout1
endif
echo "mpirun -np 2 -machinefile ~/mach $file1"
mpirun -np 2 -machinefile ~/mach $file1

if(-e $pout0) then
 echo "putting pouts into $file2"
 echo "pout.0" >  $file2
 cat $pout0 >> $file2   
 echo "pout.1" >> $file2
 cat $pout1 >> $file2   
endif

end
exit
