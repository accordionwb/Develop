#!/bin/bash
# This script make the project and execute the test
#
#
# Dynamic Library PATH Check 
#
if [ ! `echo $LD_LIBRARY_PATH | grep MATLAB` ] ; then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/MATLAB/R2017a/bin/glnxa64
fi

# Make the program
make

#
# Part of I/O file name variables
#
instoredir=/home/wangbing/Public
allruns=`ls $instoredir`
outstoredir=/home/wangbing/fdscov

for run in $allruns
do
   allfiles=`ls $instoredir/$run/*.bin`
   CHID=$run
   outfilename=$outstoredir/$CHID.mat
   for file in $allfiles
   do
      infilename=$file
bin/fds2mat << ieof
$infilename
$outfilename
ieof
   done
done



