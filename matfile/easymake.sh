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
# Single run with given CHID
CHID='Facility_008'
indir=/home/tmp/fdscov
outdir=/disk/fdsmat/Facility/

if ! test -d $outdir ; then
    mkdir -p $outdir
fi

allfiles=`ls $indir/$CHID/*.bin`
outfilename=$outdir/$CHID.mat
for file in $allfiles
do
    infilename=$file
    bin/fds2mat << ieof
$infilename
$outfilename
ieof
done



