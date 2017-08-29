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
if test -z $1 ; then
CHID="Facility_039"
else 
    CHID=$1
fi

indir=/disk/fdsrun/Facility/$CHID
matfile=/disk/fdsmat/Facility/$CHID.mat

rm $matfile


sel_quantity=1
sel_plate=3
provide_offset='y'

xoff=-60
yoff=-60
zoff=0
dspace=1

cd $indir

slcf2mat << ieof
$CHID
$sel_quantity
$sel_plate
$provide_offset
$xoff $yoff $zoff $dspace
$matfile
ieof


