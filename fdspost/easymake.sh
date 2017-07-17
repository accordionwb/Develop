#!/bin/bash
make intel_linux_64
rm -rf ~/intel

batchrun=1
binary=fds2slcf
single=1

Dspace=1.

# This file try to make fdspost usage easier
# fdspost input explanation:
#
# L1: Enter Job ID string (CHID):
# 
# L2: What type of file to parse?
#        PL3D file? Enter 1
#        SLCF file? Enter 2
#        BNDF file? Enter 3
# 
# L3: Enter Sampling Factor for Data?
#   (1 for all data, 2 for every other point, etc.)
#
# L4: Domain selection:
#   y - domain size is limited
#   n - domain size is not limited
#   z - domain size is not limited and z levels are offset
#   ya, na or za - slice files are selected based on type and location.
#       The y, n, z prefix are defined as before.
#
# L5: Enter starting and ending time for averaging (s)
# 
# ---- List of Variables and bounds ----
# 
# L6: How many variables to read:
# 
# L7: Enter index for variable 1,2,3,...
#
# L8: Enter output file name:
# 
#======================================================
#
#------------------------------------------
# Script single process
#------------------------------------------

if [ $single -eq 1 ] ; then
CHID='Spectra_085'

Workdir="/home/fdsout/Spectra/$CHID"
outdir="/home/wangbing/fdscov/$CHID/" # Must end with '/'
if ! test -d $outdir  ; then
  mkdir -p $outdir
fi

sel_plate=3
provide_offset='y'

xoffset=-60
yoffset=-60
zoffset=0
cd $Workdir

for sel_quantity in {1..5}
do
# Execute program

if [ $batchrun -eq 1 ] ; then
  $binary << ieof
$CHID
$sel_quantity
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset $Dspace
$outdir
ieof
else
  $binary
fi
done
fi

#------------------------------------------
# Script batch process
#------------------------------------------
if [ $single -eq 0 ] ; then
resultdir='/home/wangbing/FDS/Spectra'
allfolders=`ls $resultdir`

outdir='/home/wangbing/fdsout/testo'
if ! test -d $outdir  ; then
  mkdir -p $outdir
fi

for folder in $allfoders
do
  workingdir=$resultdir/$folder
  chid=$folder
  sel_plate=3
  provide_offset='y'
   xoffset=-60
   yoffset=-60
   zoffset=0
  cd $workingdir
  for sel_quantity in {1..5}
  do
    echo "quantity is $sel_quantity"

$binary << ieof
$chid
$sel_quantity
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset
$outdir
ieof

done
done
fi

