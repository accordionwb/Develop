#!/bin/bash
# This script convert FDS output .sf file to matlab .mat file based on the SLCF configuration
# 
# Two executables are involved in this script. 
# 1. fds2slcf:
#        Reads FDS .sf file based on Quantity type and SLCF surface
#        Output contains:
#           .csv sample ASCII file
#           .bin binary data file 
# 2. fds2mat:
#        Read multiple .bin data file and combines them in one .mat file
#        Variable name comes from FDS output quantity
# 
# Usage:
#  1.Define run mode:
#       debugrun=1 -> execute fds2slcf in console and enter input parameters from .stdin.
#       debugrun=0 -> execute fds2slcf in background and input parameters are given from 
#                     pre-defined script variables
#     A debug run is recommended at the first time
#
#  2.Define case coverage:
#       single=1 -> execute fds2slcf on specific run-case directory.
#       single=0 -> execute fds2slcf on whole run-case sourcediritory covering all cases
#
#  3.Define I/O directory
#
#  4.Define CHID if using single mode
#
#  5.Select SLCF plate
#     1 -> Y-Z plate
#     2 -> X-Z plate
#     3 -> X-Y plate
#
#  6.Give x,y,z offset to locate SW corner of the mesh. Dspace is important
#     For spectra simulatio where Dspace=1:
#       provide_offset='y'  |  'n'
#       xoffset=-60
#       yoffset=-60
#       zoff-set=0
#
#  7.Provide quantity selections
#     e.g.: sel_quantity={1..5}  (integer number)
#

debugrun=0  # 1=debug run | 0 no debug
exe1=fds2slcf
exe2=fds2mat
single=1  # 1=single run | 0 = batch run

if test -z $1 ; then
CHID="Facility_025"
else
    CHID=$1
fi


sourcedir="/disk/fdsrun/Facility"
tmpdir="/disk/work/fdscov"
destdir="/disk/fdsmat/Facility"


sel_plate=3  # 
provide_offset='y'

Dspace=1.0
xoffset=-60.
yoffset=-60.
zoffset=0.
sel_quantity=(1 2 3 4 5)

# ========================================================================
#
# Script pre-check working directorys
#
if ! test -d $tmpdir ; then
   mkdir -p $tmpdir
fi
#
# Check Dynamic Library PATH
#
if [ ! `echo $LD_LIBRARY_PATH | grep MATLAB` ] ; then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/MATLAB/R2017a/bin/glnxa64
fi
#
# ===================================================================
#
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

   Readdir="$sourcedir/$CHID"
   bindir="$tmpdir/$CHID/" # Must end with '/'
   if ! test -d $bindir  ; then
      mkdir -p $bindir
   fi

   cd $Readdir

      if [ $debugrun -eq 0 ] ; then # Non debug mode
   for each_quan in ${sel_quantity[@]}
   do
      echo "all quan is: ${sel_quantity[@]}"

$exe1 << ieof
$CHID
$each_quan
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset $Dspace
$bindir
ieof
done

         allfiles=`ls $bindir/*.bin`  # /tmp/fdsmisc/$CHID/*******.bin
         outfilename=$destdir/$CHID.mat
         for file in $allfiles # infilename 
         do
            infilename=$file
$exe2 << ieof
$infilename
$outfilename
ieof
rm $infilename
         done  # loop each .bin
echo "======= Generate $CHID.mat ==============="
rm -rf $bindir
      else  # debug mod
         $exe1
         $exe2
      fi
fi

#------------------------------------------
# Script batch process
#------------------------------------------
if [ $debugrun -eq 1 ] ; then
   echo "Could not run batch process in debug mode"
   exit 0
fi

if [ $single -eq 0 ] ; then  # Batch start

   allcases=`ls $sourcedir`
   for cases in $allcases
   do                    # Do loop 01: CASE
      CHID=$cases
      bindir=$tmpdir/$cases/   # Must have '/' in the end
      if ! test -d $bindir ; then
         mkdir -p $bindir
      fi
      readingdir=$sourcedir/$cases
      cd $readingdir
      for each_quan in ${sel_quantity[@]}
      do                 # Do loop 02: Quantity

         echo "quantity is $each_quan"

$exe1 << ieof
$CHID
$each_quan
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset $Dspace
$bindir
ieof
      done    # End do loop 02: quantity

      allfiles=`ls $bindir/*.bin`
      outfilename=$destdir/$cases.mat

      for file in $allfiles
      do
         infilename=$file

$exe2 << ieof
$infilename
$outfilename
ieof
rm $infilename
      done
echo "======= Generate $cases.mat ==============="
rm -rf $bindir

   done
fi
