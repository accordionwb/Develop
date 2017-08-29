#!/bin/bash
# This script convert FDS output .sf file to matlab .mat file based on the SLCF configuration
# 
# One executables are involved in this script. 
# 1. slcf2mat
#        Reads FDS .sf file based on Quantity type and SLCF surface
#        Output  matlab .mat format file
# 
# Usage:
#  1.Define case coverage:
#
#       single=1 -> execute fds2slcf on specific run-case directory.
#       single=0 -> execute fds2slcf on whole run-case sourcediritory covering all cases
#
#  2.Define I/O directory
#
#  3.Define CHID if using single mode
#
#  4.Select SLCF plate
#     1 -> Y-Z plate
#     2 -> X-Z plate
#     3 -> X-Y plate
#
#  5.Give x,y,z offset to locate SW corner of the mesh. Dspace is important
#     For spectra simulatio where Dspace=1:
#       provide_offset='y'  |  'n'
#       xoffset=-60
#       yoffset=-60
#       zoff-set=0
#
#  6.Provide quantity selections
#     e.g.: sel_quantity={1..5}  (integer number)
#

single=1  # 1=single run | 0 = batch run

if test -z $1 ; then
CHID="Facility_039"
else 
    CHID=$1
fi


sourcedir="/disk/fdsrun/Facility"
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
if ! test -d $destdir ; then
   mkdir -p $destdir
fi
#
# Check Dynamic Library PATH
#
if [ ! `echo $LD_LIBRARY_PATH | grep MATLAB` ] ; then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/MATLAB/R2017a/bin/glnxa64
fi
#
# Define execute 
exe=slcf2mat
#
# ===================================================================
#
#------------------------------------------
# Script single process
#------------------------------------------

if [ $single -eq 1 ] ; then

   Readdir="$sourcedir/$CHID"
matfile=$destdir/$CHID.mat

   cd $Readdir

   for each_quan in ${sel_quantity[@]}
   do
      echo "all quan is: ${sel_quantity[@]}"

$exe << ieof
$CHID
$each_quan
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset $Dspace
$matfile
ieof
done

fi

#mail -s "$CHID finished" wangbing_2007@tsinghua.org.cn
#.

#------------------------------------------
# Script batch process
#------------------------------------------

if [ $single -eq 0 ] ; then  # Batch start

   allcases=`ls $sourcedir`
   for chid in $allcases
   do                    # Do loop 01: CASE
      readingdir=$sourcedir/$chid
      matfile=$destdir/$chid.mat

      cd $readingdir

      for each_quan in ${sel_quantity[@]}
      do                 # Do loop 02: Quantity

         echo "quantity is $each_quan"

$exe << ieof
$chid
$each_quan
$sel_plate
$provide_offset
$xoffset  $yoffset  $zoffset $Dspace
$matfile
ieof
      done
echo "======= Generate $cases.mat ==============="
done
fi
