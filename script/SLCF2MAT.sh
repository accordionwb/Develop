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



# ==================== Input arguement defination ========================
# 1st parameter: specify the case number in integer
if test -z $1 ; then
    echo "--------------------------------------------------------"
    echo " Usage:"
    echo "   SLCF2MAT [arg_1] [arg_2] [art_3]"
    echo "      -arg_1: Start case number"
    echo "      -arg_2: Rewrite flag (true/false, default)"
    echo "      -arg_3: End case number (optional)"
    echo "   Default prefix is 'Facility_'"
    echo "   Case number should be integer between 1 and 120"
    echo " "
    echo " Modify the script to change the input&output path"
    echo "--------------------------------------------------------"
    exit
else 
    startcase=$1
fi

# 2nd parameter: If re-write the existing file
if test -z $2 ; then
    rewrite=false
else 
    rewrite=true
fi

# 3rd parameter: Indicate the end case in integer
if test -z $3 ; then
    endcase=$1
else
    endcase=$3
    if $endcase < $startcase ; then
        echo "***Error, 1st arg < 3st arg, program exit"
        exit 
    fi
fi



#================= Parameter Settings =========================


sourcedir="/disk/fdsrun/Facility"
destdir="/disk/fdsmat/Facility"

sel_plate=3  # 
provide_offset='y'

Dspace=1.0
xoffset=-60.
yoffset=-60.
zoffset=0.
sel_quantity=(1 2 3 4 5)


#================== Script pre-check working directorys ================
if ! test -d $destdir ; then
   mkdir -p $destdir
fi

#================ Check Dynamic Library PATH ============================
if [ ! `echo $LD_LIBRARY_PATH | grep MATLAB` ] ; then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/MATLAB/R2017a/bin/glnxa64
fi

# Define execute 
exe=slcf2mat

# ================ Main Loop ================================
for ((i=$startcase; i<=$endcase; i++))
do

strcase=`printf "%03d" $i`
CHID="Facility_"$strcase

   Readdir="$sourcedir/$CHID"
matfile=$destdir/$CHID.mat
if $rewrite ; then
    echo "Re-write the existing mat file"
    if test -e $matfile ; then
    rm $matfile
fi
fi
echo "*********** Program Start on case $CHID ******************"
echo " "
sleep 5

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



# Loop output
echo "=================================================="
echo "="
echo "=            Generate $CHID.mat "
echo "="
echo "=================================================="

done


#echo . | mail -s "$CHID finished" wangbing_2007@tsinghua.org.cn


