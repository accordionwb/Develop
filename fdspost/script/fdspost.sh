#!/bin/sh
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
# Script batch process

