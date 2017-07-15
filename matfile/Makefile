# This program build fortran source code to generate MATLAB '.mat' file 
# Functions defined in matlab should be included
# -------------------------------------------------------
# System PATH variables

MATLAB=/opt/MATLAB/R2017a
includedir=$(MATLAB)/extern/include

srcdir = src
target = matdemo1
bindir =bin
OBJCTS = $(srcdir)/$(target).o

FC=ifort
FFLAGS=  -O2 -m64 -fpp
INCLDS=-I$(includedir)
LIBS=-lmat -lmx -L$(MATLAB)/bin/glnxa64 
TARGS=$(bindir)/$(target)


# ------------- Linux 64 ----------------
#  64 bit linux



all : $(OBJCTS)
	$(FC) -o $(TARGS) $(FFLAGS) $(OBJCTS) $(INCLDS) $(LIBS) 


#-------------- Clean Target to remove OBJCTSect and Module files -----------

.PHONY : clean
clean:
	rm $(OBJCTS)

 
#*** General Purpose Rules ***

.SUFFIXES: .F90 .c .o .f90 .cpp .f .F

.F90.o:
	$(FC) -c -o $@ $(FFLAGS) $(INCLDS) $<
.f90.o:
	$(FC) -c -o $@ $(FFLAGS) $(INCLDS) $<
.f.o:
	$(FC) -c -o $@ $(FFLAGS) $(INCLDS) $<
.F.o:
	$(FC) -c -o $@ $(FFLAGS) $(INCLDS) $<

 
