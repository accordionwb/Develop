# This program build fortran source code to generate MATLAB '.mat' file 
# Functions defined in matlab should be included
# -------------------------------------------------------
# System PATH variables

MATLAB=/opt/MATLAB/R2016a
includedir=$(MATLAB)/extern/include

srcdir = src
target = matdemo1
bindir =bin



# Definition of the OBJCTSect variables

OBJCTS = $(srcdir)/$(target).o

#*** General Purpose Rules ***

.SUFFIXES: .F90 .c .o .f90 .cpp .f

.F90.o:
	$(FC) -c -o $@ $(FFLAGS) $<
.f90.o:
	$(FC) -c -o $@ $(FFLAGS) $<
.f.o:
	$(FC) -c -o $@ $(FFLAGS) $<


# ------------- Linux 64 ----------------
#  64 bit linux


intel_linux_64 : FFLAGS= -I$(includedir) -O2 -m64 -fpp
intel_linux_64 : LIBS=-L$(MATLAB)/bin/glnxa64 -lmat -lmx
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = $(bindir)/$(target)

intel_linux_64 : $(OBJCTS)
	$(FC) -o $(exe) $(FFLAGS) $(OBJCTS) $(LIBS)

#-------------- Clean Target to remove OBJCTSect and Module files -----------

.PHONY : clean
clean:
	rm $(OBJCTS)

 
 
