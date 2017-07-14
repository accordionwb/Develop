SRCDir = src
target = fds2slcf

VPATH = src/
BIN =bin


# Definition of the OBJCTSect variables

OBJCTS = $(SRCDir)/$(target).o

#*** General Purpose Rules ***

no_target:
	@echo \******** You did not specify a make target \********
	@echo Use one of targets found in the Makefile

.SUFFIXES: .F90 .c .o .f90 .cpp

.F90.o:
	$(FC) -c $(FFLAGS) -o $(OBJCTS) $<
.f90.o:
	$(FC) -c $(FFLAGS) -o $(OBJCTS) $<


# ------------- gnu linux 64 ----------------
#  64 bit gcc linux

gnu_linux_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  
gnu_linux_64 : FC        = gfortran
gnu_linux_64 : exe       = $(BIN)/$(target)

gnu_linux_64 : $(OBJCTS)
	$(FC) -o $(exe) $(OBJCTS)

# ------------- gnu osx 64 ----------------
#  64 bit gcc linux

gnu_osx_64 : FFLAGS    = -O0 -m64 -ffree-form -frecord-marker=4  $(GITINFOGNU)
gnu_osx_64 : FC        = gfortran
gnu_osx_64 : exe       = $(BIN)/$(target)

gnu_osx_64 : $(OBJCTS)
	$(FC) -o $(exe) $(OBJCTS)



# ------------- Linux 64 ----------------
#  64 bit linux

intel_linux_64 : FFLAGS    = -O2 -m64 
intel_linux_64 : FC        = ifort
intel_linux_64 : exe       = $(BIN)/$(target)

intel_linux_64 : $(OBJCTS)
	$(FC) -o $(exe) $(OBJCTS)

# ------------- OSX 64 ----------------
#  64 bit OSX

intel_osx_64 : FFLAGS    = -O2 -m64 $(GITINFO)
intel_osx_64 : FC        = ifort
intel_osx_64 : exe       = $(BIN)/$(target)

intel_osx_64 : $(OBJCTS)
	$(FC) -o $(exe)  -mmacosx-version-min=10.4 $(OBJCTS)

#-------------- Clean Target to remove OBJCTSect and Module files -----------

.PHONY : clean
clean:
	rm $(OBJCTS)
