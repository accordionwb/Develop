# Makefile to compile lid-driven carvity flow
.SUFFIXES: .c .f90 .o
VPATH = ../../src/lid_cf
FC = pgf90
#FLGS = -convert big_endian -Bstatic
FLGS = -Wl, -Bstatic
LIBS =

TARGT = exec

OBJCTS=	 lid_cf.o\
	 result.o\
	 streaming.o\
	 collesion.o\
	 sfbound.o\
	 rhouv.o


$(TARGT) : $(OBJCTS)
	$(FC) -o $(TARGT) $(OBJCTS)

.f90.o : 
	$(FC) -c -o $@ $(FLGS) $<

.PHONY : clean
clean : 
	rm -f $(OBJCTS) $(TARGT) 



