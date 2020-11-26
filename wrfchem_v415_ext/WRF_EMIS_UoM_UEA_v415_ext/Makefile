# Makefile for use on generic system instead of make_rossby.sh
# T. Pugh
# 25/08/10

.KEEP_STATE:

.f90.o:
	$(FC) $(FCFLAGS) -c $<

PROGRAM = main

SOURCES = readnc_2d_2.f90 BNGtolatlonconv.f90 voc_factors.f90 module.f90 module_lanuv.f90 module_naei.f90 module_ineris.f90 module_emep.f90 module_retro.f90 module_tno.f90 init.f90 strconv.f90 convert_lanuv.f90 mapcf.f90 misc.f90 regrid.f90 get_values.f90 scale_time.f90 write_files.f90 write_wrfchemi.f90 main.f90

OBJ = $(SOURCES:f90=o)

# The compiler
FC = gfortran
#FC = ifort
# Flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -ffree-form -fconvert=big-endian -ffree-line-length-none -O3
FCGLAGS= -ffree-form -Wall -fbounds-check -ffree-line-length-none -byteswapio -O3 -g
#FCFLAGS = -static -check bounds -free -convert big_endian -O3
# Flags for all
#FCFLAGS += -I/home/pught/netcdf-4.1.1-gfortran/build/include -L/home/pught/netcdf-4.1.1-gfortran/build/lib
FCFLAGS += -I/usr/include -L/usr/lib
#FCFLAGS += -I/sw/lib/netcdf-gfortran/include -L/sw/lib/netcdf-gfortran/lib
# Libraries needed for linking
#LDFLAGS = -L/home/pught/netcdf-4.1.1-gfortran/build/lib -lnetcdf
LDFLAGS = -L/usr/lib -lnetcdff -lnetcdf
#LDFLAGS = -L/sw/lib/netcdf-gfortran/lib -lnetcdf

# Redefine the implicit rule to put the object in the right place.
$(PROGRAM):
	$(FC) -o $(PROGRAM) $(FCFLAGS) $(SOURCES) $(LDFLAGS)

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD $(PROGRAM)

veryclean: clean
	rm -f *~ 
