# Use the same compiler and compatible options as the LIBSTELL build.
#FC = gfortran
#FFLAGS ?= -O0 -g -fcheck=all -fbacktrace -Wimplicit-interface -Wall -Wextra
#CPPFLAGS ?=
#LDFLAGS ?=

# Override these paths on the make command line or in config.mk.
STELLOPT ?= $(STELLOPT_PATH)
MODDIR ?= $(STELLOPT)/LIBSTELL/Release
LIBSTELL ?= $(STELLOPT)/LIBSTELL/Release/libstell.a

# Supply dependencies used by YOUR LIBSTELL build, in static-link order.
# Examples: nf-config --flibs, LAPACK/BLAS, HDF5, FFTW, OpenMP, or MPI.
LDLIBS ?= 
-include $(STELLOPT)/make.inc

FC = $(MPI_COMPILE)
FFLAGS = $(FLAGS_R) -g -fcheck=all -fbacktrace -Wimplicit-interface -Wall -Wextra
FFLAGS += $(NETCDF_INC)
LDLIBS += $(LIBS) $(NETCDF_LIB)

.PHONY: all clean dsym check
all: compare_wout

compare_wout: cylindrical_curl.o compare_wout.o $(LIBSTELL)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ cylindrical_curl.o compare_wout.o $(LIBSTELL) $(LDLIBS)

cylindrical_curl.o: cylindrical_curl.f90
	$(FC) $(CPPFLAGS) $(FFLAGS) -I$(MODDIR) -c $< -o $@

compare_wout.o: compare_wout.f90 cylindrical_curl.o
	$(FC) $(CPPFLAGS) $(FFLAGS) -I$(MODDIR) -c $< -o $@

# macOS: run after linking, while the LIBSTELL archive remains available.
dsym: compare_wout
	dsymutil ./compare_wout

# Standalone analytic test; needs only stel_kinds.mod, not libstell.a.
check: cylindrical_curl.o tests/test_curl.f90
	$(FC) $(FFLAGS) -I$(MODDIR) -o test_curl tests/test_curl.f90 cylindrical_curl.o
	./test_curl

clean:
	rm -f compare_wout test_curl *.o cylindrical_curl.mod
	rm -rf compare_wout.dSYM
