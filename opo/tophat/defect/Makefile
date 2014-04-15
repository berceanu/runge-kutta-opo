# Start of the makefile
# Defining variables
f90comp = gfortran

# when compiling for debugging
switch = -g3 -Wall -Wextra -Wconversion -pedantic -fimplicit-none -fbounds-check -fbacktrace -ffree-line-length-none -ffpe-trap=zero,overflow,underflow -std=f2003

# when compiling for production
#switch = -O3 -ffree-line-length-none

# when compiling for profiling
#switch = -O2 -ffree-line-length-none -g -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls

# List of includes for the compiler
includes = -I/usr/include/

# for the serial version
#libs = -lfftw3 -lm

# for the multi-threaded version
libs = -lfftw3_threads -lfftw3 -lm -lpthread

# Makefile

subroutines.mod: subroutines.o subroutines.f90
	$(f90comp) -c $(switch) subroutines.f90
subroutines.o: global.mod subroutines.f90
	$(f90comp) -c $(switch) subroutines.f90
rk_adaptive.mod: rk_adaptive.o rk_adaptive.f90
	$(f90comp) -c $(switch) rk_adaptive.f90 $(includes)
rk_adaptive.o: global.mod FFTW3.mod rk_adaptive.f90
	$(f90comp) -c $(switch) rk_adaptive.f90 $(includes)
global.mod: global.o global.f90
	$(f90comp) -c $(switch) global.f90
global.o: ode_path.mod global.f90
	$(f90comp) -c $(switch) global.f90
ode_path.mod: ode_path.o ode_path.f90
	$(f90comp) -c $(switch) ode_path.f90
ode_path.o: nrtype.mod ode_path.f90
	$(f90comp) -c $(switch) ode_path.f90
nrtype.mod: nrtype.o nrtype.f90
	$(f90comp) -c $(switch) nrtype.f90
nrtype.o: nrtype.f90
	$(f90comp) -c $(switch) nrtype.f90
FFTW3.mod: FFTW3.o FFTW3.f90
	$(f90comp) -c $(switch) FFTW3.f90 $(includes)
FFTW3.o: FFTW3.f90
	$(f90comp) -c $(switch) FFTW3.f90 $(includes)

# Cleaning everything
clean:
	rm *.o
	rm *.mod
	rm *.out
# End of the makefile
