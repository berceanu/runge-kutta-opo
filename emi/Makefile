
TARGET = a.out

######################################################################
#  1 Dimensional Discrete Fast Fourier Transform from intel libraries          
#  
#  30 April 2009        ifort compiler         by Emiliano Cancellieri
######################################################################

## name of the fortran compiler
F90 = ifort   

## flags to use for fortran compilation
F90FLAGS = -r8 -fast -O3

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90FLAGS)

## LLIBS are libraries to use at link time
MKL_PATH=/opt/intel/mkl/10.1.2.024/lib/em64t
LLIBS = -L$(MKL_PATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm

#################################################################
#      objects files
#################################################################

OBJS =  RK_gaussian_probe.o		\
        mkl_dfti.o			\
	mod_modules.o			\
	mod_subroutines.o		\
	mod_rk_adaptive.o

#################################################################
#      rules 
#################################################################

.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(F90FLAGS) -c $<

$(TARGET) : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $(OBJS) $(LLIBS)

clean:
	/bin/rm -f $(OBJS)
	/bin/rm -f *.mod
	/bin/rm -f $(TARGET)

#################################################################
#      dependencies
#################################################################

RK_gaussian_probe.o :		RK_gaussian_probe.f90	\
				mod_modules.o		\
				mod_subroutines.o	\
				mkl_dfti.o		\
				mod_rk_adaptive.o

mod_rk_adaptive.o :		mod_rk_adaptive.f90	\
				mod_modules.o		\
				mod_subroutines.o

mod_modules.o :			mod_modules.f90		\
				mkl_dfti.o

mod_subroutines.o :		mod_subroutines.f90	\
				mod_modules.o		\
				mkl_dfti.o

mkl_dfti.o :			mkl_dfti.f90
