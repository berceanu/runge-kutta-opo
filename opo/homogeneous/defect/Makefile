# Start of the makefile
# Defining variables
f90comp = gfortran
switch = 

# List of includes for the compiler 
NAG_PATH=/opt/NAG/fnl6i04dcl
NAG_MODULES=$(NAG_PATH)/nagfl90_modules/
NAG_LIBS=$(NAG_PATH)/lib
includes = -I$(NAG_MODULES)  
libs = $(NAG_LIBS)/libnagfl90_nag.a

# Makefile
filter_en.out: filter_en.o subroutines.o global.o ode_path.o nrtype.o
	$(f90comp) -o filter_en.out $(switch) filter_en.o subroutines.o global.o ode_path.o nrtype.o $(libs)
filter_en.o: subroutines.mod global.mod filter_en.f90
	$(f90comp) -c $(switch) filter_en.f90 $(includes)

filter_mom.out: filter_mom.o subroutines.o global.o ode_path.o nrtype.o
	$(f90comp) -o filter_mom.out $(switch) filter_mom.o subroutines.o global.o ode_path.o nrtype.o $(libs)
filter_mom.o: subroutines.mod global.mod filter_mom.f90
	$(f90comp) -c $(switch) filter_mom.f90 $(includes)

RK_OPO_sp.out: RK_OPO_sp.o subroutines.o rk_adaptive.o global.o ode_path.o nrtype.o
	$(f90comp) -o RK_OPO_sp.out $(switch) RK_OPO_sp.o subroutines.o rk_adaptive.o global.o ode_path.o nrtype.o $(libs)
RK_OPO_sp.o: subroutines.mod rk_adaptive.mod global.mod RK_OPO_sp.f90
	$(f90comp) -c $(switch) RK_OPO_sp.f90 $(includes)

RK_OPO_rk.out: RK_OPO_rk.o subroutines.o rk_adaptive.o global.o ode_path.o nrtype.o
	$(f90comp) -o RK_OPO_rk.out $(switch) RK_OPO_rk.o subroutines.o rk_adaptive.o global.o ode_path.o nrtype.o $(libs)
RK_OPO_rk.o: subroutines.mod rk_adaptive.mod global.mod RK_OPO_rk.f90
	$(f90comp) -c $(switch) RK_OPO_rk.f90 $(includes)
subroutines.mod: subroutines.o subroutines.f90
	$(f90comp) -c $(switch) subroutines.f90
subroutines.o: global.mod subroutines.f90
	$(f90comp) -c $(switch) subroutines.f90
rk_adaptive.mod: rk_adaptive.o rk_adaptive.f90
	$(f90comp) -c $(switch) rk_adaptive.f90 $(includes)
rk_adaptive.o: global.mod rk_adaptive.f90
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

# Cleaning everything
clean:
	rm *.o
	rm *.mod
	rm *.out
# End of the makefile
