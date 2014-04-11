      Program RK_OPO_rk
        USE global
        USE subroutines
        USE rk_adaptive
        IMPLICIT NONE

        REAL :: x1_r,x2_r,h1_r,hmin_r
        REAL, EXTERNAL ::  findfermpart,findcoopernum
        COMPLEX, EXTERNAL :: findfermcond1, findfermcond2

        call read

        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pot_c(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(pump(Nx,Ny))

	!initialize photon potential
        call init_pot_c

        !initialize wf
        if(run.eq.1) then
           pdb=(0.0,0.0)
        else if(run.eq.0) then
           write(label,FMT="(i3)") in_sswf_rk
           call init_pdb
        end if

        !use top hat pump
        call init_pump_th

        call setg

        x1_r=0.0
        x2_r=tot_h
        h1_r=0.001
        hmin_r=0.0
        CALL odeint_rk(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)

        deallocate(pdb, kinetic, pot_c, pump_spatial)
        deallocate(pump)
      end Program RK_OPO_rk

