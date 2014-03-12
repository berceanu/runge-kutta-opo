      Program RK_OPO_rk
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_trig
        USE global
        USE subroutines
        USE rk_adaptive
        IMPLICIT NONE

        REAL(SP)  :: x1_r,x2_r,h1_r,hmin_r    
        REAL, EXTERNAL ::  findfermpart,findcoopernum    
        COMPLEX, EXTERNAL :: findfermcond1, findfermcond2

        call read 

        allocate(trig_x(2*Nx), trig_y(2*Ny))
        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pot_c(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(pump(Nx,Ny))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)

		!initialize photon potential
        call init_pot_c
        
        !initialize wf
        if(run.eq.1) then
           pdb=(0.0,0.0)
        else if(run.eq.0) then
           write(label,FMT="(i3)") in_sswf_rk
           call init_pdb
        end if
        
        !use homogeneous pump
        call init_pump_homo
        
        call setg
  
        x1_r=0.0    
        x2_r=tot_h    
        h1_r=0.001    
        hmin_r=0.0    
        CALL odeint_rk(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)    
  
        deallocate(trig_x, trig_y)
        deallocate(pdb, kinetic, pot_c, pump_spatial)
        deallocate(pump)
      end Program RK_OPO_rk

