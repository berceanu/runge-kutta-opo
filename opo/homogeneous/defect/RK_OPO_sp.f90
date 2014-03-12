      Program RK_OPO_sp
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_3d, nag_fft_trig
        USE global
        USE subroutines
        USE rk_adaptive
        IMPLICIT NONE

        REAL(SP)  :: x1_r,x2_r,h1_r,hmin_r    
        REAL, EXTERNAL ::  findfermpart,findcoopernum    
        COMPLEX, EXTERNAL :: findfermcond1, findfermcond2

        open(unit=file_time_0, file="times_0.dat", status='replace')

        call read 

        tot_h=t_stfft+(Nt+1)*(dxsav_sp)

        allocate(trig_x(2*Nx), trig_y(2*Ny), trig_t(2*Nt))
        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pot_c(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(pump(Nx,Ny))
        allocate(y_tot_0(Nx,Ny,Nt), int_sp(Nt))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)
        call nag_fft_trig(trig_t)

		!initialize variables
		!photon potential
		call init_pot_c

        pdb=(0.0,0.0)
        write(label,FMT="(i3)") in_sswf_sp
		!initialize wavefunction
        call init_pdb
        !use homogeneous pump
        call init_pump_homo

		!calculate kinetic energy
        call setg

		!GP time evolution
        x1_r=0.0    
        x2_r=tot_h    
        h1_r=0.001    
        hmin_r=0.0    
        CALL odeint_sp(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)    

        !FFT to energy and momentum space    
        y_tot_0(:,:,:)=&    
             nag_fft_3d(y_tot_0(:,:,:),trig_1=trig_x,trig_2=trig_y,trig_3=trig_t)                 
  
		!same time evolution to file for integrating later in energy window
  	    call export_evolution
  	    
		!calculate the energy spectrum
        call eval_spectr_0

        deallocate(trig_x, trig_y, trig_t)
        deallocate(pdb, kinetic, pot_c, pump_spatial)
        deallocate(pump)
        deallocate(y_tot_0, int_sp)

        close(file_time_0)

      end Program RK_OPO_sp
