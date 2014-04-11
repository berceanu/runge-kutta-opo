      Program RK_OPO_sp
        USE global
        USE subroutines
        USE rk_adaptive
        USE FFTW3
        IMPLICIT NONE

        real(8) :: x1_r,x2_r,h1_r,hmin_r
        real(8), EXTERNAL ::  findfermpart,findcoopernum
        complex(8), EXTERNAL :: findfermcond1, findfermcond2

        ! fft stuff
        ! forward means real space to momentum space
        type(C_PTR) :: plan_forward
        complex(C_DOUBLE_COMPLEX), pointer :: in_forward(:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: out_forward(:,:,:)
        type(C_PTR) :: p, q
        ! array dimensions
        integer(C_INT) :: dimx, dimy, dimt

        dimx=size(y_tot_0,1)
        dimy=size(y_tot_0,2)
        dimt=size(y_tot_0,3)

        open(unit=file_time_0, file="times_0.dat", status='replace')
        call read
        tot_h=t_stfft+(Nt+1)*(dxsav_sp)

        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pot_c(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(pump(Nx,Ny))
        allocate(y_tot_0(Nx,Ny,Nt), int_sp(Nt))

        !initialize variables
        !photon potential
        call init_pot_c

        pdb=(0.0,0.0)
        write(label,FMT="(i3)") in_sswf_sp
        !initialize wavefunction
        call init_pdb
        !use top hat pump
        call init_pump_th

        !calculate kinetic energy
        call setg

        !GP time evolution
        x1_r=0.0
        x2_r=tot_h
        h1_r=0.001
        hmin_r=0.0
        CALL odeint_sp(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)

        ! allocating memory contiguously using C function
        p = fftw_alloc_complex(int(dimx*dimy*dimt, C_SIZE_T))
        q = fftw_alloc_complex(int(dimx*dimy*dimt, C_SIZE_T))

        ! here we use the usual fortran order
        call c_f_pointer(p, in_forward, [dimx,dimy,dimt])
        call c_f_pointer(q, out_forward, [dimx,dimy,dimt])

        ! prepare plans needed by fftw3
        ! here we must make sure we reverse the array dimensions for FFTW
        plan_forward = fftw_plan_dft_3d(dimt, dimy, dimx, in_forward, out_forward, FFTW_FORWARD, FFTW_MEASURE)

        !FFT to energy and momentum space

        in_forward = y_tot_0
        call fftw_execute_dft(plan_forward, in_forward, out_forward)
        y_tot_0 = out_forward

        !save time evolution to file for integrating later in energy window
        call export_evolution

        !calculate the energy spectrum
        call eval_spectr_0

        deallocate(pdb, kinetic, pot_c, pump_spatial)
        deallocate(pump)
        deallocate(y_tot_0, int_sp)

        close(file_time_0)

        ! avoiding any potential memory leaks
        call fftw_destroy_plan(plan_forward)
        call fftw_free(p)
        call fftw_free(q)

      end Program RK_OPO_sp
