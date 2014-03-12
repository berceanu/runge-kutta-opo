      Program RK_OPO_sp
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_3d, nag_fft_trig
        USE global
        USE subroutines
        USE rk_adaptive
        IMPLICIT NONE

        REAL(SP)  :: x1_r,x2_r,h1_r,hmin_r
        REAL, EXTERNAL ::  findfermpart,findcoopernum
        COMPLEX, EXTERNAL :: findfermcond1, findfermcond2

        !to label the time of the en filt images
        time_kount=1
        n_save=1

        open(unit=file_time, file="times.dat", status='replace')
        write(file_time,*) n_save

        call read

        !allow for one more step before end
        tot_h=t_stfft+tot_save*Nt*dxsav_sp

        allocate(trig_x(2*Nx), trig_y(2*Ny), trig_t(2*Nt))
        allocate(trig_redt(2*(2*Nredt+1)))
        allocate(pot_c_t(Nx,Ny,tot_h))
        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(pump(Nx,Ny))
        allocate(y_tot(Nx,Ny,Nt))
        allocate(int_sp(Nt), int_sp_ok(Nt), int_sp_s(Nt), int_sp_p(Nt), int_sp_i(Nt))
        allocate(y_enfilt(Nx,Ny,2*Nredt+1))
        allocate(check(tot_save))

        !to call eval_spectr only once
        check=0

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)
        call nag_fft_trig(trig_t)
        call nag_fft_trig(trig_redt)

        !initialize variables
        !photon potential
        call init_pot_c_t

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
        call odeint_sp(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)

        deallocate(trig_x, trig_y, trig_t, trig_redt)
        deallocate(pdb, kinetic, pot_c_t, pump_spatial)
        deallocate(pump)
        deallocate(y_tot, y_enfilt)
        deallocate(int_sp, int_sp_ok, int_sp_p, int_sp_s, int_sp_i)
        deallocate(check)

        close(file_time)

      end Program RK_OPO_sp
