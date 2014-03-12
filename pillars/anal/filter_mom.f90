      Program filter_mom
        use nag_fft, only: nag_fft_2d, nag_fft_trig
        use global
        USE subroutines
        implicit none

        integer :: ix, iy
        integer :: kx, ky
        integer :: tnuok
        real :: sx,sy
        real :: mom_x, mom_y
        real :: mom_x_maxs

        call read

        allocate(trig_x(2*Nx), trig_y(2*Ny))
        allocate(pdb(Nx,Ny,2), wave_f_mom(Nx,Ny,2))
        allocate(wave_f_flt(Nx,Ny,2), wave_f_fltspc(Nx,Ny,2))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)

        open(unit=29, file="max_signal.dat", status='replace')
        write(29,*) '#', 'kount', 'mom_x_maxs', 'maxval'

        do tnuok=kount_st, kount_end
           write(label,FMT="(i3)") tnuok
           call init_pdb
           wave_f_mom(:,:,1)=&
                nag_fft_2d(pdb(:,:,1),trig_m=trig_x,trig_n=trig_y)
           wave_f_mom(:,:,2)=&
                nag_fft_2d(pdb(:,:,2),trig_m=trig_x,trig_n=trig_y)
           call write_momentum
           call mom_filter

           kx_max=maxloc( abs(wave_f_flt(:,1,1)) )
           if ( kx_max(1).ge.Nx/2+2 ) mom_x_maxs=pi*(kx_max(1)-1-Nx)/Lx
           if ( kx_max(1).le.Nx/2+1 ) mom_x_maxs=pi*(kx_max(1)-1)/Lx
           write(29,*) tnuok, mom_x_maxs, maxval( abs(wave_f_flt(:,1,1)) )

           open(unit=26, file="opo_ph_spcfilter"//trim(adjustl(label))//".dat", status='replace')
           write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')
           wave_f_fltspc(:,:,1)=&
                nag_fft_2d(wave_f_flt(:,:,1),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
           do iy=1, Ny
              sy=-Ly+(iy-1)*ay
              do ix=1, Nx
                 sx=-Lx+(ix-1)*ax
                 write(26,*) sx, sy, &
                      abs(wave_f_fltspc(ix,iy,1))*abs(wave_f_fltspc(ix,iy,1))
              end do
              write(26,*)
           end do
           close(26)
        end do
        close(29)

        deallocate(trig_x, trig_y)
        deallocate(pdb, wave_f_mom)
        deallocate(wave_f_flt, wave_f_fltspc)

      end Program filter_mom
