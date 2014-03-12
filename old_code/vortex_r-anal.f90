      !-----------------------------------------------------------------!
      ! Module global                                                   !
      !-----------------------------------------------------------------!
      Module global
        implicit none
        character (len=3) :: label
        complex, allocatable :: wave_f_spc(:,:,:), wave_f_mom(:,:,:)
        complex, allocatable :: wave_f_flt(:,:,:), wave_f_fltspc(:,:,:)
        complex, allocatable :: wave_f_fltresc(:,:,:)
        real, allocatable :: re_y1(:,:), im_y1(:,:)
        real, allocatable :: re_y2(:,:), im_y2(:,:)
        real, allocatable :: trig_x(:), trig_y(:)
        real, allocatable :: modprobe(:,:)
        integer, parameter :: SP2 = KIND(1.0)
        integer, parameter :: dp2=kind(1.0d0)
        integer :: Nx, Ny
        integer :: Lz
        integer :: kount
        integer :: kx_max(1)
        real :: Lx, Ly, ax, ay
        real :: delta
        real :: kappa_C, kappa_X, sigma_p, sigma_pb, sigma_t, t_init, k_p, k_pb
        real :: f_p, f_pb, assymx,assymy, omega_p, omega_pb
        real :: norm, tot_h
        real :: eps_r
        real :: mom_cent, mom_cut
        real :: dxsav
      end Module global

      !-----------------------------------------------------------------!
      ! Module nrtype                                                   !
      !-----------------------------------------------------------------!
      Module nrtype
        !Symbolic names or kind types o 4-,2-,and 1-byte integers:
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
        INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
        !Symbolic names or kind types o single-and double-precision reals:
        INTEGER, PARAMETER :: SP = KIND(1.0)
        INTEGER, PARAMETER :: DP = KIND(1.0D0)
        !Symbolic names or kind types o single-and double-precision complex:
        INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
        INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
        !Symbolic name for kind type o de ault logical:
        INTEGER, PARAMETER :: LGT = KIND(.true.)
        !Frequently used mathematical constants (with precision to spare):
        REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
      end Module nrtype

      !-----------------------------------------------------------------!
      ! Main Program                                                    !
      !                                                                 !
      !-----------------------------------------------------------------!

! in order to run the program
! ifc -r8 -I/usr/local/shared/nag/nagf90_mod -L/usr/local/shared/nag/nagf90_mod vortex_r-anal.f90 -o analys.o -lnag -lnagfl90 -Vaxlib

      Program dynamika
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_trig
        USE global
        IMPLICIT NONE

        call read 

        allocate(trig_x(2*Nx), trig_y(2*Ny))
        allocate(re_y1(Nx,Ny), im_y1(Nx,Ny))
        allocate(re_y2(Nx,Ny), im_y2(Nx,Ny))
        allocate(wave_f_spc(Nx,Ny,2), wave_f_mom(Nx,Ny,2))
        allocate(wave_f_flt(Nx,Ny,2), wave_f_fltspc(Nx,Ny,2))
        allocate(wave_f_fltresc(Nx,Ny,2))
        allocate(modprobe(Nx,Ny))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)

        call analysis

        deallocate(trig_x, trig_y)
        deallocate(re_y1, im_y1)
        deallocate(re_y2, im_y2)
        deallocate(wave_f_spc, wave_f_mom)
        deallocate(wave_f_flt, wave_f_fltspc)
        deallocate(wave_f_fltresc)
        deallocate(modprobe)

      end Program dynamika

      !-----------------------------------------------------------------!
      ! Subroutine read                                                 !
      !-----------------------------------------------------------------!
      Subroutine read
        USE global
        IMPLICIT NONE

        OPEN(unit=22, file='INPUT', status='old')
        read(22,fmt="(10x,F10.6)") kappa_C  !decay rate for photons
        read(22,fmt="(10x,F10.6)") kappa_X  !decay rate for excitons
        read(22,fmt="(10x,F10.6)") delta    ! detuning
        read(22,fmt="(10x,I2)")    Lz       ! Lz of the vortex, Lz=1, 2, ...
        read(22,fmt="(10x,F10.6)") k_p      !pump  angle
        read(22,fmt="(10x,F10.6)") k_pb     !probe angle
        read(22,fmt="(10x,F10.6)") omega_p  !pump frequency
        read(22,fmt="(10x,F10.6)") omega_pb !probe frequency
        read(22,fmt="(10x,F10.6)") sigma_p  !pump  spatial
        read(22,fmt="(10x,F10.6)") sigma_pb !probe spatial
        read(22,fmt="(10x,F10.6)") sigma_t  !pulse 
        read(22,fmt="(10x,F10.6)") t_init   !initial time of the pulse
        read(22,fmt="(10x,F10.6)") f_p      ! pump strenght
        read(22,fmt="(10x,F10.6)") f_pb     ! probe strenght
        read(22,fmt="(10x,F10.6)") assymx   ! x position of probe
        read(22,fmt="(10x,F10.6)") assymy   ! y position of probe
        read(22,fmt="(10x,F10.6)") mom_cent ! center of filtering in momentum
        read(22,fmt="(10x,F10.6)") mom_cut  ! filtering in momentum
        read(22,fmt="(10x,F10.6)") Lx       ! x size of box (L)
        read(22,fmt="(10x,F10.6)") Ly       ! y size of box (L)
        read(22,fmt="(10x,I5)")    Nx       ! Number of points in x direction
        read(22,fmt="(10x,I5)")    Ny       ! Number of points in y direction
        read(22,fmt="(10x,F10.6)") tot_h    ! total time
        READ(22,FMT="(10x,F12.10)") dxsav
        READ(22,FMT="(10x,F12.10)") eps_r
        CLOSE(22)

        ax=2.0*Lx/Nx
        ay=2.0*Ly/Ny
        norm=ax*ay
        
      end Subroutine read
 
      !-----------------------------------------------------------------!
      ! Subroutine analysis                                             !
      !-----------------------------------------------------------------!
      Subroutine analysis
        use nag_fft, only: nag_fft_2d, nag_fft_trig
        USE global
        USE nrtype
        implicit none
        
        integer :: ix, iy, kcutx, kcuty
        integer :: kx, ky
        integer :: kount_st, kount_end
        real :: sx,sy
        real :: mom_x, mom_y
        real :: mom_x_maxs

        open(unit=29, file="max_signal.dat", status='replace')
        write(29,*) '#', 'kount', 'mom_x_maxs', 'maxval'
        
        call read

        !call read_probe_print
       
        write(*,*) "evaluate currents"
        write(*,*) "starting from frame number"
        read(*,*) kount_st

        write(*,*) "up to frame number"
        read(*,*) kount_end

        do kount=kount_st, kount_end
        !do kount=100, 200!int(tot_h/dxsav)
        !kount=25 

           call read_files
           !fft to momentum space
           wave_f_mom(:,:,1)=&
                nag_fft_2d(wave_f_spc(:,:,1),trig_m=trig_x,trig_n=trig_y)
           wave_f_mom(:,:,2)=&
                nag_fft_2d(wave_f_spc(:,:,2),trig_m=trig_x,trig_n=trig_y)

           call write_momentum
           call mom_filter

           !evaluate the signal peak in momentum
           kx_max=maxloc( abs(wave_f_flt(:,1,1)) )

           !region(2)
           if ( kx_max(1).ge.Nx/2+2 ) mom_x_maxs=pi*(kx_max(1)-1-Nx)/Lx
           !region(1)
           if ( kx_max(1).le.Nx/2+1 ) mom_x_maxs=pi*(kx_max(1)-1)/Lx

           write(29,*) kount, mom_x_maxs, maxval( abs(wave_f_flt(:,1,1)) )
           
           !one can do it if one is sure a signal does exists
           !if(Lz.eq.0) then
           !   mom_cent=mom_x_maxs
           !   call mom_filter
           !end if
           
           !write the wave function in real space and its phase
           !after momentum filtering
           open(unit=26, file="opo_ph_spcfilter"//trim(adjustl(label))//".dat", status='replace')
           write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')

           open(unit=28, file="opo_phasespc"//trim(adjustl(label))//".dat", status='replace') 
           write(28, fmt=' ("#", 1x, "x", 12x, "y", 12x, "theta") ')
           
           wave_f_fltspc(:,:,1)=&
                nag_fft_2d(wave_f_flt(:,:,1),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
           
           do iy=1, Ny
              sy=-Ly+(iy-1)*ay
              do ix=1, Nx
                 sx=-Lx+(ix-1)*ax
                 write(26,*) sx, sy, &
                      abs(wave_f_fltspc(ix,iy,1))*abs(wave_f_fltspc(ix,iy,1))
                 
                 write(28,*) sx, sy, real(pi+(0.0,-1.0)*log( -wave_f_fltspc(ix,iy,1)/( abs(wave_f_fltspc(ix,iy,1)) ) ))
              end do
              write(26,*)
              write(28,*)
           end do
           close(26)
           close(28)

!!$           !write the wave function cut y=const in real space
!!$           open(unit=24, file="y0cut_vortex"//trim(adjustl(label))//"_sq.dat", status='replace')
!!$           iy=int((0.0+Ly)/ay)+1
!!$           sy=-Ly+(iy-1)*ay
!!$           do ix=1, Nx
!!$              sx=-Lx+(ix-1)*ax
!!$              write(24, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ')&
!!$                   sx, sy, real(wave_f_fltspc(ix,iy,1))**2+aimag(wave_f_fltspc(ix,iy,1))**2
!!$           end do
!!$           close(24)

           !n.b. we also rescale down to zero mom_x_maxs
           !so that to subtract the oscillating phase
           call rescale_zero

           !call write_momflt_resc
           
!!$           !write the wave function in real space and its phase
!!$           !after momentum filtering and rescaling to zero momentum
!!$           open(unit=26, file="opo_ph_rescspcfilter"//trim(adjustl(label))//".dat", status='replace')
!!$           write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')
!!$
!!$           open(unit=28, file="opo_rescphasespc"//trim(adjustl(label))//".dat", status='replace') 
!!$           write(28, fmt=' ("#", 1x, "x", 12x, "y", 12x, "theta") ')
!!$           
!!$           wave_f_fltspc(:,:,1)=&
!!$                nag_fft_2d(wave_f_fltresc(:,:,1),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
!!$           
!!$           do iy=1, Ny
!!$              sy=-Ly+(iy-1)*ay
!!$              do ix=1, Nx
!!$                 sx=-Lx+(ix-1)*ax
!!$                 write(26,*) sx, sy, &
!!$                      abs(wave_f_fltspc(ix,iy,1))*abs(wave_f_fltspc(ix,iy,1))
!!$
!!$                 write(28,*) sx, sy, real(pi+(0.0,-1.0)*log( -wave_f_fltspc(ix,iy,1)/( abs(wave_f_fltspc(ix,iy,1)) ) ))
!!$              end do
!!$              write(26,*)
!!$              write(28,*)
!!$           end do
!!$           close(26)
!!$           close(28)
!!$
!!$           
           !open(unit=27, file="opo_ex_rescspcfilter"//trim(adjustl(label))//".dat", status='replace')
           !write(27, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(2)|^2") ')

        end do
        close(29)

      end Subroutine analysis

      !-----------------------------------------------------------------!
      ! Subroutine write_momflt_resc                                    !
      !-----------------------------------------------------------------!
      Subroutine write_momflt_resc
        USE global
        USE nrtype
        implicit none

        integer :: kx, ky
        real :: mom_x, mom_y

        !write the (photon) wave function in momentum space
        open(unit=24, file="opo_momflt_ph"//trim(adjustl(label))//".dat", status='replace')
        write(24, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')
       
        !write the wave function in momentum space
        !(n.b. one has to rearrange)
                   
        do ky=Ny/2+2, Ny
           mom_y=pi*(ky-1-Ny)/Ly
           !region(4)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_fltresc(kx,ky,1))*abs(wave_f_fltresc(kx,ky,1))
           end do
           !region(3)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_fltresc(kx,ky,1))*abs(wave_f_fltresc(kx,ky,1))
           end do
           write(24,*)
        end do
           
        do ky=1, Ny/2+1
           mom_y=pi*(ky-1)/Ly
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_fltresc(kx,ky,1))*abs(wave_f_fltresc(kx,ky,1))
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_fltresc(kx,ky,1))*abs(wave_f_fltresc(kx,ky,1))
           end do
           write(24,*)
        end do
        close(24)
        
      end Subroutine write_momflt_resc

      !-----------------------------------------------------------------!
      ! Subroutine rescale_zero                                         !
      !-----------------------------------------------------------------!
      Subroutine rescale_zero
        USE global
        USE nrtype
        implicit none

        integer :: kx, ky
        real :: mom_x, mom_y

        do ky=1, Ny
           do kx=1, Nx
              if( kx-1+kx_max(1).le.Nx ) then
                 wave_f_fltresc(kx,ky,1) = wave_f_flt(kx-1+kx_max(1),ky,1)
              else if ( kx-1+kx_max(1).gt.Nx ) then
                 wave_f_fltresc(kx,ky,1) = wave_f_flt(kx-1+kx_max(1)-Nx,ky,1)
              end if
           end do
        end do
        
      end Subroutine rescale_zero

      !-----------------------------------------------------------------!
      ! Subroutine mom_filter                                           !
      !-----------------------------------------------------------------!
      Subroutine mom_filter
        USE global
        USE nrtype
        implicit none

        integer :: kx, ky
        real :: mom_x, mom_y

        !filtering in momentum
        !define mom_cent 
        !and mom_cut in the INPUT file
        !n.b. 0.818 corresponds to 7 degrees
        
        wave_f_flt=(0.0d0,0.0d0)
           
        do ky=Ny/2+2, Ny
           mom_y=pi*(ky-1-Ny)/Ly
           !region(4)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              if( sqrt((mom_x-mom_cent)*(mom_x-mom_cent)+mom_y*mom_y) .le. mom_cut ) then
                 wave_f_flt(kx,ky,1) = wave_f_mom(kx,ky,1)
              end if
           end do
           !region(3)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              if( sqrt((mom_x-mom_cent)*(mom_x-mom_cent)+mom_y*mom_y) .le. mom_cut ) then
                 wave_f_flt(kx,ky,1) = wave_f_mom(kx,ky,1)
              end if
           end do
        end do
           
        do ky=1, Ny/2+1
           mom_y=pi*(ky-1)/Ly
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              if( sqrt((mom_x-mom_cent)*(mom_x-mom_cent)+mom_y*mom_y) .le. mom_cut ) then
                 wave_f_flt(kx,ky,1) = wave_f_mom(kx,ky,1)
              end if
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              if( sqrt((mom_x-mom_cent)*(mom_x-mom_cent)+mom_y*mom_y) .le. mom_cut ) then
                 wave_f_flt(kx,ky,1) = wave_f_mom(kx,ky,1)
              end if
           end do
        end do
        
      end Subroutine mom_filter

      !-----------------------------------------------------------------!
      ! Subroutine write_momntum                                        !
      !-----------------------------------------------------------------!
      Subroutine write_momentum
        USE global
        USE nrtype
        implicit none

        integer :: kx, ky
        real :: mom_x, mom_y

        !write the (photon) wave function in momentum space
        open(unit=24, file="opo_mom_ph"//trim(adjustl(label))//".dat", status='replace')
        write(24, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')
       
        !write the (photon) wave function in momentum space
        open(unit=25, file="cutky0opo_mom_ph"//trim(adjustl(label))//".dat", status='replace')
        write(25, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')

        !write the wave function in momentum space
        !(n.b. one has to rearrange)
                   
        do ky=Ny/2+2, Ny
           mom_y=pi*(ky-1-Ny)/Ly
           !region(4)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
           end do
           !region(3)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
           end do
           write(24,*)
        end do
           
        do ky=1, Ny/2+1
           mom_y=pi*(ky-1)/Ly
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(24,*) mom_x, mom_y,&
                   abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
           end do
           write(24,*)
        end do
        close(24)

        ky=1
        mom_y=pi*(ky-1)/Ly
        !region(2)
        do kx=Nx/2+2, Nx
           mom_x=pi*(kx-1-Nx)/Lx
           write(25,*) mom_x, mom_y,&
                abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
        end do
        !region(1)
        do kx=1, Nx/2+1
           mom_x=pi*(kx-1)/Lx
           write(25,*) mom_x, mom_y,&
                abs(wave_f_mom(kx,ky,1))*abs(wave_f_mom(kx,ky,1))
        end do
        close(25)
        
        write(*,*) 'max_val_sqrt_mom_ph='
        write(*,*) maxval( sqrt( abs(wave_f_mom(:,:,1))*abs(wave_f_mom(:,:,1)) ) )
        
        write(*,*) 'max_val_sqrt_mom_ex='
        write(*,*) maxval( sqrt( abs(wave_f_mom(:,:,2))*abs(wave_f_mom(:,:,2)) ) )

        !the exciton part
!!$        open(unit=24, file="opo_mom_ex"//trim(adjustl(label))//".dat", status='replace')
!!$        write(24, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(2)|^2") ') 
!!$        
!!$        do ky=Ny/2+2, Ny
!!$           mom_y=pi*(ky-1-Ny)/Ly
!!$           !region(4)
!!$           do kx=Nx/2+2, Nx
!!$              mom_x=pi*(kx-1-Nx)/Lx
!!$              write(24,*) mom_x, mom_y,&
!!$                   abs(wave_f_mom(kx,ky,2))*abs(wave_f_mom(kx,ky,2))
!!$           end do
!!$           !region(3)
!!$           do kx=1, Nx/2+1
!!$              mom_x=pi*(kx-1)/Lx
!!$              write(24,*) mom_x, mom_y,&
!!$                   abs(wave_f_mom(kx,ky,2))*abs(wave_f_mom(kx,ky,2))
!!$           end do
!!$           write(24,*)
!!$        end do
!!$           
!!$        do ky=1, Ny/2+1
!!$           mom_y=pi*(ky-1)/Ly
!!$           !region(2)
!!$           do kx=Nx/2+2, Nx
!!$              mom_x=pi*(kx-1-Nx)/Lx
!!$              write(24,*) mom_x, mom_y,&
!!$                   abs(wave_f_mom(kx,ky,2))*abs(wave_f_mom(kx,ky,2))
!!$           end do
!!$           !region(1)
!!$           do kx=1, Nx/2+1
!!$              mom_x=pi*(kx-1)/Lx
!!$              write(24,*) mom_x, mom_y,&
!!$                   abs(wave_f_mom(kx,ky,2))*abs(wave_f_mom(kx,ky,2))
!!$           end do
!!$           write(24,*)
!!$        end do
!!$        close(24)
           

      end Subroutine write_momentum

      !-----------------------------------------------------------------!
      ! Subroutine read_files                                           !
      !-----------------------------------------------------------------!
      Subroutine read_files
        USE global
        USE nrtype
        implicit none

        integer :: ix, iy
        real :: sx,sy

        write(label,FMT="(i3)") kount
        !read file with complex photon wf in space
        open(unit=22, file="phcplx-opo_spc"//trim(adjustl(label))//".dat", status='old')
        read(22, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(1))", 1x, "aimag(psi(1))") ')
        open(unit=23, file="excplx-opo_spc"//trim(adjustl(label))//".dat", status='old')
        read(23, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(2))", 1x, "aimag(psi(2))") ')
        do iy=1, Ny
           sy=-Ly+(iy-1)*ay
           do ix=1, Nx
              sx=-Lx+(ix-1)*ax
              read(22, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, re_y1(ix,iy), im_y1(ix,iy)
              read(23, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, re_y2(ix,iy), im_y2(ix,iy)
              
              wave_f_spc(ix,iy,1)=(1.0d0,0.0d0)*re_y1(ix,iy)+(0.0d0,1.0d0)*im_y1(ix,iy)
              wave_f_spc(ix,iy,2)=(1.0d0,0.0d0)*re_y2(ix,iy)+(0.0d0,1.0d0)*im_y2(ix,iy)
           end do
           read(22,*)
           read(23,*)
        end do
        close(22)
        close(23)

        

      end Subroutine read_files

      !-----------------------------------------------------------------!
      ! Subroutine read_probe_print                                     !
      !-----------------------------------------------------------------!
      Subroutine read_probe_print
        USE global
        USE nrtype
        implicit none

        integer :: ix, iy
        real :: sx,sy

        open(unit=23, file='probe.dat', status='old')

        do iy=1, Ny
           sy=-Ly+(iy-1)*ay
           do ix=1, Nx
              sx=-Lx+(ix-1)*ax
              read(23,*) sx, sy, modprobe(ix,iy)
           end do
           read(23,*)
        end do
        close(23)

        open(unit=23, file='y0cut_probe.dat', status='replace')

        !iy=Ny/2+1
        iy=int((-6.0+Ly)/ay)+1
        sy=-Ly+(iy-1)*ay
        do ix=1, Nx
           sx=-Lx+(ix-1)*ax
           write(23,*) sx, sy, modprobe(ix,iy)
        end do
        close(23)

      end Subroutine read_probe_print
