      !-----------------------------------------------------------------!
      ! Module global                                                   !
      !-----------------------------------------------------------------!
      Module global
        implicit none
        character (len=3) :: label
        complex, allocatable :: psi(:,:,:), pump_spatial(:,:)
        complex, allocatable :: y_tot_0(:,:,:), y_tot_1(:,:,:)
        complex, allocatable :: probe_spatial(:,:)
        complex, allocatable :: pdb(:,:,:)
        complex, allocatable :: pump(:,:), probe(:,:)
        real, allocatable :: trig_x(:), trig_y(:), trig_t(:)
        real, allocatable :: kinetic(:,:)
        real, allocatable :: int_sp(:)
        integer, parameter :: SP2 = KIND(1.0)
        integer, parameter :: dp2=kind(1.0d0)
        integer :: Nx, Ny
        integer :: Lz
        integer :: in_sswf, Nt
        integer :: file_time_0=29, file_time_1=30
        integer :: trigg
        real :: Lx, Ly, ax, ay
        real :: delta
        real :: kappa_C, kappa_X, sigma_p, sigma_pb, sigma_t, t_init, k_p, k_pb
        real :: f_p, f_pb, assymx,assymy, omega_p, omega_pb
        real :: norm, tot_h
        real :: eps_r
        real :: mom_cent, mom_cut
        real :: t_stfft, t_stfft2
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
      ! Module ode_path                                                 !
      !-----------------------------------------------------------------!
      Module ode_path 
!On output nok and nbad are the number of good and bad (but retried
!and fixed) steps taken. If save_steps is set to true in the calling program,
!then intermediate values are stored in xp and yp at in ervals greater
!than dxsav. kount is the total number of saved steps.
        USE nrtype
        INTEGER(I4B) :: nok,nbad,kount, kount_1
        !LOGICAL(LGT), SAVE :: save_steps=.false.
        LOGICAL(LGT), SAVE :: save_steps=.true.
        REAL(kind=4) :: dxsav
      end Module ode_path

      !-----------------------------------------------------------------!
      ! Main Program                                                    !
      !                                                                 !
      !-----------------------------------------------------------------!

! in order to run the program
! ifc -r8 -I/usr/local/shared/nag/nagf90_mod -L/usr/local/shared/nag/nagf90_mod spectrum_RK_OPO.f90 -o eval_spectr.o -lnag -lnagfl90 -Vaxlib

      ! on tuesday:
      ! ifort -r8 spectrum_RK_OPO.f90 -o eval_spectr.o -I/opt/NAG/fnl6i04dcl/nagfl90_modules/ /opt/NAG/fnl6i04dcl/lib/libnagfl90_nag.a /opt/NAG/fll6i22dcl/lib/libnag_nag.a

      Program dynamika
        USE nag_fft, ONLY: nag_fft_2d, nag_fft_3d, nag_fft_trig
        USE global
        USE ode_path 
        IMPLICIT NONE

        write(*,*) "evaluate spectrum after trigger as well? yes(1) no(0)"
        read(*,*) trigg

        open(unit=file_time_0, file="times_0.dat", status='replace')
        if (trigg.eq.1) then
           open(unit=file_time_1, file="times_1.dat", status='replace')
        end if

        CALL read 

        !check that t_init>=t_stfft+(Nt-1)*dxsav
        if(t_init.lt.t_stfft+(Nt-1)*dxsav) stop
        
        if(trigg.eq.0) then
           !tot_h=t_stfft+(Nt-1)*(dxsav)
           tot_h=t_stfft+(Nt+1)*(dxsav)
        else if (trigg.eq.1) then
           !tot_h=t_stfft2+(Nt-1)*(dxsav)
           tot_h=t_stfft2+(Nt+1)*(dxsav)
        end if

        allocate(trig_x(2*Nx), trig_y(2*Ny), trig_t(2*Nt))
        allocate(pdb(Nx,Ny,2), kinetic(Nx,Ny), pump_spatial(Nx,Ny))
        allocate(probe_spatial(Nx,Ny))
        allocate(pump(Nx,Ny), probe(Nx,Ny))
        allocate(y_tot_0(Nx,Ny,Nt), y_tot_1(Nx,Ny,Nt), int_sp(Nt))

        call nag_fft_trig(trig_x)
        call nag_fft_trig(trig_y)
        call nag_fft_trig(trig_t)
        
        call init
        call setg
        
        CALL dynamic
        call eval_spectr_0

        if(trigg.eq.1) call eval_spectr_1
           
        deallocate(trig_x, trig_y, trig_t)
        deallocate(pdb, kinetic, pump_spatial)
        deallocate(probe_spatial)
        deallocate(pump, probe)
        deallocate(y_tot_0, y_tot_1, int_sp)

        close(file_time_0)
        if(trigg.eq.1) then
           close(file_time_1)
        end if
           
      end Program dynamika

      !-----------------------------------------------------------------!
      ! Subroutine read                                                 !
      !-----------------------------------------------------------------!
      Subroutine read
        USE global   
        USE ode_path 
        IMPLICIT NONE

        OPEN(unit=22, file='INPUT_spectr', status='old')
        read(22,fmt="(10x,F10.6)")  kappa_C  !decay rate for photons
        read(22,fmt="(10x,F10.6)")  kappa_X  !decay rate for excitons
        read(22,fmt="(10x,F10.6)")  delta    ! detuning
        read(22,fmt="(10x,I2)")     Lz       ! Lz of the vortex, Lz=1, 2, ...
        read(22,fmt="(10x,F10.6)")  k_p      !pump  angle
        read(22,fmt="(10x,F10.6)")  k_pb     !probe angle
        read(22,fmt="(10x,F10.6)")  omega_p  !pump frequency
        read(22,fmt="(10x,F10.6)")  omega_pb !probe frequency
        read(22,fmt="(10x,F10.6)")  sigma_p  !pump  spatial
        read(22,fmt="(10x,F10.6)")  sigma_pb !probe spatial
        read(22,fmt="(10x,F10.6)")  sigma_t  !pulse 
        read(22,fmt="(10x,F10.6)")  t_init   !initial time of the pulse
        read(22,fmt="(10x,F10.6)")  f_p      ! pump strenght
        read(22,fmt="(10x,F10.6)")  f_pb     ! probe strenght
        read(22,fmt="(10x,F10.6)")  assymx   ! x position of probe
        read(22,fmt="(10x,F10.6)")  assymy   ! y position of probe
        read(22,fmt="(10x,F10.6)")  mom_cent ! center of filtering in momentum
        read(22,fmt="(10x,F10.6)")  mom_cut  ! filtering in momentum
        read(22,fmt="(10x,F10.6)")  Lx       ! x size of box (L)
        read(22,fmt="(10x,F10.6)")  Ly       ! y size of box (L)
        read(22,fmt="(10x,I5)")     Nx       ! Number of points in x direction
        read(22,fmt="(10x,I5)")     Ny       ! Number of points in y direction
        read(22,fmt="(10x,F12.10)") dxsav    ! time interval for fft
        read(22,fmt="(10x,F12.10)") eps_r
        read(22,fmt="(10x,I5)")     in_sswf  ! initial ss wf to read
        read(22,fmt="(10x,F12.10)") t_stfft  ! starting time for fft n.b. has to be multiple of dxsav
        read(22,fmt="(10x,F12.10)") t_stfft2 ! starting time for second fft n.b. has to be multiple of dxsav
        read(22,fmt="(10x,I5)")     Nt       ! resolution in energy (even)
        CLOSE(22)

        ax=2.0*Lx/Nx
        ay=2.0*Ly/Ny
        norm=ax*ay

      end Subroutine read
 
      !-----------------------------------------------------------------!
      ! Subroutine init                                                 !
      !-----------------------------------------------------------------!
      SUBROUTINE  init
        USE global
        USE nrtype
        IMPLICIT NONE

        integer :: ix, iy
        complex :: theta
        real  :: sx, sy, erre
        real :: re_y1(Nx,Ny), im_y1(Nx,Ny)
        real :: re_y2(Nx,Ny), im_y2(Nx,Ny)

        !choose the starting 
        !wavefunction psi0 in space
        pdb=(0.0,0.0)

        write(label,FMT="(i3)") in_sswf
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

              pdb(ix,iy,1)=(1.0,0.0)*re_y1(ix,iy)+(0.0,1.0)*im_y1(ix,iy)
              pdb(ix,iy,2)=(1.0,0.0)*re_y2(ix,iy)+(0.0,1.0)*im_y2(ix,iy)
           end do
           read(22,*)
           read(23,*)
        end do
        close(22)
        close(23)

        open(unit=23, file='probe.dat', status='replace')
        open(unit=24, file='probe_phase.dat', status='replace')
        open(unit=25, file='pump.dat', status='replace')

        !choose the pump profile in space
        do iy=1, Ny
           sy=-Ly+(iy-1)*ay
           do ix=1, Nx
              sx=-Lx+(ix-1)*ax
              
              !gaussian pump
              !pump_spatial(ix,iy)=f_p*exp(-0.5*( sx*sx+sy*sy )/(sigma_p*sigma_p)) +&
              !     (0.0,0.0)
              !smoothen top hat pump
              pump_spatial(ix,iy)=f_p*0.5*&
                   ( tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)+sigma_p ))-&
                   tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)-sigma_p )) ) + (0.0,0.0)

              !add a small probe with a vortex
              erre=sqrt( (sx-assymx)*(sx-assymx)+(sy-assymy)*(sy-assymy) )
              if(erre .eq. 0.0) then
                 theta=(0.0,0.0)
                 probe_spatial(ix,iy)=(0.0,0.0)
              else
                 theta=Lz*(  (0.0,-1.0)*&
                      log( ((sx-assymx)+(0.0,1.0)*(sy-assymy))/erre ) + pi  )
                 !!vortex Lz=2
                 !probe_spatial(ix,iy)=f_pb*( exp(1.0)/(2*sigma_pb*sigma_pb) )*&
                 !     erre*erre*&
                 !     exp(-0.5*(erre*erre)/(sigma_pb*sigma_pb))*&
                 !     exp((0.0,1.0)*theta) + (0.0,0.0)
                 !vortex Lz=1
                 probe_spatial(ix,iy)=f_pb*( exp(0.5)/(sigma_pb) )*&
                      erre*&
                      exp(-0.5*(erre*erre)/(sigma_pb*sigma_pb))*&
                      exp((0.0,1.0)*theta) + (0.0,0.0)
              end if
                 
              write(23,*) sx, sy, abs(probe_spatial(ix,iy))
              write(24,*) sx, sy, real(theta)
              !write(24,*) sx, sy, real(pi+(0.0,-1.0)*log( -probe_spatial(ix,iy)/( abs(probe_spatial(ix,iy)) ) ))
              write(25,*) sx, sy, abs(pump_spatial(ix,iy))
           end do
           write(23,*)
           write(24,*)
           write(25,*)
        end do
        close(23)
        close(24)
        close(25)
        
        ! multiply by k_p, k_pb part
        do iy=1, Ny
           sy=-Ly+(iy-1)*ay
           do ix=1, Nx
              sx=-Lx+(ix-1)*ax
              pump_spatial(ix,iy)= pump_spatial(ix,iy)*cos(k_p*sx)+(0.0,1.0)*pump_spatial(ix,iy)*sin(k_p*sx)
           end do
        end do
        
        do iy=1, Ny
           sy=-Ly+(iy-1)*ay
           do ix=1, Nx
              sx=-Lx+(ix-1)*ax
              probe_spatial(ix,iy)=probe_spatial(ix,iy)*cos(k_pb*sx)+(0.0,1.0)*probe_spatial(ix,iy)*sin(k_pb*sx)
           end do
        end do
        
        call setg

        RETURN
      END SUBROUTINE  init

      !-----------------------------------------------------------------!
      ! Subroutine setg to calculate k*k on FFT grid                    !
      !-----------------------------------------------------------------!
      Subroutine setg
        use global
        USE nrtype
        implicit none
        
        integer :: j,k
        
        DO j=1,(Ny/2+1)
           DO k=1,(Nx/2+1)
              kinetic(k,j)=pi*pi*(&
                   &(k-1)*(k-1)/(Lx*Lx)+(j-1)*(j-1)/(Ly*Ly))
           END DO
        END DO
        
        DO j=(Ny/2+2),Ny
           DO k=(Nx/2+2),Nx
              kinetic(k,j)=pi*pi*( &
                   & (k-1-Nx)*(k-1-Nx)/(Lx*Lx)+(j-1-Ny)*(j-1-Ny)/(Ly*Ly))
           END DO
        END DO
        
        DO j=1,(Ny/2+1)
           DO k=(Nx/2+2),Nx
              kinetic(k,j)=pi*pi*(&
                   &(k-1-Nx)*(k-1-Nx)/(Lx*Lx)+(j-1)*(j-1)/(Ly*Ly))
           END DO
        END DO
        
        DO j=(Ny/2+2),Ny
           DO k=1,(Nx/2+1)  
              kinetic(k,j)=pi*pi*(&
                   &(k-1)*(k-1)/(Lx*Lx)+(j-1-Ny)*(j-1-Ny)/(Ly*Ly))
           END DO
        END DO
        
        return
      end subroutine setg
      
      !-----------------------------------------------------------------!
      ! Subroutine dynamic                                              !
      !-----------------------------------------------------------------!
      Subroutine  dynamic
        USE global
        USE nrtype
        USE ode_path 
        IMPLICIT NONE
        
        INTERFACE
           SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
             USE nrtype
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: ystart
             REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
             INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                  USE nrtype
                  REAL(SP), INTENT(IN) :: x
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
                END SUBROUTINE derivs
                SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
                  USE nrtype
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: y
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal
                  REAL(SP), INTENT(INOUT) :: x
                  REAL(SP), INTENT(IN) :: htry,eps
                  REAL(SP), INTENT(OUT) :: hdid,hnext
                  INTERFACE
                     SUBROUTINE derivs(x,y,dydx)
                       USE nrtype
                       REAL(SP), INTENT(IN) :: x
                       COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
                       COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
                     END SUBROUTINE derivs
                  END INTERFACE
                END SUBROUTINE rkqs
             END INTERFACE
           END SUBROUTINE odeint
           SUBROUTINE derivs(x,y,dydx)
             USE nrtype
             IMPLICIT NONE
             REAL(SP), INTENT(IN) :: x
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
           END SUBROUTINE derivs
           SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
             USE nrtype
             IMPLICIT NONE
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal
             REAL(SP), INTENT(INOUT) :: x
             REAL(SP), INTENT(IN) :: htry,eps
             REAL(SP), INTENT(OUT) :: hdid,hnext
             INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                  USE nrtype
                  IMPLICIT NONE
                  REAL(SP), INTENT(IN) :: x
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
                END SUBROUTINE derivs
             END INTERFACE
           END SUBROUTINE rkqs
        END INTERFACE
        
        INTEGER :: i,j
        REAL, EXTERNAL ::  findfermpart,findcoopernum
        COMPLEX, EXTERNAL :: findfermcond1, findfermcond2
        !COMPLEX(SP), DIMENSION(:),ALLOCATABLE :: ystart_r
        REAL(SP)  :: x1_r,x2_r,h1_r,hmin_r
        
        x1_r=0.0
        x2_r=tot_h
        !eps_r=0.0001 
        h1_r=0.001
        hmin_r=0.0
        
        !ALLOCATE(ystart_r(2))
        
        !ystart_r(1)=CMPLX(1.0,1.0)
        !ystart_r(2)=CMPLX(1.0,1.0)

        CALL odeint(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r,derivs,rkqs)
        !write(*,*) ystart_r(1) 
        !write(*,*) ystart_r(2)
        
        !DO i=1, kount
        !write(71,*) xp(i), REAL(yp(Ne,i)*CONJG(yp(Ne,i)))
        !ENDDO
        
        RETURN
      END SUBROUTINE  dynamic

      !-----------------------------------------------------------------!
      ! Subroutine derivs(x,y,dydx)                                     !
      !-----------------------------------------------------------------!
      Subroutine derivs(x,y,dydx)
        use nag_fft, only: nag_fft_2d, nag_fft_trig
        USE global
        USE nrtype
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
        INTEGER :: ix,iy
        REAL :: sx, sy,part_r,part_i
        COMPLEX, DIMENSION(size(pdb,1),size(pdb,2)) :: pom
        
        pump=pump_spatial*(cos(omega_p*x)+(0.0,-1.0)*sin(omega_p*x))
        probe=probe_spatial*(cos(omega_pb*x)+(0.0,-1.0)*sin(omega_pb*x))*&
             &exp(-0.5*(x-t_init)*(x-t_init)/(sigma_t*sigma_t))
        
        !DO iy=1, Ny
        !   DO ix=1,Nx
        dydx(:,:,1)= (0.0,-1.0)*y(:,:,2)-(0.0,1.0)*delta*y(:,:,1)-&
             kappa_C*y(:,:,1)+(0.0,-1.0)*(pump(:,:)+probe(:,:))
        dydx(:,:,2)= -kappa_X*y(:,:,2)+&
             (0.0,-1.0)*ABS(y(:,:,2))*ABS(y(:,:,2))*y(:,:,2)/norm+ &
             &(0.0,-1.0)*y(:,:,1)
        
        !   ENDDO
        !ENDDO
        
        !fft to momentum space
        pom(:,:)=nag_fft_2d(y(:,:,1),trig_m=trig_x,trig_n=trig_y)

        pom=kinetic*pom
        
        !fft back to real space
        pom(:,:)=nag_fft_2d(pom(:,:),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
        
        dydx(:,:,1)=dydx(:,:,1)+pom*CMPLX(0.0,-1.0)
        
        !dydx=dydx/sqrt(norm)

      END SUBROUTINE derivs

      !-----------------------------------------------------------------!
      ! Subroutine eval_spectr_0                                        !
      !-----------------------------------------------------------------!
      Subroutine eval_spectr_0
        use nag_fft, only: nag_fft_2d, nag_fft_trig,  nag_fft_3d
        use global
        use nrtype
        USE ode_path 
        implicit none

        integer :: i_t
        integer :: ix, iy, kx, ky
        integer :: i_tmax(1), kx_max(1)
        complex :: y_enfilt(Nx,Ny)
        real :: omega
        real :: sx, sy, mom_x, mom_y
        real :: omega_pmax, omega_smax, omega_imax
        real :: int_sp_pom(Nt)
        real :: mom_x_maxs, mom_x_maxi

        !fft to momentum and energy space
        y_tot_0(:,:,:)=&
             nag_fft_3d(y_tot_0(:,:,:),trig_1=trig_x,trig_2=trig_y,trig_3=trig_t)
                
        !write the wave function in momentum and energy space for ky=0
        open(unit=23, file="spectr_om-vs-kx_no-trigg.dat", status='replace')
        write(23, fmt=' ("#", 1x, "mom_x", 19x, "omega", 19x, "abs(psi(1))**2") ')
        
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
              
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           write(23,*)
        end do
                
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           write(23,*)
        end do
        close(23)

        !evaluate and write the integrated spectrum
        open(unit=24, file="int-spectr_no-trigg.dat", status='replace')
        write(24, fmt=' ("#", 1x, "omega", 20x, "int_omega") ')

        int_sp=sum(sum(abs(y_tot_0), dim=1), dim=1)
           
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           write(24, *) -omega, int_sp(i_t)
        end do
                
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
         
           write(24, *) -omega, int_sp(i_t)
        end do
        close(24)

        !evaluate peaks at pump, signal and idler
        open(unit=25, file="en-peaks.dat", status='replace')
        
        !find location of the pump energy, omega_pmax
        i_tmax=maxloc(int_sp)
        
        if( i_tmax(1).ge.Nt/2+2 ) &
             omega_pmax=2.0*pi*(i_tmax(1)-1-Nt)/( (Nt-1)*dxsav )
        if( i_tmax(1).le.Nt/2+1) & 
             omega_pmax=2.0*pi*(i_tmax(1)-1)/( (Nt-1)*dxsav )
                
        write(25,*) "omega_pmax=", -omega_pmax
        write(25,*) "int_pmax=", maxval(int_sp)

        !find location of the next highest peak after omega_pmax (signal)
        !n.b. remember the code evaluate the energies with the wrong sign
        int_sp_pom=int_sp
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )

           if( abs(omega-omega_pmax).lt.0.4d0 ) int_sp_pom(i_t)=0.0
        end do
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
         
           if( abs(omega-omega_pmax).lt.0.4d0 ) int_sp_pom(i_t)=0.0
        end do
        
        i_tmax(1)=Nt+1
        i_tmax=maxloc(int_sp_pom)
        
        if(i_tmax(1).ne.Nt+1) then
           if( i_tmax(1).ge.Nt/2+2 ) &
                omega_smax=2.0*pi*(i_tmax(1)-1-Nt)/( (Nt-1)*dxsav )
           if( i_tmax(1).le.Nt/2+1) &
                omega_smax=2.0*pi*(i_tmax(1)-1)/( (Nt-1)*dxsav )
           
           write(25,*) "omega_smax=", -omega_smax
           write(25,*) "int_smax=", maxval(int_sp_pom)

           !write the signal en-filtered wf in momentum
           open(unit=27, file="opo_ph-mom_enfilt_signal.dat", status='replace')
           write(27, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')
           
           y_enfilt(:,:)=y_tot_0(:,:,i_tmax(1))
           
           do ky=Ny/2+2, Ny
              mom_y=pi*(ky-1-Ny)/Ly
              !region(4)
              do kx=Nx/2+2, Nx
                 mom_x=pi*(kx-1-Nx)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              !region(3)
              do kx=1, Nx/2+1
                 mom_x=pi*(kx-1)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              write(27,*)
           end do
           
           do ky=1, Ny/2+1
              mom_y=pi*(ky-1)/Ly
              !region(2)
              do kx=Nx/2+2, Nx
                 mom_x=pi*(kx-1-Nx)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              !region(1)
              do kx=1, Nx/2+1
                 mom_x=pi*(kx-1)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              write(27,*)
           end do
           close(27)
           
           !evaluate the signal peak in momentum
           kx_max=maxloc( abs(y_enfilt(:,1)) )
           !region(2)
           if ( kx_max(1).ge.Nx/2+2 ) mom_x_maxs=pi*(kx_max(1)-1-Nx)/Lx
           !region(1)
           if ( kx_max(1).le.Nx/2+1 ) mom_x_maxs=pi*(kx_max(1)-1)/Lx

           write(25,*) "mom_x_maxs=", mom_x_maxs
           
           !write the signal en-filtered wf in space
           open(unit=26, file="opo_ph-spc_enfilt_signal.dat", status='replace')
           write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')
           
           y_enfilt(:,:)=nag_fft_2d(y_tot_0(:,:,i_tmax(1)),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
           
           do iy=1, Ny
              sy=-Ly+(iy-1)*ay
              do ix=1, Nx
                 sx=-Lx+(ix-1)*ax
                 write(26,*) sx, sy, abs(y_enfilt(ix,iy))*abs(y_enfilt(ix,iy))
                 
              end do
              write(26,*)
           end do
           close(26)
        end if
        
        !find location of the next highest peak after 
        !omega_pmax & omega_smax (idler)
        !n.b. remember the code evaluate the energies with the wrong sign
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           if( abs(omega-omega_smax).lt.0.6d0 ) int_sp_pom(i_t)=0.0
        end do
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
           
           if( abs(omega-omega_smax).lt.0.6d0 ) int_sp_pom(i_t)=0.0
        end do
        
        i_tmax(1)=Nt+1
        i_tmax=maxloc(int_sp_pom)
        
        if(i_tmax(1).ne.Nt+1) then
           if( i_tmax(1).ge.Nt/2+2 ) &
                omega_imax=2.0*pi*(i_tmax(1)-1-Nt)/( (Nt-1)*dxsav )
           if( i_tmax(1).le.Nt/2+1) &
                omega_imax=2.0*pi*(i_tmax(1)-1)/( (Nt-1)*dxsav )
           
           write(25,*) "omega_imax=", -omega_imax
           write(25,*) "int_imax=", maxval(int_sp_pom)

           !write the idler en-filtered wf in momentum
           open(unit=27, file="opo_ph-mom_enfilt_idler.dat", status='replace')
           write(27, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')
           
           y_enfilt(:,:)=y_tot_0(:,:,i_tmax(1))
           
           do ky=Ny/2+2, Ny
              mom_y=pi*(ky-1-Ny)/Ly
              !region(4)
              do kx=Nx/2+2, Nx
                 mom_x=pi*(kx-1-Nx)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              !region(3)
              do kx=1, Nx/2+1
                 mom_x=pi*(kx-1)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              write(27,*)
           end do
           
           do ky=1, Ny/2+1
              mom_y=pi*(ky-1)/Ly
              !region(2)
              do kx=Nx/2+2, Nx
                 mom_x=pi*(kx-1-Nx)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              !region(1)
              do kx=1, Nx/2+1
                 mom_x=pi*(kx-1)/Lx
                 
                 write(27,*) mom_x, mom_y, &
                      abs(y_enfilt(kx,ky))*abs(y_enfilt(kx,ky))
              end do
              write(27,*)
           end do
           close(27)
           
           !evaluate the idler peak in momentum
           kx_max=maxloc( abs(y_enfilt(:,1)) )
           !region(2)
           if ( kx_max(1).ge.Nx/2+2 ) mom_x_maxi=pi*(kx_max(1)-1-Nx)/Lx
           !region(1)
           if ( kx_max(1).le.Nx/2+1 ) mom_x_maxi=pi*(kx_max(1)-1)/Lx

           write(25,*) "mom_x_maxi=", mom_x_maxi
           
           !write the signal en-filtered wf in space
           open(unit=26, file="opo_ph-spc_enfilt_idler.dat", status='replace')
           write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')
           
           y_enfilt(:,:)=nag_fft_2d(y_tot_0(:,:,i_tmax(1)),inverse=.true.,trig_m=trig_x,trig_n=trig_y)
           
           do iy=1, Ny
              sy=-Ly+(iy-1)*ay
              do ix=1, Nx
                 sx=-Lx+(ix-1)*ax
                 write(26,*) sx, sy, abs(y_enfilt(ix,iy))*abs(y_enfilt(ix,iy))
                 
              end do
              write(26,*)
           end do
           close(26)
           
           write(25,*) "2*omega_pmax=", -2.0*omega_pmax
           write(25,*) "omega_smax+omega_imax=", -omega_smax-omega_imax
        
           write(25,*) "mom_x_maxs+mom_x_maxi", mom_x_maxs+mom_x_maxi

        end if
        
        close(25)

      end Subroutine eval_spectr_0
      
      !-----------------------------------------------------------------!
      ! Subroutine eval_spectr_1                                        !
      !-----------------------------------------------------------------!
      Subroutine eval_spectr_1
        use nag_fft, only: nag_fft_2d, nag_fft_trig,  nag_fft_3d
        use global
        use nrtype
        USE ode_path 
        implicit none

        integer :: i_t
        integer :: ix, iy, kx, ky
        integer :: loc_max_0(3), loc_max_1(3)
        real :: omega
        real :: sx, sy, mom_x, mom_y
        real :: max_sqy_0, max_sqy_1
        real :: omega_pmax_0, omega_pmax_1
        real :: mom_x_max_0, mom_x_max_1
        real :: mom_y_max_0, mom_y_max_1

        !fft to momentum and energy space
        y_tot_1(:,:,:)=&
             nag_fft_3d(y_tot_1(:,:,:),trig_1=trig_x,trig_2=trig_y,trig_3=trig_t)
                
        !write the wave function in momentum and energy space for ky=0
        open(unit=23, file="spectr_om-vs-kx_full-trigg.dat", status='replace')
        write(23, fmt=' ("#", 1x, "mom_x", 19x, "omega", 19x, "abs(psi(1))**2") ')
        
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))
              
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))
           end do
           write(23,*)
        end do
                
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(23, *) mom_x, -omega, abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))
           end do
           write(23,*)
        end do
        close(23)

        !evaluate and write the integrated spectrum
        open(unit=24, file="int-spectr_full-trigg.dat", status='replace')
        write(24, fmt=' ("#", 1x, "omega", 20x, "int_omega") ')

        int_sp=sum(sum(abs(y_tot_1), dim=1), dim=1)
           
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           write(24, *) -omega, int_sp(i_t)
        end do
                
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
         
           write(24, *) -omega, int_sp(i_t)
        end do
        close(24)
        
        max_sqy_0=maxval( abs(y_tot_0)*abs(y_tot_0) )
        max_sqy_1=maxval( abs(y_tot_1)*abs(y_tot_1) )

        open(unit=25, file="spectr_om-vs-kx_subtr-trigg.dat", status='replace')
        write(25, fmt=' ("#", 1x, "mom_x", 19x, "omega", 19x, "abs(psi(1))**2") ')        
        write(25,*) "#maxval_sqy_0=", max_sqy_0
        write(25,*) "#maxval_sqy_1=", max_sqy_1
        
        loc_max_0=maxloc( abs(y_tot_0)*abs(y_tot_0) )
        loc_max_1=maxloc( abs(y_tot_1)*abs(y_tot_1) )

        if( loc_max_0(3).ge.Nt/2+2 ) &
             omega_pmax_0=2.0*pi*(loc_max_0(3)-1-Nt)/( (Nt-1)*dxsav )
        if( loc_max_0(3).le.Nt/2+1) & 
             omega_pmax_0=2.0*pi*(loc_max_0(3)-1)/( (Nt-1)*dxsav )
        if( loc_max_1(3).ge.Nt/2+2 ) &
             omega_pmax_1=2.0*pi*(loc_max_1(3)-1-Nt)/( (Nt-1)*dxsav )
        if( loc_max_1(3).le.Nt/2+1) & 
             omega_pmax_1=2.0*pi*(loc_max_1(3)-1)/( (Nt-1)*dxsav )

        write(25,*) "#omega_pmax_0=", omega_pmax_0
        write(25,*) "#omega_pmax_1=", omega_pmax_1

        if ( loc_max_0(1).ge.Nx/2+2 ) mom_x_max_0=pi*(loc_max_0(1)-1-Nx)/Lx
        if ( loc_max_0(1).le.Nx/2+1 ) mom_x_max_0=pi*(loc_max_0(1)-1)/Lx
        if ( loc_max_1(1).ge.Nx/2+2 ) mom_x_max_1=pi*(loc_max_1(1)-1-Nx)/Lx
        if ( loc_max_1(1).le.Nx/2+1 ) mom_x_max_1=pi*(loc_max_1(1)-1)/Lx
        
        write(25,*) "#mom_x_max_0=", mom_x_max_0
        write(25,*) "#mom_x_max_1=", mom_x_max_1

        if ( loc_max_0(2).ge.Nx/2+2 ) mom_y_max_0=pi*(loc_max_0(2)-1-Nx)/Lx
        if ( loc_max_0(2).le.Nx/2+1 ) mom_y_max_0=pi*(loc_max_0(2)-1)/Lx
        if ( loc_max_1(2).ge.Nx/2+2 ) mom_y_max_1=pi*(loc_max_1(2)-1-Nx)/Lx
        if ( loc_max_1(2).le.Nx/2+1 ) mom_y_max_1=pi*(loc_max_1(2)-1)/Lx

        write(25,*) "#mom_y_max_0=", mom_y_max_0
        write(25,*) "#mom_y_max_1=", mom_y_max_1
        
        do i_t=Nt/2+2, Nt
           omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(25, *) mom_x, -omega,&
                   abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))-&
                   abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
              
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(25, *) mom_x, -omega,&
                   abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))-&
                   abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           write(25,*)
        end do
                
        do i_t=1, Nt/2+1
           omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav )
           
           !region(2)
           do kx=Nx/2+2, Nx
              mom_x=pi*(kx-1-Nx)/Lx
              write(25, *) mom_x, -omega,&
                   abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))-&
                   abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           !region(1)
           do kx=1, Nx/2+1
              mom_x=pi*(kx-1)/Lx
              write(25, *) mom_x, -omega,&
                   abs(y_tot_1(kx,1,i_t))*abs(y_tot_1(kx,1,i_t))-&
                   abs(y_tot_0(kx,1,i_t))*abs(y_tot_0(kx,1,i_t))
           end do
           write(25,*)
        end do
        close(25)
        
      end Subroutine eval_spectr_1

      !-----------------------------------------------------------------!
      ! Subroutine odeint                                               !
      !-----------------------------------------------------------------!
      SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
        use global
        USE nrtype 
        !USE nrutil, ONLY : nrerror,reallocate
        USE ode_path
        IMPLICIT NONE
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: ystart
        REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
        Integer :: i
        INTERFACE
           !SUBROUTINE reallocate_rv(p,n)
           !USE nrtype 
           !REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
           !INTEGER(I4B), INTENT(IN) :: n
           ! yp=>reallocate(yp,size(yp,1),size(xp))
           SUBROUTINE derivs(x,y,dydx)
             USE nrtype
             IMPLICIT NONE
             REAL(SP), INTENT(IN) :: x
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
           END SUBROUTINE derivs
           SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
             USE nrtype
             IMPLICIT NONE
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal
             REAL(SP), INTENT(INOUT) :: x
             REAL(SP), INTENT(IN) :: htry,eps
             REAL(SP), INTENT(OUT) :: hdid,hnext
             INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                  USE nrtype
                  IMPLICIT NONE
                  REAL(SP), INTENT(IN) :: x
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
                END SUBROUTINE derivs
             END INTERFACE
           END SUBROUTINE rkqs
        END INTERFACE
        
!!$INTERFACE reallocate
!!$MODULE PROCEDURE reallocate_rv,reallocate_rm,&
!!$reallocate_iv,reallocate_im,reallocate_hv, reallocate_new
!!$END INTERFACE

        REAL(SP), PARAMETER :: TINY=1.0e-30_sp
        
        INTEGER(I4B), PARAMETER :: MAXSTP=1000000000
!Runge-Kutta driver with adaptive step size control.Integrate the array
!of starting values ystart from x1 to x2 with accuracy eps storing
!intermediate results in the module variables in ode_path. h1 should be
!set as a guessed first stepsize, hmin as the minimum allowed stepsize
!(can be zero). On output ystart is replaced by values at the end of the
!integration interval. derivs is the user-supplied subroutine for
!calculating the right-hand-side derivative, while rkqs is the name of
!he stepper routine to be used.
        INTEGER(I4B) :: nstp
        REAL(SP) :: h,hdid,hnext,x,xsav
        COMPLEX(SP), DIMENSION(size(ystart,1),size(ystart,2),size(ystart,3)) :: dydx,y,yscal
        REAL :: dphase
        !Complex :: pom_cond
        REAL, External :: findfermpart,findcoopernum
        COMPLEX, External :: findfermcond
        !REAL(SP), DIMENSION(:), POINTER, EXTERNAL:: reallocate_rv

        !xsav=0.01
        !write(*,*) xsav
        
        x=x1
        h=sign(h1,x2-x1)
        nok=0
        nbad=0
        kount=0
        kount_1=0
        y(:,:,:)=ystart(:,:,:)

!nullify(xp,yp)
!Pointers nullified here, but memory not deallocated. If odeint is
!called multiple times, calling program should deallocate xp and yp
!between calls.

        if (save_steps) then
           xsav=x-2.0_sp*dxsav
           !   allocate(xp(256))
           !   allocate(yp(size(ystart),size(xp)))
        end if
        do nstp=1,MAXSTP 
           !Take at most MAXSTP steps.
           !write(*,*) "stop1"
           !write(*,*) x, y(1),y(2)

           call derivs(x,y,dydx)

           yscal(:,:,:)=abs(y(:,:,:))+abs(h*dydx(:,:,:))+TINY

!!$           !Scaling used to monitor accuracy. This general purpose choice can be
!!$           !modified if need be.
!!$           if (save_steps .and. (abs(x-xsav) > abs(dxsav))) then
!!$                !Store intermedia e results.
!!$
!!$              !write(*,*) 'x=', x
!!$              !write(*,*) 'xsav=', xsav
!!$              call save_a_step
!!$           end if

           if (save_steps .and. x.ge.t_stfft ) then
              call save_a_step
           end if
           
           if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x 
           !If stepsize can overshoot,decrease.
           call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
           if (hdid == h) then
              nok=nok+1
           else
              nbad=nbad+1
           end if
           if ((x-x2)*(x2-x1) >= 0.0) then 
              !Are we done?
              ystart(:,:,:)=y(:,:,:)
              if (save_steps) call save_a_step 
              !Save final step.
              RETURN 
              !Normal exit.
           end if
           if (abs(hnext) < hmin) write(*,*) "stepsize smaller than minimum in odeint"
           h=hnext
        end do

        write(*,*) "too many steps in odeint"

      CONTAINS
        
        !-----------------------------------------------------------------!
        ! Subroutine  save_a_step                                         !
        !-----------------------------------------------------------------!
        SUBROUTINE save_a_step
          !use nag_fft, only: nag_fft_2d, nag_fft_trig
          USE global
          USE nrtype

          !real :: sx,sy
          !integer :: ix, iy
          integer :: i_t

          if( x.ge.t_stfft+kount*dxsav .and. kount+1.le.Nt) then
             write(file_time_0,*) kount, x
             call flush(file_time_0)
             i_t=kount+1
             
             y_tot_0(:,:,i_t)=y(:,:,1)
             
             kount=kount+1
          end if
          
          if (trigg.eq.1) then
             
             if( x.ge.t_stfft2+kount_1*dxsav .and. kount_1+1.le.Nt) then
                write(file_time_1,*) kount_1, x
                call flush(file_time_1)
                i_t=kount_1+1
             
                y_tot_1(:,:,i_t)=y(:,:,1)
             
                kount_1=kount_1+1
                
             end if
          end if

          !if (kount > size(xp)) then
          !   xp=>reallocate(xp,2*size(xp))
          !   yp=>reallocate(yp,size(yp,1),size(xp))
          !end if
          !xp(kount)=x
          !yp(:,kount)=y(:)
          xsav=x
          
!!$          write(label,FMT="(i3)") kount
!!$          
!!$          open(unit=22, file="phcplx-opo_spc"//trim(adjustl(label))//".dat", status='replace')
!!$          write(22, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(1))", 1x, "aimag(psi(1))") ')
!!$          open(unit=23, file="excplx-opo_spc"//trim(adjustl(label))//".dat", status='replace')
!!$          write(23, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(2))", 1x, "aimag(psi(2))") ')
!!$          
!!$          !If(Mod(kount,10).eq.0) then
!!$          do iy=1, Ny
!!$             sy=-Ly+(iy-1)*ay
!!$             do ix=1, Nx
!!$                sx=-Lx+(ix-1)*ax
!!$                write(22, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, real(y(ix,iy,1)), aimag(y(ix,iy,1))
!!$                write(23, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, real(y(ix,iy,2)), aimag(y(ix,iy,2))
!!$             end do
!!$             write(22,*)
!!$             write(23,*)
!!$          end do
!!$          close(22)
!!$          close(23)
          
        END SUBROUTINE save_a_step

      END SUBROUTINE odeint

      !-----------------------------------------------------------------!
      ! Subroutine  rkqs                                                !
      !-----------------------------------------------------------------!
      SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
        USE nrtype; 
        !USE nrutil, ONLY : assert_eq,nrerror
        !USE nr, ONLY : rkck
        IMPLICIT NONE
        !EXTERNAL  rkck
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(INOUT) :: y
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal
        REAL(SP), INTENT(INOUT) :: x
        REAL(SP), INTENT(IN) :: htry,eps
        REAL(SP), INTENT(OUT) :: hdid,hnext
        INTEGER, EXTERNAL :: assert_eq
        INTERFACE
           SUBROUTINE derivs(x,y,dydx)
             USE nrtype
             IMPLICIT NONE
             REAL(SP), INTENT(IN) :: x
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
           END SUBROUTINE derivs
           SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
             USE nrtype
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y,dydx
             REAL(SP), INTENT(IN) :: x,h
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: yout,yerr
             INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                  USE nrtype
                  REAL(SP), INTENT(IN) :: x
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
                  COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
                END SUBROUTINE derivs
             END INTERFACE
           END SUBROUTINE rkck
           !SUBROUTINE derivs(x,y,dydx)
        END INTERFACE

!Fifth order Runge-Kutta step with monitoring of local runca ion error
!o ensure accuracy and adjus stepsize.Input are he dependent variable
!vector y and its derivative dydx at the star ing value of the
!independent variable x Also input are he s epsize to be attempted htry
!he required accuracy eps and the vector yscal agains which the error
!is scaled.y dydx and yscal areallof hesamelengh.Onoupu,y and x are
!replaced by their new values,hdid is he s epsize ha was ac ually
!accomplished,and hnext is he estimated nex stepsize.derivs is the
!user-supplied subroutine ha computes the right-hand-side deriva ives.
        INTEGER(I4B) :: ndum
        REAL(SP) :: errmax,h,htemp,xnew
        COMPLEX(SP), DIMENSION(size(y,1),size(y,2),size(y,3)) :: yerr,ytemp
        REAL(SP), PARAMETER :: SAFETY=0.9_sp,PGROW=-0.2_sp,PSHRNK=-0.25_sp,&
             ERRCON=1.89e-4
        !The value ERRCON equals (5/SAFETY)**(1/PGROW),see use below.
        
        ndum=assert_eq(size(y),size(dydx),size(yscal), "rkqs") 
        h=htry 
        !Set step size to the initial trial value. 
        do
           call rkck(y,dydx,x,h,ytemp,yerr,derivs) 
           !Take a step.
           errmax=maxval(abs(yerr(:,:,:)/yscal(:,:,:)))/eps 
           !Evaluate accuracy.
           if (errmax <= 1.0) exit 
           !Step succeeded.
           htemp=SAFETY*h*(errmax**PSHRNK) 
           !Truncation error too large, reduce stepsize.
           h=sign(max(abs(htemp),0.1_sp*abs(h)),h) 
           !No more than a factor of 10.
           xnew=x+h
           if (xnew == x) write(*,*)  "stepsize underflow in rkqs"
        end do
        !Go back for another try.
        if (errmax > ERRCON) then 
           !Compute size of nex step.
           hnext=SAFETY*h*(errmax**PGROW) 
        else 
           !No more han a fac or of 5 increase. 
           hnext=5.0_sp*h
        end if
        hdid=h  
        x=x+h   
        y(:,:,:)=ytemp(:,:,:)  
      END SUBROUTINE rkqs

      !-----------------------------------------------------------------!
      ! Subroutine rkck                                                 !
      !-----------------------------------------------------------------!
      SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
        USE nrtype 
        !USE nrutil, ONLY : assert_eq
        IMPLICIT NONE
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y,dydx
        REAL(SP), INTENT(IN) :: x,h
        COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: yout,yerr
        INTEGER, EXTERNAL :: assert_eq
        INTERFACE
           SUBROUTINE derivs(x,y,dydx)
             USE nrtype
             IMPLICIT NONE
             REAL(SP), INTENT(IN) :: x
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(IN) :: y
             COMPLEX(SP), DIMENSION(:,:,:), INTENT(OUT) :: dydx
           END SUBROUTINE derivs
        END INTERFACE
!Given values for N variables y and their derivatives dydx known at x
!use the fifth order Cash-Karp Runge-Kutta method to advance the
!solution over an interval h and return the incremented variables as
!yout Also return an estimate of he local truncation error in yout
!using the embedded fourth order method. The user supplies the
!subroutine derivs(x,y,dydx),which returns derivatives dydx at x
        INTEGER(I4B) :: ndum
        COMPLEX(SP), DIMENSION(size(y,1),size(y,2),size(y,3)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        REAL(SP), PARAMETER :: A2=0.2_sp,A3=0.3_sp,A4=0.6_sp,A5=1.0_sp,&
             A6=0.875_sp,B21=0.2_sp,B31=3.0_sp/40.0_sp,B32=9.0_sp/40.0_sp,&
             B41=0.3_sp,B42=-0.9_sp,B43=1.2_sp,B51=-11.0_sp/54.0_sp,&
             B52=2.5_sp,B53=-70.0_sp/27.0_sp,B54=35.0_sp/27.0_sp,&
             B61=1631.0_sp/55296.0_sp,B62=175.0_sp/512.0_sp,&
             B63=575.0_sp/13824.0_sp,B64=44275.0_sp/110592.0_sp,&
             B65=253.0_sp/4096.0_sp,C1=37.0_sp/378.0_sp,&
             C3=250.0_sp/621.0_sp,C4=125.0_sp/594.0_sp,&
             C6=512.0_sp/1771.0_sp,DC1=C1-2825.0_sp/27648.0_sp,&
             DC3=C3-18575.0_sp/48384.0_sp,DC4=C4-13525.0_sp/55296.0_sp,&
             DC5=-277.0_sp/14336.0_sp,DC6=C6-0.25_sp
        
        ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),"rkck")
        ytemp=y+B21*h*dydx 
        !First step.
        call derivs(x+A2*h,ytemp,ak2) 
        !Second step.
        ytemp=y+h*(B31*dydx+B32*ak2)
        call derivs(x+A3*h,ytemp,ak3) 
        !Third step.
        ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs(x+A4*h,ytemp,ak4) 
        !Fourth step.
        ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs(x+A5*h,ytemp,ak5) 
        !Fifth step.
        ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs(x+A6*h,ytemp,ak6) 
        !Sixth step.
        yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6) 
        !Accumulate increments with proper weights.
        yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
        !Estimate error as diference between fourth and fifth order methods.
      END SUBROUTINE rkck

      !-----------------------------------------------------------------!
      ! Function assert_eq(n1,n2,n3,string)                             !
      !-----------------------------------------------------------------!
      FUNCTION assert_eq(n1,n2,n3,string)
        USE nrtype; 
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3
        INTEGER :: assert_eq
        if (n1 == n2 .and. n2 == n3) then
           assert_eq=n1
        else
           write (*,*) "nrerror: an assert_eq failed with this tag: &
                &program terminated by assert_eq3"
           STOP 
        end if
      END FUNCTION assert_eq





