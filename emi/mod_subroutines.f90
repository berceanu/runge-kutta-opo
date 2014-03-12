
!-----------------------------------------------------------------!
! Module subroutines                                              !
!-----------------------------------------------------------------!
Module mod_subroutines
  USE MKL_DFTI
  USE global
  IMPLICIT NONE
  SAVE
CONTAINS
  !-----------------------------------------------------------------!
  ! Subroutine read                                                 !
  !-----------------------------------------------------------------!
  Subroutine read
    IMPLICIT NONE
    REAL :: aux1,aux2

    OPEN (UNIT=30,FILE='INPUT',STATUS='old')
    READ(30,NML=indata)
    CLOSE(30)

    ax=2.0*Lx/Nx
    ay=2.0*Ly/Ny
    norm=ax*ay
    npointxy(1)=Nx
    npointxy(2)=Ny
    normDFTIsqrt=sqrt(1.0*Nx*Ny)

    ixbin=Nx/2+1+NINT(bin_x/ax)
    ixbfin=Nx/2+1+NINT(bfin_x/ax)

    iybin=Ny/2+1+INT(bin_y/ay)
    iybfin=Ny/2+1+INT(bfin_y/ay)

    aux1=0.5*(delta+(k_p*0.868)**2.)
    aux2=0.5*sqrt((delta/2.2+(k_p*0.868)**2.)**2.0+4.)
    thetasin=0.197*(k_p/0.868)/(1.528+0.0022*(aux1-aux2))
    write(*,*) thetasin,norm
!stop

  end Subroutine read

  !-----------------------------------------------------------------!
  ! Subroutine init                                                 !
  !-----------------------------------------------------------------!
  SUBROUTINE  init
    IMPLICIT NONE
    integer :: ix, iy
    complex :: theta
    real  :: sx, sy, erre,dist,bx,by,kpb_x,kpb_y

    bx=ixbfin-ixbin
    by=iybfin-iybin

    ALLOCATE(pdb(Nx,Ny,2),pot_x(Nx,NY),pot_c(Nx,Ny))
    ALLOCATE(kinetic(Nx,Ny),kvecx(Nx,Ny),kvecy(Ny,Ny))
    ALLOCATE(pump_spatial(Nx,Ny), probe_spatial(Nx,Ny))
    ALLOCATE(pump(Nx,Ny),probe(Nx,Ny))
    ALLOCATE(dragPOTx(Nx,Ny),dragPOTy(Nx,Ny))
    ALLOCATE(psiE(Nx,Ny,nene))
    ALLOCATE(psiref(Nx,Ny,nene))
    ALLOCATE(auxpsi(Nx,Ny,nene))
    ALLOCATE(transfE(nene))

    ti=1
    tti=1
    tavs=1
    tavm=1
    tavp=1
    psiE=(0.0,0.0)
    pdb=(0.0,0.0)
    pot_x=0.0
    pot_c=0.0

    !choose the pump profile in space
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax

          !gaussian pump
          pump_spatial(ix,iy)=f_p*exp(-0.5*( (sx-assymx_p)*(sx-assymx_p)+ &
               (sy-assymy_p)*(sy-assymy_p) )/(sigma_p*sigma_p))+&
               (0.0,0.0)

          !smoothen top hat pump with hyperbolic tangent
          !pump_spatial(ix,iy)=f_p*0.5*&
          !     ( tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)+sigma_p ))-&
          !     tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)-sigma_p )) ) + (0.0,0.0)

          !gaussian with decreasing sigma_p
          !dist=sqrt((sx-assymx_p)**2.+(sy-assymy_p)**2.)*sigma_p/(Lx*sqrt(2.))
          !pump_spatial(ix,iy)=f_p*exp(-0.5*((sx-assymx_p)/(sigma_p-dist))**2.)*&
          !        exp(-0.5*((sy-assymy_p)/(sigma_p-dist))**2.)

          !add a small Gaussian probe
          erre=sqrt((sx-assymx_pb)*(sx-assymx_pb)+(sy-assymy_pb)*(sy-assymy_pb))
          probe_spatial(ix,iy)=f_pb*&
               exp(-0.5*(erre*erre)/(sigma_pb*sigma_pb))&
               + (0.0,0.0)

          !add a Top Hat probe constant in time
          !probe_spatial(ix,iy)=f_pb*0.5*&
          !     ( tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)+sigma_pb ))-&
          !     tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)-sigma_pb )) ) + (0.0,0.0)

          !add a small probe with a vortex
          !erre=sqrt((sx-assymx_pb)*(sx-assymx_pb)+(sy-assymy_pb)*(sy-assymy_pb))
          !if(erre .eq. 0.0) then
          !   theta=(0.0,0.0)
          !   probe_spatial(ix,iy)=(0.0,0.0)
          !else
          !   theta=Lz*(  (0.0,-1.0)*&
          !        log( ((sx-assymx_pb)+(0.0,1.0)*(sy-assymy_pb))/erre ) + pi  )
          !   probe_spatial(ix,iy)=f_pb*( exp(1.0)/(2*sigma_pb*sigma_pb) )*&
          !        erre*erre*&
          !        exp(-0.5*(erre*erre)/(sigma_pb*sigma_pb))*&
          !        exp((0.0,1.0)*theta) + (0.0,0.0)
          !end if

       end do
    end do


!!$       ! multiply by k_p, k_pb part
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          pump_spatial(ix,iy)= pump_spatial(ix,iy)*cos(k_p*sx)+(0.0,1.0)*&
               pump_spatial(ix,iy)*sin(k_p*sx)
       end do
    end do

    theta_pb=theta_pb*2.0*pi/360
    kpb_x=k_pb*cos(theta_pb)
    kpb_y=k_pb*sin(theta_pb)
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
!          probe_spatial(ix,iy)=probe_spatial(ix,iy)*cos(k_pb*sx)+(0.0,1.0)*&
!               probe_spatial(ix,iy)*sin(k_pb*sx)

          probe_spatial(ix,iy)=probe_spatial(ix,iy)*&
                               (cos(kpb_x*sx)+(0.0,1.0)*sin(kpb_x*sx))*&
                               (cos(kpb_y*sy)+(0.0,1.0)*sin(kpb_y*sy))
       end do
    end do


!!$        !add a barrier with height "barrier" between bin_x and bfin_x
!!$        do ix=ixbin,ixbfin
!!$           pot_c(ix,:)=barrier
!!$        enddo
!!$

!!$    !add a square difect with height "barrier" between bin_x and bfin_x
!    do ix=ixbin,ixbfin
!       do iy=iybin,iybfin
!          pot_c(ix,iy)=barrier
!       enddo
!    enddo

    !smoothen top hat barrier with hyperbolic tangent
!    do ix=1,Nx
!          sx=-Lx+(ix-1)*ax
!       do iy=1,Ny
!       sy=-Ly+(iy-1)*ay
!          pot_c(ix,iy)=barrier*0.25*&
!               ( tanh((1.0/10.0)*( sqrt(sx*sx)+bx ))-&
!                 tanh((1.0/10.0)*( sqrt(sx*sx)-bx )) )*&
!               ( tanh((1.0/10.0)*( sqrt(sy*sy)+by ))-&
!                 tanh((1.0/10.0)*( sqrt(sy*sy)-by )) )+&
!               + (0.0,0.0)
!        enddo
!   enddo	

    !smoothen top hat barrier with hyperbolic tangent
!    do ix=1,Nx
!       sx=-Lx+(ix-1)*ax
!       do iy=1,Ny
!          sy=-Ly+(iy-1)*ay
!          pot_c(ix,iy)=barrier*0.25*&
!               ( tanh((0.25)*( sqrt(sx*sx)+bx ))-&
!                 tanh((0.25)*( sqrt(sx*sx)-bx )) )*&
!               ( tanh((0.25)*( sqrt(sy*sy)+by ))-&
!                 tanh((0.25)*( sqrt(sy*sy)-by )) )+&
!               + (0.0,0.0)
!        enddo
!   enddo
!!$
    !smoothen top hat barrier with hyperbolic tangent
!    do ix=1,Nx
!       sx=-Lx+(ix-1)*ax
!       do iy=1,Ny
!       sy=-Ly+(iy-1)*ay
!          pot_c(ix,iy)=barrier*0.25*&
!               ( tanh((1.0/10.0)*( sx-bin_x ))-&
!                 tanh((1.0/10.0)*( sx-bfin_x )) )*&
!               ( tanh((1.0/10.0)*( sy-bin_y ))-&
!                 tanh((1.0/10.0)*( sy-bfin_y )) )+&
!               + (0.0,0.0)
!        enddo
!   enddo


      do ix=1,Nx
          sx=-Lx+(ix-1)*ax
          do iy=1,Ny
             sy=-Ly+(iy-1)*ay
              erre=sx*sx+sy*sy
!             pot_c(ix,iy)=0.5*barrier*(tanh((sqrt(sx**2.+sy**2.)+br_x)/(0.1*br_x)) - tanh((sqrt(sx**2.+sy**2.)-br_x)/(0.1*br_x)))
!             pot_c(ix,iy)=barrier*exp(-((sx**2.+sy**2.)/(br_x*br_x)))
             if ( erre .le. br_x*br_x ) pot_c(ix,iy)=barrier
          enddo
      enddo

!pot_x=pot_c

    !smoothen top hat barrier with hyperbolic tangent broken in the middle
!    if (barrier .ne. 0.0) then
!       do ix=1,Nx
!          sx=-Lx+(ix-1)*ax
!          do iy=1,Ny
!             sy=-Ly+(iy-1)*ay
!             if ( ABS(sx) .le. bx .and. ABS(sy) .le. by ) then
!                pot_c(ix,iy)=0.1*&
!                     ( tanh((0.2)*( sqrt(sx*sx)+bx ))-&
!                     tanh((0.2)*( sqrt(sx*sx)-bx )) )*&
!                     ( tanh((0.2)*( sqrt(sy*sy)+by ))-&
!                     tanh((0.2)*( sqrt(sy*sy)-by )) )+&
!                     + (0.0,0.0) + barrier
!             else
!                pot_c(ix,iy)=0.25*&
!                     ( tanh((0.2)*( sqrt(sx*sx)+bx ))-&
!                     tanh((0.2)*( sqrt(sx*sx)-bx )) )*&
!                     ( tanh((0.2)*( sqrt(sy*sy)+by ))-&
!                     tanh((0.2)*( sqrt(sy*sy)-by )) )+&
!                     + (0.0,0.0)
!             endif
!          enddo
!       enddo
!    end if

!!$        !add a circular difect with height "barrier" with a radio of bfin_x lenght
!!$        do ix=1,Nx
!!$           sx=-Lx+(ix-1)*ax
!!$           do iy=1,Ny
!!$              sy=-Ly+(iy-1)*ay
!!$              if ( sqrt( sx*sx+sy*sy ) .lt. ABS( bfin_x ) ) then
!!$                 pot_c(ix,iy)=barrier
!!$              endif
!!$           enddo
!!$        enddo

!!$        !add a random potential to the excitonic part
!!$        do iy=1, Ny
!!$           sy=-Ly+(iy-1)*ay
!!$           do ix=1, Nx
!!$              sx=-Lx+(ix-1)*ax
!!$              CALL RANDOM_NUMBER(rnd)
!!$              pot(ix,iy)=0.1*rnd
!!$           enddo
!!$        enddo


    ! write the pump profile
    open(unit=25, file='pump.dat', status='replace')
    WRITE(25,*)'variables = "x","y","pump"'
    WRITE(25,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(25,*) sx, sy, abs(pump_spatial(ix,iy))
       enddo
    enddo
    close(25)

    ! write the probe profile
    open(unit=25, file='probe.dat', status='replace')
    WRITE(25,*)'variables = "x","y","probe"'
    WRITE(25,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(25,*) sx, sy, abs(probe_spatial(ix,iy))
       enddo
    enddo
    close(25)
    ! write the probe phase
    open(unit=25, file='phase.dat', status='replace')
    WRITE(25,*)'variables = "x","y","phase"'
    WRITE(25,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(25,*) sx, sy, real(theta)
       enddo
    enddo
    close(25)

    ! write the potential profile added to the X part
    open(unit=25, file='pot_x.dat', status='replace')
    WRITE(25,*)'variables = "x","y","pot_x"'
    WRITE(25,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(25,*) sx, sy, pot_x(ix,iy)
       enddo
    enddo
    close(25)

    ! write the potential profile added to the C part
    open(unit=25, file='pot_c.dat', status='replace')
    WRITE(25,*)'variables = "x","y","pot_c"'
    WRITE(25,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(25,*) sx, sy, pot_c(ix,iy)
       enddo
    enddo
    close(25)

    !    call writemom(pump_spatial,100)
    RETURN
  END SUBROUTINE  init

  !-----------------------------------------------------------------!
  ! Subroutine setg to calculate k*k on FFT grid                    !
  !-----------------------------------------------------------------!
  Subroutine setg
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
  ! Subroutine setgk to calculate k on FFT grid                      !
  !-----------------------------------------------------------------!
  Subroutine setgk
    implicit none
    integer :: j,k

    DO j=1,(Ny/2+1)
       DO k=1,(Nx/2+1)
          kvecx(k,j)=pi*(k-1)/Lx
          kvecy(k,j)=pi*(j-1)/Ly
       END DO
    END DO
    DO j=(Ny/2+2),Ny
       DO k=(Nx/2+2),Nx
          kvecx(k,j)=pi*(k-1-Nx)/Lx
          kvecy(k,j)=pi*(j-1-Ny)/Ly
       END DO
    END DO
    DO j=1,(Ny/2+1)
       DO k=(Nx/2+2),Nx
          kvecx(k,j)=pi*(k-1-Nx)/Lx
          kvecy(k,j)=pi*(j-1)/Ly
       END DO
    END DO
    DO j=(Ny/2+2),Ny
       DO k=1,(Nx/2+1)
          kvecx(k,j)=pi*(k-1)/Lx
          kvecy(k,j)=pi*(j-1-Ny)/Ly
       END DO
    END DO

!    call writemom2(REAL(kvecx),4)
!    call writemom2(REAL(kvecy),5)

    kvecx=(0.0,1.0)*kvecx
    kvecy=(0.0,1.0)*kvecy

    return
  end subroutine setgk

  Subroutine dragPOT
    IMPLICIT NONE
    integer :: ix,iy
    
    !fft to momentum space
    transformingF(:,:)=pot_c(:,:)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    transformingF=kvecx*transformingF

    !fft back to real space
    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    dragPOTx=REAL(transformingF)

    !fft to momentum space
    transformingF(:,:)=pot_c(:,:)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    transformingF=kvecy*transformingF

    !fft back to real space
    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    dragPOTy=REAL(transformingF)

    call writepot(dragPOTx,1)
    call writepot(dragPOTy,2)

  endsubroutine dragPOT

  !-----------------------------------------------------------------!
  ! Subroutine derivs(x,y,dydx)                                     !
  !-----------------------------------------------------------------!
  Subroutine derivs(x,y,dydx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    COMPLEX, DIMENSION(:,:,:), INTENT(IN) :: y
    COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: dydx
    INTEGER :: ix,iy
    REAL :: sx, sy

!!$    !generate the energy dependence of the pump
    pump=pump_spatial*exp((0.0,-1.0)*omega_p*x)
!!$    !generate the energy and time dependence of the probe
!!$    probe=probe_spatial*(cos(omega_pb*x)+(0.0,-1.0)*sin(omega_pb*x))*&
!!$         &exp(-0.5*(x-t_init-nene*dxene-5.*sigma_tpb)**2.0/sigma_tpb**2.0)
!!$   !generate the energy dependence of a time independent probe
!!$    if ( x .lt. t_init+nene*dxene ) then
!!$       probe=(0.0,0.0)
!!$    else
       probe=probe_spatial*exp((0.0,-1.0)*omega_pb*x)
!!$    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !right hand side of the equation of motion for the photon population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !photon part
    dydx(:,:,1)= (0.0,-1.0)*y(:,:,2)-(0.0,1.0)*delta*y(:,:,1)-&
         kappa_C*y(:,:,1)+(0.0,-1.0)*(pump(:,:)+probe(:,:))
    !add a potential to the photon part
    dydx(:,:,1)=dydx(:,:,1)+(0.0,1.0)*pot_c(:,:)*y(:,:,1)

    !exciton part
    dydx(:,:,2)= -kappa_X*y(:,:,2)+&
         (0.0,-1.0)*ABS(y(:,:,2))*ABS(y(:,:,2))*y(:,:,2)+ &
         &(0.0,-1.0)*y(:,:,1)

    !add a potential in the exciton part
    dydx(:,:,2)=dydx(:,:,2)+(0.0,1.0)*pot_x(:,:)*y(:,:,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !adding the kinetic energy by means of the FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !fft to momentum space
    transformingF(:,:)=y(:,:,1)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    transformingF=kinetic*transformingF

    !fft back to real space
    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    dydx(:,:,1)=dydx(:,:,1)+transformingF*CMPLX(0.0,-1.0)

  END SUBROUTINE derivs

  !-----------------------------------------------------------------!
  ! Subroutine energytransform                                      !
  !-----------------------------------------------------------------!
  Subroutine energytransform
    IMPLICIT NONE
    INTEGER :: i,ix,iy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! to evaluate integrals in energy-k space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$ Transforming from x to k
!!$ transforming psiE
    do i=1,nene
       transformingF(:,:)=psiE(:,:,i)
       Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
       psiE(:,:,i)=transformingF(:,:)/normDFTIsqrt
    enddo
!!$ transforming psiref
    do i=1,nene
       transformingF(:,:)=psiref(:,:,i)
       Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
       psiref(:,:,i)=transformingF(:,:)/normDFTIsqrt
    enddo

!!$ Preparing the DFT for t-energy transform 
    Status = DftiFreeDescriptor( My_Desc1_Handle )
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, &
         DFTI_COMPLEX, 1, nene )
    Status = DftiCommitDescriptor( My_Desc1_Handle )

    ! FFT from t to e of the psiE
    do ix=1,Nx
       do iy=1,Ny
          do i=1,nene
             transfE(i)=psiE(ix,iy,i)
          enddo
          Status = DftiComputeForward( My_Desc1_Handle, transfE )
          do i=1,nene
             psiE(ix,iy,i)=transfE(i)/sqrt(REAL(nene))
          enddo
       enddo
    enddo

    ! FFT from t to e of the psiref
    do ix=1,Nx
       do iy=1,Ny
          do i=1,nene
             transfE(i)=psiref(ix,iy,i)
          enddo
          Status = DftiComputeForward( My_Desc1_Handle, transfE )
          do i=1,nene
             psiref(ix,iy,i)=transfE(i)/sqrt(REAL(nene))
          enddo
       enddo
    enddo

!!$ writing psiE_energy and it's integrals
    CALL writeene(psiref(:,:,:))
    CALL intene(psiref)
    CALL system('mv opo_ene.dat opo_eneTOT.dat')
    CALL system('mv int_ene.dat int_eneTOT.dat')

!!$ writing psiref_energy and it's integrals
!    CALL writeene(psiref(:,:,:))
!    CALL writeeneim(psiref(:,:,:))
!    CALL system('mv opo_ene.dat opo_eneREF.dat')

!    auxpsi=ABS(psiE)-ABS(psiref)
!    CALL writeene(auxpsi(:,:,:))
!    CALL system('mv opo_ene.dat opo_eneDIF.dat')
!    CALL energy_filter(auxpsi(:,:,:))
!    CALL energy_filter2(auxpsi(:,:,:))
!    CALL energy_filter3(auxpsi(:,:,:))
    CALL writeene(auxpsi(:,:,:))
    CALL intene(auxpsi)
    CALL intint(auxpsi)

!    CALL energy_filter(psiE(:,:,:))
    CALL energy_filter(psiref(:,:,:))
!    CALL energy_filter2(psiE(:,:,:))
!    CALL energy_filter2(psiref(:,:,:))
!    CALL energy_filter3(psiE(:,:,:))
!    CALL energy_filter3(psiref(:,:,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! going back to r-t space to write psiE and psiref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FFT from e to t
    do ix=1,Nx
       do iy=1,Ny
          do i=1,nene
             transfE(i)=psiE(ix,iy,i)
          enddo
          Status = DftiComputeBackward( My_Desc1_Handle, transfE )
          do i=1,nene
             psiE(ix,iy,i)=transfE(i)/sqrt(REAL(nene))
          enddo
       enddo
    enddo

    ! FFT from e to t
    do ix=1,Nx
       do iy=1,Ny
          do i=1,nene
             transfE(i)=psiref(ix,iy,i)
          enddo
          Status = DftiComputeBackward( My_Desc1_Handle, transfE )
          do i=1,nene
             psiref(ix,iy,i)=transfE(i)/sqrt(REAL(nene))
          enddo
       enddo
    enddo

    Status = DftiFreeDescriptor( My_Desc1_Handle )
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, &
         DFTI_COMPLEX, 2, npointxy )
    Status = DftiCommitDescriptor( My_Desc1_Handle )

    ! FFT from k to x
    do i=1,nene
       transformingF(:,:)=psiE(:,:,i)
       Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
       psiE(:,:,i)=transformingF(:,:)/normDFTIsqrt
    enddo

    ! FFT from k to x
    do i=1,nene
       transformingF(:,:)=psiref(:,:,i)
       Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
       psiref(:,:,i)=transformingF(:,:)/normDFTIsqrt
    enddo

    do i=10,nene,10
      CALL writespace_filt_total(psiE(:,:,i),i/10)
      CALL writespace_filt_ref(psiref(:,:,i),i/10)
    enddo

    do i=10,nene,10
       psiE(:,:,i)=(0.0,0.0)+ABS(psiE(:,:,i))-ABS(psiref(:,:,i))
       CALL writespace_filt_ene(psiE(:,:,i),i/10)
    enddo

  END SUBROUTINE energytransform

  !-----------------------------------------------------------------!
  ! Subroutine energy filter                                        !
  !-----------------------------------------------------------------!
  Subroutine energy_filter(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(INOUT) :: psi(:,:,:)
    INTEGER :: ie,kx,ky
    REAL :: dene,ene,mom_x,mom_y

    dene=2*pi/(nene*dxene)
    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       IF ( -ene .LT. ene_cent-ene_cut .or. -ene .GT. ene_cent+ene_cut) THEN
          psi(:,:,ie)=0.0
       ENDIF
    enddo
    do ie=1,(nene)/2+1
       ene=dene*ie
       IF ( -ene .LT. ene_cent-ene_cut .or. -ene .GT. ene_cent+ene_cut) THEN
          psi(:,:,ie)=0.0
       END IF
    enddo

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if

       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if

       end do
    end do

  END SUBROUTINE energy_filter

  !-----------------------------------------------------------------!
  ! Subroutine energy filter2                                        !
  !-----------------------------------------------------------------!
  Subroutine energy_filter2(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(INOUT) :: psi(:,:,:)
    INTEGER :: ie,kx,ky
    REAL :: dene,ene,mom_x,mom_y

    dene=2*pi/(nene*dxene)
    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       IF ( -ene .LT. ene_cent-ene_cut .or. -ene .GT. ene_cent+ene_cut) THEN
          psi(:,:,ie)=0.0
       ENDIF
    enddo
    do ie=1,(nene)/2+1
       ene=dene*ie
       IF ( -ene .LT. ene_cent-ene_cut .or. -ene .GT. ene_cent+ene_cut) THEN
          psi(:,:,ie)=0.0
       END IF
    enddo

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent2)*(mom_x-mom_cent2)+mom_y*mom_y) .GT. mom_cut2 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
          end if
       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent2)*(mom_x-mom_cent2)+mom_y*mom_y) .GT. mom_cut2 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
          end if
       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent2)*(mom_x-mom_cent2)+mom_y*mom_y) .GT. mom_cut2 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
          end if
       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .GT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent2)*(mom_x-mom_cent2)+mom_y*mom_y) .GT. mom_cut2 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
          end if
       end do
    end do

  END SUBROUTINE energy_filter2

  !-----------------------------------------------------------------!
  ! Subroutine energy filter3                                       !
  !-----------------------------------------------------------------!
  Subroutine energy_filter3(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(INOUT) :: psi(:,:,:)
    INTEGER :: ie,kx,ky
    REAL :: dene,ene,mom_x,mom_y

    dene=2*pi/(nene*dxene)
    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       IF ( -ene .GT. ene_cent-ene_cut .and. -ene .LT. ene_cent+ene_cut) THEN
!          write(*,*) '1',-ene
          psi(:,:,ie)=0.0
       ENDIF
    enddo
    do ie=1,(nene)/2+1
       ene=dene*ie
       IF ( -ene .GT. ene_cent-ene_cut .and. -ene .LT. ene_cent+ene_cut) THEN
!          write(*,*) '2',-ene
          psi(:,:,ie)=0.0
       END IF
    enddo

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly
!          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .LT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)) .LT. mom_cut1 ) then
!             write(*,*) 'mom',mom_x,mom_y
             psi(kx,ky,:) = (0.0,0.0)
          end if
       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly
!          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .LT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)) .LT. mom_cut1 ) then
!             write(*,*) 'mom',mom_x,mom_y
             psi(kx,ky,:) = (0.0,0.0)
          end if
       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly
!          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .LT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)) .LT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly
!          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .LT. mom_cut1 ) then
          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)) .LT. mom_cut1 ) then
             psi(kx,ky,:) = (0.0,0.0)
          end if
       end do
    end do

  END SUBROUTINE energy_filter3

  !-----------------------------------------------------------------!
  ! Subroutine writeAV                                              !
  !-----------------------------------------------------------------!
  Subroutine writeAV(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy,nc,vsound,vfluid,nclinear
    INTEGER :: ix,iy

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_AV"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) sx, sy, (ABS(psi(ix,iy))**2.)
       end do
    end do
    close(22)

    nc=SUM(ABS(psi)**2.)*norm
    nc=nc*(2*Lx*2*Ly)/(pi*(sigma_p*sigma_p-br_x*br_x))
    write(400,*) i,nc
    vsound=sqrt(nc)*thetasin*thetasin*2.81
    vfluid=k_p*14.5*(1-thetasin*thetasin)
    write(401,*) i,vsound,vfluid
    write(402,*) i,vsound/vfluid

    nclinear=(SUM(ABS(psi)**2.)*norm)/(pi*(sigma_p*sigma_p-br_x*br_x))
    write(500,*) i,nclinear

  END SUBROUTINE writeAV

  !-----------------------------------------------------------------!
  ! Subroutine writeAVvel                                           !
  !-----------------------------------------------------------------!
  Subroutine writeAVvel(psix,psiy,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psix(:,:),psiy(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy,avvel
    INTEGER :: ix,iy
    CHARACTER (LEN=80):: string

    string="(2X,F10.3,2X,F10.3,2X,E10.3,2X,E10.3)"

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_AVvel"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","vx","vy"'
    WRITE(22,*)'zone I=',ny/4,',J=',nx/4,',F=POINT'
    do iy=1, Ny,4
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx,4
          sx=-Lx+(ix-1)*ax
          write(22,string) sx, sy, psix(ix,iy), psiy(ix,iy)
       end do
    end do
    close(22)

!!$    write(label,FMT="(i3)") i
!!$    open(unit=22, file="opo_vel"//trim(adjustl(label))//".dat", status='replace')
!!$    WRITE(22,*)'variables = "x","y","vx","vy"'
!!$    WRITE(22,*)'zone I=',ny/4,',J=',nx/4,',F=POINT'
!!$    do iy=1, Ny,4
!!$       sy=-Ly+(iy-1)*ay
!!$       do ix=1, Nx,4
!!$          sx=-Lx+(ix-1)*ax
!!$          if (sqrt(sx**2.+sy**2.) .lt. sigma_p-30.) then
!!$          write(22,string) sx, sy, psix(ix,iy), psiy(ix,iy)
!!$       else
!!$          write(22,string) sx, sy, 0.0, 0.0
!!$       endif
!!$       end do
!!$    end do
!!$    close(22)

  END SUBROUTINE writeAVvel

  !-----------------------------------------------------------------!
  ! Subroutine xderiv                                               !
  !-----------------------------------------------------------------!
  Subroutine xderiv(psi,velx)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    REAL, INTENT(OUT) ::velx(:,:)
    COMPLEX :: psimf(Nx,Ny),aux(Nx,Ny)
    INTEGER :: ix,kx,ky
    REAL :: sx,mom_x,mom_y,max,kxmax

    max=0.0
    kxmax=0.0

    transformingF(:,:)=psi(:,:)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    do ky=Ny/2+2, Ny
       mom_y=pi*(ky-1-Ny)/Ly
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          if ( ABS(transformingF(kx,ky))**2. .gt. max ) then
             max=ABS(transformingF(kx,ky))**2.
             kxmax=mom_x
          endif
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          if ( ABS(transformingF(kx,ky))**2. .gt. max ) then
             max=ABS(transformingF(kx,ky))**2.
             kxmax=mom_x
          endif
       end do
    end do
    do ky=1, Ny/2+1
       mom_y=pi*(ky-1)/Ly
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          if ( ABS(transformingF(kx,ky))**2. .gt. max ) then
             max=ABS(transformingF(kx,ky))**2.
             kxmax=mom_x
          endif
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          if ( ABS(transformingF(kx,ky))**2. .gt. max ) then
             max=ABS(transformingF(kx,ky))**2.
             kxmax=mom_x
          endif
       end do
    end do

    do ix=1,Nx
       sx=-Lx+(ix-1)*ax
       psimf(ix,:)=psi(ix,:)*exp(-(0.0,1.0)*kxmax*sx)
    enddo

    !fft to momentum space
    transformingF(:,:)=psimf(:,:)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    transformingF=kvecx*transformingF

    !fft back to real space
    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    aux=CONJG(psimf)*transformingF
!    call writespace2(AIMAG(aux),4)

    velx=AIMAG(aux)/(ABS(psi)**2.)

  END SUBROUTINE xderiv

  !-----------------------------------------------------------------!
  ! Subroutine yderiv                                               !
  !-----------------------------------------------------------------!
  Subroutine yderiv(psi,vely)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    REAL, INTENT(OUT) ::vely(:,:)
    COMPLEX :: aux(Nx,Ny)

    !fft to momentum space
    transformingF(:,:)=psi(:,:)
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    transformingF=kvecy*transformingF

    !fft back to real space
    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    aux=CONJG(psi)*transformingF
!    call writespace2(AIMAG(aux),8)

    vely=AIMAG(aux)/(ABS(psi)**2.)

  END SUBROUTINE yderiv

  !-----------------------------------------------------------------!
  ! Subroutine dragforcex                                           !
  !-----------------------------------------------------------------!
  Subroutine dragforcex(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: aux(Nx,Ny)
    REAL :: force

    aux=(ABS(psi)**2.)*dragPOTx
!    force=SUM(aux)/SUM(ABS(psi)**2.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    force=-SUM(aux)/SUM(ABS(psi)**2.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(222,*) i,force

  endsubroutine dragforcex

  !-----------------------------------------------------------------!
  ! Subroutine dragforcey                                           !
  !-----------------------------------------------------------------!
  Subroutine dragforcey(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: aux(Nx,Ny)
    REAL :: force

    aux=(ABS(psi)**2.)*dragPOTy
    force=SUM(aux)/SUM(ABS(psi)**2.)
    write(333,*) i,force

  endsubroutine dragforcey

  !-----------------------------------------------------------------!
  ! Subroutine writecomplex                                         !
  !-----------------------------------------------------------------!
  Subroutine writecomplex(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy
    INTEGER :: ix,iy

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_COMPLEX"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) sx, sy, psi(ix,iy)
       end do
    end do
    close(22)
  END SUBROUTINE writecomplex

  !-----------------------------------------------------------------!
  ! Subroutine writeAVmom                                           !
  !-----------------------------------------------------------------!
  Subroutine writeAVmom(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: mom_x, mom_y
    INTEGER :: kx, ky

    write(label,FMT="(i3)") i
    open(unit=23, file="opo_AVmom"//trim(adjustl(label))//".dat", status='replace')
    WRITE(23,*)'variables = "x","y","opo_mom"'
    WRITE(23,*)'zone I=',ny,',J=',nx,',F=POINT'

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          write(23,*) mom_x, mom_y, abs(psi(kx,ky))**2.
!          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

!          write(23,*) mom_x, mom_y, abs(psi(kx,ky))**2.
          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
    end do

   do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly

!          write(23,*) mom_x, mom_y, abs(psi(kx,ky))**2.
          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly

!          write(23,*) mom_x, mom_y, abs(psi(kx,ky))**2.
          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
    end do
    close(23)

  END SUBROUTINE writeAVmom

  !-----------------------------------------------------------------!
  ! Subroutine transfmom                                            !
  !-----------------------------------------------------------------!
  Subroutine transfmom(psiin,psiout)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psiin(:,:)
    REAL, INTENT(OUT) :: psiout(:,:)

    !fft to momentum space
    transformingf=psiin
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    psiout=ABS(transformingF)

  END SUBROUTINE transfmom

  !-----------------------------------------------------------------!
  ! Subroutine writemom                                             !
  !-----------------------------------------------------------------!
  Subroutine writemom(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    real :: mom_x, mom_y
    integer :: kx, ky
    INTEGER, EXTERNAL :: assert_eq

    !fft to momentum space
    transformingf=psi
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    !write the wave function in momentum space psi
    !(n.b. one has to rearrange)
    write(label,FMT="(i3)") i
    open(unit=23, file="opo_mom"//trim(adjustl(label))//".dat", status='replace')
    WRITE(23,*)'variables = "x","y","opo_mom"'
    WRITE(23,*)'zone I=',ny,',J=',nx,',F=POINT'

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
    end do
    close(23)

  END SUBROUTINE writemom

  !-----------------------------------------------------------------!
  ! Subroutine writemom2                                            !
  !-----------------------------------------------------------------!
  Subroutine writemom2(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    real :: mom_x, mom_y
    integer :: kx, ky
    INTEGER, EXTERNAL :: assert_eq

    !write the wave function in momentum space psi
    !(n.b. one has to rearrange)
    write(label,FMT="(i3)") i
    open(unit=23, file="opo_mom"//trim(adjustl(label))//".dat", status='replace')
    WRITE(23,*)'variables = "x","y","opo_mom"'
    WRITE(23,*)'zone I=',ny,',J=',nx,',F=POINT'

    do ky=Ny/2+2, Ny
       mom_y=pi*(ky-1-Ny)/Ly
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx

          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx

          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
    end do
    do ky=1, Ny/2+1
       mom_y=pi*(ky-1)/Ly
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx

          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx

          write(23,*) mom_x, mom_y, psi(kx,ky)

       end do
    end do
    close(23)

  END SUBROUTINE writemom2

  !-----------------------------------------------------------------!
  ! Subroutine writemom_filter                                      !
  !-----------------------------------------------------------------!
  Subroutine writemom_filter(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    real :: mom_x, mom_y, kcutx, kcuty
    integer :: ix, iy, sx, sy
    integer :: kx, ky
    INTEGER, EXTERNAL :: assert_eq

    !fft to momentum space
    transformingf=psi
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt
    !write the wave function in momentum space psi
    !(n.b. one has to rearrange)

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .le. mom_cut1 ) then
             transformingF(kx,ky) = (0.0,0.0)
          end if

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .le. mom_cut1 ) then
             transformingF(kx,ky) = (0.0,0.0)
          end if

       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .le. mom_cut1 ) then
             transformingF(kx,ky) = (0.0,0.0)
          end if

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly

          if( sqrt((mom_x-mom_cent1)*(mom_x-mom_cent1)+mom_y*mom_y) .le. mom_cut1 ) then
             transformingF(kx,ky) = (0.0,0.0)
          end if

       end do
    end do
    close(23)

    Status = DftiComputeBackward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    write(label,FMT="(i3)") i
    !write the wave function in real space after momentum filtering
    !(n.b. one has to rearrange)
    open(unit=24, file="opo_spcfilter"//trim(adjustl(label))//".dat", status='replace')
    WRITE(24,*)'variables = "x","y","opo_spcfilter"'
    WRITE(24,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(24,*) sx, sy, ABS(transformingF(ix,iy))*ABS(transformingF(ix,iy))
       end do
    end do
    close(24)

  END SUBROUTINE writemom_filter

  !-----------------------------------------------------------------!
  ! Subroutine writeene                                             !
  !-----------------------------------------------------------------!
  Subroutine writeene(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:,:)
    INTEGER :: ie,kx
    REAL :: mom_x,dene,ene,aux

    dene=2.*pi/(nene*dxene)

    open(unit=23, file="opo_ene.dat", status='replace')
    WRITE(23,*)'variables = "kx","e","opo_ene"'
    WRITE(23,*)'zone I=',nx,',J=',nene,',F=POINT'

    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx

          write(23,*) mom_x, -ene, SUM(ABS(psi(kx,:,ie)))**2.0

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx

          write(23,*) mom_x, -ene, SUM(ABS(psi(kx,:,ie)))**2.0

       end do
    enddo

    do ie=1,(nene)/2+1
       ene=dene*ie
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx

          write(23,*) mom_x, -ene, SUM(ABS(psi(kx,:,ie)))**2.0

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx

          write(23,*) mom_x, -ene, SUM(ABS(psi(kx,:,ie)))**2.0

       end do
    enddo

    close(23)

  END SUBROUTINE writeene

  !-----------------------------------------------------------------!
  ! Subroutine writeeneim                                           !
  !-----------------------------------------------------------------!
  Subroutine writeeneim(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:,:)
    INTEGER :: ie,kx
    REAL :: mom_x,dene,ene,aux

    dene=2.*pi/(nene*dxene)
    open(unit=23, file="opo_eneIM.dat", status='replace')
    WRITE(23,*)'variables = "kx","e","opo_ene"'
    WRITE(23,*)'zone I=',nx,',J=',nene,',F=POINT'

    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          write(23,*) mom_x, -ene, SUM(AIMAG(psi(kx,:,ie)))
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          write(23,*) mom_x, -ene, SUM(AIMAG(psi(kx,:,ie)))
       end do
    enddo
    do ie=1,(nene)/2+1
       ene=dene*ie
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          write(23,*) mom_x, -ene, SUM(AIMAG(psi(kx,:,ie)))
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          write(23,*) mom_x, -ene, SUM(AIMAG(psi(kx,:,ie)))
       end do
    enddo
    close(23)

  END SUBROUTINE writeeneim

  !-----------------------------------------------------------------!
  ! Subroutine intene                                               !
  !-----------------------------------------------------------------!
  Subroutine intene(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:,:)
    INTEGER :: ie,kx,ky
    REAL :: int_ene,dene,ene

    dene=2*pi/(nene*dxene)

    open(unit=23, file="int_ene.dat", status='replace')
    do ie=(nene)/2+2,nene
       ene=dene*(ie-(nene))
       int_ene=SUM(SUM(ABS(psi(:,:,ie)),1),1)
       write(23,*) -ene, int_ene
    enddo

    do ie=1,(nene)/2+1
       ene=dene*ie
       int_ene=SUM(SUM(ABS(psi(:,:,ie)),1),1)
       write(23,*) -ene, int_ene
    enddo
    close(23)

  END SUBROUTINE intene

  !-----------------------------------------------------------------!
  ! Subroutine intint                                               !
  !-----------------------------------------------------------------!
  Subroutine intint(psi)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:,:)
    INTEGER :: ie
    REAL :: int_ene,ene

    int_ene=0.
    open(unit=23, file="int_int.dat", status='replace')
    do ie=1,nene
       int_ene=int_ene+SUM(SUM(ABS(psi(:,:,ie)),1),1)
    enddo
    write(23,*) int_ene
    close(23)

  END SUBROUTINE intint


  !-----------------------------------------------------------------!
  ! Subroutine writespace_filt_ene                                  !
  !-----------------------------------------------------------------!
  Subroutine writespace_filt_ene(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    INTEGER :: ix,iy,sx,sy
    REAL :: trans

    write(label,FMT="(i3)") i
    open(unit=22, file="filt_ene_spc"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) INT(sx), INT(sy), ABS(psi(ix,iy))*ABS(psi(ix,iy))
       end do
    end do
    close(22)

    trans=0.0
    do ix=ixbfin,Nx
       trans=trans+ABS(SUM(ABS(psi(ix,:))))*ABS(SUM(ABS(psi(ix,:))))
    enddo
    write(200,*) i,trans

    trans=0.0
    do ix=1,ixbfin-1
       trans=trans+ABS(SUM(ABS(psi(ix,:))))*ABS(SUM(ABS(psi(ix,:))))
    enddo
    write(300,*) i,trans

  END SUBROUTINE writespace_filt_ene

  !-----------------------------------------------------------------!
  ! Subroutine writespace_filt_total                                !
  !-----------------------------------------------------------------!
  Subroutine writespace_filt_total(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    INTEGER :: ix,iy,sx,sy
    REAL :: trans

    write(label,FMT="(i3)") i
    open(unit=22, file="filt_pulse"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) INT(sx), INT(sy), ABS(psi(ix,iy))*ABS(psi(ix,iy))
       end do
    end do
    close(22)

  END SUBROUTINE writespace_filt_total

  !-----------------------------------------------------------------!
  ! Subroutine writespace_filt_ref                                  !
  !-----------------------------------------------------------------!
  Subroutine writespace_filt_ref(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    INTEGER :: ix,iy,sx,sy
    REAL :: trans

    write(label,FMT="(i3)") i
    open(unit=22, file="filt_steady"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) INT(sx), INT(sy), ABS(psi(ix,iy))*ABS(psi(ix,iy))
       end do
    end do
    close(22)

  END SUBROUTINE writespace_filt_ref

  !-----------------------------------------------------------------!
  ! Subroutine writemom_filt_ene                                    !
  !-----------------------------------------------------------------!
  Subroutine writemom_filt_ene(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    real :: mom_x, mom_y
    integer :: kx, ky
    INTEGER, EXTERNAL :: assert_eq

    !fft to momentum space
    transformingf=psi
    Status = DftiComputeForward( My_Desc1_Handle, DFTItransformingF )
    transformingF=transformingF/normDFTIsqrt

    !write the wave function in momentum space psi
    !(n.b. one has to rearrange)
    write(label,FMT="(i3)") i
    open(unit=23, file="opo_filt_ene_mom"//trim(adjustl(label))//".dat", status='replace')
    WRITE(23,*)'variables = "x","y","opo_mom"'
    WRITE(23,*)'zone I=',ny,',J=',nx,',F=POINT'

    do ky=Ny/2+2, Ny
       !region(4)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
       !region(3)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1-Ny)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
    end do

    do ky=1, Ny/2+1
       !region(2)
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          mom_y=pi*(ky-1)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
       !region(1)
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          mom_y=pi*(ky-1)/Ly

          write(23,*) mom_x, mom_y, abs(transformingF(kx,ky))*abs(transformingF(kx,ky))

       end do
    end do
    close(23)

  END SUBROUTINE writemom_filt_ene

  !-----------------------------------------------------------------!
  ! Subroutine  save_a_step                                         !
  !-----------------------------------------------------------------!
  SUBROUTINE save_a_step(y,kount)
    IMPLICIT NONE
    COMPLEX, DIMENSION(:,:,:), INTENT(IN) :: y
    INTEGER, INTENT(IN) :: kount
    !      CALL writecomplex(y(:,:,1),kount)
          CALL writespace(y(:,:,1),kount)
    !      CALL writemom(y(:,:,1),kount)
    !      CALL writemom_filter(y(:,:,1),kount)

    !write(111,*) kount,sum(ABS(y(:,:,1)))

  END SUBROUTINE save_a_step

  !-----------------------------------------------------------------!
  ! Subroutine writespace                                           !
  !-----------------------------------------------------------------!
  Subroutine writespace(psi,i)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy
    INTEGER :: ix,iy

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_spc"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) sx, sy, ABS(psi(ix,iy))*ABS(psi(ix,iy))
       end do
    end do
    close(22)
  END SUBROUTINE writespace

  !-----------------------------------------------------------------!
  ! Subroutine writespace2                                          !
  !-----------------------------------------------------------------!
  Subroutine writespace2(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy
    INTEGER :: ix,iy

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_Dspc"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) sx, sy, psi(ix,iy)
       end do
    end do
    close(22)
  END SUBROUTINE writespace2

  Subroutine writepot(psi,i)
    IMPLICIT NONE
    REAL, INTENT(IN) :: psi(:,:)
    INTEGER, INTENT(IN) :: i
    REAL :: sx,sy
    INTEGER :: ix,iy

    write(label,FMT="(i3)") i
    open(unit=22, file="opo_Dpot"//trim(adjustl(label))//".dat", status='replace')
    WRITE(22,*)'variables = "x","y","opo_spc"'
    WRITE(22,*)'zone I=',ny,',J=',nx,',F=POINT'
    do iy=1, Ny
       sy=-Ly+(iy-1)*ay
       do ix=1, Nx
          sx=-Lx+(ix-1)*ax
          write(22,*) sx, sy, psi(ix,iy)
       end do
    end do
    close(22)
  END SUBROUTINE writepot

end Module mod_subroutines

!-----------------------------------------------------------------!
! Function assert_eq(n1,n2,n3,string)                             !
!-----------------------------------------------------------------!
FUNCTION assert_eq(n1,n2,n3,string)
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
