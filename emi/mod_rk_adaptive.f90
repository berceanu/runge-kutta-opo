MODULE mod_rk_adaptive
  USE MKL_DFTI
  USE global
  USE mod_subroutines
  IMPLICIT NONE
  SAVE

CONTAINS
  !-----------------------------------------------------------------!
  ! Subroutine odeint                                               !
  !-----------------------------------------------------------------!
  SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin)
    IMPLICIT NONE
    COMPLEX, DIMENSION(:,:,:), INTENT(INOUT) :: ystart
    REAL, DIMENSION(Nx,Ny) :: auxreal,auxmom
    REAL, DIMENSION(Nx,Ny) :: psireal,psimom
    REAL, DIMENSION(Nx,Ny) :: phase,auxvel_x,auxvel_y
    REAL, DIMENSION(Nx,Ny) :: velocity_x,velocity_y
    REAL, INTENT(IN) :: x1,x2,eps,h1,hmin
    INTEGER, EXTERNAL :: assert_eq
    REAL, PARAMETER :: TINY=1.0e-30
    INTEGER, PARAMETER :: MAXSTP=1000000000
!Runge-Kutta driver with adaptive step size control.Integrate the array
!of starting values ystart from x1 to x2 with accuracy eps storing
!intermediate results in the module variables in ode_path. h1 should be
!set as a guessed first stepsize, hmin as the minimum allowed stepsize
!(can be zero). On output ystart is replaced by values at the end of the
!integration interval. derivs is the user-supplied subroutine for
!calculating the right-hand-side derivative, while rkqs is the name of
!he stepper routine to be used.
    INTEGER :: nstp
    REAL :: h,hdid,hnext,x,xsav,xene
    COMPLEX, DIMENSION(size(ystart,1),size(ystart,2),size(ystart,3)) :: dydx,y,yscal
    REAL :: dphase
    Complex :: pom_cond
    REAL, External :: findfermpart,findcoopernum
    COMPLEX, External :: findfermcond
        
    psireal=0.0
    psimom=0.0
    auxreal=0.0
    auxmom=0.0
    velocity_x=0.0
    velocity_y=0.0
    phase=0.0
    auxvel_x=0.0
    auxvel_y=0.0
    x=x1
    h=sign(h1,x2-x1)
    nok=0
    nbad=0
    kount=0
    y(:,:,:)=ystart(:,:,:)

!nullify(xp,yp)
!Pointers nullified here, but memory not deallocated. If odeint is
!called multiple times, calling program should deallocate xp and yp
!between calls.

    xsav=x-2.0*dxsav
    xene=x-2.0*dxene

    write(100,*) abs(x-xsav) , abs(dxsav)
    do nstp=1,MAXSTP 
       !Take at most MAXSTP steps.
       !write(*,*) "stop1"
       !write(*,*) x, y(1),y(2)

       call derivs(x,y,dydx)
       yscal(:,:,:)=abs(y(:,:,:))+abs(h*dydx(:,:,:))+TINY

       !Scaling used to monitor accuracy. This general purpose choice can be
       !modified if need be.

       if ( (abs(x-xene) .gt. abs(dxene)) ) then
          !archivio la funzione d'onda in spazio reale per calcolare lo spettro
          if ( x .gt. t_init .and. x .lt. t_init+nene*dxene ) then
             psiref(:,:,ti)=y(:,:,1)
             ti=ti+1
          endif
          if ( x .gt. t_init+nene*dxene ) then
             psiE(:,:,tti)=y(:,:,1)
             tti=tti+1
          endif
          if ( x .gt. t_init ) then
             if ( (abs(x-xene) .gt. abs(dxene/20.)) ) then
                auxreal=ABS(y(:,:,1))
                psireal=(psireal*(tavs-1)+auxreal)/REAL(tavs)
                tavs=tavs+1

                call transfmom(y(:,:,1),auxmom(:,:))
                psimom=(psimom*(tavm-1)+auxmom)/REAL(tavm)
                tavm=tavm+1

                call xderiv(y(:,:,1),auxvel_x(:,:))
                call yderiv(y(:,:,1),auxvel_y(:,:))
                velocity_x=(velocity_x*(tavp-1)+auxvel_x)/REAL(tavp)
                velocity_y=(velocity_y*(tavp-1)+auxvel_y)/REAL(tavp)
!                velocity_x=auxvel_x
!                velocity_y=auxvel_y
                tavp=tavp+1
             endif
          endif
          xene=x
       endif
       write(100,*) abs(x-xsav) , kount

       if ( (abs(x-xsav) .gt. abs(dxsav)) ) then 
          !Store intermedia e results.
          kount=kount+1
          xsav=x
          if (x .gt. t_init) then
             call save_a_step(y,kount)
             call writeAV(psireal(:,:),kount)
!             call writeAV(ABS(y(:,:,1)),kount)
!             call writeAVmom(psimom(:,:),kount)
             call writeAVvel(velocity_x(:,:),velocity_y(:,:),kount)
             call dragforcex(psireal(:,:),kount)
!             call dragforcey(psireal(:,:),kount)
           endif
          write(100,*) 
       endif

       if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x 
       !If stepsize can overshoot,decrease.
       call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext)
       if (hdid == h) then
          nok=nok+1
       else
          nbad=nbad+1
       end if
       if ((x-x2)*(x2-x1) >= 0.0) then 
          !Are we done?
          ystart(:,:,:)=y(:,:,:)
          call save_a_step(y,kount) 
          !Save final step.
          RETURN 
          !Normal exit.
       end if
       if (abs(hnext) < hmin) write(*,*) "step smaller than min in odeint"
       h=hnext
    end do
    write(*,*) "too many steps in odeint"

  END SUBROUTINE odeint

  !-----------------------------------------------------------------!
  ! Subroutine  rkqs                                                !
  !-----------------------------------------------------------------!
  SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)
    IMPLICIT NONE
    COMPLEX, DIMENSION(:,:,:), INTENT(INOUT) :: y
    COMPLEX, DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal
    REAL, INTENT(INOUT) :: x
    REAL, INTENT(IN) :: htry,eps
    REAL, INTENT(OUT) :: hdid,hnext
    INTEGER, EXTERNAL :: assert_eq

!Fifth order Runge-Kutta step with monitoring of local runca ion error
!o ensure accuracy and adjus stepsize.Input are he dependent variable
!vector y and its derivative dydx at the star ing value of the
!independent variable x Also input are he s epsize to be attempted htry
!he required accuracy eps and the vector yscal agains which the error
!is scaled.y dydx and yscal areallof hesamelengh.Onoupu,y and x are
!replaced by their new values,hdid is he s epsize ha was ac ually
!accomplished,and hnext is he estimated nex stepsize.derivs is the
!user-supplied subroutine ha computes the right-hand-side deriva ives.
    INTEGER :: ndum
    REAL :: errmax,h,htemp,xnew
    COMPLEX, DIMENSION(size(y,1),size(y,2),size(y,3)) :: yerr,ytemp
    REAL, PARAMETER :: SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89e-4
    !The value ERRCON equals (5/SAFETY)**(1/PGROW),see use below.
        
    ndum=assert_eq(size(y),size(dydx),size(yscal), "rkqs") 
    h=htry 
    !Set step size to the initial trial value. 
    do
       call rkck(y,dydx,x,h,ytemp,yerr) 
       !Take a step.
       errmax=maxval(abs(yerr(:,:,:)/yscal(:,:,:)))/eps 
       !Evaluate accuracy.
       if (errmax <= 1.0) exit 
       !Step succeeded.
       htemp=SAFETY*h*(errmax**PSHRNK) 
       !Truncation error too large, reduce stepsize.
       h=sign(max(abs(htemp),0.1*abs(h)),h) 
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
       hnext=5.0*h
    end if
    hdid=h  
    x=x+h   
    y(:,:,:)=ytemp(:,:,:)  
  END SUBROUTINE rkqs

  !-----------------------------------------------------------------!
  ! Subroutine rkck                                                 !
  !-----------------------------------------------------------------!
  SUBROUTINE rkck(y,dydx,x,h,yout,yerr)
    IMPLICIT NONE
    COMPLEX, DIMENSION(:,:,:), INTENT(IN) :: y,dydx
    REAL, INTENT(IN) :: x,h
    COMPLEX, DIMENSION(:,:,:), INTENT(OUT) :: yout,yerr
    INTEGER, EXTERNAL :: assert_eq

!Given values for N variables y and their derivatives dydx known at x
!use the fifth order Cash-Karp Runge-Kutta method to advance the
!solution over an interval h and return the incremented variables as
!yout Also return an estimate of he local truncation error in yout
!using the embedded fourth order method. The user supplies the
!subroutine derivs(x,y,dydx),which returns derivatives dydx at x
    INTEGER :: ndum
    COMPLEX, DIMENSION(size(y,1),size(y,2),size(y,3)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    REAL, PARAMETER :: A2=0.2,A3=0.3,A4=0.6,A5=1.0,&
         A6=0.875,B21=0.2,B31=3.0/40.0,B32=9.0/40.0,&
         B41=0.3,B42=-0.9,B43=1.2,B51=-11.0/54.0,&
         B52=2.5,B53=-70.0/27.0,B54=35.0/27.0,&
         B61=1631.0/55296.0,B62=175.0/512.0,&
         B63=575.0/13824.0,B64=44275.0/110592.0,&
         B65=253.0/4096.0,C1=37.0/378.0,&
         C3=250.0/621.0,C4=125.0/594.0,&
         C6=512.0/1771.0,DC1=C1-2825.0/27648.0,&
         DC3=C3-18575.0/48384.0,DC4=C4-13525.0/55296.0,&
         DC5=-277.0/14336.0,DC6=C6-0.25

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

end MODULE mod_rk_adaptive
