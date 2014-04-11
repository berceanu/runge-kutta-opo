MODULE rk_adaptive
  USE global    
  IMPLICIT NONE
  SAVE

CONTAINS

  Subroutine derivs(x,y,dydx)

    USE FFTW3
    IMPLICIT NONE

    ! fft stuff
    ! forward means real space to momentum space, backward the opposite
    type(C_PTR) :: plan_forward, plan_backward
    complex(C_DOUBLE_COMPLEX), pointer :: in_forward(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: out_forward(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: in_backward(:,:)
    complex(C_DOUBLE_COMPLEX), pointer :: out_backward(:,:)
    type(C_PTR) :: p, q, r, s
    integer(C_INT) :: dimx, dimy  ! array dimensions

    real(8), INTENT(IN) :: x
    complex(8), DIMENSION(:,:,:), INTENT(IN) :: y
    complex(8), DIMENSION(:,:,:), INTENT(OUT) :: dydx

    dimx=size(pdb,1)
    dimy=size(pdb,2)

    ! allocating memory contiguously using C function
    p = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    q = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    r = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    s = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
 
    ! here we use the usual fortran order
    call c_f_pointer(p, in_forward, [dimx,dimy])
    call c_f_pointer(q, out_forward, [dimx,dimy])
    call c_f_pointer(r, in_backward, [dimx,dimy])
    call c_f_pointer(s, out_backward, [dimx,dimy])
 
    ! prepare plans needed by fftw3
    ! here we must make sure we reverse the array dimensions for FFTW
    plan_forward = fftw_plan_dft_2d(dimy, dimx, in_forward, out_forward, FFTW_FORWARD, FFTW_ESTIMATE)
    plan_backward = fftw_plan_dft_2d(dimy, dimx, in_backward, out_backward, FFTW_BACKWARD, FFTW_ESTIMATE)


    !generate the energy dependence of the pump
    pump=pump_spatial*(cos(omega_p*x)+(0.0,-1.0)*sin(omega_p*x))

    !right hand side of the equation of motion for the photon population

    !photon part
    dydx(:,:,1)= (0.0,-1.0)*y(:,:,2)-(0.0,1.0)*delta*y(:,:,1)-&
         kappa_C*y(:,:,1)+(0.0,-1.0)*pump(:,:)
    !add a potential to the photon part
    dydx(:,:,1) = dydx(:,:,1)+(0.0,1.0)*pot_c(:,:)*y(:,:,1)
    !exciton part
    dydx(:,:,2)= -kappa_X*y(:,:,2)+&
         (0.0,-1.0)*ABS(y(:,:,2))*ABS(y(:,:,2))*y(:,:,2)/norm+ &
         &(0.0,-1.0)*y(:,:,1)

    !adding the kinetic energy by means of the FFT

    !fft to momentum space
    in_forward(:,:)=y(:,:,1)
    call fftw_execute_dft(plan_forward, in_forward, out_forward)
    in_backward=kinetic*out_forward
    !fft back to real space
    call fftw_execute_dft(plan_backward, in_backward, out_backward)
    dydx(:,:,1) = dydx(:,:,1) + out_backward * cmplx(0.0,-1.0)
 
    ! avoiding any potential memory leaks
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
 
    call fftw_free(p)
    call fftw_free(q)
    call fftw_free(r)
    call fftw_free(s)

  END SUBROUTINE derivs

  SUBROUTINE odeint_rk(ystart,x1,x2,eps,h1,hmin)
    IMPLICIT NONE    

    complex(8), DIMENSION(:,:,:), INTENT(INOUT) :: ystart    
    real(8), INTENT(IN) :: x1,x2,eps,h1,hmin    
  
    real(8), PARAMETER :: TINY=1.0e-30    
    INTEGER(I4B), PARAMETER :: MAXSTP=1000000
    !Runge-Kutta driver with adaptive step size control.Integrate the array
    !of starting values ystart from x1 to x2 with accuracy eps storing
    !intermediate results in the module variables in ode_path. h1 should be
    !set as a guessed first stepsize, hmin as the minimum allowed stepsize
    !(can be zero). On output ystart is replaced by values at the end of the
    !integration interval. derivs is the user-supplied subroutine for
    !calculating the right-hand-side derivative, while rkqs is the name of
    !he stepper routine to be used.
    
    INTEGER(I4B) :: nstp    
    real(8) :: h,hdid,hnext,x,xsav    
    complex(8), DIMENSION(size(ystart,1),size(ystart,2),size(ystart,3)) :: dydx,y,yscal    
    real(8), External :: findfermpart,findcoopernum    
    complex(8), External :: findfermcond    


    x=x1    
    h=sign(h1,x2-x1)    
    nok=0    
    nbad=0    
    kount=0    
    y(:,:,:)=ystart(:,:,:)    
    if (save_steps) then    
       xsav=x-2.0*dxsav_rk    
    end if    
    do nstp=1,MAXSTP     
    !Take at most MAXSTP steps
       call derivs(x,y,dydx)    
       yscal(:,:,:)=abs(y(:,:,:))+abs(h*dydx(:,:,:))+TINY    
       !Scaling used to monitor accuracy. This general purpose choice can be
       !modified if need be.
       if (save_steps .and. (abs(x-xsav) > abs(dxsav_rk))) &     
       !Store intermediate results.
            call save_a_step_rk    
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
          if (save_steps) call save_a_step_rk     
          !Save final step.
          RETURN     
       !Normal exit.
       end if    
       if (abs(hnext) < hmin) write(*,*) "stepsize smaller than minimum in odeint"    
       h=hnext    
    end do    
    write(*,*) "too many steps in odeint"    
  CONTAINS    
    SUBROUTINE save_a_step_rk    
      real(8) :: sx,sy    
      integer :: ix, iy    
      write(label,FMT="(i3)") kount    
      open(unit=22, file="phcplx-opo_spc"//trim(adjustl(label))//".dat", status='replace')    
      write(22, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(1))", 1x, "aimag(psi(1))") ')    
      open(unit=23, file="excplx-opo_spc"//trim(adjustl(label))//".dat", status='replace')    
      write(23, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(2))", 1x, "aimag(psi(2))") ')    
      do iy=1, Ny    
         sy=-Ly+(iy-1)*ay    
         do ix=1, Nx
            sx=-Lx+(ix-1)*ax
            write(22, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy,&
                   &real(y(ix,iy,1))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70, aimag(y(ix,iy,1))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70
            write(23, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy,&
                   &real(y(ix,iy,2))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70, aimag(y(ix,iy,2))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70
         end do
         write(22,*)    
         write(23,*)    
      end do    
      close(22)    
      close(23)    
      kount=kount+1    
      xsav=x    
    END SUBROUTINE save_a_step_rk    
  END SUBROUTINE odeint_rk    

  SUBROUTINE odeint_sp(ystart,x1,x2,eps,h1,hmin)
    IMPLICIT NONE    
    complex(8), DIMENSION(:,:,:), INTENT(INOUT) :: ystart    
    real(8), INTENT(IN) :: x1,x2,eps,h1,hmin    
  
    real(8), PARAMETER :: TINY=1.0e-30    
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
    real(8) :: h,hdid,hnext,x,xsav    
    complex(8), DIMENSION(size(ystart,1),size(ystart,2),size(ystart,3)) :: dydx,y,yscal    
    real(8), External :: findfermpart,findcoopernum    
    complex(8), External :: findfermcond    
    x=x1    
    h=sign(h1,x2-x1)    
    nok=0    
    nbad=0    
    kount=0    
    y(:,:,:)=ystart(:,:,:)    
    if (save_steps) then    
       xsav=x-2.0*dxsav_sp    
    end if    
    do nstp=1,MAXSTP
    !Take at most MAXSTP steps.
       call derivs(x,y,dydx)    
       yscal(:,:,:)=abs(y(:,:,:))+abs(h*dydx(:,:,:))+TINY
       !Scaling used to monitor accuracy. This general purpose choice can be
       !modified if need be.
       if (save_steps .and. x.ge.t_stfft ) then    
       !Store intermediate results.
          call save_a_step_sp    
       end if    
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
          if (save_steps) call save_a_step_sp     
          !Save final step.
          RETURN     
       !Normal exit.
       end if    
       if (abs(hnext) < hmin) write(*,*) "stepsize smaller than minimum in odeint"    
       h=hnext    
    end do    
    write(*,*) "too many steps in odeint"    
  CONTAINS    
    SUBROUTINE save_a_step_sp    
      integer :: i_t    
      if( x.ge.t_stfft+kount*dxsav_sp .and. kount+1.le.Nt) then    
         write(file_time_0,*) kount, x    
         call flush(file_time_0)    
         i_t=kount+1    
         y_tot_0(:,:,i_t)=y(:,:,1)    
         kount=kount+1    
      end if    
      xsav=x    
    END SUBROUTINE save_a_step_sp    
  END SUBROUTINE odeint_sp    


  SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)    
    IMPLICIT NONE    
    complex(8), DIMENSION(:,:,:), INTENT(INOUT) :: y    
    complex(8), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal    
    real(8), INTENT(INOUT) :: x    
    real(8), INTENT(IN) :: htry,eps    
    real(8), INTENT(OUT) :: hdid,hnext    
    INTEGER, EXTERNAL :: assert_eq    
    !Fifth order Runge-Kutta step with monitoring of local truncation error    
    !to ensure accuracy and adjust stepsize. Input are the dependent variable    
    !vector y and its derivative dydx at the starting value of the    
    !independent variable x. Also input are the stepsize to be attempted htry    
    !the required accuracy eps and the vector yscal against which the error    
    !is scaled. y dydx and yscal are all of the same lengh. On output, y and x are    
    !replaced by their new values, hdid is the stepsize that was actually    
    !accomplished, and hnext is the estimated next stepsize. derivs is the    
    !user-supplied subroutine that computes the right-hand-side derivatives.    
    INTEGER(I4B) :: ndum    
    real(8) :: errmax,h,htemp,xnew    
    complex(8), DIMENSION(size(y,1),size(y,2),size(y,3)) :: yerr,ytemp    
    real(8), PARAMETER :: SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,&    
    ERRCON=1.89e-4    
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

  SUBROUTINE rkck(y,dydx,x,h,yout,yerr)
    IMPLICIT NONE
    complex(8), DIMENSION(:,:,:), INTENT(IN) :: y,dydx
    real(8), INTENT(IN) :: x,h
    complex(8), DIMENSION(:,:,:), INTENT(OUT) :: yout,yerr
    INTEGER, EXTERNAL :: assert_eq
    !Given values for N variables y and their derivatives dydx known at x
    !use the fifth order Cash-Karp Runge-Kutta method to advance the
    !solution over an interval h and return the incremented variables as
    !yout Also return an estimate of he local truncation error in yout
    !using the embedded fourth order method. The user supplies the
    !subroutine derivs(x,y,dydx),which returns derivatives dydx at x

    INTEGER(I4B) :: ndum
    complex(8), DIMENSION(size(y,1),size(y,2),size(y,3)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    real(8), PARAMETER :: A2=0.2,A3=0.3,A4=0.6,A5=1.0,&
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
  
end MODULE rk_adaptive

FUNCTION assert_eq(n1,n2,n3,string)  
  CHARACTER(LEN=*), INTENT(IN) :: string  
  INTEGER, INTENT(IN) :: n1,n2,n3  
  INTEGER :: assert_eq  
  if (n1 == n2 .and. n2 == n3) then  
     assert_eq=n1  
  else  
     write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
     STOP   
  end if  
END FUNCTION assert_eq
