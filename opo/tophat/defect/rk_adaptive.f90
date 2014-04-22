module rk_adaptive_module
  use FFTW3
  use global
  implicit none

  private ! make everything private by default

  ! fft stuff
  ! forward means real space to momentum space, backward the opposite
  type(C_PTR) :: plan_forward, plan_backward
  complex(C_DOUBLE_COMPLEX), pointer :: in_forward(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: out_forward(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: in_backward(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: out_backward(:,:)
  type(C_PTR) :: p, q, r, s
  integer(C_INT) :: dimx, dimy

  public :: odeint_rk, setg

contains

  subroutine create_fftw
    integer(C_INT) iret ! fftw multi-thread initialization return code

    iret = fftw_init_threads()
    write(*,*) iret
    call fftw_plan_with_nthreads(8)

    ! fft stuff
    dimx=size(pdb,1)
    dimy=size(pdb,2)

    ! allocating memory contiguously using C function
    p = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    q = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    r = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
    s = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))

    ! make pointers from C to FORTRAN
    ! here we use the usual fortran order
    call c_f_pointer(p, in_forward, [dimx,dimy])
    call c_f_pointer(q, out_forward, [dimx,dimy])
    call c_f_pointer(r, in_backward, [dimx,dimy])
    call c_f_pointer(s, out_backward, [dimx,dimy])
 
    ! prepare plans needed by fftw3
    ! here we must make sure we reverse the array dimensions for FFTW
    plan_forward = fftw_plan_dft_2d(dimy, dimx, in_forward, out_forward, FFTW_FORWARD, FFTW_PATIENT)
    plan_backward = fftw_plan_dft_2d(dimy, dimx, in_backward, out_backward, FFTW_BACKWARD, FFTW_PATIENT)

  end subroutine create_fftw


  subroutine destroy_fftw
    ! avoiding any potential memory leaks
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
 
    call fftw_free(p)
    call fftw_free(q)
    call fftw_free(r)
    call fftw_free(s)

    call fftw_cleanup_threads()

  end subroutine destroy_fftw


  subroutine derivs(x,y,dydx)
    real(dp), INTENT(IN) :: x
    complex(dpc), DIMENSION(:,:,:), INTENT(IN) :: y
    complex(dpc), DIMENSION(:,:,:), INTENT(OUT) :: dydx
 
    !generate the energy dependence of the pump
    pump=pump_spatial*(cos(omega_p*x)-I*sin(omega_p*x))

    !right hand side of the equation of motion

    !photon part
    dydx(:,:,1)= -I*y(:,:,2)-I*delta*y(:,:,1)-&
         kappa_C*y(:,:,1)-I*pump(:,:)
    !add a potential to the photon part
    dydx(:,:,1) = dydx(:,:,1)+I*pot_c(:,:)*y(:,:,1)
    !exciton part
    dydx(:,:,2)= -kappa_X*y(:,:,2)-&
         I*ABS(y(:,:,2))**2*y(:,:,2)/norm- &
         &I*y(:,:,1)

    !adding the kinetic energy by means of the FFT

    !fft to momentum space
    in_forward(:,:)=y(:,:,1)
    call fftw_execute_dft(plan_forward, in_forward, out_forward)
    out_forward = out_forward/sqrt(real(dimx*dimy, dp)) !normalization
    in_backward=kinetic*out_forward
    !fft back to real space
    call fftw_execute_dft(plan_backward, in_backward, out_backward)
    out_backward = out_backward/sqrt(real(dimx*dimy, dp)) !normalization
    dydx(:,:,1) = dydx(:,:,1) - I * out_backward
 
  end subroutine derivs

  subroutine odeint_rk(ystart,x1,x2,eps,h1,hmin)
  !Runge-Kutta driver with adaptive step size control.Integrate the array
  !of starting values ystart from x1 to x2 with accuracy eps storing
  !intermediate results in the module variables in ode_path. h1 should be
  !set as a guessed first stepsize, hmin as the minimum allowed stepsize
  !(can be zero). On output ystart is replaced by values at the end of the
  !integration interval. derivs is the user-supplied subroutine for
  !calculating the right-hand-side derivative, while rkqs is the name of
  !he stepper routine to be used.

    complex(dpc), DIMENSION(:,:,:), INTENT(INOUT) :: ystart    
    real(dp), INTENT(IN) :: x1,x2,eps,h1,hmin    
  
    real(dp), PARAMETER :: TINY=1.0e-30_dp
    integer, PARAMETER :: MAXSTP=1000000000

    integer :: nstp    
    real(dp) :: h,hdid,hnext,x,xsav    
    complex(dpc), DIMENSION(size(ystart,1),size(ystart,2),size(ystart,3)) :: dydx,y,yscal    

    x=x1    
    h=sign(h1,x2-x1)    
    nok=0    
    nbad=0    
    kount=0    
    y(:,:,:)=ystart(:,:,:)    
    if (save_steps) then    
       xsav=x-2*dxsav_rk    
    end if    
    call create_fftw
    do nstp=1,MAXSTP     
    !Take at most MAXSTP steps
       call derivs(x,y,dydx)    
       yscal(:,:,:)=abs(y(:,:,:))+abs(h*dydx(:,:,:))+TINY    
       !Scaling used to monitor accuracy. This general purpose choice can be
       !modified if need be.
       if (save_steps .and. (abs(x-xsav) > abs(dxsav_rk))) &     
       !Store intermediate results.
            call save_a_step_rk    
       if ((x+h-x2)*(x+h-x1) > 0.0_dp) h=x2-x     
       !If stepsize can overshoot,decrease.
       call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext)    
       if (hdid == h) then    
          nok=nok+1    
       else    
          nbad=nbad+1    
       end if    
       if ((x-x2)*(x2-x1) >= 0.0_dp) then !Are we done?
          ystart(:,:,:)=y(:,:,:)    
          if (save_steps) call save_a_step_rk !Save final step.
          call destroy_fftw ! clear fft memory
          RETURN !Normal exit.
       end if    
       if (abs(hnext) < hmin) write(*,*) "stepsize smaller than minimum in odeint"    
       h=hnext    
    end do    
    call destroy_fftw ! clear fft memory
    write(*,*) "too many steps in odeint"    
  CONTAINS    
    subroutine save_a_step_rk    
      real(dp) :: sx,sy    
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
                   &real(y(ix,iy,1))*norm_c, aimag(y(ix,iy,1))*norm_c
            write(23, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy,&
                   &real(y(ix,iy,2))*norm_c, aimag(y(ix,iy,2))*norm_c
         end do
         write(22,*)    
         write(23,*)    
      end do    
      close(22)    
      close(23)    
      kount=kount+1    
      xsav=x    
    end subroutine save_a_step_rk    
  end subroutine odeint_rk    


  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)    
    complex(dpc), DIMENSION(:,:,:), INTENT(INOUT) :: y    
    complex(dpc), DIMENSION(:,:,:), INTENT(IN) :: dydx,yscal    
    real(dp), INTENT(INOUT) :: x    
    real(dp), INTENT(IN) :: htry,eps    
    real(dp), INTENT(OUT) :: hdid,hnext    
    INTEGER :: assert_eq    
    !Fifth order Runge-Kutta step with monitoring of local truncation error    
    !to ensure accuracy and adjust stepsize. Input are the dependent variable    
    !vector y and its derivative dydx at the starting value of the    
    !independent variable x. Also input are the stepsize to be attempted htry    
    !the required accuracy eps and the vector yscal against which the error    
    !is scaled. y dydx and yscal are all of the same lengh. On output, y and x are    
    !replaced by their new values, hdid is the stepsize that was actually    
    !accomplished, and hnext is the estimated next stepsize. derivs is the    
    !user-supplied subroutine that computes the right-hand-side derivatives.    
    integer :: ndum    
    real(dp) :: errmax,h,htemp,xnew    
    complex(dpc), DIMENSION(size(y,1),size(y,2),size(y,3)) :: yerr,ytemp    
    real(dp), PARAMETER :: SAFETY=0.9_dp, PGROW=-0.2_dp, PSHRNK=-0.25_dp,&    
    ERRCON=1.89e-4_dp    
    !The value ERRCON equals (5/SAFETY)**(1/PGROW),see use below.    
    ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')     
    h=htry     
    !Set step size to the initial trial value.     
    do    
       call rkck(y,dydx,x,h,ytemp,yerr)     
       !Take a step.    
       errmax=maxval(abs(yerr(:,:,:)/yscal(:,:,:)))/eps     
       !Evaluate accuracy.    
       if (errmax <= 1.0_dp) exit     
       !Step succeeded.    
       htemp=SAFETY*h*(errmax**PSHRNK)     
       !Truncation error too large, reduce stepsize.    
       h=sign(max(abs(htemp),0.1_dp*abs(h)),h)     
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
       hnext=5.0_dp*h    
    end if    
    hdid=h      
    x=x+h       
    y(:,:,:)=ytemp(:,:,:)      
  end subroutine rkqs    

  subroutine rkck(y,dydx,x,h,yout,yerr)
    complex(dpc), DIMENSION(:,:,:), INTENT(IN) :: y,dydx
    real(dp), INTENT(IN) :: x,h
    complex(dpc), DIMENSION(:,:,:), INTENT(OUT) :: yout,yerr
    INTEGER :: assert_eq
    !Given values for N variables y and their derivatives dydx known at x
    !use the fifth order Cash-Karp Runge-Kutta method to advance the
    !solution over an interval h and return the incremented variables as
    !yout Also return an estimate of he local truncation error in yout
    !using the embedded fourth order method. The user supplies the
    !subroutine derivs(x,y,dydx),which returns derivatives dydx at x

    integer :: ndum
    complex(dpc), DIMENSION(size(y,1),size(y,2),size(y,3)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    real(dp), parameter :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
             A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,&
             B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
             B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
             B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
             B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
             B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
             C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,&
             C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
             DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
             DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp
    ndum=assert_eq(size(y),size(dydx),size(yout),'rkck')
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
  end subroutine rkck
  
  subroutine setg(Nx, Ny, Lx, Ly, kinetic)
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: Lx, Ly
    real(dp), dimension(Nx, Ny), intent(out) :: kinetic

    integer :: j,k    
    real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp 

    do j=1,(Ny/2+1)    
       do k=1,(Nx/2+1)    
          kinetic(k,j)=pi**2*(&    
               &(k-1)**2/(Lx**2)+(j-1)**2/(Ly**2))    
       end do    
    end do    
    do j=(Ny/2+2),Ny    
       do k=(Nx/2+2),Nx    
          kinetic(k,j)=pi**2*( &    
               & (k-1-Nx)**2/(Lx**2)+(j-1-Ny)**2/(Ly**2))    
       end do    
    end do    
    do j=1,(Ny/2+1)    
       do k=(Nx/2+2),Nx    
          kinetic(k,j)=pi**2*(&    
               &(k-1-Nx)**2/(Lx**2)+(j-1)**2/(Ly**2))    
       end do    
    end do    
    do j=(Ny/2+2),Ny    
       do k=1,(Nx/2+1)      
          kinetic(k,j)=pi**2*(&    
               &(k-1)**2/(Lx**2)+(j-1-Ny)**2/(Ly**2))    
       end do    
    end do    
  end subroutine

  function assert_eq(n1, n2, n3, string) result(res)
    integer, intent(in) :: n1,n2,n3  
    character(len=4), intent(in) :: string  
    integer :: res 
  
    if ((n1 == n2).and.(n2 == n3)) then  
       res = n1  
    else  
       write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
       stop   
    end if  
  end function assert_eq

end module
