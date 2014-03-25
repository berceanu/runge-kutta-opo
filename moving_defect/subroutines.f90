Module subroutines
  USE global
  IMPLICIT NONE
  SAVE
CONTAINS

  Subroutine read
    IMPLICIT NONE

    OPEN(UNIT=22,FILE='INPUT_timesp',STATUS='old')
    READ(22,NML=indata)
    CLOSE(22)
    
    ax=2.0*Lx/Nx    
    ay=2.0*Ly/Ny    
    norm=ax*ay    
    f_p=f_p*256/sqrt(Nx*Ny*1.0)*sqrt(Lx*Ly*1.0)/70    

  end Subroutine read

  SUBROUTINE  init_pdb
    IMPLICIT NONE
    
    integer :: ix, iy
    real  :: sx, sy  
    real :: re_y1, im_y1, re_y2, im_y2
  
    open(unit=22, file="phcplx-opo_spc"//trim(adjustl(label))//".dat", status='old')    
    read(22, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(1))", 1x, "aimag(psi(1))") ')    
    open(unit=23, file="excplx-opo_spc"//trim(adjustl(label))//".dat", status='old')    
    read(23, fmt=' ("#", 1x, "x", 12x, "y", 12x, "real(psi(2))", 1x, "aimag(psi(2))") ')    
    do iy=1, Ny    
       do ix=1, Nx    
          read(22, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, re_y1, im_y1
          read(23, fmt=' (1x, d12.5, 1x, d12.5, 1x, d12.5, 1x, d12.5) ') sx, sy, re_y2, im_y2
        
          pdb(ix,iy,1)=&    
               (1.0,0.0)*re_y1*256/sqrt(1.0*Nx*Ny)*sqrt(Lx*Ly*1.0)/70+&    
               (0.0,1.0)*im_y1*256/sqrt(1.0*Nx*Ny)*sqrt(Lx*Ly*1.0)/70    
          pdb(ix,iy,2)=&    
               (1.0,0.0)*re_y2*256/sqrt(1.0*Nx*Ny)*sqrt(Lx*Ly*1.0)/70+&    
               (0.0,1.0)*im_y2*256/sqrt(1.0*Nx*Ny)*sqrt(Lx*Ly*1.0)/70    
      end do    
       read(22,*)    
       read(23,*)    
    end do    
    close(22)    
    close(23)

  END SUBROUTINE  init_pdb

  Subroutine init_pot_c_t
    implicit none

    integer :: indx, indy, t

    pot_c_t = 0.0

    indx = (def_x_pos+Lx)/ax + 1
    indy = (def_y_pos+Ly)/ay + 1
    !for old coordinates of (Nx/2,Ny/2) set def_x_pos = -ax, def_y_pos=-ay

    do t=1, tot_h
        pot_c_t(indx,indy,t)=gv
    end do

  end Subroutine init_pot_c_t


  SUBROUTINE  init_pump_th
    IMPLICIT NONE
    
    integer :: ix, iy
    real  :: sx, sy  
		
    !top hat pump    
    open(unit=25, file='pump.dat', status='replace')    
    do iy=1, Ny    
       sy=-Ly+(iy-1)*ay    
       do ix=1, Nx    
          sx=-Lx+(ix-1)*ax    
          pump_spatial(ix,iy)=f_p*0.5*&    
               ( tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)+sigma_p ))-&    
               tanh((1.0/10.0)*( sqrt(sx*sx+sy*sy)-sigma_p )) ) + (0.0,0.0)    
          write(25,*) sx, sy, abs(pump_spatial(ix,iy))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70    
       end do    
       write(25,*)    
    end do    
    close(25)
        
    do iy=1, Ny    
       !sy=-Ly+(iy-1)*ay    
       do ix=1, Nx    
          sx=-Lx+(ix-1)*ax    
          pump_spatial(ix,iy)= pump_spatial(ix,iy)*cos(k_p*sx)+(0.0,1.0)*pump_spatial(ix,iy)*sin(k_p*sx)    
       end do    
    end do

  END SUBROUTINE  init_pump_th

  SUBROUTINE  init_pump_homo
    IMPLICIT NONE
    
    integer :: ix, iy
    real  :: sx, sy  

    !homogeneous pumping		
    open(unit=25, file='pump.dat', status='replace')    
    do iy=1, Ny    
       sy=-Ly+(iy-1)*ay    
       do ix=1, Nx    
          sx=-Lx+(ix-1)*ax    
          pump_spatial(ix,iy)=f_p + (0.0,0.0)    
          write(25,*) sx, sy, abs(pump_spatial(ix,iy))*sqrt(1.0*Nx*Ny)/256/sqrt(Lx*Ly*1.0)*70    
       end do    
       write(25,*)    
    end do    
    close(25)    
        
    do iy=1, Ny    
       !sy=-Ly+(iy-1)*ay    
       do ix=1, Nx    
          sx=-Lx+(ix-1)*ax    
          pump_spatial(ix,iy)= pump_spatial(ix,iy)*cos(k_p*sx)+(0.0,1.0)*pump_spatial(ix,iy)*sin(k_p*sx)    
       end do    
    end do

  END SUBROUTINE  init_pump_homo

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

  end subroutine setg


  Subroutine eval_spectr
    USE nag_fft, ONLY: nag_fft_2d, nag_fft_3d, nag_fft_trig
    implicit none
    integer :: i_t, idt
    integer :: kx
    integer :: i_t_1, i_t_2, i_t_3, i_t_4
    real :: omega_1, omega_2, omega_3, omega_4
    real :: omega, width_p
    real :: mom_x

    !fft to momentum and energy space
    y_tot(:,:,:)=&
         nag_fft_3d(y_tot(:,:,:),trig_1=trig_x,trig_2=trig_y,trig_3=trig_t)

    write(label_n,FMT="(i3)") n_save
    !full spectum
    open(unit=23, file="spectr_om-vs-kx"//trim(adjustl(label_n))//".dat", status='replace')
    write(23, fmt=' ("#", 1x, "mom_x", 19x, "omega", 19x, "abs(psi(1))**2") ')
    do i_t=Nt/2+2, Nt
       omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav_sp )
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          write(23, *) mom_x, -omega, abs(y_tot(kx,1,i_t))*abs(y_tot(kx,1,i_t))
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          write(23, *) mom_x, -omega, abs(y_tot(kx,1,i_t))*abs(y_tot(kx,1,i_t))
       end do
       write(23,*)
    end do
    do i_t=1, Nt/2+1
       omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav_sp )
       do kx=Nx/2+2, Nx
          mom_x=pi*(kx-1-Nx)/Lx
          write(23, *) mom_x, -omega, abs(y_tot(kx,1,i_t))*abs(y_tot(kx,1,i_t))
       end do
       do kx=1, Nx/2+1
          mom_x=pi*(kx-1)/Lx
          write(23, *) mom_x, -omega, abs(y_tot(kx,1,i_t))*abs(y_tot(kx,1,i_t))
       end do
       write(23,*)
    end do
    close(23)

    !integrated spectrum
    write(label_n,FMT="(i3)") n_save
    open(unit=24, file="int-spectr"//trim(adjustl(label_n))//".dat", status='replace')
    write(24, fmt=' ("#", 1x, "omega", 20x, "int_omega") ')

    !write(label_n,FMT="(i3)") n_save
    !open(unit=34, file="int-spectr_ok"//trim(adjustl(label_n))//".dat", status='replace')
    !write(34, fmt=' ("#", 1x, "omega", 20x, "int_omega") ')

    int_sp=sum(sum(abs(y_tot), dim=1), dim=1)
    int_sp_ok=0.0
    int_sp_s=0.0
    int_sp_p=0.0
    int_sp_i=0.0

    idt=Nt+1
    do i_t=Nt/2+2, Nt
       omega=2.0*pi*(i_t-1-Nt)/( (Nt-1)*dxsav_sp )
       write(24, *) -omega, int_sp(i_t)
       !reorder integrated spectrum properly
       idt=idt-1
       int_sp_ok(idt)=int_sp(i_t)
    end do
    do i_t=1, Nt/2+1
       omega=2.0*pi*(i_t-1)/( (Nt-1)*dxsav_sp )
       write(24, *) -omega, int_sp(i_t)
       idt=idt-1
       int_sp_ok(idt)=int_sp(i_t)
    end do
    close(24)

    !do i_t=1, Nt
        !omega = 2.0*pi*(i_t-1-Nt/2)/( (Nt-1)*dxsav_sp )
        !write(34,*) omega, int_sp_ok(i_t)
    !end do
    !close(34)

    !split spectrum into 3 parts, for p,s and i
    !calculate the 4 breaking points (energies)
    width_p=0.3
    omega_1=-1
    omega_2=omega_p-width_p
    omega_3=omega_p+width_p
    omega_4=0.5

    !find indices corresponding to these energies
    i_t_1 = 1 + Nt/2 + omega_1/(2.0*pi) * (Nt-1)*dxsav_sp
    i_t_2 = 1 + Nt/2 + omega_2/(2.0*pi) * (Nt-1)*dxsav_sp
    i_t_3 = 1 + Nt/2 + omega_3/(2.0*pi) * (Nt-1)*dxsav_sp
    i_t_4 = 1 + Nt/2 + omega_4/(2.0*pi) * (Nt-1)*dxsav_sp

    !split the spectrum in 3
    int_sp_s(i_t_1:i_t_2)=int_sp_ok(i_t_1:i_t_2)
    int_sp_p(i_t_2+1:i_t_3-1)=int_sp_ok(i_t_2+1:i_t_3-1)
    int_sp_i(i_t_3:i_t_4)=int_sp_ok(i_t_3:i_t_4)

    !find location of the pump energy, omega_pmax
    i_tmax_p=maxloc(int_sp_p)

    omega_pmax = 2.0*pi*(i_tmax_p(1)-1-Nt/2)/( (Nt-1)*dxsav_sp )

    !find location of the signal
    i_tmax_s=maxloc(int_sp_s)

    omega_smax = 2.0*pi*(i_tmax_s(1)-1-Nt/2)/( (Nt-1)*dxsav_sp )

    !find location of idler 
    i_tmax_i=maxloc(int_sp_i)

    omega_imax = 2.0*pi*(i_tmax_i(1)-1-Nt/2)/( (Nt-1)*dxsav_sp )

    write(*,*) omega_1, omega_2, omega_3, omega_4, i_t_1, i_t_2, i_t_3, i_t_4,&
               i_tmax_p, i_tmax_s, i_tmax_i, omega_pmax, omega_smax, omega_imax

    !call resolve_signal
    call resolve_signal_in_mom

  end Subroutine eval_spectr


  Subroutine resolve_signal
    use nag_fft, only: nag_fft_2d, nag_fft_trig,  nag_fft_3d
    implicit none

    integer :: i_t, j_t
    integer :: ix, iy
    real :: sx, sy

    j_t=0
    do i_t=i_tmax_s(1)-Nredt, i_tmax_s(1)+Nredt
       if(( i_t.gt.Nt ).or.( i_t.lt.1 ))&
         write(*,*) 'signal filter index out of range!'
       j_t=j_t+1
       y_enfilt(:,:,j_t)=y_tot(:,:,i_t)
    end do

    !fft back in time and space
    y_enfilt(:,:,:)=&
         nag_fft_3d(y_enfilt(:,:,:),&
         inverse=.true.,trig_1=trig_x,trig_2=trig_y,trig_3=trig_redt)

    !time_kount=1 defined at the beginning of the programme
    do j_t=1, 2*Nredt+1

       write(label,FMT="(i3)") time_kount
       !write the signal en-filtered wf in space
       open(unit=26, file="opo_ph-spcenfilt_signal"//trim(adjustl(label))//".dat", status='replace')
       write(26, fmt=' ("#", 1x, "x", 12x, "y", 12x, "|psi(1)|^2") ')
       do iy=1, Ny
          sy=-Ly+(iy-1)*ay
          do ix=1, Nx
             sx=-Lx+(ix-1)*ax
             write(26,*) sx, sy,&
                  abs(y_enfilt(ix,iy,j_t))**2
          end do
          write(26,*)
       end do
       close(26)
       time_kount=time_kount+1
    end do

  end Subroutine resolve_signal

  Subroutine resolve_signal_in_mom
    use nag_fft, only: nag_fft_1d, nag_fft_trig
    implicit none

    integer :: i_t, j_t, i, j
    integer :: kx, ky
    real :: mom_x, mom_y

    j_t=0
    do i_t=i_tmax_s(1)-Nredt, i_tmax_s(1)+Nredt
       if(( i_t.gt.Nt ).or.( i_t.lt.1 ))&
         write(*,*) 'signal filter index out of range!'
       j_t=j_t+1
       y_enfilt(:,:,j_t)=y_tot(:,:,i_t)
    end do

    do i=1, Nx
        do j=1, Ny
                y_enfilt(i,j,:)=nag_fft_1d(y_enfilt(i,j,:), inverse=.true., trig = trig_redt)    
        end do
    end do


    !time_kount=1 defined at the beginning of the programme
    do j_t=1, 2*Nredt+1

       write(label,FMT="(i3)") time_kount
       open(unit=27, file="opo_ph-momenfilt_signal"//trim(adjustl(label))//".dat", status='replace')    
       write(27, fmt=' ("#", 1x, "kx", 12x, "ky", 12x, "|psi(1)|^2") ')    

       do ky=Ny/2+2, Ny    
          mom_y=pi*(ky-1-Ny)/Ly    
          do kx=Nx/2+2, Nx    
             mom_x=pi*(kx-1-Nx)/Lx    
             write(27,*) mom_x, mom_y, &    
                  abs(y_enfilt(kx,ky,j_t))**2
          end do    
          do kx=1, Nx/2+1    
             mom_x=pi*(kx-1)/Lx    
             write(27,*) mom_x, mom_y, &    
                  abs(y_enfilt(kx,ky,j_t))**2
          end do    
          write(27,*)    
       end do    
       do ky=1, Ny/2+1    
          mom_y=pi*(ky-1)/Ly    
          do kx=Nx/2+2, Nx    
             mom_x=pi*(kx-1-Nx)/Lx    
             write(27,*) mom_x, mom_y, &    
                  abs(y_enfilt(kx,ky,j_t))**2
          end do    
          do kx=1, Nx/2+1    
             mom_x=pi*(kx-1)/Lx    
             write(27,*) mom_x, mom_y, &    
                  abs(y_enfilt(kx,ky,j_t))**2
          end do    
          write(27,*)    
       end do    
       close(27)    

       time_kount=time_kount+1
    end do

  end Subroutine resolve_signal_in_mom

end Module subroutines
