Module global
  USE ode_path
  implicit none
  character (len=3) :: label, label_n
  integer, allocatable :: check(:)
  complex, allocatable :: y_enfilt(:,:,:)
  complex, allocatable :: pump_spatial(:,:)
  complex, allocatable :: y_tot(:,:,:)
  complex, allocatable :: pdb(:,:,:)
  complex, allocatable :: pump(:,:)
  real, allocatable :: trig_x(:), trig_y(:), trig_t(:), trig_redt(:)
  real, allocatable :: kinetic(:,:), pot_c_t(:,:,:)
  real, allocatable :: int_sp(:), int_sp_ok(:), int_sp_s(:), int_sp_p(:), int_sp_i(:)
  integer, parameter :: SP2 = KIND(1.0)
  integer, parameter :: dp2=kind(1.0d0)
  integer :: in_sswf_sp
  integer :: i_tmax_p(1), i_tmax_s(1), i_tmax_i(1)
  integer :: Nx, Ny, Nt, Nredt
  integer :: file_time=29
  integer :: n_save, time_kount, tot_save
  real, parameter :: delta=0.0
  real :: Lx, Ly, ax, ay
  real :: kappa_C, kappa_X
  real :: f_p, omega_p, sigma_p, k_p
  real :: omega_pmax, omega_smax, omega_imax
  real :: norm
  real :: eps_r
  real :: startx, endx, v, gv
  integer :: tot_h, t_stfft
  
  NAMELIST /indata/                &
       &     in_sswf_sp,           &   ! startingn wf for spectrum
       &     kappa_C,              &   ! decay rate for photons
       &     kappa_X,              &   ! decay rate for excitons
       &     k_p,                  &   ! pump  angle
       &     omega_p,              &   ! pump frequency
       &     sigma_p,              &   ! pump  spatial
       &     f_p,                  &   ! pump strenght
       &     Lx,                   &   ! x size of box (L)
       &     Ly,                   &   ! y size of box (L)
       &     Nx,                   &   ! Number of points in x direction
       &     Ny,                   &   ! Number of points in y direction
       &     Nt,                   &   ! point in spectra
       &     Nredt,                &   ! en points around signal
       &     tot_save,             &   ! how many intervals of Nt
       &     startx,               &   ! starting position of defect
       &     endx,                 &   ! final coordinate of defect
       &     v,                    &   ! defect velocity in FM units
       &     t_stfft,              &   ! time steps after which start recording
       &     dxsav_sp,             &   ! time interval for spectra evaluation
       &     eps_r,                &   ! max error for odeint subroutine
       &     gv                        ! strength of defect potential
  
end Module global
