Module global
  USE ode_path
  implicit none
  character (len=3) :: label
  complex, allocatable :: wave_f_mom(:,:,:)
  complex, allocatable :: y_enfilt(:,:)
  complex, allocatable :: wave_f_flt(:,:,:), wave_f_fltspc(:,:,:)
  complex, allocatable :: psi(:,:,:), pump_spatial(:,:)
  complex, allocatable :: y_tot_0(:,:,:)
  complex, allocatable :: pdb(:,:,:)
  complex, allocatable :: pump(:,:)
  real(8), allocatable :: kinetic(:,:), pot_c(:,:)
  real(8), allocatable :: int_sp(:)
  integer :: run, in_sswf_rk, in_sswf_sp
  integer :: kx_max(1)
  integer :: Nx, Ny, Nt
  integer :: file_time_0=29
  integer :: kount_st, kount_end
  real(8) :: Lx, Ly, ax, ay
  real(8) :: delta
  real(8) :: kappa_C, kappa_X, sigma_p, k_p
  real(8) :: f_p, omega_p
  real(8) :: omega_cut_p, omega_s, omega_cut_s, omega_id, omega_cut_id
  real(8) :: norm, tot_h
  real(8) :: eps_r, gv, def_x_pos, def_y_pos
  real(8) :: t_stfft
  real(8) :: mom_cent, mom_cut

  NAMELIST /indata/               &
       &     run,                  &   ! run from scratch?
       &     in_sswf_rk,           &   ! if not, starting wf
       &     in_sswf_sp,           &   ! startingn wf for spectrum
       &     kount_st,             &   ! starting from frame number
       &     kount_end,            &   ! up to frame number
       &     kappa_C,              &   ! decay rate for photons
       &     kappa_X,              &   ! decay rate for excitons
       &     delta,                &   ! detuning
       &     k_p,                  &   ! pump  angle
       &     omega_p,              &   ! pump frequency
       &     sigma_p,              &   ! pump  spatial
       &     f_p,                  &   ! pump strenght
       &     mom_cent,             &   ! center of momentum
       &     mom_cut,              &   ! momentum window
       &     omega_cut_p,          &   ! pump energy window
       &     omega_s,              &   ! center of signal in en
       &     omega_cut_s,          &   ! signal energy window
       &     omega_id,             &   ! center of idler in en
       &     omega_cut_id,         &   ! idler energy window
       &     Lx,                   &   ! x size of box (L)
       &     Ly,                   &   ! y size of box (L)
       &     Nx,                   &   ! Number of points in x direction
       &     Ny,                   &   ! Number of points in y direction
       &     Nt,                   &   ! point in spectra
       &     tot_h,                &   !
       &     t_stfft,              &   !
       &     dxsav_rk,             &   ! time interval for saving function
       &     dxsav_sp,             &   ! time interval for spectra evaluation
       &     eps_r,                &   ! max error for odeint subroutine
       &     def_x_pos,            &   ! x coordinate of defect
       &     def_y_pos,            &   ! y coordinate of defect
       &     gv                        ! strength of defect potential

end Module global
