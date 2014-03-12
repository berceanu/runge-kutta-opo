!-----------------------------------------------------------------!
! Module global                                                   !
!-----------------------------------------------------------------!
Module global
  USE MKL_DFTI
  implicit none
  character (len=3) :: label
  complex, allocatable :: pump_spatial(:,:)
  complex, allocatable :: probe_spatial(:,:)
  complex, allocatable :: pdb(:,:,:),psiref(:,:,:)
  complex, allocatable :: pump(:,:), probe(:,:)
  complex, allocatable ::  psiE(:,:,:),auxpsi(:,:,:),transfE(:)
  complex, allocatable :: kvecx(:,:),kvecy(:,:)
  real, allocatable :: kinetic(:,:),pot_x(:,:),pot_c(:,:)
  real, allocatable :: dragPOTx(:,:),dragPOTy(:,:)
  integer :: Nx, Ny, ixbin, ixbfin,iybin,iybfin
  integer :: Lz,nene,ti,tti,tavs,tavm,tavp
  real  :: Lx, Ly, ax, ay
  real  :: delta,rnd
  real  :: kappa_C, kappa_X, sigma_p, sigma_pb, sigma_tp, sigma_tpb, k_p, k_pb
  real  :: theta_pb
  real  :: f_p, f_pb, assymx_p,assymy_p,assymx_pb,assymy_pb, omega_p, omega_pb
  real  :: norm,t_init
  real  :: eps_r
  real  :: mom_cent1, mom_cut1, mom_cent2, mom_cut2
  real  :: ene_cent, ene_cut, bc_x,bc_y,br_x,br_y
  real  :: barrier,bin_x,bfin_x,bin_y,bfin_y,normDFTIsqrt
  real  :: dxsav,dxene,thetasin
  integer :: nok,nbad,kount
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  integer :: npointxy(1:2),status
  complex :: transformingF(200,200)
  complex :: DFTItransformingF(40000)
  equivalence(transformingF,DFTItransformingF)
  REAL, PARAMETER :: PI=3.141592653589793238462643383279502884197

  NAMELIST /indata/                &
       &     kappa_C,              &   ! decay rate for photons
       &     kappa_X,              &   ! decay rate for excitons
       &     delta,                &   ! detuning
       &     Lz,                   &   ! Lz of the vortex, Lz=1, 2, ...
       &     k_p,                  &   ! pump  angle
       &     k_pb,                 &   ! probe angle
       &     theta_pb,             &   ! probe angle
       &     omega_p,              &   ! pump frequency
       &     omega_pb,             &   ! probe frequency
       &     sigma_p,              &   ! pump  spatial
       &     sigma_pb,             &   ! probe spatial
       &     sigma_tp,             &   ! pump pulse
       &     sigma_tpb,            &   ! probe pulse
       &     t_init,               &   ! initial time of the pulse
       &     f_p,                  &   ! pump strenght
       &     f_pb,                 &   ! probe strenght
       &     assymx_p,             &   ! x position of pump
       &     assymy_p,             &   ! y position of pump
       &     assymx_pb,            &   ! x position of probe
       &     assymy_pb,            &   ! y position of probe
       &     mom_cent1,            &
       &     mom_cut1,             &
       &     mom_cent2,            &
       &     mom_cut2,             &
       &     ene_cent,             &
       &     ene_cut,              &
       &     Lx,                   &   ! x size of box (L)
       &     Ly,                   &   ! y size of box (L)
       &     Nx,                   &   ! Number of points in x direction
       &     Ny,                   &   ! Number of points in y direction
       &     barrier,              &   ! barrier height
       &     bc_x,                &   ! barrier start
       &     br_x,               &   ! barrier end
       &     bc_y,                &   ! barrier start
       &     br_y,               &   ! barrier end
       &     dxsav,                &   ! time interval for saving function
       &     dxene,                &   ! time interval for spectra evaluation
       &     nene,                 &   ! point in spectra
       &     eps_r                     ! max error for odeint subroutine

end Module global
