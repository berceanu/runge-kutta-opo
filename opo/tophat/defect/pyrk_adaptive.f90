module rk_adaptive_interface
  use iso_c_binding, only: c_double, c_int
  use rk_adaptive_module, only: setg
  implicit none

contains

  subroutine c_setg(Nx, Ny, Lx, Ly, kinetic) bind(c)
    integer(c_int), intent(in) :: Nx, Ny
    real(c_double), intent(in) :: Lx, Ly
    real(c_double), dimension(Nx, Ny), intent(out) :: kinetic

    call setg(Nx, Ny, Lx, Ly, kinetic)
  end subroutine

end module
