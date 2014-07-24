Program filter_en
  use FFTW3
  use global
  use subroutines
  implicit none

  ! fft stuff
  ! forward means real space to momentum space, backward the opposite
  type(C_PTR) :: plan_backward
  complex(C_DOUBLE_COMPLEX), pointer :: in_backward(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: out_backward(:,:)
  type(C_PTR) :: p, q
  integer(C_INT) :: dimx, dimy
  integer(C_INT) iret ! fftw multi-thread initialization return code

  call read_input 

  allocate(y_tot_0(Nx,Ny,Nt))
  allocate(y_enfilt(Nx,Ny))

  call import_evolution


  iret = fftw_init_threads()
  call fftw_plan_with_nthreads(4)

  ! fft stuff
  dimx=size(y_enfilt,1)
  dimy=size(y_enfilt,2)

  ! allocating memory contiguously using C function
  p = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))
  q = fftw_alloc_complex(int(dimx*dimy, C_SIZE_T))

  ! make pointers from C to FORTRAN
  ! here we use the usual fortran order
  call c_f_pointer(p, in_backward, [dimx,dimy])
  call c_f_pointer(q, out_backward, [dimx,dimy])

  ! prepare plans needed by fftw3
  ! here we must make sure we reverse the array dimensions for FFTW
  plan_backward = fftw_plan_dft_2d(dimy, dimx, in_backward, out_backward, FFTW_BACKWARD, FFTW_PATIENT)

  !file for checking momentum and energy conservation
  open(unit=25, file="en-peaks.dat", status='replace')
  
  !filter the 3 peaks in energy (integrating)
  call filter_peak(omega_p, omega_cut_p, 'pump')
  call write_kx_max('pump')
  call write_peak_mom('pump')
  in_backward=y_enfilt
  !fft back to real space
  call fftw_execute_dft(plan_backward, in_backward, out_backward)
  out_backward = out_backward/sqrt(real(dimx*dimy, dp)) !normalization
  y_enfilt = out_backward
  call write_peak_spc('pump')

  call filter_peak(omega_s, omega_cut_s, 'signal')
  call write_kx_max('signal')        
  call write_peak_mom('signal')
  in_backward=y_enfilt
  !fft back to real space
  call fftw_execute_dft(plan_backward, in_backward, out_backward)
  out_backward = out_backward/sqrt(real(dimx*dimy, dp)) !normalization
  y_enfilt = out_backward
  call write_peak_spc('signal')

  call filter_peak(omega_id, omega_cut_id, 'idler')
  call write_kx_max('idler')        
  call write_peak_mom('idler')
  in_backward=y_enfilt
  !fft back to real space
  call fftw_execute_dft(plan_backward, in_backward, out_backward)
  out_backward = out_backward/sqrt(real(dimx*dimy, dp)) !normalization
  y_enfilt = out_backward
  call write_peak_spc('idler')

  close(25)

  deallocate(y_tot_0)
  deallocate(y_enfilt)

  ! avoiding any potential memory leaks
  call fftw_destroy_plan(plan_backward)

  call fftw_free(p)
  call fftw_free(q)

  call fftw_cleanup_threads()

end Program filter_en
