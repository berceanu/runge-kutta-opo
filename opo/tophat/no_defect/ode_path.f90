Module ode_path 
  USE nrtype
  INTEGER(I4B) :: nok,nbad,kount
  LOGICAL(LGT), SAVE :: save_steps=.true.
  REAL(kind=4) :: dxsav_rk, dxsav_sp
end Module ode_path

