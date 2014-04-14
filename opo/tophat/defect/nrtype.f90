module nrtype
  integer, parameter:: dp=kind(0.d0)             ! double
  integer, parameter :: dpc = kind((0.d0, 0.d0)) ! double complex
  complex(dpc), parameter :: one = (1.d0, 0.d0)  ! complex 1
  complex(dpc), parameter :: I = (0.d0, 1.d0)    ! complex i
end module nrtype
