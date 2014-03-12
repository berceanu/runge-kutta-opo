      !-----------------------------------------------------------------!
      ! Main Program                                                    !
      !                                                                 !
      !-----------------------------------------------------------------!
      Program dynamika
        USE MKL_DFTI
        USE global
        USE mod_subroutines
        USE mod_rk_adaptive
        IMPLICIT NONE
        REAL  :: x1_r,x2_r,h1_r,hmin_r

        CALL read
        CALL init
        Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, &
             DFTI_COMPLEX, 2, npointxy )
        Status = DftiCommitDescriptor( My_Desc1_Handle )
        CALL setg
        CALL setgk
        CALL dragPOT

        x1_r=0.0
!        x2_r=t_init+2.0*nene*dxene
        x2_r=t_init+1.*nene*dxene
        h1_r=0.00000001
        hmin_r=0.0

        CALL odeint(pdb,x1_r,x2_r,eps_r,h1_r,hmin_r)
        CALL energytransform

      end Program dynamika

