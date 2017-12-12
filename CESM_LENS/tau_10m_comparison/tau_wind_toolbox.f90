
  SUBROUTINE U_V_FROM_TAUX_TAUY_U10(TAUX,TAUY,U10,nx,ny,Uwind,Vwind)
  IMPLICIT NONE

  INTEGER,INTENT(in) :: nx,ny
  REAL(4),DIMENSION(ny,nx),INTENT(in)  :: TAUX,TAUY,U10
  REAL(4),DIMENSION(ny,nx),INTENT(out) :: Uwind,Vwind

  INTEGER :: ji,jj
  REAL(8), PARAMETER :: rho_air = 1.2       ! kg m-3
  REAL(8), PARAMETER :: c4      = 0.0027    ! m s-1   
  REAL(8), PARAMETER :: c5      = 0.000142  ! (dimensionless)
  REAL(8), PARAMETER :: c6      = 0.0000764 ! m-1 s
  REAL(8) :: Cd

  DO ji=1,nx   ! outer loop
    DO jj=1,ny ! inner loop

       Cd = c4/U10(jj,ji) + c5 + c6*U10(jj,ji)
       
       Uwind(jj,ji) = Taux(jj,ji)/(rho_air*Cd*U10(jj,ji))
       Vwind(jj,ji) = Tauy(jj,ji)/(rho_air*Cd*U10(jj,ji))

    ENDDO
  ENDDO

  END SUBROUTINE

