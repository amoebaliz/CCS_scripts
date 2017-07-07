
  SUBROUTINE Q2_FROM_D2_AND_MSL(D2,MSL,nx,ny,Q2)

  IMPLICIT NONE

  INTEGER,INTENT(in) :: nx,ny
  REAL(4),DIMENSION(ny,nx),INTENT(in)  :: D2,MSL
  REAL(4),DIMENSION(ny,nx),INTENT(out) :: Q2

  INTEGER :: ji,jj
  REAL(8), PARAMETER :: reps = 0.62197
  REAL(8) :: psat


  DO ji=1,nx   ! outer loop
    DO jj=1,ny ! inner loop

      psat = 100*( 10**(10.79574*(1 - 273.16/D2(jj,ji)) - 5.028*LOG10(D2(jj,ji)/273.16) &
           &       + 1.50475*10**(-4)*(1 - 10**(-8.2969*(D2(jj,ji)/273.16 - 1)))        &
           &       + 0.42873*10**(-3)*(10**(4.76955*(1 - 273.16/D2(jj,ji))) - 1) + 0.78614) )

      Q2(jj,ji) = reps*psat / ( MSL(jj,ji) - (1. - reps)*psat )

    ENDDO
  ENDDO

  END SUBROUTINE

