SUBROUTINE regrid_atmos(Nx, Ny, Xinp, Yinp, Finp, Amin, Amax,          &
&                   Imax, Jmax, Xout, Yout, Jout, Iout, Fout)
!
! Regird atmospheric forcing fields using ROMS 2D regridding algorithm
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates gridded data, Finp, to model locations    !
!  Xout and Yout.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
!     gtype      C-grid type (integer)                                 !
!     iflag      Interpolation flag (integer, 0: linear, 1: cubic)     !
!     Nx         X-dimension size for gridded data, Finp (integer)     !
!     Ny         Y-dimension size for gridded data, Finp (integer)     !
!     Finp       Gridded data to interpolate from (real)               !
!     Amin       Gridded data minimum value (integer)                  !
!     Amax       Gridded data maximum value (integer)                  !
!     LBi        Fout I-dimension Lower bound (integer)                !
!     UBi        Fout I-dimension Upper bound (integer)                !
!     LBj        Fout J-dimension Lower bound (integer)                !
!     UBj        Fout J-dimension Upper bound (integer)                !
!     Imin       Fout starting data I-index (integer)                  !
!     Imax       Fout ending   data I-index (integer)                  !
!     Jmin       Fout starting data J-index (integer)                  !
!     Jmax       Fout ending   data J-index (integer)                  !
!     Xout       X-locations (longitude) to interpolate (real)         !
!     Yout       Y-locations (latitude)  to interpolate (real)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated field (real)                             !
!                                                                      !
!=======================================================================
!
IMPLICIT NONE
!
!  constants from mod_kinds and mod_scalars
!
integer, parameter  :: r8 = selected_real_kind(12,300)  ! 64-bit
real(r8)            :: Eradius = 6371315.0_r8           ! m
real(r8), parameter :: pi = 3.14159265358979323846_r8   ! pi
real(r8), parameter :: deg2rad = pi/180.0_r8            ! degrees to radians conversion
!
!  Imported variable declarations:
!
!  ERAinterim lon, lat, data, dimensions
!
real(r8), DIMENSION(Ny,Nx), intent(in)    :: Xinp
real(r8), DIMENSION(Ny,Nx), intent(in)    :: Yinp
real(r8), DIMENSION(Ny,Nx), intent(inout) :: Finp
integer, intent(in) :: Ny, Nx
real(r8), intent(inout) :: Amin, Amax
!
!  CCS lon, lat, dimensions
!
real(r8), DIMENSION(Jmax,Imax), intent(in) :: Xout
real(r8), DIMENSION(Jmax,Imax), intent(in) :: Yout
integer, intent(in) :: Jmax, Imax
!
!  Interpolated Output
!
real(r8), DIMENSION(Jmax,Imax), intent(out) :: Fout
real(r8), DIMENSION(Jmax,Imax), intent(out) :: Jout
real(r8), DIMENSION(Jmax,Imax), intent(out) :: Iout
!
!  Local variable declarations:
!  
!  regrid (and others)
integer :: i, j
real(r8), dimension(Ny,Nx) :: angle
!real(r8), dimension(Jmax,Imax) :: Iout
!real(r8), dimension(Jmax,Imax) :: Jout
!  hindices
integer  :: np, mp, Imi, Jmi
real(r8) :: aa2, ang, bb2, diag2, dx, dy, phi
real(r8) :: xfac, xpp, yfac, ypp
!  linterp2d
integer  :: i1, i2, j1, j2
real(r8) :: p1, p2, q1, q2
!
!  Set input gridded data rotation angle.
!
DO i=1,Nx
   DO j=1,Ny
      angle(j,i)=0.0_r8
   END DO
END DO
!
!  Initialize local fractional coordinates arrays to avoid
!  deframentation.
!
Iout=0.0_r8
Jout=0.0_r8
!
!  CALL hindices:
!  Find fractional indices (Iout,Jout) of the grid cells in Finp
!  containing positions to intepolate.
!
!=======================================================================
!                                                                      !
!  Given any geographical locations Xpos and Ypos, this routine finds  !
!  the corresponding array cell indices (Ipos, Jpos) of gridded  data  !
!  Xgrd and Ygrd containing each requested location. This indices are  !
!  usually used for interpolation.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     LBi         I-dimension Lower bound of gridded data.             !
!     UBi         I-dimension Upper bound of gridded data.             !
!     LBj         J-dimension Lower bound of gridded data.             !
!     UBj         J-dimension Upper bound of gridded data.             !
!     Is          Starting gridded data I-index to search.             !
!     Ie          Ending   gridded data I-index to search.             !
!     Js          Starting gridded data J-index to search.             !
!     Je          Ending   gridded data J-index to search.             !
!     angler      Gridded data angle between X-axis and true EAST      !
!                   (radians).                                         !
!     Xgrd        Gridded data X-locations (usually, longitude).       !
!     Ygrd        Gridded data Y-locations (usually, latitude).        !
!     LBm         I-dimension Lower bound of requested locations.      !
!     UBm         I-dimension Upper bound of requested locations.      !
!     LBn         J-dimension Lower bound of requested locations.      !
!     UBn         J-dimension Upper bound of requested locations.      !
!     Ms          Starting requested locations I-index to search.      !
!     Me          Ending   requested locations I-index to search.      !
!     Ns          Starting requested locations J-index to search.      !
!     Ne          Ending   requested locations J-index to search.      !
!     Xpos        Requested X-locations to process (usually longitude).!
!     Ypos        Requested Y-locations to process (usually latitude). !
!     IJspv       Unbounded special value to assign.                   !
!     rectangular Logical switch indicating that gridded data has a    !
!                   plaid distribution.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-cell index containing locations in data. !
!     Jpos       Fractional J-cell index containing locations in data. !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
DO mp=1,Imax
   DO np=1,Jmax
!
!  The gridded data has a plaid distribution so the search is trivial.
!
      I_LOOP : DO i=1,Nx-1
        IF ((Xinp(1  ,i).le.Xout(np,mp)).and.                     &
        &            (Xinp(1,i+1).gt.Xout(np,mp))) THEN
           Imi=i
           EXIT I_LOOP
        END IF
      END DO I_LOOP

      J_LOOP : DO j=1,Ny-1
        IF ((Yinp(j,1  ).le.Yout(np,mp)).and.                     &
        &            (Yinp(j+1,1).gt.Yout(np,mp))) THEN
           Jmi=j
           EXIT J_LOOP
        END IF
      END DO J_LOOP
!
!  Knowing the correct cell, calculate the exact indices, accounting
!  for a possibly rotated grid.  If spherical, convert all positions
!  to meters first.
!
      yfac=Eradius*deg2rad
      xfac=yfac*COS(Yout(np,mp)*deg2rad)
      xpp=(Xout(np,mp)-Xinp(Jmi,Imi))*xfac
      ypp=(Yout(np,mp)-Yinp(Jmi,Imi))*yfac
!
!  Use Law of Cosines to get cell parallelogram "shear" angle.
!
      diag2=((Xinp(Jmi,Imi+1)-Xinp(Jmi+1,Imi))*xfac)**2+      &
      &            ((Yinp(Jmi,Imi+1)-Yinp(Jmi+1,Imi))*yfac)**2
      aa2=((Xinp(Jmi,Imi)-Xinp(Jmi,Imi+1))*xfac)**2+          &
      &          ((Yinp(Jmi,Imi)-Yinp(Jmi,Imi+1))*yfac)**2
      bb2=((Xinp(Jmi,Imi)-Xinp(Jmi+1,Imi))*xfac)**2+          &
      &          ((Yinp(Jmi,Imi)-Yinp(Jmi+1,Imi))*yfac)**2
      phi=ASIN((diag2-aa2-bb2)/(2.0_r8*SQRT(aa2*bb2)))
!
!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.
!
      ang=angle(Jmi,Imi)
      dx=xpp*COS(ang)+ypp*SIN(ang)
      dy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.
!
      dx=dx+dy*TAN(phi)
      dy=dy/COS(phi)
!
!  Scale with cell side lengths to translate into cell indices.
!
      dx=MIN(MAX(0.0_r8,dx/SQRT(aa2)),1.0_r8)
      dy=MIN(MAX(0.0_r8,dy/SQRT(bb2)),1.0_r8)
      Iout(np,mp)=REAL(Imi,r8)+dx ! Switching Imi and Jmi smooths field
      Jout(np,mp)=REAL(Jmi,r8)+dy

   END DO
END DO
!
!  CALL linterp2d:
!  Compute global interpolated field minimum and maximum values.
!  Notice that gridded data values are overwritten.
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field, Finp, this routine linearly interpolate !
!  to locations (Xout,Yout).  To facilitate the  interpolation  within !
!  any irregularly gridded 2D field,  the fractional grid cell indices !
!  (Iout,Jout) with respect Finp are needed at input.  Notice that the !
!  routine "hindices" can be used to compute these indices.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     LBi        I-dimension Lower bound of data to interpolate, Fout. !
!     UBi        I-dimension Upper bound of data to interpolate, Fout. !
!     LBj        J-dimension Lower bound of data to interpolate, Fout. !
!     UBj        J-dimension Upper bound of data to interpolate, Fout. !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Iout       I-fractional Xinp grid cell containing Xout.          !
!     Jout       J-fractional Yinp grid cell containing Yout.          !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!     Fmin       Interpolated field minimum value.                     !
!     Fmax       Interpolated field maximum value.                     !
!                                                                      !
!=======================================================================
!
! From mod_scalars
!
!-----------------------------------------------------------------------
!  Linearly interpolate requested field
!-----------------------------------------------------------------------
!
Amin=1.0E+35_r8
Amax=-1.0E+35_r8
! Itterates over CCS Grid dimensions
DO i=1,Imax
   DO j=1,Jmax

      i1=INT(Iout(j,i))
      i2=i1+1
      j1=INT(Jout(j,i))
      j2=j1+1
      IF (((1.le.i1).and.(i1.le.Nx)).and.                        &
      &        ((1.le.j1).and.(j1.le.Ny))) THEN
         p2=REAL(i2-i1,r8)*(Iout(j,i)-REAL(i1,r8))
         q2=REAL(j2-j1,r8)*(Jout(j,i)-REAL(j1,r8))
         p1=1.0_r8-p2
         q1=1.0_r8-q2
         Fout(j,i)=p1*q1*Finp(j1,i1)+                            &
         &                p1*q2*Finp(j2,i1)+                     & ! changed from p2*q1
         &                p2*q2*Finp(j2,i2)+                     &
         &                p2*q1*Finp(j1,i2)                        ! changed from p1*q2
         Amin=MIN(Amin,Fout(j,i))
         Amax=MAX(Amax,Fout(j,i))
      END IF
   END DO
END DO
RETURN
END SUBROUTINE

SUBROUTINE regrid_winds(Nx, Ny, Xinp, Yinp, FUinp, FVinp,       &
&                       AUmin, AUmax, AVmin, Avmax, Imax, Jmax, &
&                       angler, Xout, Yout, FUout, FVout)

! ---------------------------------------
! FOR ROTATING WINDS TO CURVILINEAR GRID 
! ---------------------------------------
IMPLICIT NONE
!
!  constants from mod_kinds and mod_scalars
!
integer, parameter  :: r8 = selected_real_kind(12,300)  ! 64-bit
!
!  Imported variable declarations:
!
!  ERAinterim lon, lat, data, dimensions
!
real(r8), DIMENSION(Ny,Nx), intent(in)    :: Xinp
real(r8), DIMENSION(Ny,Nx), intent(in)    :: Yinp
real(r8), DIMENSION(Ny,Nx), intent(inout) :: FUinp
real(r8), DIMENSION(Ny,Nx), intent(inout) :: FVinp
integer, intent(in) :: Nx, Ny
real(r8), intent(inout) :: AUmin, AUmax, AVmin, AVmax
!
!  CCS lon, lat, dimensions
!
real(r8), DIMENSION(Jmax,Imax), intent(in) :: angler
real(r8), DIMENSION(Jmax,Imax), intent(in) :: Xout
real(r8), DIMENSION(Jmax,Imax), intent(in) :: Yout
integer, intent(in) :: Imax, Jmax
!
!  Interpolated Output
!
real(r8), DIMENSION(Jmax,Imax) :: Jout
real(r8), DIMENSION(Jmax,Imax) :: Iout
real(r8), DIMENSION(Jmax,Imax), intent(out) :: FUout
real(r8), DIMENSION(Jmax,Imax), intent(out) :: FVout
!
!  Local variable declarations:
!  
integer  :: i, j
real(r8) :: cff1, cff2
!
! ---------------------------------------

CALL regrid_atmos(Nx, Ny, Xinp, Yinp, FUinp, AUmin, AUmax, Imax, Jmax, Xout, Yout, Jout, Iout, FUout)

CALL regrid_atmos(Nx, Ny, Xinp, Yinp, FVinp, AVmin, AVmax, Imax, Jmax, Xout, Yout, Jout, Iout, FVout)

DO i=1,Imax
   DO j=1,Jmax
      cff1=FUout(j,i)*COS(angler(j,i))+    &        
&          FVout(j,i)*SIN(angler(j,i))
      cff2=FVout(j,i)*COS(angler(j,i))-    &
&          FUout(j,i)*SIN(angler(j,i))
      FUout(j,i)=cff1
      FVout(j,i)=cff2
   END DO
END DO

RETURN
END SUBROUTINE
