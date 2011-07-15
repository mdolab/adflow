!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
!
!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidDissFluxScalar.f90                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-24-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE INVISCIDDISSFLUXSCALAR_CD()
  USE INPUTPHYSICS_SPATIAL_D
  USE ITERATION_SPATIAL_D
  USE CGNSGRID_SPATIAL_D
  USE INPUTDISCRETIZATION_SPATIAL_D
  USE CONSTANTS_SPATIAL_D
  USE BLOCKPOINTERS_SPATIAL_D
  USE FLOWVARREFSTATE_SPATIAL_D
  IMPLICIT NONE
!
!      ******************************************************************
!      *                                                                *
!      * inviscidDissFluxScalar computes the scalar artificial          *
!      * dissipation, see AIAA paper 81-1259, for a given block.        *
!      * Therefore it is assumed that the pointers in  blockPointers    *
!      * already point to the correct block.                            *
!      *                                                                *
!      ******************************************************************
!
!
!      Local parameter.
!
  REAL(kind=realtype), PARAMETER :: dssmax=0.25_realType
!
!      Local variables.
!
  INTEGER(kind=inttype) :: i, j, k, ind
  REAL(kind=realtype) :: sslim, rhoi
  REAL(kind=realtype) :: sfil, fis2, fis4
  REAL(kind=realtype) :: ppor, rrad, dis2, dis4
  REAL(kind=realtype) :: dss1, dss2, ddw, fs
  REAL(kind=realtype), DIMENSION(0:ib, 0:jb, 0:kb) :: ss
  REAL(kind=realtype) :: pwr1
  REAL(kind=realtype) :: min6
  REAL(kind=realtype) :: min5
  REAL(kind=realtype) :: min4
  REAL(kind=realtype) :: min3
  REAL(kind=realtype) :: min2
  REAL(kind=realtype) :: min1
  INTRINSIC MAX
  REAL(kind=realtype) :: x6
  REAL(kind=realtype) :: x5
  REAL(kind=realtype) :: x4
  REAL(kind=realtype) :: x3
  INTRINSIC ABS
  REAL(kind=realtype) :: x2
  REAL(kind=realtype) :: x1
  INTRINSIC MIN
  REAL(kind=realtype) :: y6
  REAL(kind=realtype) :: y5
  REAL(kind=realtype) :: y4
  REAL(kind=realtype) :: y3
  REAL(kind=realtype) :: y2
  REAL(kind=realtype) :: y1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! Check if rFil == 0. If so, the dissipative flux needs not to
! be computed.
  IF (rfil .EQ. zero) THEN
    RETURN
  ELSE
! Determine the variables used to compute the switch.
! For the inviscid case this is the pressure; for the viscous
! case it is the entropy.
    SELECT CASE  (equations) 
    CASE (eulerequations) 
! Inviscid case. Pressure switch is based on the pressure.
! Also set the value of sslim. To be fully consistent this
! must have the dimension of pressure and it is therefore
! set to a fraction of the free stream value.
      sslim = 0.001_realType*pinfcorr
! Copy the pressure in ss. Only fill the entries used in
! the discretization, i.e. ignore the corner halo's.
      DO k=0,kb
        DO j=2,jl
          DO i=2,il
            ss(i, j, k) = p(i, j, k)
          END DO
        END DO
      END DO
      DO k=2,kl
        DO j=2,jl
          ss(0, j, k) = p(0, j, k)
          ss(1, j, k) = p(1, j, k)
          ss(ie, j, k) = p(ie, j, k)
          ss(ib, j, k) = p(ib, j, k)
        END DO
      END DO
      DO k=2,kl
        DO i=2,il
          ss(i, 0, k) = p(i, 0, k)
          ss(i, 1, k) = p(i, 1, k)
          ss(i, je, k) = p(i, je, k)
          ss(i, jb, k) = p(i, jb, k)
        END DO
      END DO
    CASE (nsequations, ransequations) 
!===============================================================
! Viscous case. Pressure switch is based on the entropy.
! Also set the value of sslim. To be fully consistent this
! must have the dimension of entropy and it is therefore
! set to a fraction of the free stream value.
      pwr1 = rhoinf**gammainf
      sslim = 0.001_realType*pinfcorr/pwr1
! Store the entropy in ss. Only fill the entries used in
! the discretization, i.e. ignore the corner halo's.
      DO k=0,kb
        DO j=2,jl
          DO i=2,il
            pwr1 = w(i, j, k, irho)**gamma(i, j, k)
            ss(i, j, k) = p(i, j, k)/pwr1
          END DO
        END DO
      END DO
      DO k=2,kl
        DO j=2,jl
          pwr1 = w(0, j, k, irho)**gamma(0, j, k)
          ss(0, j, k) = p(0, j, k)/pwr1
          pwr1 = w(1, j, k, irho)**gamma(1, j, k)
          ss(1, j, k) = p(1, j, k)/pwr1
          pwr1 = w(ie, j, k, irho)**gamma(ie, j, k)
          ss(ie, j, k) = p(ie, j, k)/pwr1
          pwr1 = w(ib, j, k, irho)**gamma(ib, j, k)
          ss(ib, j, k) = p(ib, j, k)/pwr1
        END DO
      END DO
      DO k=2,kl
        DO i=2,il
          pwr1 = w(i, 0, k, irho)**gamma(i, 0, k)
          ss(i, 0, k) = p(i, 0, k)/pwr1
          pwr1 = w(i, 1, k, irho)**gamma(i, 1, k)
          ss(i, 1, k) = p(i, 1, k)/pwr1
          pwr1 = w(i, je, k, irho)**gamma(i, je, k)
          ss(i, je, k) = p(i, je, k)/pwr1
          pwr1 = w(i, jb, k, irho)**gamma(i, jb, k)
          ss(i, jb, k) = p(i, jb, k)/pwr1
        END DO
      END DO
    END SELECT
! Set a couple of constants for the scheme.
    fis2 = rfil*vis2
    fis4 = rfil*vis4
    sfil = one - rfil
! Replace the total energy by rho times the total enthalpy.
! In this way the numerical solution is total enthalpy preserving
! for the steady Euler equations. Also replace the velocities by
! the momentum. Only done for the entries used in the
! discretization, i.e. ignore the corner halo's.
    DO k=0,kb
      DO j=2,jl
        DO i=2,il
          w(i, j, k, ivx) = w(i, j, k, irho)*w(i, j, k, ivx)
          w(i, j, k, ivy) = w(i, j, k, irho)*w(i, j, k, ivy)
          w(i, j, k, ivz) = w(i, j, k, irho)*w(i, j, k, ivz)
          w(i, j, k, irhoe) = w(i, j, k, irhoe) + p(i, j, k)
        END DO
      END DO
    END DO
    DO k=2,kl
      DO j=2,jl
        w(0, j, k, ivx) = w(0, j, k, irho)*w(0, j, k, ivx)
        w(0, j, k, ivy) = w(0, j, k, irho)*w(0, j, k, ivy)
        w(0, j, k, ivz) = w(0, j, k, irho)*w(0, j, k, ivz)
        w(0, j, k, irhoe) = w(0, j, k, irhoe) + p(0, j, k)
        w(1, j, k, ivx) = w(1, j, k, irho)*w(1, j, k, ivx)
        w(1, j, k, ivy) = w(1, j, k, irho)*w(1, j, k, ivy)
        w(1, j, k, ivz) = w(1, j, k, irho)*w(1, j, k, ivz)
        w(1, j, k, irhoe) = w(1, j, k, irhoe) + p(1, j, k)
        w(ie, j, k, ivx) = w(ie, j, k, irho)*w(ie, j, k, ivx)
        w(ie, j, k, ivy) = w(ie, j, k, irho)*w(ie, j, k, ivy)
        w(ie, j, k, ivz) = w(ie, j, k, irho)*w(ie, j, k, ivz)
        w(ie, j, k, irhoe) = w(ie, j, k, irhoe) + p(ie, j, k)
        w(ib, j, k, ivx) = w(ib, j, k, irho)*w(ib, j, k, ivx)
        w(ib, j, k, ivy) = w(ib, j, k, irho)*w(ib, j, k, ivy)
        w(ib, j, k, ivz) = w(ib, j, k, irho)*w(ib, j, k, ivz)
        w(ib, j, k, irhoe) = w(ib, j, k, irhoe) + p(ib, j, k)
      END DO
    END DO
    DO k=2,kl
      DO i=2,il
        w(i, 0, k, ivx) = w(i, 0, k, irho)*w(i, 0, k, ivx)
        w(i, 0, k, ivy) = w(i, 0, k, irho)*w(i, 0, k, ivy)
        w(i, 0, k, ivz) = w(i, 0, k, irho)*w(i, 0, k, ivz)
        w(i, 0, k, irhoe) = w(i, 0, k, irhoe) + p(i, 0, k)
        w(i, 1, k, ivx) = w(i, 1, k, irho)*w(i, 1, k, ivx)
        w(i, 1, k, ivy) = w(i, 1, k, irho)*w(i, 1, k, ivy)
        w(i, 1, k, ivz) = w(i, 1, k, irho)*w(i, 1, k, ivz)
        w(i, 1, k, irhoe) = w(i, 1, k, irhoe) + p(i, 1, k)
        w(i, je, k, ivx) = w(i, je, k, irho)*w(i, je, k, ivx)
        w(i, je, k, ivy) = w(i, je, k, irho)*w(i, je, k, ivy)
        w(i, je, k, ivz) = w(i, je, k, irho)*w(i, je, k, ivz)
        w(i, je, k, irhoe) = w(i, je, k, irhoe) + p(i, je, k)
        w(i, jb, k, ivx) = w(i, jb, k, irho)*w(i, jb, k, ivx)
        w(i, jb, k, ivy) = w(i, jb, k, irho)*w(i, jb, k, ivy)
        w(i, jb, k, ivz) = w(i, jb, k, irho)*w(i, jb, k, ivz)
        w(i, jb, k, irhoe) = w(i, jb, k, irhoe) + p(i, jb, k)
      END DO
    END DO
! Initialize the dissipative residual to a certain times,
! possibly zero, the previously stored value. Owned cells
! only, because the halo values do not matter.
    DO k=2,kl
      DO j=2,jl
        DO i=2,il
          fw(i, j, k, irho) = sfil*fw(i, j, k, irho)
          fw(i, j, k, imx) = sfil*fw(i, j, k, imx)
          fw(i, j, k, imy) = sfil*fw(i, j, k, imy)
          fw(i, j, k, imz) = sfil*fw(i, j, k, imz)
          fw(i, j, k, irhoe) = sfil*fw(i, j, k, irhoe)
        END DO
      END DO
    END DO
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the i-direction.                         *
!      *                                                                *
!      ******************************************************************
!
    DO k=2,kl
      DO j=2,jl
        x1 = (ss(2, j, k)-two*ss(1, j, k)+ss(0, j, k))/(ss(2, j, k)+two*&
&          ss(1, j, k)+ss(0, j, k)+sslim)
        IF (x1 .GE. 0.) THEN
          dss1 = x1
        ELSE
          dss1 = -x1
        END IF
! Loop in i-direction.
        DO i=1,il
          x2 = (ss(i+2, j, k)-two*ss(i+1, j, k)+ss(i, j, k))/(ss(i+2, j&
&            , k)+two*ss(i+1, j, k)+ss(i, j, k)+sslim)
          IF (x2 .GE. 0.) THEN
            dss2 = x2
          ELSE
            dss2 = -x2
          END IF
! Compute the dissipation coefficients for this face.
          ppor = zero
          IF (pori(i, j, k) .EQ. normalflux) ppor = half
          rrad = ppor*(radi(i, j, k)+radi(i+1, j, k))
! Modification for FD Preconditioner Note: This lumping
! actually still results in a greater than 3 cell stencil
! in any direction. Since this seems to work slightly
! better than the dis2=sigma*fis4*rrad, we will just use
! a 5-cell stencil for doing the PC
          IF (lumpeddiss) THEN
            IF (dss1 .LT. dss2) THEN
              y1 = dss2
            ELSE
              y1 = dss1
            END IF
            IF (dssmax .GT. y1) THEN
              min1 = y1
            ELSE
              min1 = dssmax
            END IF
            dis2 = fis2*rrad*min1 + sigma*fis4*rrad
!dis2 = sigma*fis4*rrad 
            dis4 = 0.0
          ELSE
            IF (dss1 .LT. dss2) THEN
              y2 = dss2
            ELSE
              y2 = dss1
            END IF
            IF (dssmax .GT. y2) THEN
              min2 = y2
            ELSE
              min2 = dssmax
            END IF
            dis2 = fis2*rrad*min2
            dis4 = DIM_CD(fis4*rrad, dis2)
          END IF
! Compute and scatter the dissipative flux.
! Density. Store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw = w(i+1, j, k, irho) - w(i, j, k, irho)
          fs = dis2*ddw - dis4*(w(i+2, j, k, irho)-w(i-1, j, k, irho)-&
&            three*ddw)
          fw(i+1, j, k, irho) = fw(i+1, j, k, irho) + fs
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
          ind = indfamilyi(i, j, k)
          massflowfamilydiss(ind, spectralsol) = massflowfamilydiss(ind&
&            , spectralsol) - factfamilyi(i, j, k)*fs
! X-momentum.
          ddw = w(i+1, j, k, ivx) - w(i, j, k, ivx)
          fs = dis2*ddw - dis4*(w(i+2, j, k, ivx)-w(i-1, j, k, ivx)-&
&            three*ddw)
          fw(i+1, j, k, imx) = fw(i+1, j, k, imx) + fs
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! Y-momentum.
          ddw = w(i+1, j, k, ivy) - w(i, j, k, ivy)
          fs = dis2*ddw - dis4*(w(i+2, j, k, ivy)-w(i-1, j, k, ivy)-&
&            three*ddw)
          fw(i+1, j, k, imy) = fw(i+1, j, k, imy) + fs
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! Z-momentum.
          ddw = w(i+1, j, k, ivz) - w(i, j, k, ivz)
          fs = dis2*ddw - dis4*(w(i+2, j, k, ivz)-w(i-1, j, k, ivz)-&
&            three*ddw)
          fw(i+1, j, k, imz) = fw(i+1, j, k, imz) + fs
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! Energy.
          ddw = w(i+1, j, k, irhoe) - w(i, j, k, irhoe)
          fs = dis2*ddw - dis4*(w(i+2, j, k, irhoe)-w(i-1, j, k, irhoe)-&
&            three*ddw)
          fw(i+1, j, k, irhoe) = fw(i+1, j, k, irhoe) + fs
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
! Set dss1 to dss2 for the next face.
          dss1 = dss2
        END DO
      END DO
    END DO
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the j-direction.                         *
!      *                                                                *
!      ******************************************************************
!
    DO k=2,kl
      DO i=2,il
        x3 = (ss(i, 2, k)-two*ss(i, 1, k)+ss(i, 0, k))/(ss(i, 2, k)+two*&
&          ss(i, 1, k)+ss(i, 0, k)+sslim)
        IF (x3 .GE. 0.) THEN
          dss1 = x3
        ELSE
          dss1 = -x3
        END IF
! Loop in j-direction.
        DO j=1,jl
          x4 = (ss(i, j+2, k)-two*ss(i, j+1, k)+ss(i, j, k))/(ss(i, j+2&
&            , k)+two*ss(i, j+1, k)+ss(i, j, k)+sslim)
          IF (x4 .GE. 0.) THEN
            dss2 = x4
          ELSE
            dss2 = -x4
          END IF
! Compute the dissipation coefficients for this face.
          ppor = zero
          IF (porj(i, j, k) .EQ. normalflux) ppor = half
          rrad = ppor*(radj(i, j, k)+radj(i, j+1, k))
! Modification for FD Preconditioner
          IF (lumpeddiss) THEN
            IF (dss1 .LT. dss2) THEN
              y3 = dss2
            ELSE
              y3 = dss1
            END IF
            IF (dssmax .GT. y3) THEN
              min3 = y3
            ELSE
              min3 = dssmax
            END IF
            dis2 = fis2*rrad*min3 + sigma*fis4*rrad
!dis2 = sigma*fis4*rrad 
            dis4 = 0.0
          ELSE
            IF (dss1 .LT. dss2) THEN
              y4 = dss2
            ELSE
              y4 = dss1
            END IF
            IF (dssmax .GT. y4) THEN
              min4 = y4
            ELSE
              min4 = dssmax
            END IF
            dis2 = fis2*rrad*min4
            dis4 = DIM_CD(fis4*rrad, dis2)
          END IF
! Compute and scatter the dissipative flux.
! Density. Store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw = w(i, j+1, k, irho) - w(i, j, k, irho)
          fs = dis2*ddw - dis4*(w(i, j+2, k, irho)-w(i, j-1, k, irho)-&
&            three*ddw)
          fw(i, j+1, k, irho) = fw(i, j+1, k, irho) + fs
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
          ind = indfamilyj(i, j, k)
          massflowfamilydiss(ind, spectralsol) = massflowfamilydiss(ind&
&            , spectralsol) - factfamilyj(i, j, k)*fs
! X-momentum.
          ddw = w(i, j+1, k, ivx) - w(i, j, k, ivx)
          fs = dis2*ddw - dis4*(w(i, j+2, k, ivx)-w(i, j-1, k, ivx)-&
&            three*ddw)
          fw(i, j+1, k, imx) = fw(i, j+1, k, imx) + fs
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! Y-momentum.
          ddw = w(i, j+1, k, ivy) - w(i, j, k, ivy)
          fs = dis2*ddw - dis4*(w(i, j+2, k, ivy)-w(i, j-1, k, ivy)-&
&            three*ddw)
          fw(i, j+1, k, imy) = fw(i, j+1, k, imy) + fs
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! Z-momentum.
          ddw = w(i, j+1, k, ivz) - w(i, j, k, ivz)
          fs = dis2*ddw - dis4*(w(i, j+2, k, ivz)-w(i, j-1, k, ivz)-&
&            three*ddw)
          fw(i, j+1, k, imz) = fw(i, j+1, k, imz) + fs
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! Energy.
          ddw = w(i, j+1, k, irhoe) - w(i, j, k, irhoe)
          fs = dis2*ddw - dis4*(w(i, j+2, k, irhoe)-w(i, j-1, k, irhoe)-&
&            three*ddw)
          fw(i, j+1, k, irhoe) = fw(i, j+1, k, irhoe) + fs
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
! Set dss1 to dss2 for the next face.
          dss1 = dss2
        END DO
      END DO
    END DO
!
!      ******************************************************************
!      *                                                                *
!      * Dissipative fluxes in the k-direction.                         *
!      *                                                                *
!      ******************************************************************
!
    DO j=2,jl
      DO i=2,il
        x5 = (ss(i, j, 2)-two*ss(i, j, 1)+ss(i, j, 0))/(ss(i, j, 2)+two*&
&          ss(i, j, 1)+ss(i, j, 0)+sslim)
        IF (x5 .GE. 0.) THEN
          dss1 = x5
        ELSE
          dss1 = -x5
        END IF
! Loop in k-direction.
        DO k=1,kl
          x6 = (ss(i, j, k+2)-two*ss(i, j, k+1)+ss(i, j, k))/(ss(i, j, k&
&            +2)+two*ss(i, j, k+1)+ss(i, j, k)+sslim)
          IF (x6 .GE. 0.) THEN
            dss2 = x6
          ELSE
            dss2 = -x6
          END IF
! Compute the dissipation coefficients for this face.
          ppor = zero
          IF (pork(i, j, k) .EQ. normalflux) ppor = half
          rrad = ppor*(radk(i, j, k)+radk(i, j, k+1))
! Modification for FD Preconditioner
          IF (lumpeddiss) THEN
            IF (dss1 .LT. dss2) THEN
              y5 = dss2
            ELSE
              y5 = dss1
            END IF
            IF (dssmax .GT. y5) THEN
              min5 = y5
            ELSE
              min5 = dssmax
            END IF
            dis2 = fis2*rrad*min5 + sigma*fis4*rrad
!dis2 = sigma*fis4*rrad 
            dis4 = 0.0
          ELSE
            IF (dss1 .LT. dss2) THEN
              y6 = dss2
            ELSE
              y6 = dss1
            END IF
            IF (dssmax .GT. y6) THEN
              min6 = y6
            ELSE
              min6 = dssmax
            END IF
            dis2 = fis2*rrad*min6
            dis4 = DIM_CD(fis4*rrad, dis2)
          END IF
! Compute and scatter the dissipative flux.
! Density. Store it in the mass flow of the
! appropriate sliding mesh interface.
          ddw = w(i, j, k+1, irho) - w(i, j, k, irho)
          fs = dis2*ddw - dis4*(w(i, j, k+2, irho)-w(i, j, k-1, irho)-&
&            three*ddw)
          fw(i, j, k+1, irho) = fw(i, j, k+1, irho) + fs
          fw(i, j, k, irho) = fw(i, j, k, irho) - fs
          ind = indfamilyk(i, j, k)
          massflowfamilydiss(ind, spectralsol) = massflowfamilydiss(ind&
&            , spectralsol) - factfamilyk(i, j, k)*fs
! X-momentum.
          ddw = w(i, j, k+1, ivx) - w(i, j, k, ivx)
          fs = dis2*ddw - dis4*(w(i, j, k+2, ivx)-w(i, j, k-1, ivx)-&
&            three*ddw)
          fw(i, j, k+1, imx) = fw(i, j, k+1, imx) + fs
          fw(i, j, k, imx) = fw(i, j, k, imx) - fs
! Y-momentum.
          ddw = w(i, j, k+1, ivy) - w(i, j, k, ivy)
          fs = dis2*ddw - dis4*(w(i, j, k+2, ivy)-w(i, j, k-1, ivy)-&
&            three*ddw)
          fw(i, j, k+1, imy) = fw(i, j, k+1, imy) + fs
          fw(i, j, k, imy) = fw(i, j, k, imy) - fs
! Z-momentum.
          ddw = w(i, j, k+1, ivz) - w(i, j, k, ivz)
          fs = dis2*ddw - dis4*(w(i, j, k+2, ivz)-w(i, j, k-1, ivz)-&
&            three*ddw)
          fw(i, j, k+1, imz) = fw(i, j, k+1, imz) + fs
          fw(i, j, k, imz) = fw(i, j, k, imz) - fs
! Energy.
          ddw = w(i, j, k+1, irhoe) - w(i, j, k, irhoe)
          fs = dis2*ddw - dis4*(w(i, j, k+2, irhoe)-w(i, j, k-1, irhoe)-&
&            three*ddw)
          fw(i, j, k+1, irhoe) = fw(i, j, k+1, irhoe) + fs
          fw(i, j, k, irhoe) = fw(i, j, k, irhoe) - fs
! Set dss1 to dss2 for the next face.
          dss1 = dss2
        END DO
      END DO
    END DO
! Replace rho times the total enthalpy by the total energy and
! store the velocities again instead of the momentum. Only for
! those entries that have been altered, i.e. ignore the
! corner halo's.
    DO k=0,kb
      DO j=2,jl
        DO i=2,il
          rhoi = one/w(i, j, k, irho)
          w(i, j, k, ivx) = w(i, j, k, ivx)*rhoi
          w(i, j, k, ivy) = w(i, j, k, ivy)*rhoi
          w(i, j, k, ivz) = w(i, j, k, ivz)*rhoi
          w(i, j, k, irhoe) = w(i, j, k, irhoe) - p(i, j, k)
        END DO
      END DO
    END DO
    DO k=2,kl
      DO j=2,jl
        rhoi = one/w(0, j, k, irho)
        w(0, j, k, ivx) = w(0, j, k, ivx)*rhoi
        w(0, j, k, ivy) = w(0, j, k, ivy)*rhoi
        w(0, j, k, ivz) = w(0, j, k, ivz)*rhoi
        w(0, j, k, irhoe) = w(0, j, k, irhoe) - p(0, j, k)
        rhoi = one/w(1, j, k, irho)
        w(1, j, k, ivx) = w(1, j, k, ivx)*rhoi
        w(1, j, k, ivy) = w(1, j, k, ivy)*rhoi
        w(1, j, k, ivz) = w(1, j, k, ivz)*rhoi
        w(1, j, k, irhoe) = w(1, j, k, irhoe) - p(1, j, k)
        rhoi = one/w(ie, j, k, irho)
        w(ie, j, k, ivx) = w(ie, j, k, ivx)*rhoi
        w(ie, j, k, ivy) = w(ie, j, k, ivy)*rhoi
        w(ie, j, k, ivz) = w(ie, j, k, ivz)*rhoi
        w(ie, j, k, irhoe) = w(ie, j, k, irhoe) - p(ie, j, k)
        rhoi = one/w(ib, j, k, irho)
        w(ib, j, k, ivx) = w(ib, j, k, ivx)*rhoi
        w(ib, j, k, ivy) = w(ib, j, k, ivy)*rhoi
        w(ib, j, k, ivz) = w(ib, j, k, ivz)*rhoi
        w(ib, j, k, irhoe) = w(ib, j, k, irhoe) - p(ib, j, k)
      END DO
    END DO
    DO k=2,kl
      DO i=2,il
        rhoi = one/w(i, 0, k, irho)
        w(i, 0, k, ivx) = w(i, 0, k, ivx)*rhoi
        w(i, 0, k, ivy) = w(i, 0, k, ivy)*rhoi
        w(i, 0, k, ivz) = w(i, 0, k, ivz)*rhoi
        w(i, 0, k, irhoe) = w(i, 0, k, irhoe) - p(i, 0, k)
        rhoi = one/w(i, 1, k, irho)
        w(i, 1, k, ivx) = w(i, 1, k, ivx)*rhoi
        w(i, 1, k, ivy) = w(i, 1, k, ivy)*rhoi
        w(i, 1, k, ivz) = w(i, 1, k, ivz)*rhoi
        w(i, 1, k, irhoe) = w(i, 1, k, irhoe) - p(i, 1, k)
        rhoi = one/w(i, je, k, irho)
        w(i, je, k, ivx) = w(i, je, k, ivx)*rhoi
        w(i, je, k, ivy) = w(i, je, k, ivy)*rhoi
        w(i, je, k, ivz) = w(i, je, k, ivz)*rhoi
        w(i, je, k, irhoe) = w(i, je, k, irhoe) - p(i, je, k)
        rhoi = one/w(i, jb, k, irho)
        w(i, jb, k, ivx) = w(i, jb, k, ivx)*rhoi
        w(i, jb, k, ivy) = w(i, jb, k, ivy)*rhoi
        w(i, jb, k, ivz) = w(i, jb, k, ivz)*rhoi
        w(i, jb, k, irhoe) = w(i, jb, k, irhoe) - p(i, jb, k)
      END DO
    END DO
  END IF
END SUBROUTINE INVISCIDDISSFLUXSCALAR_CD
