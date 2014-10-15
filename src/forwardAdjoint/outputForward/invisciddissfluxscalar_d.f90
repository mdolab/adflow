   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of invisciddissfluxscalar in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *w
   !   with respect to varying inputs: *p *w
   !   Plus diff mem management of: p:in w:in
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
   SUBROUTINE INVISCIDDISSFLUXSCALAR_D()
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
   USE BLOCKPOINTERS_D
   USE CGNSGRID
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE INPUTDISCRETIZATION
   USE INPUTPHYSICS
   USE ITERATION
   IMPLICIT NONE
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: dssmax=0.25_realType
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k, ind
   REAL(kind=realtype) :: sslim, rhoi
   REAL(kind=realtype) :: rhoid
   REAL(kind=realtype) :: sfil, fis2, fis4
   REAL(kind=realtype) :: ppor, rrad, dis2, dis4
   REAL(kind=realtype) :: dss1, dss2, ddw, fs
   REAL(kind=realtype), DIMENSION(0:ib, 0:jb, 0:kb) :: ss
   INTRINSIC ABS
   INTRINSIC MAX
   INTRINSIC MIN
   REAL(kind=realtype) :: DIM
   REAL(kind=realtype) :: DIM_D
   REAL(kind=realtype) :: pwr1
   REAL(kind=realtype) :: min6
   REAL(kind=realtype) :: min5
   REAL(kind=realtype) :: min4
   REAL(kind=realtype) :: min3
   REAL(kind=realtype) :: min2
   REAL(kind=realtype) :: min1
   REAL(kind=realtype) :: x6
   REAL(kind=realtype) :: x5
   REAL(kind=realtype) :: x4
   REAL(kind=realtype) :: x3
   REAL(kind=realtype) :: x2
   REAL(kind=realtype) :: x1
   REAL(kind=realtype) :: abs0
   REAL(kind=realtype) :: y6
   REAL(kind=realtype) :: y5
   REAL(kind=realtype) :: y4
   REAL(kind=realtype) :: y3
   REAL(kind=realtype) :: y2
   REAL(kind=realtype) :: y1
   IF (rfil .GE. 0.) THEN
   abs0 = rfil
   ELSE
   abs0 = -rfil
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Check if rFil == 0. If so, the dissipative flux needs not to
   ! be computed.
   IF (abs0 .LT. thresholdreal) THEN
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
   wd(i, j, k, ivx) = wd(i, j, k, irho)*w(i, j, k, ivx) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivx)
   w(i, j, k, ivx) = w(i, j, k, irho)*w(i, j, k, ivx)
   wd(i, j, k, ivy) = wd(i, j, k, irho)*w(i, j, k, ivy) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivy)
   w(i, j, k, ivy) = w(i, j, k, irho)*w(i, j, k, ivy)
   wd(i, j, k, ivz) = wd(i, j, k, irho)*w(i, j, k, ivz) + w(i, j&
   &           , k, irho)*wd(i, j, k, ivz)
   w(i, j, k, ivz) = w(i, j, k, irho)*w(i, j, k, ivz)
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) + pd(i, j, k)
   w(i, j, k, irhoe) = w(i, j, k, irhoe) + p(i, j, k)
   END DO
   END DO
   END DO
   DO k=2,kl
   DO j=2,jl
   wd(0, j, k, ivx) = wd(0, j, k, irho)*w(0, j, k, ivx) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivx)
   w(0, j, k, ivx) = w(0, j, k, irho)*w(0, j, k, ivx)
   wd(0, j, k, ivy) = wd(0, j, k, irho)*w(0, j, k, ivy) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivy)
   w(0, j, k, ivy) = w(0, j, k, irho)*w(0, j, k, ivy)
   wd(0, j, k, ivz) = wd(0, j, k, irho)*w(0, j, k, ivz) + w(0, j, k&
   &         , irho)*wd(0, j, k, ivz)
   w(0, j, k, ivz) = w(0, j, k, irho)*w(0, j, k, ivz)
   wd(0, j, k, irhoe) = wd(0, j, k, irhoe) + pd(0, j, k)
   w(0, j, k, irhoe) = w(0, j, k, irhoe) + p(0, j, k)
   wd(1, j, k, ivx) = wd(1, j, k, irho)*w(1, j, k, ivx) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivx)
   w(1, j, k, ivx) = w(1, j, k, irho)*w(1, j, k, ivx)
   wd(1, j, k, ivy) = wd(1, j, k, irho)*w(1, j, k, ivy) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivy)
   w(1, j, k, ivy) = w(1, j, k, irho)*w(1, j, k, ivy)
   wd(1, j, k, ivz) = wd(1, j, k, irho)*w(1, j, k, ivz) + w(1, j, k&
   &         , irho)*wd(1, j, k, ivz)
   w(1, j, k, ivz) = w(1, j, k, irho)*w(1, j, k, ivz)
   wd(1, j, k, irhoe) = wd(1, j, k, irhoe) + pd(1, j, k)
   w(1, j, k, irhoe) = w(1, j, k, irhoe) + p(1, j, k)
   wd(ie, j, k, ivx) = wd(ie, j, k, irho)*w(ie, j, k, ivx) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivx)
   w(ie, j, k, ivx) = w(ie, j, k, irho)*w(ie, j, k, ivx)
   wd(ie, j, k, ivy) = wd(ie, j, k, irho)*w(ie, j, k, ivy) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivy)
   w(ie, j, k, ivy) = w(ie, j, k, irho)*w(ie, j, k, ivy)
   wd(ie, j, k, ivz) = wd(ie, j, k, irho)*w(ie, j, k, ivz) + w(ie, &
   &         j, k, irho)*wd(ie, j, k, ivz)
   w(ie, j, k, ivz) = w(ie, j, k, irho)*w(ie, j, k, ivz)
   wd(ie, j, k, irhoe) = wd(ie, j, k, irhoe) + pd(ie, j, k)
   w(ie, j, k, irhoe) = w(ie, j, k, irhoe) + p(ie, j, k)
   wd(ib, j, k, ivx) = wd(ib, j, k, irho)*w(ib, j, k, ivx) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivx)
   w(ib, j, k, ivx) = w(ib, j, k, irho)*w(ib, j, k, ivx)
   wd(ib, j, k, ivy) = wd(ib, j, k, irho)*w(ib, j, k, ivy) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivy)
   w(ib, j, k, ivy) = w(ib, j, k, irho)*w(ib, j, k, ivy)
   wd(ib, j, k, ivz) = wd(ib, j, k, irho)*w(ib, j, k, ivz) + w(ib, &
   &         j, k, irho)*wd(ib, j, k, ivz)
   w(ib, j, k, ivz) = w(ib, j, k, irho)*w(ib, j, k, ivz)
   wd(ib, j, k, irhoe) = wd(ib, j, k, irhoe) + pd(ib, j, k)
   w(ib, j, k, irhoe) = w(ib, j, k, irhoe) + p(ib, j, k)
   END DO
   END DO
   DO k=2,kl
   DO i=2,il
   wd(i, 0, k, ivx) = wd(i, 0, k, irho)*w(i, 0, k, ivx) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivx)
   w(i, 0, k, ivx) = w(i, 0, k, irho)*w(i, 0, k, ivx)
   wd(i, 0, k, ivy) = wd(i, 0, k, irho)*w(i, 0, k, ivy) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivy)
   w(i, 0, k, ivy) = w(i, 0, k, irho)*w(i, 0, k, ivy)
   wd(i, 0, k, ivz) = wd(i, 0, k, irho)*w(i, 0, k, ivz) + w(i, 0, k&
   &         , irho)*wd(i, 0, k, ivz)
   w(i, 0, k, ivz) = w(i, 0, k, irho)*w(i, 0, k, ivz)
   wd(i, 0, k, irhoe) = wd(i, 0, k, irhoe) + pd(i, 0, k)
   w(i, 0, k, irhoe) = w(i, 0, k, irhoe) + p(i, 0, k)
   wd(i, 1, k, ivx) = wd(i, 1, k, irho)*w(i, 1, k, ivx) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivx)
   w(i, 1, k, ivx) = w(i, 1, k, irho)*w(i, 1, k, ivx)
   wd(i, 1, k, ivy) = wd(i, 1, k, irho)*w(i, 1, k, ivy) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivy)
   w(i, 1, k, ivy) = w(i, 1, k, irho)*w(i, 1, k, ivy)
   wd(i, 1, k, ivz) = wd(i, 1, k, irho)*w(i, 1, k, ivz) + w(i, 1, k&
   &         , irho)*wd(i, 1, k, ivz)
   w(i, 1, k, ivz) = w(i, 1, k, irho)*w(i, 1, k, ivz)
   wd(i, 1, k, irhoe) = wd(i, 1, k, irhoe) + pd(i, 1, k)
   w(i, 1, k, irhoe) = w(i, 1, k, irhoe) + p(i, 1, k)
   wd(i, je, k, ivx) = wd(i, je, k, irho)*w(i, je, k, ivx) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivx)
   w(i, je, k, ivx) = w(i, je, k, irho)*w(i, je, k, ivx)
   wd(i, je, k, ivy) = wd(i, je, k, irho)*w(i, je, k, ivy) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivy)
   w(i, je, k, ivy) = w(i, je, k, irho)*w(i, je, k, ivy)
   wd(i, je, k, ivz) = wd(i, je, k, irho)*w(i, je, k, ivz) + w(i, &
   &         je, k, irho)*wd(i, je, k, ivz)
   w(i, je, k, ivz) = w(i, je, k, irho)*w(i, je, k, ivz)
   wd(i, je, k, irhoe) = wd(i, je, k, irhoe) + pd(i, je, k)
   w(i, je, k, irhoe) = w(i, je, k, irhoe) + p(i, je, k)
   wd(i, jb, k, ivx) = wd(i, jb, k, irho)*w(i, jb, k, ivx) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivx)
   w(i, jb, k, ivx) = w(i, jb, k, irho)*w(i, jb, k, ivx)
   wd(i, jb, k, ivy) = wd(i, jb, k, irho)*w(i, jb, k, ivy) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivy)
   w(i, jb, k, ivy) = w(i, jb, k, irho)*w(i, jb, k, ivy)
   wd(i, jb, k, ivz) = wd(i, jb, k, irho)*w(i, jb, k, ivz) + w(i, &
   &         jb, k, irho)*wd(i, jb, k, ivz)
   w(i, jb, k, ivz) = w(i, jb, k, irho)*w(i, jb, k, ivz)
   wd(i, jb, k, irhoe) = wd(i, jb, k, irhoe) + pd(i, jb, k)
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
   &         ss(1, j, k)+ss(0, j, k)+sslim)
   IF (x1 .GE. 0.) THEN
   dss1 = x1
   ELSE
   dss1 = -x1
   END IF
   ! Loop in i-direction.
   DO i=1,il
   x2 = (ss(i+2, j, k)-two*ss(i+1, j, k)+ss(i, j, k))/(ss(i+2, j&
   &           , k)+two*ss(i+1, j, k)+ss(i, j, k)+sslim)
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
   dis4 = DIM(fis4*rrad, dis2)
   END IF
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddw = w(i+1, j, k, irho) - w(i, j, k, irho)
   fs = dis2*ddw - dis4*(w(i+2, j, k, irho)-w(i-1, j, k, irho)-&
   &           three*ddw)
   fw(i+1, j, k, irho) = fw(i+1, j, k, irho) + fs
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddw = w(i+1, j, k, ivx) - w(i, j, k, ivx)
   fs = dis2*ddw - dis4*(w(i+2, j, k, ivx)-w(i-1, j, k, ivx)-&
   &           three*ddw)
   fw(i+1, j, k, imx) = fw(i+1, j, k, imx) + fs
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddw = w(i+1, j, k, ivy) - w(i, j, k, ivy)
   fs = dis2*ddw - dis4*(w(i+2, j, k, ivy)-w(i-1, j, k, ivy)-&
   &           three*ddw)
   fw(i+1, j, k, imy) = fw(i+1, j, k, imy) + fs
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddw = w(i+1, j, k, ivz) - w(i, j, k, ivz)
   fs = dis2*ddw - dis4*(w(i+2, j, k, ivz)-w(i-1, j, k, ivz)-&
   &           three*ddw)
   fw(i+1, j, k, imz) = fw(i+1, j, k, imz) + fs
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddw = w(i+1, j, k, irhoe) - w(i, j, k, irhoe)
   fs = dis2*ddw - dis4*(w(i+2, j, k, irhoe)-w(i-1, j, k, irhoe)-&
   &           three*ddw)
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
   &         ss(i, 1, k)+ss(i, 0, k)+sslim)
   IF (x3 .GE. 0.) THEN
   dss1 = x3
   ELSE
   dss1 = -x3
   END IF
   ! Loop in j-direction.
   DO j=1,jl
   x4 = (ss(i, j+2, k)-two*ss(i, j+1, k)+ss(i, j, k))/(ss(i, j+2&
   &           , k)+two*ss(i, j+1, k)+ss(i, j, k)+sslim)
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
   dis4 = DIM(fis4*rrad, dis2)
   END IF
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddw = w(i, j+1, k, irho) - w(i, j, k, irho)
   fs = dis2*ddw - dis4*(w(i, j+2, k, irho)-w(i, j-1, k, irho)-&
   &           three*ddw)
   fw(i, j+1, k, irho) = fw(i, j+1, k, irho) + fs
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddw = w(i, j+1, k, ivx) - w(i, j, k, ivx)
   fs = dis2*ddw - dis4*(w(i, j+2, k, ivx)-w(i, j-1, k, ivx)-&
   &           three*ddw)
   fw(i, j+1, k, imx) = fw(i, j+1, k, imx) + fs
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddw = w(i, j+1, k, ivy) - w(i, j, k, ivy)
   fs = dis2*ddw - dis4*(w(i, j+2, k, ivy)-w(i, j-1, k, ivy)-&
   &           three*ddw)
   fw(i, j+1, k, imy) = fw(i, j+1, k, imy) + fs
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddw = w(i, j+1, k, ivz) - w(i, j, k, ivz)
   fs = dis2*ddw - dis4*(w(i, j+2, k, ivz)-w(i, j-1, k, ivz)-&
   &           three*ddw)
   fw(i, j+1, k, imz) = fw(i, j+1, k, imz) + fs
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddw = w(i, j+1, k, irhoe) - w(i, j, k, irhoe)
   fs = dis2*ddw - dis4*(w(i, j+2, k, irhoe)-w(i, j-1, k, irhoe)-&
   &           three*ddw)
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
   &         ss(i, j, 1)+ss(i, j, 0)+sslim)
   IF (x5 .GE. 0.) THEN
   dss1 = x5
   ELSE
   dss1 = -x5
   END IF
   ! Loop in k-direction.
   DO k=1,kl
   x6 = (ss(i, j, k+2)-two*ss(i, j, k+1)+ss(i, j, k))/(ss(i, j, k&
   &           +2)+two*ss(i, j, k+1)+ss(i, j, k)+sslim)
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
   dis4 = DIM(fis4*rrad, dis2)
   END IF
   ! Compute and scatter the dissipative flux.
   ! Density. Store it in the mass flow of the
   ! appropriate sliding mesh interface.
   ddw = w(i, j, k+1, irho) - w(i, j, k, irho)
   fs = dis2*ddw - dis4*(w(i, j, k+2, irho)-w(i, j, k-1, irho)-&
   &           three*ddw)
   fw(i, j, k+1, irho) = fw(i, j, k+1, irho) + fs
   fw(i, j, k, irho) = fw(i, j, k, irho) - fs
   ! X-momentum.
   ddw = w(i, j, k+1, ivx) - w(i, j, k, ivx)
   fs = dis2*ddw - dis4*(w(i, j, k+2, ivx)-w(i, j, k-1, ivx)-&
   &           three*ddw)
   fw(i, j, k+1, imx) = fw(i, j, k+1, imx) + fs
   fw(i, j, k, imx) = fw(i, j, k, imx) - fs
   ! Y-momentum.
   ddw = w(i, j, k+1, ivy) - w(i, j, k, ivy)
   fs = dis2*ddw - dis4*(w(i, j, k+2, ivy)-w(i, j, k-1, ivy)-&
   &           three*ddw)
   fw(i, j, k+1, imy) = fw(i, j, k+1, imy) + fs
   fw(i, j, k, imy) = fw(i, j, k, imy) - fs
   ! Z-momentum.
   ddw = w(i, j, k+1, ivz) - w(i, j, k, ivz)
   fs = dis2*ddw - dis4*(w(i, j, k+2, ivz)-w(i, j, k-1, ivz)-&
   &           three*ddw)
   fw(i, j, k+1, imz) = fw(i, j, k+1, imz) + fs
   fw(i, j, k, imz) = fw(i, j, k, imz) - fs
   ! Energy.
   ddw = w(i, j, k+1, irhoe) - w(i, j, k, irhoe)
   fs = dis2*ddw - dis4*(w(i, j, k+2, irhoe)-w(i, j, k-1, irhoe)-&
   &           three*ddw)
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
   rhoid = -(one*wd(i, j, k, irho)/w(i, j, k, irho)**2)
   rhoi = one/w(i, j, k, irho)
   wd(i, j, k, ivx) = wd(i, j, k, ivx)*rhoi + w(i, j, k, ivx)*&
   &           rhoid
   w(i, j, k, ivx) = w(i, j, k, ivx)*rhoi
   wd(i, j, k, ivy) = wd(i, j, k, ivy)*rhoi + w(i, j, k, ivy)*&
   &           rhoid
   w(i, j, k, ivy) = w(i, j, k, ivy)*rhoi
   wd(i, j, k, ivz) = wd(i, j, k, ivz)*rhoi + w(i, j, k, ivz)*&
   &           rhoid
   w(i, j, k, ivz) = w(i, j, k, ivz)*rhoi
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - pd(i, j, k)
   w(i, j, k, irhoe) = w(i, j, k, irhoe) - p(i, j, k)
   END DO
   END DO
   END DO
   DO k=2,kl
   DO j=2,jl
   rhoid = -(one*wd(0, j, k, irho)/w(0, j, k, irho)**2)
   rhoi = one/w(0, j, k, irho)
   wd(0, j, k, ivx) = wd(0, j, k, ivx)*rhoi + w(0, j, k, ivx)*rhoid
   w(0, j, k, ivx) = w(0, j, k, ivx)*rhoi
   wd(0, j, k, ivy) = wd(0, j, k, ivy)*rhoi + w(0, j, k, ivy)*rhoid
   w(0, j, k, ivy) = w(0, j, k, ivy)*rhoi
   wd(0, j, k, ivz) = wd(0, j, k, ivz)*rhoi + w(0, j, k, ivz)*rhoid
   w(0, j, k, ivz) = w(0, j, k, ivz)*rhoi
   wd(0, j, k, irhoe) = wd(0, j, k, irhoe) - pd(0, j, k)
   w(0, j, k, irhoe) = w(0, j, k, irhoe) - p(0, j, k)
   rhoid = -(one*wd(1, j, k, irho)/w(1, j, k, irho)**2)
   rhoi = one/w(1, j, k, irho)
   wd(1, j, k, ivx) = wd(1, j, k, ivx)*rhoi + w(1, j, k, ivx)*rhoid
   w(1, j, k, ivx) = w(1, j, k, ivx)*rhoi
   wd(1, j, k, ivy) = wd(1, j, k, ivy)*rhoi + w(1, j, k, ivy)*rhoid
   w(1, j, k, ivy) = w(1, j, k, ivy)*rhoi
   wd(1, j, k, ivz) = wd(1, j, k, ivz)*rhoi + w(1, j, k, ivz)*rhoid
   w(1, j, k, ivz) = w(1, j, k, ivz)*rhoi
   wd(1, j, k, irhoe) = wd(1, j, k, irhoe) - pd(1, j, k)
   w(1, j, k, irhoe) = w(1, j, k, irhoe) - p(1, j, k)
   rhoid = -(one*wd(ie, j, k, irho)/w(ie, j, k, irho)**2)
   rhoi = one/w(ie, j, k, irho)
   wd(ie, j, k, ivx) = wd(ie, j, k, ivx)*rhoi + w(ie, j, k, ivx)*&
   &         rhoid
   w(ie, j, k, ivx) = w(ie, j, k, ivx)*rhoi
   wd(ie, j, k, ivy) = wd(ie, j, k, ivy)*rhoi + w(ie, j, k, ivy)*&
   &         rhoid
   w(ie, j, k, ivy) = w(ie, j, k, ivy)*rhoi
   wd(ie, j, k, ivz) = wd(ie, j, k, ivz)*rhoi + w(ie, j, k, ivz)*&
   &         rhoid
   w(ie, j, k, ivz) = w(ie, j, k, ivz)*rhoi
   wd(ie, j, k, irhoe) = wd(ie, j, k, irhoe) - pd(ie, j, k)
   w(ie, j, k, irhoe) = w(ie, j, k, irhoe) - p(ie, j, k)
   rhoid = -(one*wd(ib, j, k, irho)/w(ib, j, k, irho)**2)
   rhoi = one/w(ib, j, k, irho)
   wd(ib, j, k, ivx) = wd(ib, j, k, ivx)*rhoi + w(ib, j, k, ivx)*&
   &         rhoid
   w(ib, j, k, ivx) = w(ib, j, k, ivx)*rhoi
   wd(ib, j, k, ivy) = wd(ib, j, k, ivy)*rhoi + w(ib, j, k, ivy)*&
   &         rhoid
   w(ib, j, k, ivy) = w(ib, j, k, ivy)*rhoi
   wd(ib, j, k, ivz) = wd(ib, j, k, ivz)*rhoi + w(ib, j, k, ivz)*&
   &         rhoid
   w(ib, j, k, ivz) = w(ib, j, k, ivz)*rhoi
   wd(ib, j, k, irhoe) = wd(ib, j, k, irhoe) - pd(ib, j, k)
   w(ib, j, k, irhoe) = w(ib, j, k, irhoe) - p(ib, j, k)
   END DO
   END DO
   DO k=2,kl
   DO i=2,il
   rhoid = -(one*wd(i, 0, k, irho)/w(i, 0, k, irho)**2)
   rhoi = one/w(i, 0, k, irho)
   wd(i, 0, k, ivx) = wd(i, 0, k, ivx)*rhoi + w(i, 0, k, ivx)*rhoid
   w(i, 0, k, ivx) = w(i, 0, k, ivx)*rhoi
   wd(i, 0, k, ivy) = wd(i, 0, k, ivy)*rhoi + w(i, 0, k, ivy)*rhoid
   w(i, 0, k, ivy) = w(i, 0, k, ivy)*rhoi
   wd(i, 0, k, ivz) = wd(i, 0, k, ivz)*rhoi + w(i, 0, k, ivz)*rhoid
   w(i, 0, k, ivz) = w(i, 0, k, ivz)*rhoi
   wd(i, 0, k, irhoe) = wd(i, 0, k, irhoe) - pd(i, 0, k)
   w(i, 0, k, irhoe) = w(i, 0, k, irhoe) - p(i, 0, k)
   rhoid = -(one*wd(i, 1, k, irho)/w(i, 1, k, irho)**2)
   rhoi = one/w(i, 1, k, irho)
   wd(i, 1, k, ivx) = wd(i, 1, k, ivx)*rhoi + w(i, 1, k, ivx)*rhoid
   w(i, 1, k, ivx) = w(i, 1, k, ivx)*rhoi
   wd(i, 1, k, ivy) = wd(i, 1, k, ivy)*rhoi + w(i, 1, k, ivy)*rhoid
   w(i, 1, k, ivy) = w(i, 1, k, ivy)*rhoi
   wd(i, 1, k, ivz) = wd(i, 1, k, ivz)*rhoi + w(i, 1, k, ivz)*rhoid
   w(i, 1, k, ivz) = w(i, 1, k, ivz)*rhoi
   wd(i, 1, k, irhoe) = wd(i, 1, k, irhoe) - pd(i, 1, k)
   w(i, 1, k, irhoe) = w(i, 1, k, irhoe) - p(i, 1, k)
   rhoid = -(one*wd(i, je, k, irho)/w(i, je, k, irho)**2)
   rhoi = one/w(i, je, k, irho)
   wd(i, je, k, ivx) = wd(i, je, k, ivx)*rhoi + w(i, je, k, ivx)*&
   &         rhoid
   w(i, je, k, ivx) = w(i, je, k, ivx)*rhoi
   wd(i, je, k, ivy) = wd(i, je, k, ivy)*rhoi + w(i, je, k, ivy)*&
   &         rhoid
   w(i, je, k, ivy) = w(i, je, k, ivy)*rhoi
   wd(i, je, k, ivz) = wd(i, je, k, ivz)*rhoi + w(i, je, k, ivz)*&
   &         rhoid
   w(i, je, k, ivz) = w(i, je, k, ivz)*rhoi
   wd(i, je, k, irhoe) = wd(i, je, k, irhoe) - pd(i, je, k)
   w(i, je, k, irhoe) = w(i, je, k, irhoe) - p(i, je, k)
   rhoid = -(one*wd(i, jb, k, irho)/w(i, jb, k, irho)**2)
   rhoi = one/w(i, jb, k, irho)
   wd(i, jb, k, ivx) = wd(i, jb, k, ivx)*rhoi + w(i, jb, k, ivx)*&
   &         rhoid
   w(i, jb, k, ivx) = w(i, jb, k, ivx)*rhoi
   wd(i, jb, k, ivy) = wd(i, jb, k, ivy)*rhoi + w(i, jb, k, ivy)*&
   &         rhoid
   w(i, jb, k, ivy) = w(i, jb, k, ivy)*rhoi
   wd(i, jb, k, ivz) = wd(i, jb, k, ivz)*rhoi + w(i, jb, k, ivz)*&
   &         rhoid
   w(i, jb, k, ivz) = w(i, jb, k, ivz)*rhoi
   wd(i, jb, k, irhoe) = wd(i, jb, k, irhoe) - pd(i, jb, k)
   w(i, jb, k, irhoe) = w(i, jb, k, irhoe) - p(i, jb, k)
   END DO
   END DO
   END IF
   END SUBROUTINE INVISCIDDISSFLUXSCALAR_D
