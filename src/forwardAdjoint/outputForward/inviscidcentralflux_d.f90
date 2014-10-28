   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of inviscidcentralflux in forward (tangent) mode (with options i4 dr8 r8):
   !   variations   of useful results: *dw
   !   with respect to varying inputs: *p *dw *w
   !   Plus diff mem management of: p:in dw:in w:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          inviscidCentralFlux.f90                         *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 03-24-2003                                      *
   !      * Last modified: 10-29-2007                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE INVISCIDCENTRALFLUX_D()
   !
   !      ******************************************************************
   !      *                                                                *
   !      * inviscidCentralFlux computes the Euler fluxes using a central  *
   !      * discretization for a given block. Therefore it is assumed that *
   !      * the pointers in block pointer already point to the correct     *
   !      * block on the correct multigrid level.                          *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_D
   USE CGNSGRID
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE INPUTPHYSICS
   IMPLICIT NONE
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k, ind
   REAL(kind=realtype) :: qsp, qsm, rqsp, rqsm, porvel, porflux
   REAL(kind=realtype) :: qspd, qsmd, rqspd, rqsmd
   REAL(kind=realtype) :: pa, fs, sface, vnp, vnm
   REAL(kind=realtype) :: pad, fsd, vnpd, vnmd
   REAL(kind=realtype) :: wx, wy, wz, rvol
   REAL(kind=realtype) :: rvold
   sface = zero
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Advective fluxes in the i-direction.                           *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=2,jl
   DO i=1,il
   ! Set the dot product of the grid velocity and the
   ! normal in i-direction for a moving face.
   IF (addgridvelocities) sface = sfacei(i, j, k)
   ! Compute the normal velocities of the left and right state.
   vnpd = si(i, j, k, 1)*wd(i+1, j, k, ivx) + si(i, j, k, 2)*wd(i+1&
   &         , j, k, ivy) + si(i, j, k, 3)*wd(i+1, j, k, ivz)
   vnp = w(i+1, j, k, ivx)*si(i, j, k, 1) + w(i+1, j, k, ivy)*si(i&
   &         , j, k, 2) + w(i+1, j, k, ivz)*si(i, j, k, 3)
   vnmd = si(i, j, k, 1)*wd(i, j, k, ivx) + si(i, j, k, 2)*wd(i, j&
   &         , k, ivy) + si(i, j, k, 3)*wd(i, j, k, ivz)
   vnm = w(i, j, k, ivx)*si(i, j, k, 1) + w(i, j, k, ivy)*si(i, j, &
   &         k, 2) + w(i, j, k, ivz)*si(i, j, k, 3)
   ! Set the values of the porosities for this face.
   ! porVel defines the porosity w.r.t. velocity;
   ! porFlux defines the porosity w.r.t. the entire flux.
   ! The latter is only zero for a discontinuous block
   ! boundary that must be treated conservatively.
   ! The default value of porFlux is 0.5, such that the
   ! correct central flux is scattered to both cells.
   ! In case of a boundFlux the normal velocity is set
   ! to sFace.
   porvel = one
   porflux = half
   IF (pori(i, j, k) .EQ. noflux) porflux = zero
   IF (pori(i, j, k) .EQ. boundflux) THEN
   porvel = zero
   vnp = sface
   vnm = sface
   vnmd = 0.0_8
   vnpd = 0.0_8
   END IF
   ! Incorporate porFlux in porVel.
   porvel = porvel*porflux
   ! Compute the normal velocities relative to the grid for
   ! the face as well as the mass fluxes.
   qspd = porvel*vnpd
   qsp = (vnp-sface)*porvel
   qsmd = porvel*vnmd
   qsm = (vnm-sface)*porvel
   rqspd = qspd*w(i+1, j, k, irho) + qsp*wd(i+1, j, k, irho)
   rqsp = qsp*w(i+1, j, k, irho)
   rqsmd = qsmd*w(i, j, k, irho) + qsm*wd(i, j, k, irho)
   rqsm = qsm*w(i, j, k, irho)
   ! Compute the sum of the pressure multiplied by porFlux.
   ! For the default value of porFlux, 0.5, this leads to
   ! the average pressure.
   pad = porflux*(pd(i+1, j, k)+pd(i, j, k))
   pa = porflux*(p(i+1, j, k)+p(i, j, k))
   ! Compute the fluxes and scatter them to the cells
   ! i,j,k and i+1,j,k. Store the density flux in the
   ! mass flow of the appropriate sliding mesh interface.
   fsd = rqspd + rqsmd
   fs = rqsp + rqsm
   dwd(i+1, j, k, irho) = dwd(i+1, j, k, irho) - fsd
   dw(i+1, j, k, irho) = dw(i+1, j, k, irho) - fs
   dwd(i, j, k, irho) = dwd(i, j, k, irho) + fsd
   dw(i, j, k, irho) = dw(i, j, k, irho) + fs
   fsd = rqspd*w(i+1, j, k, ivx) + rqsp*wd(i+1, j, k, ivx) + rqsmd*&
   &         w(i, j, k, ivx) + rqsm*wd(i, j, k, ivx) + si(i, j, k, 1)*pad
   fs = rqsp*w(i+1, j, k, ivx) + rqsm*w(i, j, k, ivx) + pa*si(i, j&
   &         , k, 1)
   dwd(i+1, j, k, imx) = dwd(i+1, j, k, imx) - fsd
   dw(i+1, j, k, imx) = dw(i+1, j, k, imx) - fs
   dwd(i, j, k, imx) = dwd(i, j, k, imx) + fsd
   dw(i, j, k, imx) = dw(i, j, k, imx) + fs
   fsd = rqspd*w(i+1, j, k, ivy) + rqsp*wd(i+1, j, k, ivy) + rqsmd*&
   &         w(i, j, k, ivy) + rqsm*wd(i, j, k, ivy) + si(i, j, k, 2)*pad
   fs = rqsp*w(i+1, j, k, ivy) + rqsm*w(i, j, k, ivy) + pa*si(i, j&
   &         , k, 2)
   dwd(i+1, j, k, imy) = dwd(i+1, j, k, imy) - fsd
   dw(i+1, j, k, imy) = dw(i+1, j, k, imy) - fs
   dwd(i, j, k, imy) = dwd(i, j, k, imy) + fsd
   dw(i, j, k, imy) = dw(i, j, k, imy) + fs
   fsd = rqspd*w(i+1, j, k, ivz) + rqsp*wd(i+1, j, k, ivz) + rqsmd*&
   &         w(i, j, k, ivz) + rqsm*wd(i, j, k, ivz) + si(i, j, k, 3)*pad
   fs = rqsp*w(i+1, j, k, ivz) + rqsm*w(i, j, k, ivz) + pa*si(i, j&
   &         , k, 3)
   dwd(i+1, j, k, imz) = dwd(i+1, j, k, imz) - fsd
   dw(i+1, j, k, imz) = dw(i+1, j, k, imz) - fs
   dwd(i, j, k, imz) = dwd(i, j, k, imz) + fsd
   dw(i, j, k, imz) = dw(i, j, k, imz) + fs
   fsd = qspd*w(i+1, j, k, irhoe) + qsp*wd(i+1, j, k, irhoe) + qsmd&
   &         *w(i, j, k, irhoe) + qsm*wd(i, j, k, irhoe) + porflux*(vnpd*p(&
   &         i+1, j, k)+vnp*pd(i+1, j, k)+vnmd*p(i, j, k)+vnm*pd(i, j, k))
   fs = qsp*w(i+1, j, k, irhoe) + qsm*w(i, j, k, irhoe) + porflux*(&
   &         vnp*p(i+1, j, k)+vnm*p(i, j, k))
   dwd(i+1, j, k, irhoe) = dwd(i+1, j, k, irhoe) - fsd
   dw(i+1, j, k, irhoe) = dw(i+1, j, k, irhoe) - fs
   dwd(i, j, k, irhoe) = dwd(i, j, k, irhoe) + fsd
   dw(i, j, k, irhoe) = dw(i, j, k, irhoe) + fs
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Advective fluxes in the j-direction.                           *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=1,jl
   DO i=2,il
   ! Set the dot product of the grid velocity and the
   ! normal in j-direction for a moving face.
   IF (addgridvelocities) sface = sfacej(i, j, k)
   ! Compute the normal velocities of the left and right state.
   vnpd = sj(i, j, k, 1)*wd(i, j+1, k, ivx) + sj(i, j, k, 2)*wd(i, &
   &         j+1, k, ivy) + sj(i, j, k, 3)*wd(i, j+1, k, ivz)
   vnp = w(i, j+1, k, ivx)*sj(i, j, k, 1) + w(i, j+1, k, ivy)*sj(i&
   &         , j, k, 2) + w(i, j+1, k, ivz)*sj(i, j, k, 3)
   vnmd = sj(i, j, k, 1)*wd(i, j, k, ivx) + sj(i, j, k, 2)*wd(i, j&
   &         , k, ivy) + sj(i, j, k, 3)*wd(i, j, k, ivz)
   vnm = w(i, j, k, ivx)*sj(i, j, k, 1) + w(i, j, k, ivy)*sj(i, j, &
   &         k, 2) + w(i, j, k, ivz)*sj(i, j, k, 3)
   ! Set the values of the porosities for this face.
   ! porVel defines the porosity w.r.t. velocity;
   ! porFlux defines the porosity w.r.t. the entire flux.
   ! The latter is only zero for a discontinuous block
   ! boundary that must be treated conservatively.
   ! The default value of porFlux is 0.5, such that the
   ! correct central flux is scattered to both cells.
   ! In case of a boundFlux the normal velocity is set
   ! to sFace.
   porvel = one
   porflux = half
   IF (porj(i, j, k) .EQ. noflux) porflux = zero
   IF (porj(i, j, k) .EQ. boundflux) THEN
   porvel = zero
   vnp = sface
   vnm = sface
   vnmd = 0.0_8
   vnpd = 0.0_8
   END IF
   ! Incorporate porFlux in porVel.
   porvel = porvel*porflux
   ! Compute the normal velocities for the face as well as the
   ! mass fluxes.
   qspd = porvel*vnpd
   qsp = (vnp-sface)*porvel
   qsmd = porvel*vnmd
   qsm = (vnm-sface)*porvel
   rqspd = qspd*w(i, j+1, k, irho) + qsp*wd(i, j+1, k, irho)
   rqsp = qsp*w(i, j+1, k, irho)
   rqsmd = qsmd*w(i, j, k, irho) + qsm*wd(i, j, k, irho)
   rqsm = qsm*w(i, j, k, irho)
   ! Compute the sum of the pressure multiplied by porFlux.
   ! For the default value of porFlux, 0.5, this leads to
   ! the average pressure.
   pad = porflux*(pd(i, j+1, k)+pd(i, j, k))
   pa = porflux*(p(i, j+1, k)+p(i, j, k))
   ! Compute the fluxes and scatter them to the cells
   ! i,j,k and i,j+1,k. Store the density flux in the
   ! mass flow of the appropriate sliding mesh interface.
   fsd = rqspd + rqsmd
   fs = rqsp + rqsm
   dwd(i, j+1, k, irho) = dwd(i, j+1, k, irho) - fsd
   dw(i, j+1, k, irho) = dw(i, j+1, k, irho) - fs
   dwd(i, j, k, irho) = dwd(i, j, k, irho) + fsd
   dw(i, j, k, irho) = dw(i, j, k, irho) + fs
   fsd = rqspd*w(i, j+1, k, ivx) + rqsp*wd(i, j+1, k, ivx) + rqsmd*&
   &         w(i, j, k, ivx) + rqsm*wd(i, j, k, ivx) + sj(i, j, k, 1)*pad
   fs = rqsp*w(i, j+1, k, ivx) + rqsm*w(i, j, k, ivx) + pa*sj(i, j&
   &         , k, 1)
   dwd(i, j+1, k, imx) = dwd(i, j+1, k, imx) - fsd
   dw(i, j+1, k, imx) = dw(i, j+1, k, imx) - fs
   dwd(i, j, k, imx) = dwd(i, j, k, imx) + fsd
   dw(i, j, k, imx) = dw(i, j, k, imx) + fs
   fsd = rqspd*w(i, j+1, k, ivy) + rqsp*wd(i, j+1, k, ivy) + rqsmd*&
   &         w(i, j, k, ivy) + rqsm*wd(i, j, k, ivy) + sj(i, j, k, 2)*pad
   fs = rqsp*w(i, j+1, k, ivy) + rqsm*w(i, j, k, ivy) + pa*sj(i, j&
   &         , k, 2)
   dwd(i, j+1, k, imy) = dwd(i, j+1, k, imy) - fsd
   dw(i, j+1, k, imy) = dw(i, j+1, k, imy) - fs
   dwd(i, j, k, imy) = dwd(i, j, k, imy) + fsd
   dw(i, j, k, imy) = dw(i, j, k, imy) + fs
   fsd = rqspd*w(i, j+1, k, ivz) + rqsp*wd(i, j+1, k, ivz) + rqsmd*&
   &         w(i, j, k, ivz) + rqsm*wd(i, j, k, ivz) + sj(i, j, k, 3)*pad
   fs = rqsp*w(i, j+1, k, ivz) + rqsm*w(i, j, k, ivz) + pa*sj(i, j&
   &         , k, 3)
   dwd(i, j+1, k, imz) = dwd(i, j+1, k, imz) - fsd
   dw(i, j+1, k, imz) = dw(i, j+1, k, imz) - fs
   dwd(i, j, k, imz) = dwd(i, j, k, imz) + fsd
   dw(i, j, k, imz) = dw(i, j, k, imz) + fs
   fsd = qspd*w(i, j+1, k, irhoe) + qsp*wd(i, j+1, k, irhoe) + qsmd&
   &         *w(i, j, k, irhoe) + qsm*wd(i, j, k, irhoe) + porflux*(vnpd*p(&
   &         i, j+1, k)+vnp*pd(i, j+1, k)+vnmd*p(i, j, k)+vnm*pd(i, j, k))
   fs = qsp*w(i, j+1, k, irhoe) + qsm*w(i, j, k, irhoe) + porflux*(&
   &         vnp*p(i, j+1, k)+vnm*p(i, j, k))
   dwd(i, j+1, k, irhoe) = dwd(i, j+1, k, irhoe) - fsd
   dw(i, j+1, k, irhoe) = dw(i, j+1, k, irhoe) - fs
   dwd(i, j, k, irhoe) = dwd(i, j, k, irhoe) + fsd
   dw(i, j, k, irhoe) = dw(i, j, k, irhoe) + fs
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Advective fluxes in the k-direction.                           *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=1,kl
   DO j=2,jl
   DO i=2,il
   ! Set the dot product of the grid velocity and the
   ! normal in k-direction for a moving face.
   IF (addgridvelocities) sface = sfacek(i, j, k)
   ! Compute the normal velocities of the left and right state.
   vnpd = sk(i, j, k, 1)*wd(i, j, k+1, ivx) + sk(i, j, k, 2)*wd(i, &
   &         j, k+1, ivy) + sk(i, j, k, 3)*wd(i, j, k+1, ivz)
   vnp = w(i, j, k+1, ivx)*sk(i, j, k, 1) + w(i, j, k+1, ivy)*sk(i&
   &         , j, k, 2) + w(i, j, k+1, ivz)*sk(i, j, k, 3)
   vnmd = sk(i, j, k, 1)*wd(i, j, k, ivx) + sk(i, j, k, 2)*wd(i, j&
   &         , k, ivy) + sk(i, j, k, 3)*wd(i, j, k, ivz)
   vnm = w(i, j, k, ivx)*sk(i, j, k, 1) + w(i, j, k, ivy)*sk(i, j, &
   &         k, 2) + w(i, j, k, ivz)*sk(i, j, k, 3)
   ! Set the values of the porosities for this face.
   ! porVel defines the porosity w.r.t. velocity;
   ! porFlux defines the porosity w.r.t. the entire flux.
   ! The latter is only zero for a discontinuous block
   ! block boundary that must be treated conservatively.
   ! The default value of porFlux is 0.5, such that the
   ! correct central flux is scattered to both cells.
   ! In case of a boundFlux the normal velocity is set
   ! to sFace.
   porvel = one
   porflux = half
   IF (pork(i, j, k) .EQ. noflux) porflux = zero
   IF (pork(i, j, k) .EQ. boundflux) THEN
   porvel = zero
   vnp = sface
   vnm = sface
   vnmd = 0.0_8
   vnpd = 0.0_8
   END IF
   ! Incorporate porFlux in porVel.
   porvel = porvel*porflux
   ! Compute the normal velocities for the face as well as the
   ! mass fluxes.
   qspd = porvel*vnpd
   qsp = (vnp-sface)*porvel
   qsmd = porvel*vnmd
   qsm = (vnm-sface)*porvel
   rqspd = qspd*w(i, j, k+1, irho) + qsp*wd(i, j, k+1, irho)
   rqsp = qsp*w(i, j, k+1, irho)
   rqsmd = qsmd*w(i, j, k, irho) + qsm*wd(i, j, k, irho)
   rqsm = qsm*w(i, j, k, irho)
   ! Compute the sum of the pressure multiplied by porFlux.
   ! For the default value of porFlux, 0.5, this leads to
   ! the average pressure.
   pad = porflux*(pd(i, j, k+1)+pd(i, j, k))
   pa = porflux*(p(i, j, k+1)+p(i, j, k))
   ! Compute the fluxes and scatter them to the cells
   ! i,j,k and i,j,k+1. Store the density flux in the
   ! mass flow of the appropriate sliding mesh interface.
   fsd = rqspd + rqsmd
   fs = rqsp + rqsm
   dwd(i, j, k+1, irho) = dwd(i, j, k+1, irho) - fsd
   dw(i, j, k+1, irho) = dw(i, j, k+1, irho) - fs
   dwd(i, j, k, irho) = dwd(i, j, k, irho) + fsd
   dw(i, j, k, irho) = dw(i, j, k, irho) + fs
   fsd = rqspd*w(i, j, k+1, ivx) + rqsp*wd(i, j, k+1, ivx) + rqsmd*&
   &         w(i, j, k, ivx) + rqsm*wd(i, j, k, ivx) + sk(i, j, k, 1)*pad
   fs = rqsp*w(i, j, k+1, ivx) + rqsm*w(i, j, k, ivx) + pa*sk(i, j&
   &         , k, 1)
   dwd(i, j, k+1, imx) = dwd(i, j, k+1, imx) - fsd
   dw(i, j, k+1, imx) = dw(i, j, k+1, imx) - fs
   dwd(i, j, k, imx) = dwd(i, j, k, imx) + fsd
   dw(i, j, k, imx) = dw(i, j, k, imx) + fs
   fsd = rqspd*w(i, j, k+1, ivy) + rqsp*wd(i, j, k+1, ivy) + rqsmd*&
   &         w(i, j, k, ivy) + rqsm*wd(i, j, k, ivy) + sk(i, j, k, 2)*pad
   fs = rqsp*w(i, j, k+1, ivy) + rqsm*w(i, j, k, ivy) + pa*sk(i, j&
   &         , k, 2)
   dwd(i, j, k+1, imy) = dwd(i, j, k+1, imy) - fsd
   dw(i, j, k+1, imy) = dw(i, j, k+1, imy) - fs
   dwd(i, j, k, imy) = dwd(i, j, k, imy) + fsd
   dw(i, j, k, imy) = dw(i, j, k, imy) + fs
   fsd = rqspd*w(i, j, k+1, ivz) + rqsp*wd(i, j, k+1, ivz) + rqsmd*&
   &         w(i, j, k, ivz) + rqsm*wd(i, j, k, ivz) + sk(i, j, k, 3)*pad
   fs = rqsp*w(i, j, k+1, ivz) + rqsm*w(i, j, k, ivz) + pa*sk(i, j&
   &         , k, 3)
   dwd(i, j, k+1, imz) = dwd(i, j, k+1, imz) - fsd
   dw(i, j, k+1, imz) = dw(i, j, k+1, imz) - fs
   dwd(i, j, k, imz) = dwd(i, j, k, imz) + fsd
   dw(i, j, k, imz) = dw(i, j, k, imz) + fs
   fsd = qspd*w(i, j, k+1, irhoe) + qsp*wd(i, j, k+1, irhoe) + qsmd&
   &         *w(i, j, k, irhoe) + qsm*wd(i, j, k, irhoe) + porflux*(vnpd*p(&
   &         i, j, k+1)+vnp*pd(i, j, k+1)+vnmd*p(i, j, k)+vnm*pd(i, j, k))
   fs = qsp*w(i, j, k+1, irhoe) + qsm*w(i, j, k, irhoe) + porflux*(&
   &         vnp*p(i, j, k+1)+vnm*p(i, j, k))
   dwd(i, j, k+1, irhoe) = dwd(i, j, k+1, irhoe) - fsd
   dw(i, j, k+1, irhoe) = dw(i, j, k+1, irhoe) - fs
   dwd(i, j, k, irhoe) = dwd(i, j, k, irhoe) + fsd
   dw(i, j, k, irhoe) = dw(i, j, k, irhoe) + fs
   END DO
   END DO
   END DO
   ! Add the rotational source terms for a moving block in a
   ! steady state computation. These source terms account for the
   ! centrifugal acceleration and the coriolis term. However, as
   ! the the equations are solved in the inertial frame and not
   ! in the moving frame, the form is different than what you
   ! normally find in a text book.
   IF (blockismoving .AND. equationmode .EQ. steady) THEN
   ! Compute the three nonDimensional angular velocities.
   wx = timeref*cgnsdoms(nbkglobal)%rotrate(1)
   wy = timeref*cgnsdoms(nbkglobal)%rotrate(2)
   wz = timeref*cgnsdoms(nbkglobal)%rotrate(3)
   ! Loop over the internal cells of this block to compute the
   ! rotational terms for the momentum equations.
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   rvold = vol(i, j, k)*wd(i, j, k, irho)
   rvol = w(i, j, k, irho)*vol(i, j, k)
   dwd(i, j, k, imx) = dwd(i, j, k, imx) + rvold*(wy*w(i, j, k, &
   &           ivz)-wz*w(i, j, k, ivy)) + rvol*(wy*wd(i, j, k, ivz)-wz*wd(i&
   &           , j, k, ivy))
   dw(i, j, k, imx) = dw(i, j, k, imx) + rvol*(wy*w(i, j, k, ivz)&
   &           -wz*w(i, j, k, ivy))
   dwd(i, j, k, imy) = dwd(i, j, k, imy) + rvold*(wz*w(i, j, k, &
   &           ivx)-wx*w(i, j, k, ivz)) + rvol*(wz*wd(i, j, k, ivx)-wx*wd(i&
   &           , j, k, ivz))
   dw(i, j, k, imy) = dw(i, j, k, imy) + rvol*(wz*w(i, j, k, ivx)&
   &           -wx*w(i, j, k, ivz))
   dwd(i, j, k, imz) = dwd(i, j, k, imz) + rvold*(wx*w(i, j, k, &
   &           ivy)-wy*w(i, j, k, ivx)) + rvol*(wx*wd(i, j, k, ivy)-wy*wd(i&
   &           , j, k, ivx))
   dw(i, j, k, imz) = dw(i, j, k, imz) + rvol*(wx*w(i, j, k, ivy)&
   &           -wy*w(i, j, k, ivx))
   END DO
   END DO
   END DO
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Initialize sFace to zero. This value will be used if the
   ! block is not moving.
   40 FORMAT(1x,i4,i4,i4,e20.6)
   END SUBROUTINE INVISCIDCENTRALFLUX_D
