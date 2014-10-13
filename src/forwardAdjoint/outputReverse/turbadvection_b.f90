   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !  Differentiation of turbadvection in reverse (adjoint) mode (with options i4 dr8 r8 noISIZE):
   !   gradient     of useful results: *dw *w qq
   !   with respect to varying inputs: *dw *w qq
   !   Plus diff mem management of: dw:in w:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          turbAdvection.f90                               *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 09-01-2003                                      *
   !      * Last modified: 04-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE TURBADVECTION_B(madv, nadv, offset, qq, qqb)
   !
   !      ******************************************************************
   !      *                                                                *
   !      * turbAdvection discretizes the advection part of the turbulent  *
   !      * transport equations. As the advection part is the same for all *
   !      * models, this generic routine can be used. Both the             *
   !      * discretization and the central jacobian are computed in this   *
   !      * subroutine. The former can either be 1st or 2nd order          *
   !      * accurate; the latter is always based on the 1st order upwind   *
   !      * discretization. When the discretization must be second order   *
   !      * accurate, the fully upwind (kappa = -1) scheme in combination  *
   !      * with the minmod limiter is used.                               *
   !      *                                                                *
   !      * Only nAdv equations are treated, while the actual system has   *
   !      * size mAdv. The reason is that some equations for some          *
   !      * turbulence equations do not have an advection part, e.g. the   *
   !      * f equation in the v2-f model. The argument offset indicates    *
   !      * the offset in the w vector where this subsystem starts. As a   *
   !      * consequence it is assumed that the indices of the current      *
   !      * subsystem are contiguous, e.g. if a 2*2 system is solved the   *
   !      * Last index in w is offset+1 and offset+2 respectively.         *
   !      *                                                                *
   !      ******************************************************************
   !
   USE BLOCKPOINTERS_B
   USE TURBMOD
   IMPLICIT NONE
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nadv, madv, offset
   REAL(kind=realtype), DIMENSION(2:il, 2:jl, 2:kl, madv, madv), INTENT(&
   & INOUT) :: qq
   REAL(kind=realtype), DIMENSION(2:il, 2:jl, 2:kl, madv, madv), INTENT(&
   & INOUT) :: qqb
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k, ii, jj, kk
   REAL(kind=realtype) :: qs, voli, xa, ya, za
   REAL(kind=realtype) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk
   REAL(kind=realtype) :: uub, dwtb, dwtm1b, dwtp1b, dwtib, dwtjb, dwtkb
   REAL(kind=realtype), DIMENSION(madv) :: impl
   INTRINSIC ABS
   INTRINSIC MAX
   INTEGER :: branch
   REAL(kind=realtype) :: abs23
   REAL(kind=realtype) :: abs22
   REAL(kind=realtype) :: abs21
   REAL(kind=realtype) :: abs20
   REAL(kind=realtype) :: abs19
   REAL(kind=realtype) :: abs18
   REAL(kind=realtype) :: abs17
   REAL(kind=realtype) :: abs16
   REAL(kind=realtype) :: abs15
   REAL(kind=realtype) :: abs14
   REAL(kind=realtype) :: abs13
   REAL(kind=realtype) :: abs12
   REAL(kind=realtype) :: abs11
   REAL(kind=realtype) :: abs10
   REAL(kind=realtype) :: abs9
   REAL(kind=realtype) :: abs8
   REAL(kind=realtype) :: abs7
   REAL(kind=realtype) :: abs6
   REAL(kind=realtype) :: abs5
   REAL(kind=realtype) :: abs4
   REAL(kind=realtype) :: abs3
   REAL(kind=realtype) :: abs2
   REAL(kind=realtype) :: abs1
   REAL(kind=realtype) :: abs0
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Initialize the grid velocity to zero. This value will be used
   ! if the block is not moving.
   qs = zero
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Upwind discretization of the convective term in k (zeta)       *
   !      * direction. Either the 1st order upwind or the second order     *
   !      * fully upwind interpolation scheme, kappa = -1, is used in      *
   !      * combination with the minmod limiter.                           *
   !      * The possible grid velocity must be taken into account.         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   ! Compute the grid velocity if present.
   ! It is taken as the average of k and k-1,
   voli = half/vol(i, j, k)
   IF (addgridvelocities) qs = (sfacek(i, j, k)+sfacek(i, j, k-1))*&
   &           voli
   ! Compute the normal velocity, where the normal direction
   ! is taken as the average of faces k and k-1.
   xa = (sk(i, j, k, 1)+sk(i, j, k-1, 1))*voli
   ya = (sk(i, j, k, 2)+sk(i, j, k-1, 2))*voli
   za = (sk(i, j, k, 3)+sk(i, j, k-1, 3))*voli
   CALL PUSHREAL8(uu)
   uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz&
   &         ) - qs
   ! Determine the situation we are having here, i.e. positive
   ! or negative normal velocity.
   IF (uu .GT. zero) THEN
   ! Velocity has a component in positive k-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Second order; store the three differences for the
   ! discretization of the derivative in k-direction.
   dwtm1 = w(i, j, k-1, jj) - w(i, j, k-2, jj)
   dwt = w(i, j, k, jj) - w(i, j, k-1, jj)
   dwtp1 = w(i, j, k+1, jj) - w(i, j, k, jj)
   ! Construct the derivative in this cell center. This
   ! is the first order upwind derivative with two
   ! nonlinear corrections.
   CALL PUSHREAL8(dwtk)
   dwtk = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs0 = dwt
   ELSE
   abs0 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs12 = dwtp1
   ELSE
   abs12 = -dwtp1
   END IF
   IF (abs0 .LT. abs12) THEN
   dwtk = dwtk + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtk = dwtk + half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs1 = dwt
   ELSE
   abs1 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs13 = dwtm1
   ELSE
   abs13 = -dwtm1
   END IF
   IF (abs1 .LT. abs13) THEN
   dwtk = dwtk - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtk = dwtk - half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwtk)
   dwtk = w(i, j, k, jj) - w(i, j, k-1, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side of
   ! the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. uu.
   ! For boundary cells k == 2, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (k .EQ. 2) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmtk1(i, j, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(1)
   ELSE
   ! Velocity has a component in negative k-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Store the three differences for the discretization of
   ! the derivative in k-direction.
   dwtm1 = w(i, j, k, jj) - w(i, j, k-1, jj)
   dwt = w(i, j, k+1, jj) - w(i, j, k, jj)
   dwtp1 = w(i, j, k+2, jj) - w(i, j, k+1, jj)
   ! Construct the derivative in this cell center. This is
   ! the first order upwind derivative with two nonlinear
   ! corrections.
   CALL PUSHREAL8(dwtk)
   dwtk = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs2 = dwt
   ELSE
   abs2 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs14 = dwtp1
   ELSE
   abs14 = -dwtp1
   END IF
   IF (abs2 .LT. abs14) THEN
   dwtk = dwtk - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtk = dwtk - half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs3 = dwt
   ELSE
   abs3 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs15 = dwtm1
   ELSE
   abs15 = -dwtm1
   END IF
   IF (abs3 .LT. abs15) THEN
   dwtk = dwtk + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtk = dwtk + half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwtk)
   dwtk = w(i, j, k+1, jj) - w(i, j, k, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side
   ! of the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. -uu.
   ! For boundary cells k == kl, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (k .EQ. kl) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmtk2(i, j, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Upwind discretization of the convective term in j (eta)        *
   !      * direction. Either the 1st order upwind or the second order     *
   !      * fully upwind interpolation scheme, kappa = -1, is used in      *
   !      * combination with the minmod limiter.                           *
   !      * The possible grid velocity must be taken into account.         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   ! Compute the grid velocity if present.
   ! It is taken as the average of j and j-1,
   voli = half/vol(i, j, k)
   IF (addgridvelocities) qs = (sfacej(i, j, k)+sfacej(i, j-1, k))*&
   &           voli
   ! Compute the normal velocity, where the normal direction
   ! is taken as the average of faces j and j-1.
   xa = (sj(i, j, k, 1)+sj(i, j-1, k, 1))*voli
   ya = (sj(i, j, k, 2)+sj(i, j-1, k, 2))*voli
   za = (sj(i, j, k, 3)+sj(i, j-1, k, 3))*voli
   CALL PUSHREAL8(uu)
   uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz&
   &         ) - qs
   ! Determine the situation we are having here, i.e. positive
   ! or negative normal velocity.
   IF (uu .GT. zero) THEN
   ! Velocity has a component in positive j-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Second order; store the three differences for the
   ! discretization of the derivative in j-direction.
   dwtm1 = w(i, j-1, k, jj) - w(i, j-2, k, jj)
   dwt = w(i, j, k, jj) - w(i, j-1, k, jj)
   dwtp1 = w(i, j+1, k, jj) - w(i, j, k, jj)
   ! Construct the derivative in this cell center. This is
   ! the first order upwind derivative with two nonlinear
   ! corrections.
   CALL PUSHREAL8(dwtj)
   dwtj = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs4 = dwt
   ELSE
   abs4 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs16 = dwtp1
   ELSE
   abs16 = -dwtp1
   END IF
   IF (abs4 .LT. abs16) THEN
   dwtj = dwtj + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtj = dwtj + half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs5 = dwt
   ELSE
   abs5 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs17 = dwtm1
   ELSE
   abs17 = -dwtm1
   END IF
   IF (abs5 .LT. abs17) THEN
   dwtj = dwtj - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtj = dwtj - half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwtj)
   dwtj = w(i, j, k, jj) - w(i, j-1, k, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side of
   ! the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. uu.
   ! For boundary cells j == 2, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (j .EQ. 2) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmtj1(i, k, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(1)
   ELSE
   ! Velocity has a component in negative j-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Store the three differences for the discretization of
   ! the derivative in j-direction.
   dwtm1 = w(i, j, k, jj) - w(i, j-1, k, jj)
   dwt = w(i, j+1, k, jj) - w(i, j, k, jj)
   dwtp1 = w(i, j+2, k, jj) - w(i, j+1, k, jj)
   ! Construct the derivative in this cell center. This is
   ! the first order upwind derivative with two nonlinear
   ! corrections.
   CALL PUSHREAL8(dwtj)
   dwtj = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs6 = dwt
   ELSE
   abs6 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs18 = dwtp1
   ELSE
   abs18 = -dwtp1
   END IF
   IF (abs6 .LT. abs18) THEN
   dwtj = dwtj - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtj = dwtj - half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs7 = dwt
   ELSE
   abs7 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs19 = dwtm1
   ELSE
   abs19 = -dwtm1
   END IF
   IF (abs7 .LT. abs19) THEN
   dwtj = dwtj + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwtj = dwtj + half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwtj)
   dwtj = w(i, j+1, k, jj) - w(i, j, k, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side
   ! of the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. -uu.
   ! For boundary cells j == jl, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (j .EQ. jl) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmtj2(i, k, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   END DO
   END DO
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Upwind discretization of the convective term in i (xi)         *
   !      * direction. Either the 1st order upwind or the second order     *
   !      * fully upwind interpolation scheme, kappa = -1, is used in      *
   !      * combination with the minmod limiter.                           *
   !      * The possible grid velocity must be taken into account.         *
   !      *                                                                *
   !      ******************************************************************
   !
   DO k=2,kl
   DO j=2,jl
   DO i=2,il
   ! Compute the grid velocity if present.
   ! It is taken as the average of i and i-1,
   voli = half/vol(i, j, k)
   IF (addgridvelocities) qs = (sfacei(i, j, k)+sfacei(i-1, j, k))*&
   &           voli
   ! Compute the normal velocity, where the normal direction
   ! is taken as the average of faces i and i-1.
   xa = (si(i, j, k, 1)+si(i-1, j, k, 1))*voli
   ya = (si(i, j, k, 2)+si(i-1, j, k, 2))*voli
   za = (si(i, j, k, 3)+si(i-1, j, k, 3))*voli
   CALL PUSHREAL8(uu)
   uu = xa*w(i, j, k, ivx) + ya*w(i, j, k, ivy) + za*w(i, j, k, ivz&
   &         ) - qs
   ! Determine the situation we are having here, i.e. positive
   ! or negative normal velocity.
   IF (uu .GT. zero) THEN
   ! Velocity has a component in positive i-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Second order; store the three differences for the
   ! discretization of the derivative in i-direction.
   dwtm1 = w(i-1, j, k, jj) - w(i-2, j, k, jj)
   dwt = w(i, j, k, jj) - w(i-1, j, k, jj)
   dwtp1 = w(i+1, j, k, jj) - w(i, j, k, jj)
   ! Construct the derivative in this cell center. This is
   ! the first order upwind derivative with two nonlinear
   ! corrections.
   CALL PUSHREAL8(dwti)
   dwti = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs8 = dwt
   ELSE
   abs8 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs20 = dwtp1
   ELSE
   abs20 = -dwtp1
   END IF
   IF (abs8 .LT. abs20) THEN
   dwti = dwti + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwti = dwti + half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs9 = dwt
   ELSE
   abs9 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs21 = dwtm1
   ELSE
   abs21 = -dwtm1
   END IF
   IF (abs9 .LT. abs21) THEN
   dwti = dwti - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwti = dwti - half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwti)
   dwti = w(i, j, k, jj) - w(i-1, j, k, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side of
   ! the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. uu.
   ! For boundary cells i == 2, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (i .EQ. 2) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmti1(j, k, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(1)
   ELSE
   ! Velocity has a component in negative i-direction.
   ! Loop over the number of advection equations.
   DO ii=1,nadv
   ! Set the value of jj such that it corresponds to the
   ! turbulent entry in w.
   CALL PUSHINTEGER4(jj)
   jj = ii + offset
   ! Check whether a first or a second order discretization
   ! must be used.
   IF (secondord) THEN
   ! Second order; store the three differences for the
   ! discretization of the derivative in i-direction.
   dwtm1 = w(i, j, k, jj) - w(i-1, j, k, jj)
   dwt = w(i+1, j, k, jj) - w(i, j, k, jj)
   dwtp1 = w(i+2, j, k, jj) - w(i+1, j, k, jj)
   ! Construct the derivative in this cell center. This is
   ! the first order upwind derivative with two nonlinear
   ! corrections.
   CALL PUSHREAL8(dwti)
   dwti = dwt
   IF (dwt*dwtp1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs10 = dwt
   ELSE
   abs10 = -dwt
   END IF
   IF (dwtp1 .GE. 0.) THEN
   abs22 = dwtp1
   ELSE
   abs22 = -dwtp1
   END IF
   IF (abs10 .LT. abs22) THEN
   dwti = dwti - half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwti = dwti - half*dwtp1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   IF (dwt*dwtm1 .GT. zero) THEN
   IF (dwt .GE. 0.) THEN
   abs11 = dwt
   ELSE
   abs11 = -dwt
   END IF
   IF (dwtm1 .GE. 0.) THEN
   abs23 = dwtm1
   ELSE
   abs23 = -dwtm1
   END IF
   IF (abs11 .LT. abs23) THEN
   dwti = dwti + half*dwt
   CALL PUSHCONTROL2B(0)
   ELSE
   dwti = dwti + half*dwtm1
   CALL PUSHCONTROL2B(1)
   END IF
   ELSE
   CALL PUSHCONTROL2B(2)
   END IF
   ELSE
   ! 1st order upwind scheme.
   CALL PUSHREAL8(dwti)
   dwti = w(i+1, j, k, jj) - w(i, j, k, jj)
   CALL PUSHCONTROL2B(3)
   END IF
   ! Update the residual. The convective term must be
   ! substracted, because it appears on the other side
   ! of the equation as the source and viscous terms.
   ! Update the central jacobian. First the term which is
   ! always present, i.e. -uu.
   ! For boundary cells i == il, the implicit treatment must
   ! be taken into account. Note that the implicit part
   ! is only based on the 1st order discretization.
   ! To improve stability the diagonal term is only taken
   ! into account when it improves stability, i.e. when
   ! it is positive.
   IF (i .EQ. il) THEN
   DO kk=1,madv
   CALL PUSHREAL8(impl(kk))
   impl(kk) = bmti2(j, k, jj, kk+offset)
   END DO
   IF (impl(ii) .LT. zero) THEN
   CALL PUSHREAL8(impl(ii))
   impl(ii) = zero
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHREAL8(impl(ii))
   impl(ii) = impl(ii)
   CALL PUSHCONTROL1B(1)
   END IF
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   END DO
   END DO
   DO k=kl,2,-1
   DO j=jl,2,-1
   DO i=il,2,-1
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub - impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub - dwti*dwb(i, j, k, idvt+ii) - qqb(i, j, k, ii, ii&
   &             )
   dwtib = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = half*dwtib
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = half*dwtib
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwti)
   wb(i+1, j, k, jj) = wb(i+1, j, k, jj) + dwtib
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtib
   GOTO 100
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb - half*dwtib
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = -(half*dwtib)
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwti)
   dwtb = dwtb + dwtib
   wb(i+2, j, k, jj) = wb(i+2, j, k, jj) + dwtp1b
   wb(i+1, j, k, jj) = wb(i+1, j, k, jj) - dwtp1b
   wb(i+1, j, k, jj) = wb(i+1, j, k, jj) + dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtm1b
   wb(i-1, j, k, jj) = wb(i-1, j, k, jj) - dwtm1b
   100        CALL POPINTEGER4(jj)
   END DO
   ELSE
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub + impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub + qqb(i, j, k, ii, ii) - dwti*dwb(i, j, k, idvt+ii&
   &             )
   dwtib = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = -(half*dwtib)
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = -(half*dwtib)
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwti)
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtib
   wb(i-1, j, k, jj) = wb(i-1, j, k, jj) - dwtib
   GOTO 110
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb + half*dwtib
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = half*dwtib
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwti)
   dwtb = dwtb + dwtib
   wb(i+1, j, k, jj) = wb(i+1, j, k, jj) + dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtb
   wb(i-1, j, k, jj) = wb(i-1, j, k, jj) - dwtb
   wb(i-1, j, k, jj) = wb(i-1, j, k, jj) + dwtm1b
   wb(i-2, j, k, jj) = wb(i-2, j, k, jj) - dwtm1b
   110        CALL POPINTEGER4(jj)
   END DO
   END IF
   voli = half/vol(i, j, k)
   xa = (si(i, j, k, 1)+si(i-1, j, k, 1))*voli
   ya = (si(i, j, k, 2)+si(i-1, j, k, 2))*voli
   za = (si(i, j, k, 3)+si(i-1, j, k, 3))*voli
   CALL POPREAL8(uu)
   wb(i, j, k, ivx) = wb(i, j, k, ivx) + xa*uub
   wb(i, j, k, ivy) = wb(i, j, k, ivy) + ya*uub
   wb(i, j, k, ivz) = wb(i, j, k, ivz) + za*uub
   END DO
   END DO
   END DO
   DO k=kl,2,-1
   DO j=jl,2,-1
   DO i=il,2,-1
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub - impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub - dwtj*dwb(i, j, k, idvt+ii) - qqb(i, j, k, ii, ii&
   &             )
   dwtjb = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = half*dwtjb
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = half*dwtjb
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwtj)
   wb(i, j+1, k, jj) = wb(i, j+1, k, jj) + dwtjb
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtjb
   GOTO 120
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb - half*dwtjb
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = -(half*dwtjb)
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwtj)
   dwtb = dwtb + dwtjb
   wb(i, j+2, k, jj) = wb(i, j+2, k, jj) + dwtp1b
   wb(i, j+1, k, jj) = wb(i, j+1, k, jj) - dwtp1b
   wb(i, j+1, k, jj) = wb(i, j+1, k, jj) + dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtm1b
   wb(i, j-1, k, jj) = wb(i, j-1, k, jj) - dwtm1b
   120        CALL POPINTEGER4(jj)
   END DO
   ELSE
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub + impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub + qqb(i, j, k, ii, ii) - dwtj*dwb(i, j, k, idvt+ii&
   &             )
   dwtjb = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = -(half*dwtjb)
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = -(half*dwtjb)
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwtj)
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtjb
   wb(i, j-1, k, jj) = wb(i, j-1, k, jj) - dwtjb
   GOTO 130
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb + half*dwtjb
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = half*dwtjb
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwtj)
   dwtb = dwtb + dwtjb
   wb(i, j+1, k, jj) = wb(i, j+1, k, jj) + dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtb
   wb(i, j-1, k, jj) = wb(i, j-1, k, jj) - dwtb
   wb(i, j-1, k, jj) = wb(i, j-1, k, jj) + dwtm1b
   wb(i, j-2, k, jj) = wb(i, j-2, k, jj) - dwtm1b
   130        CALL POPINTEGER4(jj)
   END DO
   END IF
   voli = half/vol(i, j, k)
   xa = (sj(i, j, k, 1)+sj(i, j-1, k, 1))*voli
   ya = (sj(i, j, k, 2)+sj(i, j-1, k, 2))*voli
   za = (sj(i, j, k, 3)+sj(i, j-1, k, 3))*voli
   CALL POPREAL8(uu)
   wb(i, j, k, ivx) = wb(i, j, k, ivx) + xa*uub
   wb(i, j, k, ivy) = wb(i, j, k, ivy) + ya*uub
   wb(i, j, k, ivz) = wb(i, j, k, ivz) + za*uub
   END DO
   END DO
   END DO
   DO k=kl,2,-1
   DO j=jl,2,-1
   DO i=il,2,-1
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub - impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub - dwtk*dwb(i, j, k, idvt+ii) - qqb(i, j, k, ii, ii&
   &             )
   dwtkb = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = half*dwtkb
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = half*dwtkb
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwtk)
   wb(i, j, k+1, jj) = wb(i, j, k+1, jj) + dwtkb
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtkb
   GOTO 140
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb - half*dwtkb
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = -(half*dwtkb)
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwtk)
   dwtb = dwtb + dwtkb
   wb(i, j, k+2, jj) = wb(i, j, k+2, jj) + dwtp1b
   wb(i, j, k+1, jj) = wb(i, j, k+1, jj) - dwtp1b
   wb(i, j, k+1, jj) = wb(i, j, k+1, jj) + dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtb
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtm1b
   wb(i, j, k-1, jj) = wb(i, j, k-1, jj) - dwtm1b
   140        CALL POPINTEGER4(jj)
   END DO
   ELSE
   uub = 0.0_8
   DO ii=nadv,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) THEN
   DO kk=madv,1,-1
   uub = uub + impl(kk)*qqb(i, j, k, ii, kk)
   END DO
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   CALL POPREAL8(impl(ii))
   ELSE
   CALL POPREAL8(impl(ii))
   END IF
   jj = ii + offset
   DO kk=madv,1,-1
   CALL POPREAL8(impl(kk))
   END DO
   END IF
   uub = uub + qqb(i, j, k, ii, ii) - dwtk*dwb(i, j, k, idvt+ii&
   &             )
   dwtkb = -(uu*dwb(i, j, k, idvt+ii))
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   dwtb = -(half*dwtkb)
   dwtm1b = 0.0_8
   ELSE
   dwtm1b = -(half*dwtkb)
   dwtb = 0.0_8
   END IF
   ELSE IF (branch .EQ. 2) THEN
   dwtb = 0.0_8
   dwtm1b = 0.0_8
   ELSE
   CALL POPREAL8(dwtk)
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtkb
   wb(i, j, k-1, jj) = wb(i, j, k-1, jj) - dwtkb
   GOTO 150
   END IF
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   dwtb = dwtb + half*dwtkb
   dwtp1b = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   dwtp1b = half*dwtkb
   ELSE
   dwtp1b = 0.0_8
   END IF
   CALL POPREAL8(dwtk)
   dwtb = dwtb + dwtkb
   wb(i, j, k+1, jj) = wb(i, j, k+1, jj) + dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) - dwtp1b
   wb(i, j, k, jj) = wb(i, j, k, jj) + dwtb
   wb(i, j, k-1, jj) = wb(i, j, k-1, jj) - dwtb
   wb(i, j, k-1, jj) = wb(i, j, k-1, jj) + dwtm1b
   wb(i, j, k-2, jj) = wb(i, j, k-2, jj) - dwtm1b
   150        CALL POPINTEGER4(jj)
   END DO
   END IF
   voli = half/vol(i, j, k)
   xa = (sk(i, j, k, 1)+sk(i, j, k-1, 1))*voli
   ya = (sk(i, j, k, 2)+sk(i, j, k-1, 2))*voli
   za = (sk(i, j, k, 3)+sk(i, j, k-1, 3))*voli
   CALL POPREAL8(uu)
   wb(i, j, k, ivx) = wb(i, j, k, ivx) + xa*uub
   wb(i, j, k, ivy) = wb(i, j, k, ivy) + ya*uub
   wb(i, j, k, ivz) = wb(i, j, k, ivz) + za*uub
   END DO
   END DO
   END DO
   END SUBROUTINE TURBADVECTION_B
