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
       subroutine turbAdvection(mAdv, nAdv, offset, qq)
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
       use blockPointers
       use turbMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nAdv, mAdv, offset

       real(kind=realType), dimension(2:il,2:jl,2:kl,mAdv,mAdv), &
                                                      intent(inout) :: qq
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, kk

       real(kind=realType) :: qs, voli, xa, ya, za
       real(kind=realType) :: uu, dwt, dwtm1, dwtp1, dwti, dwtj, dwtk

       real(kind=realType), dimension(mAdv) :: impl
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
       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the grid velocity if present.
             ! It is taken as the average of k and k-1,

             voli = half/vol(i,j,k)
             if( addGridVelocities ) &
               qs = (sFaceK(i,j,k) + sFaceK(i,j,k-1))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces k and k-1.

             xa = (sk(i,j,k,1) + sk(i,j,k-1,1))*voli
             ya = (sk(i,j,k,2) + sk(i,j,k-1,2))*voli
             za = (sk(i,j,k,3) + sk(i,j,k-1,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velKdir: if(uu > zero) then

               ! Velocity has a component in positive k-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Second order; store the three differences for the
                   ! discretization of the derivative in k-direction.

                   dwtm1 = w(i,j,k-1,jj) - w(i,j,k-2,jj)
                   dwt   = w(i,j,k,  jj) - w(i,j,k-1,jj)
                   dwtp1 = w(i,j,k+1,jj) - w(i,j,k,  jj)

                   ! Construct the derivative in this cell center. This
                   ! is the first order upwind derivative with two
                   ! nonlinear corrections.

                   dwtk = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwtk = dwtk + half*dwt
                     else
                       dwtk = dwtk + half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwtk = dwtk - half*dwt
                     else
                       dwtk = dwtk - half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwtk = w(i,j,k,jj) - w(i,j,k-1,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side of
                 ! the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwtk

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                 ! For boundary cells k == 2, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(k == 2) then
                   do kk=1,mAdv
                     impl(kk) = bmtk1(i,j,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                   enddo
                 endif

               enddo

             else velKdir

               ! Velocity has a component in negative k-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Store the three differences for the discretization of
                   ! the derivative in k-direction.

                   dwtm1 = w(i,j,k,  jj) - w(i,j,k-1,jj)
                   dwt   = w(i,j,k+1,jj) - w(i,j,k,  jj)
                   dwtp1 = w(i,j,k+2,jj) - w(i,j,k+1,jj)

                   ! Construct the derivative in this cell center. This is
                   ! the first order upwind derivative with two nonlinear
                   ! corrections.

                   dwtk = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwtk = dwtk - half*dwt
                     else
                       dwtk = dwtk - half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwtk = dwtk + half*dwt
                     else
                       dwtk = dwtk + half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwtk = w(i,j,k+1,jj) - w(i,j,k,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side
                 ! of the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwtk

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. -uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                 ! For boundary cells k == kl, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(k == kl) then
                   do kk=1,mAdv
                     impl(kk) = bmtk2(i,j,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                   enddo
                 endif

               enddo

             endif velKdir

           enddo
         enddo
       enddo
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
       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the grid velocity if present.
             ! It is taken as the average of j and j-1,

             voli = half/vol(i,j,k)
             if( addGridVelocities ) &
               qs = (sFaceJ(i,j,k) + sFaceJ(i,j-1,k))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces j and j-1.

             xa = (sj(i,j,k,1) + sj(i,j-1,k,1))*voli
             ya = (sj(i,j,k,2) + sj(i,j-1,k,2))*voli
             za = (sj(i,j,k,3) + sj(i,j-1,k,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velJdir: if(uu > zero) then

               ! Velocity has a component in positive j-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Second order; store the three differences for the
                   ! discretization of the derivative in j-direction.

                   dwtm1 = w(i,j-1,k,jj) - w(i,j-2,k,jj)
                   dwt   = w(i,j,  k,jj) - w(i,j-1,k,jj)
                   dwtp1 = w(i,j+1,k,jj) - w(i,j,  k,jj)

                   ! Construct the derivative in this cell center. This is
                   ! the first order upwind derivative with two nonlinear
                   ! corrections.

                   dwtj = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwtj = dwtj + half*dwt
                     else
                       dwtj = dwtj + half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwtj = dwtj - half*dwt
                     else
                       dwtj = dwtj - half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwtj = w(i,j,k,jj) - w(i,j-1,k,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side of
                 ! the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwtj

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                 ! For boundary cells j == 2, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(j == 2) then
                   do kk=1,mAdv
                     impl(kk) = bmtj1(i,k,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                   enddo
                 endif

               enddo

             else velJdir

               ! Velocity has a component in negative j-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Store the three differences for the discretization of
                   ! the derivative in j-direction.

                   dwtm1 = w(i,j,  k,jj) - w(i,j-1,k,jj)
                   dwt   = w(i,j+1,k,jj) - w(i,j,  k,jj)
                   dwtp1 = w(i,j+2,k,jj) - w(i,j+1,k,jj)

                   ! Construct the derivative in this cell center. This is
                   ! the first order upwind derivative with two nonlinear
                   ! corrections.

                   dwtj = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwtj = dwtj - half*dwt
                     else
                       dwtj = dwtj - half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwtj = dwtj + half*dwt
                     else
                       dwtj = dwtj + half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwtj = w(i,j+1,k,jj) - w(i,j,k,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side
                 ! of the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwtj

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. -uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                 ! For boundary cells j == jl, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(j == jl) then
                   do kk=1,mAdv
                     impl(kk) = bmtj2(i,k,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                   enddo
                 endif

               enddo

             endif velJdir

           enddo
         enddo
       enddo
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
       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the grid velocity if present.
             ! It is taken as the average of i and i-1,

             voli = half/vol(i,j,k)
             if( addGridVelocities ) &
               qs = (sFaceI(i,j,k) + sFaceI(i-1,j,k))*voli

             ! Compute the normal velocity, where the normal direction
             ! is taken as the average of faces i and i-1.

             xa = (si(i,j,k,1) + si(i-1,j,k,1))*voli
             ya = (si(i,j,k,2) + si(i-1,j,k,2))*voli
             za = (si(i,j,k,3) + si(i-1,j,k,3))*voli

             uu = xa*w(i,j,k,ivx) + ya*w(i,j,k,ivy) + za*w(i,j,k,ivz) - qs

             ! Determine the situation we are having here, i.e. positive
             ! or negative normal velocity.

             velIdir: if(uu > zero) then

               ! Velocity has a component in positive i-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Second order; store the three differences for the
                   ! discretization of the derivative in i-direction.

                   dwtm1 = w(i-1,j,k,jj) - w(i-2,j,k,jj)
                   dwt   = w(i,  j,k,jj) - w(i-1,j,k,jj)
                   dwtp1 = w(i+1,j,k,jj) - w(i,  j,k,jj)

                   ! Construct the derivative in this cell center. This is
                   ! the first order upwind derivative with two nonlinear
                   ! corrections.

                   dwti = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwti = dwti + half*dwt
                     else
                       dwti = dwti + half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwti = dwti - half*dwt
                     else
                       dwti = dwti - half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwti = w(i,j,k,jj) - w(i-1,j,k,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side of
                 ! the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwti

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) + uu

                 ! For boundary cells i == 2, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(i == 2) then
                   do kk=1,mAdv
                     impl(kk) = bmti1(j,k,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) + uu*impl(kk)
                   enddo
                 endif

               enddo

             else velIdir

               ! Velocity has a component in negative i-direction.
               ! Loop over the number of advection equations.

               do ii=1,nAdv

                 ! Set the value of jj such that it corresponds to the
                 ! turbulent entry in w.

                 jj = ii + offset

                 ! Check whether a first or a second order discretization
                 ! must be used.

                 if( secondOrd ) then

                   ! Second order; store the three differences for the
                   ! discretization of the derivative in i-direction.

                   dwtm1 = w(i,  j,k,jj) - w(i-1,j,k,jj)
                   dwt   = w(i+1,j,k,jj) - w(i,  j,k,jj)
                   dwtp1 = w(i+2,j,k,jj) - w(i+1,j,k,jj)

                   ! Construct the derivative in this cell center. This is
                   ! the first order upwind derivative with two nonlinear
                   ! corrections.

                   dwti = dwt

                   if(dwt*dwtp1 > zero) then
                     if(abs(dwt) < abs(dwtp1)) then
                       dwti = dwti - half*dwt
                     else
                       dwti = dwti - half*dwtp1
                     endif
                   endif

                   if(dwt*dwtm1 > zero) then
                     if(abs(dwt) < abs(dwtm1)) then
                       dwti = dwti + half*dwt
                     else
                       dwti = dwti + half*dwtm1
                     endif
                   endif

                 else

                   ! 1st order upwind scheme.

                   dwti = w(i+1,j,k,jj) - w(i,j,k,jj)

                 endif

                 ! Update the residual. The convective term must be
                 ! substracted, because it appears on the other side
                 ! of the equation as the source and viscous terms.

                 dvt(i,j,k,ii) = dvt(i,j,k,ii) - uu*dwti

                 ! Update the central jacobian. First the term which is
                 ! always present, i.e. -uu.

                 qq(i,j,k,ii,ii) = qq(i,j,k,ii,ii) - uu

                 ! For boundary cells i == il, the implicit treatment must
                 ! be taken into account. Note that the implicit part
                 ! is only based on the 1st order discretization.
                 ! To improve stability the diagonal term is only taken
                 ! into account when it improves stability, i.e. when
                 ! it is positive.

                 if(i == il) then
                   do kk=1,mAdv
                     impl(kk) = bmti2(j,k,jj,kk+offset)
                   enddo

                   impl(ii) = max(impl(ii),zero)

                   do kk=1,mAdv
                     qq(i,j,k,ii,kk) = qq(i,j,k,ii,kk) - uu*impl(kk)
                   enddo
                 endif

               enddo

             endif velIdir

           enddo
         enddo
       enddo

       end subroutine turbAdvection
