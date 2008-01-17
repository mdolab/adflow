!
!      ******************************************************************
!      *                                                                *
!      * File:          timeStep.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-17-2003                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine timeStep(onlyRadii)
!
!      ******************************************************************
!      *                                                                *
!      * timeStep computes the time step, or more precisely the time    *
!      * step divided by the volume per unit CFL, in the owned cells.   *
!      * However, for the artificial dissipation schemes, the spectral  *
!      * radIi in the halo's are needed. Therefore the loop is taken    *
!      * over the the first level of halo cells. The spectral radIi are *
!      * stored and possibly modified for high aspect ratio cells.      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use flowVarRefState
       use inputDiscretization
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use iteration
       use section
       implicit none
!
!      Subroutine argument.
!
       logical, intent(in) :: onlyRadii
!
!      Local parameters.
!
       real(kind=realType), parameter :: b = 2.0_realType
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, i, j, k

       real(kind=realType) :: plim, rlim, clim2
       real(kind=realType) :: ux, uy, uz, cc2, qs, sx, sy, sz, rmu
       real(kind=realType) :: ri, rj, rk, rij, rjk, rki
       real(kind=realType) :: vsi, vsj, vsk, rfl, dpi, dpj, dpk
       real(kind=realType) :: sFace, tmp

       logical :: radiiNeeded
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not the spectral radii are needed for the
       ! flux computation.

       radiiNeeded = radiiNeededCoarse
       if(currentLevel <= groundLevel) radiiNeeded = radiiNeededFine

       ! Return immediately if only the spectral radii must be computed
       ! and these are not needed for the flux computation.

       if(onlyRadii .and. (.not. radiiNeeded)) return

       ! Set the value of plim. To be fully consistent this must have
       ! the dimension of a pressure. Therefore a fraction of pInfCorr
       ! is used. Idem for rlim; compute clim2 as well.
 
       plim  = 0.001_realType*pInfCorr
       rlim  = 0.001_realType*rhoInf
       clim2 = 0.000001_realType*gammaInf*pInfCorr/rhoInf

       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, currentLevel, sps)

           ! Initialize sFace to zero. This value will be used if the
           ! block is not moving.

           sFace = zero
!
!          **************************************************************
!          *                                                            *
!          * Inviscid contribution, depending on the preconditioner.    *
!          * Compute the cell centered values of the spectral radii.    *
!          *                                                            *
!          **************************************************************
!
           select case (precond)

             case (noPrecond)

               ! No preconditioner. Simply the standard spectral radius.
               ! Loop over the cells, including the first level halo.

               do k=1,ke
                 do j=1,je
                   do i=1,ie

                     ! Compute the velocities and speed of sound squared.

                     ux  = w(i,j,k,ivx)
                     uy  = w(i,j,k,ivy)
                     uz  = w(i,j,k,ivz)
                     cc2 = gamma(i,j,k)*p(i,j,k)/w(i,j,k,irho)
                     cc2 = max(cc2,clim2)

                     ! Set the dot product of the grid velocity and the
                     ! normal in i-direction for a moving face. To avoid
                     ! a number of multiplications by 0.5 simply the sum
                     ! is taken.

                     if( addGridVelocities ) &
                       sFace = sFaceI(i-1,j,k) + sFaceI(i,j,k)

                     ! Spectral radius in i-direction.

                     sx = si(i-1,j,k,1) + si(i,j,k,1)
                     sy = si(i-1,j,k,2) + si(i,j,k,2)
                     sz = si(i-1,j,k,3) + si(i,j,k,3)

                     qs = ux*sx + uy*sy + uz*sz - sFace

                     radi(i,j,k) = half*(abs(qs) &
                                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                     ! The grid velocity in j-direction.

                     if( addGridVelocities ) &
                       sFace = sFaceJ(i,j-1,k) + sFaceJ(i,j,k)

                     ! Spectral radius in j-direction.

                     sx = sj(i,j-1,k,1) + sj(i,j,k,1)
                     sy = sj(i,j-1,k,2) + sj(i,j,k,2)
                     sz = sj(i,j-1,k,3) + sj(i,j,k,3)

                     qs = ux*sx + uy*sy + uz*sz - sFace

                     radJ(i,j,k) = half*(abs(qs) &
                                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                     ! The grid velocity in k-direction.

                     if( addGridVelocities ) &
                       sFace = sFaceK(i,j,k-1) + sFaceK(i,j,k)

                     ! Spectral radius in k-direction.

                     sx = sk(i,j,k-1,1) + sk(i,j,k,1)
                     sy = sk(i,j,k-1,2) + sk(i,j,k,2)
                     sz = sk(i,j,k-1,3) + sk(i,j,k,3)

                     qs = ux*sx + uy*sy + uz*sz - sFace

                     radK(i,j,k) = half*(abs(qs) &
                                 +       sqrt(cc2*(sx**2 + sy**2 + sz**2)))

                     ! Compute the inviscid contribution to the time step.

                     dtl(i,j,k) = radi(i,j,k) + radJ(i,j,k) + radK(i,j,k)

                   enddo
                 enddo
               enddo

             case (Turkel)
               call terminate("timeStep", &
                              "Turkel preconditioner not &
                              &implemented yet")


             case (ChoiMerkle)
               call terminate("timeStep", &
                              "choi merkle preconditioner not &
                              &implemented yet")
           end select
!
!          **************************************************************
!          *                                                            *
!          * Adapt the spectral radii if directional scaling must be    *
!          * applied.                                                   *
!          *                                                            *
!          **************************************************************
!
           if(dirScaling .and. currentLevel <= groundLevel) then
           ! if( dirScaling ) then

             do k=1,ke
               do j=1,je
                 do i=1,ie

                   ! Avoid division by zero by clipping radi, radJ and
                   ! radK.

                   ri = max(radi(i,j,k),eps)
                   rj = max(radJ(i,j,k),eps)
                   rk = max(radK(i,j,k),eps)

                   ! Compute the scaling in the three coordinate
                   ! directions.

                   rij = (ri/rj)**adis
                   rjk = (rj/rk)**adis
                   rki = (rk/ri)**adis

                   ! Create the scaled versions of the aspect ratios.
                   ! Note that the multiplication is done with radi, radJ
                   ! and radK, such that the influence of the clipping
                   ! is negligible.

               !   radi(i,j,k) = third*radi(i,j,k)*(one + one/rij + rki)
               !   radJ(i,j,k) = third*radJ(i,j,k)*(one + one/rjk + rij)
               !   radK(i,j,k) = third*radK(i,j,k)*(one + one/rki + rjk)

                   radi(i,j,k) = radi(i,j,k)*(one + one/rij + rki)
                   radJ(i,j,k) = radJ(i,j,k)*(one + one/rjk + rij)
                   radK(i,j,k) = radK(i,j,k)*(one + one/rki + rjk)

                 enddo
               enddo
             enddo

           endif

           ! The rest of this file can be skipped if only the spectral
           ! radii need to be computed.

           testRadiiOnly: if(.not. onlyRadii) then

             ! The viscous contribution, if needed.

             viscousTerm: if( viscous ) then

               ! Loop over the owned cell centers.

               do k=2,kl
                 do j=2,jl
                   do i=2,il

                     ! Compute the effective viscosity coefficient. The
                     ! factor 0.5 is a combination of two things. In the
                     ! standard central discretization of a second
                     ! derivative there is a factor 2 multiplying the
                     ! central node. However in the code below not the
                     ! average but the sum of the left and the right face
                     ! is taken and squared. This leads to a factor 4.
                     ! Combining both effects leads to 0.5. Furthermore,
                     ! it is divided by the volume and density to obtain
                     ! the correct dimensions and multiplied by the
                     ! non-dimensional factor factVis.

                     rmu = rlv(i,j,k)
                     if( eddyModel ) rmu = rmu + rev(i,j,k)
                     rmu = half*rmu/(w(i,j,k,irho)*vol(i,j,k))

                     ! Add the viscous contribution in i-direction to the
                     ! (inverse) of the time step.

                     sx = si(i,j,k,1) + si(i-1,j,k,1)
                     sy = si(i,j,k,2) + si(i-1,j,k,2)
                     sz = si(i,j,k,3) + si(i-1,j,k,3)

                     vsi        = rmu*(sx*sx + sy*sy + sz*sz)
                     dtl(i,j,k) = dtl(i,j,k) + vsi

                     ! Add the viscous contribution in j-direction to the
                     ! (inverse) of the time step.

                     sx = sj(i,j,k,1) + sj(i,j-1,k,1)
                     sy = sj(i,j,k,2) + sj(i,j-1,k,2)
                     sz = sj(i,j,k,3) + sj(i,j-1,k,3)

                     vsj        = rmu*(sx*sx + sy*sy + sz*sz)
                     dtl(i,j,k) = dtl(i,j,k) + vsj

                     ! Add the viscous contribution in k-direction to the
                     ! (inverse) of the time step.

                     sx = sk(i,j,k,1) + sk(i,j,k-1,1)
                     sy = sk(i,j,k,2) + sk(i,j,k-1,2)
                     sz = sk(i,j,k,3) + sk(i,j,k-1,3)

                     vsk        = rmu*(sx*sx + sy*sy + sz*sz)
                     dtl(i,j,k) = dtl(i,j,k) + vsk

                   enddo
                 enddo
               enddo

             endif viscousTerm

             ! For the spectral mode an additional term term must be
             ! taken into account, which corresponds to the contribution
             ! of the highest frequency.

             if(equationMode == timeSpectral) then

               tmp = nTimeIntervalsSpectral*pi*timeRef &
                   / sections(sectionID)%timePeriod

               ! Loop over the owned cell centers and add the term.

               do k=2,kl
                 do j=2,jl
                   do i=2,il
                     dtl(i,j,k) = dtl(i,j,k) + tmp*vol(i,j,k)
                   enddo
                 enddo
               enddo

             endif

             ! Currently the inverse of dt/vol is stored in dtl. Invert
             ! this value such that the time step per unit cfl number is
             ! stored and correct in cases of high gradients.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   dpi = abs(p(i+1,j,k) - two*p(i,j,k) + p(i-1,j,k)) &
                       /    (p(i+1,j,k) + two*p(i,j,k) + p(i-1,j,k) + plim)
                   dpj = abs(p(i,j+1,k) - two*p(i,j,k) + p(i,j-1,k)) &
                       /    (p(i,j+1,k) + two*p(i,j,k) + p(i,j-1,k) + plim)
                   dpk = abs(p(i,j,k+1) - two*p(i,j,k) + p(i,j,k-1)) &
                       /    (p(i,j,k+1) + two*p(i,j,k) + p(i,j,k-1) + plim)
                   rfl = one/(one + b*(dpi  +dpj  +dpk))
 
                   dtl(i,j,k) = rfl/dtl(i,j,k)
                 enddo
               enddo
             enddo

           endif testRadiiOnly

         enddo domains
       enddo spectralLoop

       end subroutine timeStep
