!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidUpwindFlux.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-25-2003                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidUpwindFlux(fineGrid)
!
!      ******************************************************************
!      *                                                                *
!      * inviscidUpwindFlux computes the artificial dissipation part of *
!      * the Euler fluxes by means of an approximate solution of the 1D *
!      * Riemann problem on the face. For first order schemes,          *
!      * fineGrid == .false., the states in the cells are assumed to    *
!      * be constant; for the second order schemes on the fine grid a   *
!      * nonlinear reconstruction of the left and right state is done   *
!      * for which several options exist.                               *
!      * It is assumed that the pointers in blockPointers already       *
!      * point to the correct block.                                    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use constants
       use inputDiscretization
       use inputPhysics
       use flowVarRefState
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: fineGrid
!
!      Local variables.
!
       integer(kind=porType) :: por

       integer(kind=intType) :: nwInt
       integer(kind=intType) :: i, j, k, ind
       integer(kind=intType) :: limUsed, riemannUsed

       real(kind=realType) :: sx, sy, sz, omk, opk, sFil, gammaFace
       real(kind=realType) :: factMinmod, sFace

       real(kind=realType), dimension(nw)  :: left, right
       real(kind=realType), dimension(nw)  :: du1, du2, du3
       real(kind=realType), dimension(nwf) :: flux

       logical :: firstOrderK, correctForK, rotationalPeriodic
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if rFil == 0. If so, the dissipative flux needs not to
       ! be computed.

       if(rFil == zero) return

       ! Check if the formulation for rotational periodic problems
       ! must be used.

       if( associated(rotMatrixI) ) then
         rotationalPeriodic = .true.
       else
         rotationalPeriodic = .false.
       endif

       ! Initialize the dissipative residual to a certain times,
       ! possibly zero, the previously stored value. Owned cells
       ! only, because the halo values do not matter.

       sFil = one - rFil

       do k=2,kl
         do j=2,jl
           do i=2,il
             fw(i,j,k,irho)  = sFil*fw(i,j,k,irho)
             fw(i,j,k,imx)   = sFil*fw(i,j,k,imx)
             fw(i,j,k,imy)   = sFil*fw(i,j,k,imy)
             fw(i,j,k,imz)   = sFil*fw(i,j,k,imz)
             fw(i,j,k,irhoE) = sFil*fw(i,j,k,irhoE)
           enddo
         enddo
       enddo

       ! Determine whether or not the total energy must be corrected
       ! for the presence of the turbulent kinetic energy.

       if( kPresent ) then
         if((currentLevel == groundLevel) .or. turbCoupled) then
           correctForK = .true.
         else
           correctForK = .false.
         endif
       else
         correctForK = .false.
       endif

       ! Compute the factor used in the minmod limiter.

       factMinmod = (three-kappaCoef) &
                   / max(1.e-10_realType, one-kappaCoef)

       ! Determine the limiter scheme to be used. On the fine grid the
       ! user specified scheme is used; on the coarse grid a first order
       ! scheme is computed.

       limUsed = firstOrder
       if( fineGrid ) limUsed = limiter

       ! Determine the riemann solver which must be used.

       riemannUsed = riemannCoarse
       if( fineGrid ) riemannUsed = riemann

       ! Store 1-kappa and 1+kappa a bit easier and multiply it by 0.25.

       omk = fourth*(one - kappaCoef)
       opk = fourth*(one + kappaCoef)

       ! Initialize sFace to zero. This value will be used if the
       ! block is not moving.

       sFace = zero

       ! Set the number of variables to be interpolated depending
       ! whether or not a k-equation is present. If a k-equation is
       ! present also set the logical firstOrderK. This indicates
       ! whether or not only a first order approximation is to be used
       ! for the turbulent kinetic energy.

       if( correctForK ) then
         if(orderTurb == firstOrder) then
           nwInt       = nwf
           firstOrderK = .true.
         else
           nwInt       = itu1
           firstOrderK = .false.
         endif
       else
         nwInt       = nwf
         firstOrderK = .false.
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Flux computation. A distinction is made between first and      *
!      * second order schemes to avoid the overhead for the first order *
!      * scheme.                                                        *
!      *                                                                *
!      ******************************************************************
!
       orderTest: if(limUsed == firstOrder) then
!
!        ****************************************************************
!        *                                                              *
!        * First order reconstruction. The states in the cells are      *
!        * constant. The left and right states are constructed easily.  *
!        *                                                              *
!        ****************************************************************
!
         ! Fluxes in the i-direction.

         do k=2,kl
           do j=2,jl
             do i=1,il

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
               por = porI(i,j,k)
               if( addGridVelocities ) sFace = sFaceI(i,j,k)

               ! Determine the left and right state.

               left(irho)  = w(i,j,k,irho)
               left(ivx)   = w(i,j,k,ivx)
               left(ivy)   = w(i,j,k,ivy)
               left(ivz)   = w(i,j,k,ivz)
               left(irhoE) = p(i,j,k)
               if( correctForK ) left(itu1) = w(i,j,k,itu1)

               right(irho)  = w(i+1,j,k,irho)
               right(ivx)   = w(i+1,j,k,ivx)
               right(ivy)   = w(i+1,j,k,ivy)
               right(ivz)   = w(i+1,j,k,ivz)
               right(irhoE) = p(i+1,j,k)
               if( correctForK ) right(itu1) = w(i+1,j,k,itu1)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it to the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
               fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
               fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
               fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
               fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyI(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyI(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

         ! Fluxes in j-direction.

         do k=2,kl
           do j=1,jl
             do i=2,il

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
               por = porJ(i,j,k)
               if( addGridVelocities ) sFace = sFaceJ(i,j,k)

               ! Determine the left and right state.

               left(irho)  = w(i,j,k,irho)
               left(ivx)   = w(i,j,k,ivx)
               left(ivy)   = w(i,j,k,ivy)
               left(ivz)   = w(i,j,k,ivz)
               left(irhoE) = p(i,j,k)
               if( correctForK ) left(itu1) = w(i,j,k,itu1)

               right(irho)  = w(i,j+1,k,irho)
               right(ivx)   = w(i,j+1,k,ivx)
               right(ivy)   = w(i,j+1,k,ivy)
               right(ivz)   = w(i,j+1,k,ivz)
               right(irhoE) = p(i,j+1,k)
               if( correctForK ) right(itu1) = w(i,j+1,k,itu1)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it to the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
               fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
               fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
               fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
               fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyJ(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyJ(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

         ! Fluxes in k-direction.

         do k=1,kl
           do j=2,jl
             do i=2,il

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
               por = porK(i,j,k)
               if( addGridVelocities ) sFace = sFaceK(i,j,k)

               ! Determine the left and right state.

               left(irho)  = w(i,j,k,irho)
               left(ivx)   = w(i,j,k,ivx)
               left(ivy)   = w(i,j,k,ivy)
               left(ivz)   = w(i,j,k,ivz)
               left(irhoE) = p(i,j,k)
               if( correctForK ) left(itu1) = w(i,j,k,itu1)

               right(irho)  = w(i,j,k+1,irho)
               right(ivx)   = w(i,j,k+1,ivx)
               right(ivy)   = w(i,j,k+1,ivy)
               right(ivz)   = w(i,j,k+1,ivz)
               right(irhoE) = p(i,j,k+1)
               if( correctForK ) right(itu1) = w(i,j,k+1,itu1)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
               fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
               fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
               fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
               fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyK(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyK(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

!      ==================================================================

       else orderTest

!      ==================================================================
!
!        ****************************************************************
!        *                                                              *
!        * Second order reconstruction of the left and right state.     *
!        * The three differences used in the, possibly nonlinear,       *
!        * interpolation are constructed here; the actual left and      *
!        * right states, or at least the differences from the first     *
!        * order interpolation, are computed in the subroutine          *
!        * leftRightState.                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Fluxes in the i-direction.

         do k=2,kl
           do j=2,jl
             do i=1,il

               ! Store the three differences used in the interpolation
               ! in du1, du2, du3.

               du1(irho) = w(i,  j,k,irho) - w(i-1,j,k,irho)
               du2(irho) = w(i+1,j,k,irho) - w(i,  j,k,irho)
               du3(irho) = w(i+2,j,k,irho) - w(i+1,j,k,irho)

               du1(ivx) = w(i,  j,k,ivx) - w(i-1,j,k,ivx)
               du2(ivx) = w(i+1,j,k,ivx) - w(i,  j,k,ivx)
               du3(ivx) = w(i+2,j,k,ivx) - w(i+1,j,k,ivx)

               du1(ivy) = w(i,  j,k,ivy) - w(i-1,j,k,ivy)
               du2(ivy) = w(i+1,j,k,ivy) - w(i,  j,k,ivy)
               du3(ivy) = w(i+2,j,k,ivy) - w(i+1,j,k,ivy)

               du1(ivz) = w(i,  j,k,ivz) - w(i-1,j,k,ivz)
               du2(ivz) = w(i+1,j,k,ivz) - w(i,  j,k,ivz)
               du3(ivz) = w(i+2,j,k,ivz) - w(i+1,j,k,ivz)

               du1(irhoE) = p(i,  j,k) - p(i-1,j,k)
               du2(irhoE) = p(i+1,j,k) - p(i,  j,k)
               du3(irhoE) = p(i+2,j,k) - p(i+1,j,k)

               if( correctForK ) then
                 du1(itu1) = w(i,  j,k,itu1) - w(i-1,j,k,itu1)
                 du2(itu1) = w(i+1,j,k,itu1) - w(i,  j,k,itu1)
                 du3(itu1) = w(i+2,j,k,itu1) - w(i+1,j,k,itu1)
               endif

               ! Compute the differences from the first order scheme.

               call leftRightState(du1, du2, du3, rotMatrixI, &
                                   left, right)

               ! Add the first order part to the currently stored
               ! differences, such that the correct state vector
               ! is stored.

               left(irho)  = left(irho)  + w(i,j,k,irho)
               left(ivx)   = left(ivx)   + w(i,j,k,ivx)
               left(ivy)   = left(ivy)   + w(i,j,k,ivy)
               left(ivz)   = left(ivz)   + w(i,j,k,ivz)
               left(irhoE) = left(irhoE) + p(i,j,k)

               right(irho)  = right(irho)  + w(i+1,j,k,irho)
               right(ivx)   = right(ivx)   + w(i+1,j,k,ivx)
               right(ivy)   = right(ivy)   + w(i+1,j,k,ivy)
               right(ivz)   = right(ivz)   + w(i+1,j,k,ivz)
               right(irhoE) = right(irhoE) + p(i+1,j,k)

               if( correctForK ) then
                 left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                 right(itu1) = right(itu1) + w(i+1,j,k,itu1)
               endif

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = si(i,j,k,1); sy = si(i,j,k,2); sz = si(i,j,k,3)
               por = porI(i,j,k)
               if( addGridVelocities ) sFace = sFaceI(i,j,k)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i+1,j,k))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it to the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i+1,j,k,irho)  = fw(i+1,j,k,irho)  - flux(irho)
               fw(i+1,j,k,imx)   = fw(i+1,j,k,imx)   - flux(imx)
               fw(i+1,j,k,imy)   = fw(i+1,j,k,imy)   - flux(imy)
               fw(i+1,j,k,imz)   = fw(i+1,j,k,imz)   - flux(imz)
               fw(i+1,j,k,irhoE) = fw(i+1,j,k,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyI(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyI(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

         ! Fluxes in the j-direction.

         do k=2,kl
           do j=1,jl
             do i=2,il

               ! Store the three differences used in the interpolation
               ! in du1, du2, du3.

               du1(irho) = w(i,  j,k,irho) - w(i,j-1,k,irho)
               du2(irho) = w(i,j+1,k,irho) - w(i,  j,k,irho)
               du3(irho) = w(i,j+2,k,irho) - w(i,j+1,k,irho)

               du1(ivx) = w(i,  j,k,ivx) - w(i,j-1,k,ivx)
               du2(ivx) = w(i,j+1,k,ivx) - w(i,  j,k,ivx)
               du3(ivx) = w(i,j+2,k,ivx) - w(i,j+1,k,ivx)

               du1(ivy) = w(i,  j,k,ivy) - w(i,j-1,k,ivy)
               du2(ivy) = w(i,j+1,k,ivy) - w(i,  j,k,ivy)
               du3(ivy) = w(i,j+2,k,ivy) - w(i,j+1,k,ivy)

               du1(ivz) = w(i,  j,k,ivz) - w(i,j-1,k,ivz)
               du2(ivz) = w(i,j+1,k,ivz) - w(i,  j,k,ivz)
               du3(ivz) = w(i,j+2,k,ivz) - w(i,j+1,k,ivz)

               du1(irhoE) = p(i,  j,k) - p(i,j-1,k)
               du2(irhoE) = p(i,j+1,k) - p(i,  j,k)
               du3(irhoE) = p(i,j+2,k) - p(i,j+1,k)

               if( correctForK ) then
                 du1(itu1) = w(i,  j,k,itu1) - w(i,j-1,k,itu1)
                 du2(itu1) = w(i,j+1,k,itu1) - w(i,  j,k,itu1)
                 du3(itu1) = w(i,j+2,k,itu1) - w(i,j+1,k,itu1)
               endif

               ! Compute the differences from the first order scheme.

               call leftRightState(du1, du2, du3, rotMatrixJ, &
                                   left, right)

               ! Add the first order part to the currently stored
               ! differences, such that the correct state vector
               ! is stored.

               left(irho)  = left(irho)  + w(i,j,k,irho)
               left(ivx)   = left(ivx)   + w(i,j,k,ivx)
               left(ivy)   = left(ivy)   + w(i,j,k,ivy)
               left(ivz)   = left(ivz)   + w(i,j,k,ivz)
               left(irhoE) = left(irhoE) + p(i,j,k)

               right(irho)  = right(irho)  + w(i,j+1,k,irho)
               right(ivx)   = right(ivx)   + w(i,j+1,k,ivx)
               right(ivy)   = right(ivy)   + w(i,j+1,k,ivy)
               right(ivz)   = right(ivz)   + w(i,j+1,k,ivz)
               right(irhoE) = right(irhoE) + p(i,j+1,k)

               if( correctForK ) then
                 left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                 right(itu1) = right(itu1) + w(i,j+1,k,itu1)
               endif

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = sj(i,j,k,1); sy = sj(i,j,k,2); sz = sj(i,j,k,3)
               por = porJ(i,j,k)
               if( addGridVelocities ) sFace = sFaceJ(i,j,k)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i,j+1,k))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it to the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i,j+1,k,irho)  = fw(i,j+1,k,irho)  - flux(irho)
               fw(i,j+1,k,imx)   = fw(i,j+1,k,imx)   - flux(imx)
               fw(i,j+1,k,imy)   = fw(i,j+1,k,imy)   - flux(imy)
               fw(i,j+1,k,imz)   = fw(i,j+1,k,imz)   - flux(imz)
               fw(i,j+1,k,irhoE) = fw(i,j+1,k,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyJ(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyJ(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

         ! Fluxes in the k-direction.

         do k=1,kl
           do j=2,jl
             do i=2,il

               ! Store the three differences used in the interpolation
               ! in du1, du2, du3.

               du1(irho) = w(i,  j,k,irho) - w(i,j,k-1,irho)
               du2(irho) = w(i,j,k+1,irho) - w(i,  j,k,irho)
               du3(irho) = w(i,j,k+2,irho) - w(i,j,k+1,irho)

               du1(ivx) = w(i,  j,k,ivx) - w(i,j,k-1,ivx)
               du2(ivx) = w(i,j,k+1,ivx) - w(i,  j,k,ivx)
               du3(ivx) = w(i,j,k+2,ivx) - w(i,j,k+1,ivx)

               du1(ivy) = w(i,  j,k,ivy) - w(i,j,k-1,ivy)
               du2(ivy) = w(i,j,k+1,ivy) - w(i,  j,k,ivy)
               du3(ivy) = w(i,j,k+2,ivy) - w(i,j,k+1,ivy)

               du1(ivz) = w(i,  j,k,ivz) - w(i,j,k-1,ivz)
               du2(ivz) = w(i,j,k+1,ivz) - w(i,  j,k,ivz)
               du3(ivz) = w(i,j,k+2,ivz) - w(i,j,k+1,ivz)

               du1(irhoE) = p(i,  j,k) - p(i,j,k-1)
               du2(irhoE) = p(i,j,k+1) - p(i,  j,k)
               du3(irhoE) = p(i,j,k+2) - p(i,j,k+1)

               if( correctForK ) then
                 du1(itu1) = w(i,  j,k,itu1) - w(i,j,k-1,itu1)
                 du2(itu1) = w(i,j,k+1,itu1) - w(i,  j,k,itu1)
                 du3(itu1) = w(i,j,k+2,itu1) - w(i,j,k+1,itu1)
               endif

               ! Compute the differences from the first order scheme.

               call leftRightState(du1, du2, du3, rotMatrixK, &
                                   left, right)

               ! Add the first order part to the currently stored
               ! differences, such that the correct state vector
               ! is stored.

               left(irho)  = left(irho)  + w(i,j,k,irho)
               left(ivx)   = left(ivx)   + w(i,j,k,ivx)
               left(ivy)   = left(ivy)   + w(i,j,k,ivy)
               left(ivz)   = left(ivz)   + w(i,j,k,ivz)
               left(irhoE) = left(irhoE) + p(i,j,k)

               right(irho)  = right(irho)  + w(i,j,k+1,irho)
               right(ivx)   = right(ivx)   + w(i,j,k+1,ivx)
               right(ivy)   = right(ivy)   + w(i,j,k+1,ivy)
               right(ivz)   = right(ivz)   + w(i,j,k+1,ivz)
               right(irhoE) = right(irhoE) + p(i,j,k+1)

               if( correctForK ) then
                 left(itu1)  = left(itu1)  + w(i,j,k,itu1)
                 right(itu1) = right(itu1) + w(i,j,k+1,itu1)
               endif

               ! Store the normal vector, the porosity and the
               ! mesh velocity if present.

               sx = sk(i,j,k,1); sy = sk(i,j,k,2); sz = sk(i,j,k,3)
               por = porK(i,j,k)
               if( addGridVelocities ) sFace = sFaceK(i,j,k)

               ! Compute the value of gamma on the face. Take an
               ! arithmetic average of the two states.

               gammaFace = half*(gamma(i,j,k) + gamma(i,j,k+1))

               ! Compute the dissipative flux across the interface.

               call riemannFlux(left, right, flux)

               ! And scatter it to the left and right.

               fw(i,j,k,irho)  = fw(i,j,k,irho)  + flux(irho)
               fw(i,j,k,imx)   = fw(i,j,k,imx)   + flux(imx)
               fw(i,j,k,imy)   = fw(i,j,k,imy)   + flux(imy)
               fw(i,j,k,imz)   = fw(i,j,k,imz)   + flux(imz)
               fw(i,j,k,irhoE) = fw(i,j,k,irhoE) + flux(irhoE)

               fw(i,j,k+1,irho)  = fw(i,j,k+1,irho)  - flux(irho)
               fw(i,j,k+1,imx)   = fw(i,j,k+1,imx)   - flux(imx)
               fw(i,j,k+1,imy)   = fw(i,j,k+1,imy)   - flux(imy)
               fw(i,j,k+1,imz)   = fw(i,j,k+1,imz)   - flux(imz)
               fw(i,j,k+1,irhoE) = fw(i,j,k+1,irhoE) - flux(irhoE)

               ! Store the density flux in the mass flow of the
               ! appropriate sliding mesh interface.

               ind = indFamilyK(i,j,k)
               massFlowFamilyDiss(ind,spectralSol) =       &
                       massFlowFamilyDiss(ind,spectralSol) &
                                        + factFamilyK(i,j,k)*flux(irho)
             enddo
           enddo
         enddo

       endif orderTest

!      ==================================================================

       contains

         subroutine leftRightState(du1, du2, du3, rotMatrix, left, right)
!
!        ****************************************************************
!        *                                                              *
!        * leftRightState computes the differences in the left and      *
!        * right state compared to the first order interpolation. For a *
!        * monotonic second order discretization the interpolations     *
!        * need to be nonlinear. The linear second order scheme can be  *
!        * stable (depending on the value of kappa), but it will have   *
!        * oscillations near discontinuities.                           *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local parameter.
!
         real(kind=realType), parameter :: epsLim = 1.e-10_realType
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(:), intent(inout) :: du1, du2, du3
         real(kind=realType), dimension(:), intent(out)   :: left, right

         real(kind=realType), dimension(:,:,:,:,:), pointer :: rotMatrix
!
!        Local variables.
!
         integer(kind=intType) :: l

         real(kind=realType) :: rl1, rl2, rr1, rr2, tmp, dvx, dvy, dvz

         real(kind=realType), dimension(3,3) :: rot
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Check if the velocity components should be transformed to
         ! the cylindrical frame.

         if( rotationalPeriodic ) then

           ! Store the rotation matrix a bit easier. Note that the i,j,k
           ! come from the main subroutine.

           rot(1,1) = rotMatrix(i,j,k,1,1)
           rot(1,2) = rotMatrix(i,j,k,1,2)
           rot(1,3) = rotMatrix(i,j,k,1,3)

           rot(2,1) = rotMatrix(i,j,k,2,1)
           rot(2,2) = rotMatrix(i,j,k,2,2)
           rot(2,3) = rotMatrix(i,j,k,2,3)

           rot(3,1) = rotMatrix(i,j,k,3,1)
           rot(3,2) = rotMatrix(i,j,k,3,2)
           rot(3,3) = rotMatrix(i,j,k,3,3)

           ! Apply the transformation to the velocity components
           ! of du1, du2 and du3.

           dvx = du1(ivx); dvy = du1(ivy); dvz = du1(ivz)
           du1(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
           du1(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
           du1(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

           dvx = du2(ivx); dvy = du2(ivy); dvz = du2(ivz)
           du2(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
           du2(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
           du2(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

           dvx = du3(ivx); dvy = du3(ivy); dvz = du3(ivz)
           du3(ivx) = rot(1,1)*dvx + rot(1,2)*dvy + rot(1,3)*dvz
           du3(ivy) = rot(2,1)*dvx + rot(2,2)*dvy + rot(2,3)*dvz
           du3(ivz) = rot(3,1)*dvx + rot(3,2)*dvy + rot(3,3)*dvz

         endif

         ! Determine the limiter used.

         select case (limUsed)

           case (noLimiter)

             ! Linear interpolation; no limiter.
             ! Loop over the number of variables to be interpolated.

             do l=1,nwInt
               left(l)  =  omk*du1(l) + opk*du2(l)
               right(l) = -omk*du3(l) - opk*du2(l)
             enddo

!          ==============================================================

           case (vanAlbeda)

             ! Nonlinear interpolation using the van albeda limiter.
             ! Loop over the number of variables to be interpolated.

             do l=1,nwInt

               ! Compute the limiter argument rl1, rl2, rr1 and rr2.
               ! Note the cut off to 0.0.

               tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
               rl1 = max(zero, &
                         du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
               rl2 = max(zero,du1(l)*tmp)

               rr1 = max(zero,du3(l)*tmp)
               rr2 = max(zero, &
                         du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

               ! Compute the corresponding limiter values.

               rl1 = rl1*(rl1 + one)/(rl1*rl1 + one)
               rl2 = rl2*(rl2 + one)/(rl2*rl2 + one)
               rr1 = rr1*(rr1 + one)/(rr1*rr1 + one)
               rr2 = rr2*(rr2 + one)/(rr2*rr2 + one)

               ! Compute the nonlinear corrections to the first order
               ! scheme.

               left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
               right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

             enddo

!          ==============================================================

           case (minmod)

             ! Nonlinear interpolation using the minmod limiter.
             ! Loop over the number of variables to be interpolated.

             do l=1,nwInt

               ! Compute the limiter argument rl1, rl2, rr1 and rr2.
               ! Note the cut off to 0.0.

               tmp = one/sign(max(abs(du2(l)),epsLim),du2(l))
               rl1 = max(zero, &
                         du2(l)/sign(max(abs(du1(l)),epsLim),du1(l)))
               rl2 = max(zero,du1(l)*tmp)

               rr1 = max(zero,du3(l)*tmp)
               rr2 = max(zero, &
                         du2(l)/sign(max(abs(du3(l)),epsLim),du3(l)))

               ! Compute the corresponding limiter values.

               rl1 = min(one, factMinmod*rl1)
               rl2 = min(one, factMinmod*rl2)
               rr1 = min(one, factMinmod*rr1)
               rr2 = min(one, factMinmod*rr2)

               ! Compute the nonlinear corrections to the first order
               ! scheme.

               left(l)  =  omk*rl1*du1(l) + opk*rl2*du2(l)
               right(l) = -opk*rr1*du2(l) - omk*rr2*du3(l)

             enddo

         end select

         ! In case only a first order scheme must be used for the
         ! turbulent transport equations, set the correction for the
         ! turbulent kinetic energy to 0.

         if( firstOrderK ) then
           left(itu1)  = zero
           right(itu1) = zero
         endif

         ! For rotational periodic problems transform the velocity
         ! differences back to Cartesian again. Note that now the
         ! transpose of the rotation matrix must be used.

         if( rotationalPeriodic ) then

           ! Left state.

           dvx = left(ivx); dvy = left(ivy); dvz = left(ivz)
           left(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
           left(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
           left(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

           ! Right state.

           dvx = right(ivx); dvy = right(ivy); dvz = right(ivz)
           right(ivx) = rot(1,1)*dvx + rot(2,1)*dvy + rot(3,1)*dvz
           right(ivy) = rot(1,2)*dvx + rot(2,2)*dvy + rot(3,2)*dvz
           right(ivz) = rot(1,3)*dvx + rot(2,3)*dvy + rot(3,3)*dvz

         endif

         end subroutine leftRightState

!        ================================================================

         subroutine riemannFlux(left, right, flux)
!
!        ****************************************************************
!        *                                                              *
!        * riemannFlux computes the flux for the given face and left    *
!        * and right states.                                            *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         real(kind=realType), dimension(*), intent(in)  :: left, right
         real(kind=realType), dimension(*), intent(out) :: flux
!
!        Local variables.
!
         real(kind=realType) :: porFlux, rFace
         real(kind=realType) :: Etl, Etr, z1l, z1r, tmp
         real(kind=realType) :: dr, dru, drv, drw, drE, drk
         real(kind=realType) :: rAvg, uAvg, vAvg, wAvg, hAvg, kAvg
         real(kind=realType) :: alphaAvg, a2Avg, aAvg, unAvg
         real(kind=realType) :: ovaAvg, ova2Avg, area, Eta
         real(kind=realType) :: gm1, gm53
         real(kind=realType) :: lam1, lam2, lam3
         real(kind=realType) :: abv1, abv2, abv3, abv4, abv5, abv6, abv7

         real(kind=realType), dimension(2) :: rhotmp, utmp, vtmp, wtmp
         real(kind=realType), dimension(2) :: ptmp, ktmp, Etmp
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Set the porosity for the flux. The default value, 0.5*rFil, is
         ! a scaling factor where an rFil != 1 is taken into account.

         porFlux = half*rFil
         if(por == noFlux .or. por == boundFlux)    porFlux = zero

         ! Abbreviate some expressions in which gamma occurs.

         gm1  = gammaFace - one
         gm53 = gammaFace - five*third

         ! Determine which riemann solver must be solved.

         select case (riemannUsed)

           case (Roe)

             ! Determine the preconditioner used.

             select case (precond)

               case (noPrecond)

                 ! No preconditioner used. Use the Roe scheme of the
                 ! standard equations.

                 ! Compute the square root of the left and right densities
                 ! and the inverse of the sum.

                 z1l = sqrt(left(irho))
                 z1r = sqrt(right(irho))
                 tmp = one/(z1l + z1r)

                 ! Compute some variables depending whether or not a
                 ! k-equation is present.

                 if( correctForK ) then

                   ! Store the left and right kinetic energy in ktmp,
                   ! which is needed to compute the total energy.

                   ktmp(1) = left(itu1)
                   ktmp(2) = right(itu1)

                   ! Store the difference of the turbulent kinetic energy
                   ! per unit volume, i.e. the conserved variable.

                   drk = right(irho)*right(itu1) - left(irho)*left(itu1)

                   ! Compute the average turbulent energy per unit mass
                   ! using Roe averages.

                   kAvg = tmp*(z1l*left(itu1) + z1r*right(itu1))

                 else

                   ! Set the difference of the turbulent kinetic energy
                   ! per unit volume and the averaged kinetic energy per
                   ! unit mass to zero.

                   drk  = 0.0
                   kAvg = 0.0

                 endif

                 ! Compute the total energy of the left and right state.

                 rhotmp(1) = left(irho)
                 rhotmp(2) = right(irho)

                 utmp(1) = left(ivx)
                 utmp(2) = right(ivx)

                 vtmp(1) = left(ivy)
                 vtmp(2) = right(ivy)

                 wtmp(1) = left(ivz)
                 wtmp(2) = right(ivz)

                 ptmp(1) = left(irhoE)
                 ptmp(2) = right(irhoE)

                 call etotArray(rhotmp, utmp, vtmp, wtmp, ptmp, ktmp, &
                                Etmp, correctForK, 2_intType)

                 Etl = Etmp(1)
                 Etr = Etmp(2)

                 ! Compute the difference of the conservative mean
                 ! flow variables.

                 dr  = right(irho) - left(irho)
                 dru = right(irho)*right(ivx) - left(irho)*left(ivx)
                 drv = right(irho)*right(ivy) - left(irho)*left(ivy)
                 drw = right(irho)*right(ivz) - left(irho)*left(ivz)
                 drE = Etr - Etl

                 ! Compute the Roe average variables, which can be
                 ! computed directly from the average Roe vector.

                 rAvg = fourth*(z1r + z1l)**2
                 uAvg = tmp*(z1l*left(ivx) + z1r*right(ivx))
                 vAvg = tmp*(z1l*left(ivy) + z1r*right(ivy))
                 wAvg = tmp*(z1l*left(ivz) + z1r*right(ivz))
                 hAvg = tmp*((Etl+left(irhoE)) /z1l &
                      +      (Etr+right(irhoE))/z1r)

                 ! Compute the unit vector and store the area of the
                 ! normal. Also compute the unit normal velocity of the face.

                 area  = sqrt(sx**2 + sy**2 + sz**2)
                 tmp   = one/max(1.e-25_realType,area)
                 sx    = sx*tmp
                 sy    = sy*tmp
                 sz    = sz*tmp
                 rFace = sFace*tmp

                 ! Compute some dependent variables at the Roe
                 ! average state.

                 alphaAvg = half*(uAvg**2 + vAvg**2 + wAvg**2)
                 a2Avg    = abs(gm1*(hAvg - alphaAvg) - gm53*kAvg)
                 aAvg     = sqrt(a2Avg)
                 unAvg    = uAvg*sx + vAvg*sy + wAvg*sz

                 ovaAvg  = one/aAvg
                 ova2Avg = one/a2Avg

                 ! Set for a boundary the normal velocity to rFace, the
                 ! normal velocity of the boundary.

                 if(por == boundFlux) unAvg = rFace

                 ! Compute the coefficient eta for the entropy correction.
                 ! At the moment a 1D entropy correction is used, which
                 ! removes expansion shocks. Although it also reduces the
                 ! carbuncle phenomenon, it does not remove it completely.
                 ! In other to do that a multi-dimensional entropy fix is
                 ! needed, see Sanders et. al, JCP, vol. 145, 1998,
                 ! pp. 511 - 537. Although relatively easy to implement,
                 ! an efficient implementation requires the storage of
                 ! all the left and right states, which is rather
                 ! expensive in terms of memory.

                 eta = half*(abs((left(ivx) - right(ivx))*sx        &
                     +           (left(ivy) - right(ivy))*sy        &
                     +           (left(ivz) - right(ivz))*sz)       &
                     +       abs(sqrt(gammaFace*left(irhoE)/left(irho)) &
                     -           sqrt(gammaFace*right(irhoE)/right(irho))))

                 ! Compute the absolute values of the three eigenvalues.

                 lam1 = abs(unAvg - rFace + aAvg)
                 lam2 = abs(unAvg - rFace - aAvg)
                 lam3 = abs(unAvg - rFace)

                 ! Apply the entropy correction to the eigenvalues.

                 tmp = two*eta
                 if(lam1 < tmp) lam1 = eta + fourth*lam1*lam1/eta
                 if(lam2 < tmp) lam2 = eta + fourth*lam2*lam2/eta
                 if(lam3 < tmp) lam3 = eta + fourth*lam3*lam3/eta

                 ! Multiply the eigenvalues by the area to obtain
                 ! the correct values for the dissipation term.

                 lam1 = lam1*area
                 lam2 = lam2*area
                 lam3 = lam3*area

                 ! Some abbreviations, which occur quite often in the
                 ! dissipation terms.

                 abv1 = half*(lam1 + lam2)
                 abv2 = half*(lam1 - lam2)
                 abv3 = abv1 - lam3

                 abv4 = gm1*(alphaAvg*dr - uAvg*dru -vAvg*drv &
                      -      wAvg*drw + drE) - gm53*drk
                 abv5 = sx*dru + sy*drv + sz*drw - unAvg*dr

                 abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg
                 abv7 = abv2*abv4*ovaAvg  + abv3*abv5

                 ! Compute the dissipation term, -|a| (wr - wl), which is
                 ! multiplied by porFlux. Note that porFlux is either
                 ! 0.0 or 0.5*rFil.

                 flux(irho)  = -porFlux*(lam3*dr  + abv6)
                 flux(imx)   = -porFlux*(lam3*dru + uAvg*abv6 &
                                                  + sx*abv7)
                 flux(imy)   = -porFlux*(lam3*drv + vAvg*abv6 &
                                                  + sy*abv7)
                 flux(imz)   = -porFlux*(lam3*drw + wAvg*abv6 &
                                                  + sz*abv7)
                 flux(irhoE) = -porFlux*(lam3*drE + hAvg*abv6 &
                                                  + unAvg*abv7)

      !          tmp = max(lam1,lam2,lam3)

      !          flux(irho)  = -porFlux*(tmp*dr)
      !          flux(imx)   = -porFlux*(tmp*dru)
      !          flux(imy)   = -porFlux*(tmp*drv)
      !          flux(imz)   = -porFlux*(tmp*drw)
      !          flux(irhoE) = -porFlux*(tmp*drE)

               case (Turkel)
                 call terminate("riemannFlux", &
                                "Turkel preconditioner not &
                                &implemented yet")

               case (ChoiMerkle)
                 call terminate("riemannFlux", &
                                "choi merkle preconditioner not &
                                &implemented yet")

             end select

           case (vanLeer)
             call terminate("riemannFlux", &
                            "van leer fvs not implemented yet")

           case (ausmdv)
             call terminate("riemannFlux", &
                            "ausmdv fvs not implemented yet")

         end select

         end subroutine riemannFlux

       end subroutine inviscidUpwindFlux
