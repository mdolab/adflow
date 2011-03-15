!
!      ******************************************************************
!      *                                                                *
!      * File:          inviscidUpwindFluxAdj.f90                       *
!      * Author:        Edwin van der Weide                             *
!      *                Seongim Choi,C.A.(Sandy)Mader                   *
!      * Starting date: 03-20-2006                                      *
!      * Last modified: 04-25-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine inviscidUpwindFluxAdj(wAdj,  pAdj,  dwAdj, &
                                        siAdj, sjAdj, skAdj, &
                                        sFaceIAdj,sFaceJAdj,sFaceKAdj,&
                                        iCell, jCell, kCell,finegrid,nn,level,sps)
!
!      ******************************************************************
!      *                                                                *
!      * inviscidUpwindFluxAdj computes the artificial dissipation part *
!      * the Euler fluxes by means of an approximate solution of the 1D *
!      * Riemann problem on the face. The fluxes are computed for the   *
!      * given cell of the block to which the variables in              *
!      * blockPointers currently point to.                              *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers        ! sI,sJ,sK
       use inputDiscretization  ! limiter, firstOrder
       use inputPhysics
       use flowVarRefState
       use inputTimeSpectral !nTimeIntervalsSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType) :: iCell, jCell, kCell,nn,level,sps

       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
                                                      intent(in) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral),    &
                                                      intent(in) :: pAdj
       real(kind=realType), dimension(nw,nTimeIntervalsSpectral), intent(inout) :: dwAdj
       real(kind=realType), dimension(-3:2,-3:2,-3:2,3,nTimeIntervalsSpectral), intent(in) :: siAdj, sjAdj, skAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral), intent(in) ::sFaceIAdj,sFaceJAdj,sFaceKAdj
!
!      Local variables.
!
       integer(kind=porType) :: por

       integer(kind=intType) :: nwInt
       integer(kind=intType) :: i, j, k, ii, jj, kk

       real(kind=realType) :: sx, sy, sz, omk, opk, gammaFace,gammaface2
       real(kind=realType) :: factMinmod, sFace, fact

       real(kind=realType), dimension(nw)  :: left, right
       real(kind=realType), dimension(nw)  :: du1, du2, du3
       real(kind=realType), dimension(nwf) :: flux

       logical :: firstOrderK, correctForK,finegrid

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine whether or not the total energy must be corrected
       ! for the presence of the turbulent kinetic energy.

       correctForK = kPresent

       ! Compute the factor used in the minmod limiter.

       factMinmod = (three-kappaCoef) &
                   / max(1.e-10_realType, one-kappaCoef)

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
       
       orderTest: if(limiter == firstOrder) then
!
!        ****************************************************************
!        *                                                              *
!        * First order reconstruction. The states in the cells are      *
!        * constant. The left and right states are constructed easily.  *
!        *                                                              *
!        ****************************************************************
!
         ! Fluxes in the i-direction.

         i    = iCell-1; j = jCell; k = kCell
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do ii=-1,0

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = siAdj(ii,0,0,1,sps); sy = siAdj(ii,0,0,2,sps); sz = siAdj(ii,0,0,3,sps)
           por = porI(i,j,k)
           if( addGridVelocities ) sFace = sFaceIAdj(ii,0,0,sps)

           ! Determine the left and right state.

           left(irho)  = wAdj(ii,0,0,irho,sps)
           left(ivx)   = wAdj(ii,0,0,ivx,sps)
           left(ivy)   = wAdj(ii,0,0,ivy,sps)
           left(ivz)   = wAdj(ii,0,0,ivz,sps)
           left(irhoE) = pAdj(ii,0,0,sps)
           if( correctForK ) left(itu1) = wAdj(ii,0,0,itu1,sps)

           right(irho)  = wAdj(ii+1,0,0,irho,sps)
           right(ivx)   = wAdj(ii+1,0,0,ivx,sps)
           right(ivy)   = wAdj(ii+1,0,0,ivy,sps)
           right(ivz)   = wAdj(ii+1,0,0,ivz,sps)
           right(irhoE) = pAdj(ii+1,0,0,sps)
           if( correctForK ) right(itu1) = wAdj(ii+1,0,0,itu1,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant
           gammaFace2 = half*(gamma(i,j,k) + gamma(i+1,j,k))
           !print *,'gammaface',gammaface,gammaface2

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.

!           call riemannFluxAdj(left,right,flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
            call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)

           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           ! Update i and set fact to 1 for the second face.

           i    = i + 1
           fact = one

         enddo

         ! Fluxes in j-direction.

         i    = iCell; j = jCell-1; k = kCell
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do jj=-1,0

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = sjAdj(0,jj,0,1,sps); sy = sjAdj(0,jj,0,2,sps); sz = sjAdj(0,jj,0,3,sps)
           por = porJ(i,j,k)
           if( addGridVelocities ) sFace = sFaceJAdj(0,jj,0,sps)

           ! Determine the left and right state.

           left(irho)  = wAdj(0,jj,0,irho,sps)
           left(ivx)   = wAdj(0,jj,0,ivx,sps)
           left(ivy)   = wAdj(0,jj,0,ivy,sps)
           left(ivz)   = wAdj(0,jj,0,ivz,sps)
           left(irhoE) = pAdj(0,jj,0,sps)
           if( correctForK ) left(itu1) = wAdj(0,jj,0,itu1,sps)

           right(irho)  = wAdj(0,jj+1,0,irho,sps)
           right(ivx)   = wAdj(0,jj+1,0,ivx,sps)
           right(ivy)   = wAdj(0,jj+1,0,ivy,sps)
           right(ivz)   = wAdj(0,jj+1,0,ivz,sps)
           right(irhoE) = pAdj(0,jj+1,0,sps)
           if( correctForK ) right(itu1) = wAdj(0,jj+1,0,itu1,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.

!           call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
            call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)

           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           ! Update j and set fact to 1 for the second face.

           j    = j + 1
           fact = one

         enddo

         ! Fluxes in k-direction.

         i    = iCell; j = jCell; k = kCell-1
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do kk=-1,0

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = skAdj(0,0,kk,1,sps); sy = skAdj(0,0,kk,2,sps); sz = skAdj(0,0,kk,3,sps)
           por = porK(i,j,k)
           if( addGridVelocities ) sFace = sFaceKAdj(0,0,kk,sps)

           ! Determine the left and right state.

           left(irho)  = wAdj(0,0,kk,irho,sps)
           left(ivx)   = wAdj(0,0,kk,ivx,sps)
           left(ivy)   = wAdj(0,0,kk,ivy,sps)
           left(ivz)   = wAdj(0,0,kk,ivz,sps)
           left(irhoE) = pAdj(0,0,kk,sps)
           if( correctForK ) left(itu1) = wAdj(0,0,kk,itu1,sps)

           right(irho)  = wAdj(0,0,kk+1,irho,sps)
           right(ivx)   = wAdj(0,0,kk+1,ivx,sps)
           right(ivy)   = wAdj(0,0,kk+1,ivy,sps)
           right(ivz)   = wAdj(0,0,kk+1,ivz,sps)
           right(irhoE) = pAdj(0,0,kk+1,sps)
           if( correctForK ) right(itu1) = wAdj(0,0,kk+1,itu1,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.
           
           !call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
            call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)

           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           ! Update k and set fact to 1 for the second face.

           k    = k + 1
           fact = one

         enddo

!      ==================================================================

       else orderTest
          !PRINT *,'limiter',limiter
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

         i    = iCell-1; j = jCell; k = kCell
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do ii=-1,0

           ! Store the three differences used in the interpolation
           ! in du1, du2, du3.

           du1(irho) = wAdj(ii,  0,0,irho,sps) - wAdj(ii-1,0,0,irho,sps)
           du2(irho) = wAdj(ii+1,0,0,irho,sps) - wAdj(ii,  0,0,irho,sps)
           du3(irho) = wAdj(ii+2,0,0,irho,sps) - wAdj(ii+1,0,0,irho,sps)

           du1(ivx) = wAdj(ii,  0,0,ivx,sps) - wAdj(ii-1,0,0,ivx,sps)
           du2(ivx) = wAdj(ii+1,0,0,ivx,sps) - wAdj(ii,  0,0,ivx,sps)
           du3(ivx) = wAdj(ii+2,0,0,ivx,sps) - wAdj(ii+1,0,0,ivx,sps)

           du1(ivy) = wAdj(ii,  0,0,ivy,sps) - wAdj(ii-1,0,0,ivy,sps)
           du2(ivy) = wAdj(ii+1,0,0,ivy,sps) - wAdj(ii,  0,0,ivy,sps)
           du3(ivy) = wAdj(ii+2,0,0,ivy,sps) - wAdj(ii+1,0,0,ivy,sps)

           du1(ivz) = wAdj(ii,  0,0,ivz,sps) - wAdj(ii-1,0,0,ivz,sps)
           du2(ivz) = wAdj(ii+1,0,0,ivz,sps) - wAdj(ii,  0,0,ivz,sps)
           du3(ivz) = wAdj(ii+2,0,0,ivz,sps) - wAdj(ii+1,0,0,ivz,sps)

           du1(irhoE) = pAdj(ii,  0,0,sps) - pAdj(ii-1,0,0,sps)
           du2(irhoE) = pAdj(ii+1,0,0,sps) - pAdj(ii,  0,0,sps)
           du3(irhoE) = pAdj(ii+2,0,0,sps) - pAdj(ii+1,0,0,sps)
           
!!$           print *,'pAdj',p(i,  j,k) - p(i-1,j,k),pAdj(ii,  0,0) - pAdj(ii-1,0,0),p(i,  j,k),p(i-1,j,k),pAdj(ii,  0,0), pAdj(ii-1,0,0),p(i,  j,k) -pAdj(ii,  0,0),p(i,  j,k) - p(i-1,j,k)-(pAdj(ii,  0,0) - pAdj(ii-1,0,0))
           if( correctForK ) then
             du1(itu1) = wAdj(ii,  0,0,itu1,sps) - wAdj(ii-1,0,0,itu1,sps)
             du2(itu1) = wAdj(ii+1,0,0,itu1,sps) - wAdj(ii,  0,0,itu1,sps)
             du3(itu1) = wAdj(ii+2,0,0,itu1,sps) - wAdj(ii+1,0,0,itu1,sps)
           endif

           ! Compute the differences from the first order scheme.

           call leftRightStateAdj(du1, du2, du3, left, right,nwInt,omk,opk,factminmod,firstOrderK)
           !print *,'leftrightadj',left,right,icell,jcell,kcell
           ! Add the first order part to the currently stored
           ! differences, such that the correct state vector
           ! is stored.

           left(irho)  = left(irho)  + wAdj(ii,0,0,irho,sps)
           !print *,'left',left(irho),wAdj(ii,0,0,irho)
           left(ivx)   = left(ivx)   + wAdj(ii,0,0,ivx,sps)
           left(ivy)   = left(ivy)   + wAdj(ii,0,0,ivy,sps)
           left(ivz)   = left(ivz)   + wAdj(ii,0,0,ivz,sps)
           left(irhoE) = left(irhoE) + pAdj(ii,0,0,sps)

           right(irho)  = right(irho)  + wAdj(ii+1,0,0,irho,sps)
           right(ivx)   = right(ivx)   + wAdj(ii+1,0,0,ivx,sps)
           right(ivy)   = right(ivy)   + wAdj(ii+1,0,0,ivy,sps)
           right(ivz)   = right(ivz)   + wAdj(ii+1,0,0,ivz,sps)
           right(irhoE) = right(irhoE) + pAdj(ii+1,0,0,sps)

           if( correctForK ) then
             left(itu1)  = left(itu1)  + wAdj(ii,0,0,itu1,sps)
             right(itu1) = right(itu1) + wAdj(ii+1,0,0,itu1,sps)
           endif

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = siAdj(ii,0,0,1,sps); sy = siAdj(ii,0,0,2,sps); sz = siAdj(ii,0,0,3,sps)
           por = porI(i,j,k)
           if( addGridVelocities ) sFace = sFaceIAdj(ii,0,0,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant
           gammaFace2 = half*(gamma(i,j,k) + gamma(i+1,j,k))
           !print *,'gammaface',gammaface,gammaface2

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.
!           print *,'leftrightadj',left,right
!stop
!           call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
 !          print *,'riemanninputI',left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid
           call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)
  !         print *,'fluxadjI',flux,icell,jcell,kcell
           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           ! Update i and set fact to 1 for the second face.

           i    = i + 1
           fact = one

         enddo

         ! Fluxes in the j-direction.

         i    = iCell; j = jCell-1; k = kCell
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do jj=-1,0

           ! Store the three differences used in the interpolation
           ! in du1, du2, du3.

           du1(irho) = wAdj(0,jj,  0,irho,sps) - wAdj(0,jj-1,0,irho,sps)
           du2(irho) = wAdj(0,jj+1,0,irho,sps) - wAdj(0,jj,  0,irho,sps)
           du3(irho) = wAdj(0,jj+2,0,irho,sps) - wAdj(0,jj+1,0,irho,sps)

           du1(ivx) = wAdj(0,jj,  0,ivx,sps) - wAdj(0,jj-1,0,ivx,sps)
           du2(ivx) = wAdj(0,jj+1,0,ivx,sps) - wAdj(0,jj,  0,ivx,sps)
           du3(ivx) = wAdj(0,jj+2,0,ivx,sps) - wAdj(0,jj+1,0,ivx,sps)

           du1(ivy) = wAdj(0,jj,  0,ivy,sps) - wAdj(0,jj-1,0,ivy,sps)
           du2(ivy) = wAdj(0,jj+1,0,ivy,sps) - wAdj(0,jj,  0,ivy,sps)
           du3(ivy) = wAdj(0,jj+2,0,ivy,sps) - wAdj(0,jj+1,0,ivy,sps)

           du1(ivz) = wAdj(0,jj,  0,ivz,sps) - wAdj(0,jj-1,0,ivz,sps)
           du2(ivz) = wAdj(0,jj+1,0,ivz,sps) - wAdj(0,jj,  0,ivz,sps)
           du3(ivz) = wAdj(0,jj+2,0,ivz,sps) - wAdj(0,jj+1,0,ivz,sps)

           du1(irhoE) = pAdj(0,jj,  0,sps) - pAdj(0,jj-1,0,sps)
           du2(irhoE) = pAdj(0,jj+1,0,sps) - pAdj(0,jj,  0,sps)
           du3(irhoE) = pAdj(0,jj+2,0,sps) - pAdj(0,jj+1,0,sps)

           if( correctForK ) then
             du1(itu1) = wAdj(0,jj,  0,itu1,sps) - wAdj(0,jj-1,0,itu1,sps)
             du2(itu1) = wAdj(0,jj+1,0,itu1,sps) - wAdj(0,jj,  0,itu1,sps)
             du3(itu1) = wAdj(0,jj+2,0,itu1,sps) - wAdj(0,jj+1,0,itu1,sps)
           endif

           ! Compute the differences from the first order scheme.

           call leftRightStateAdj(du1, du2, du3, left, right,nwInt,omk,opk,factminmod,firstOrderK)

           ! Add the first order part to the currently stored
           ! differences, such that the correct state vector
           ! is stored.

           left(irho)  = left(irho)  + wAdj(0,jj,0,irho,sps)
           left(ivx)   = left(ivx)   + wAdj(0,jj,0,ivx,sps)
           left(ivy)   = left(ivy)   + wAdj(0,jj,0,ivy,sps)
           left(ivz)   = left(ivz)   + wAdj(0,jj,0,ivz,sps)
           left(irhoE) = left(irhoE) + pAdj(0,jj,0,sps)

           right(irho)  = right(irho)  + wAdj(0,jj+1,0,irho,sps)
           right(ivx)   = right(ivx)   + wAdj(0,jj+1,0,ivx,sps)
           right(ivy)   = right(ivy)   + wAdj(0,jj+1,0,ivy,sps)
           right(ivz)   = right(ivz)   + wAdj(0,jj+1,0,ivz,sps)
           right(irhoE) = right(irhoE) + pAdj(0,jj+1,0,sps)

           if( correctForK ) then
             left(itu1)  = left(itu1)  + wAdj(0,jj,0,itu1,sps)
             right(itu1) = right(itu1) + wAdj(0,jj+1,0,itu1,sps)
           endif

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = sjAdj(0,jj,0,1,sps); sy = sjAdj(0,jj,0,2,sps); sz = sjAdj(0,jj,0,3,sps)
           por = porJ(i,j,k)
           if( addGridVelocities ) sFace = sFaceJAdj(0,jj,0,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.

           !call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
            call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)

           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           ! Update j and set fact to 1 for the second face.

           j    = j + 1
           fact = one

         enddo

         ! Fluxes in the k-direction.

         i    = iCell; j = jCell; k = kCell-1
         fact = -one

         ! Loop over the two faces which contribute to the residual of
         ! the cell considered.

         do kk=-1,0

           ! Store the three differences used in the interpolation
           ! in du1, du2, du3.

           du1(irho) = wAdj(0,0,kk,  irho,sps) - wAdj(0,0,kk-1,irho,sps)
           du2(irho) = wAdj(0,0,kk+1,irho,sps) - wAdj(0,0,kk,  irho,sps)
           du3(irho) = wAdj(0,0,kk+2,irho,sps) - wAdj(0,0,kk+1,irho,sps)

           du1(ivx) = wAdj(0,0,kk,  ivx,sps) - wAdj(0,0,kk-1,ivx,sps)
           du2(ivx) = wAdj(0,0,kk+1,ivx,sps) - wAdj(0,0,kk,  ivx,sps)
           du3(ivx) = wAdj(0,0,kk+2,ivx,sps) - wAdj(0,0,kk+1,ivx,sps)

           du1(ivy) = wAdj(0,0,kk,  ivy,sps) - wAdj(0,0,kk-1,ivy,sps)
           du2(ivy) = wAdj(0,0,kk+1,ivy,sps) - wAdj(0,0,kk,  ivy,sps)
           du3(ivy) = wAdj(0,0,kk+2,ivy,sps) - wAdj(0,0,kk+1,ivy,sps)

           du1(ivz) = wAdj(0,0,kk,  ivz,sps) - wAdj(0,0,kk-1,ivz,sps)
           du2(ivz) = wAdj(0,0,kk+1,ivz,sps) - wAdj(0,0,kk,  ivz,sps)
           du3(ivz) = wAdj(0,0,kk+2,ivz,sps) - wAdj(0,0,kk+1,ivz,sps)

           du1(irhoE) = pAdj(0,0,kk,sps)   - pAdj(0,0,kk-1,sps)
           du2(irhoE) = pAdj(0,0,kk+1,sps) - pAdj(0,0,kk,sps)
           du3(irhoE) = pAdj(0,0,kk+2,sps) - pAdj(0,0,kk+1,sps)

           if( correctForK ) then
             du1(itu1) = wAdj(0,0,kk,  itu1,sps) - wAdj(0,0,kk-1,itu1,sps)
             du2(itu1) = wAdj(0,0,kk+1,itu1,sps) - wAdj(0,0,kk,  itu1,sps)
             du3(itu1) = wAdj(0,0,kk+2,itu1,sps) - wAdj(0,0,kk+1,itu1,sps)
           endif

           ! Compute the differences from the first order scheme.

           call leftRightStateAdj(du1, du2, du3, left, right,nwInt,omk,opk,factminmod,firstOrderK)

           ! Add the first order part to the currently stored
           ! differences, such that the correct state vector
           ! is stored.

           left(irho)  = left(irho)  + wAdj(0,0,kk,irho,sps)
           left(ivx)   = left(ivx)   + wAdj(0,0,kk,ivx,sps)
           left(ivy)   = left(ivy)   + wAdj(0,0,kk,ivy,sps)
           left(ivz)   = left(ivz)   + wAdj(0,0,kk,ivz,sps)
           left(irhoE) = left(irhoE) + pAdj(0,0,kk,sps)

           right(irho)  = right(irho)  + wAdj(0,0,kk+1,irho,sps)
           right(ivx)   = right(ivx)   + wAdj(0,0,kk+1,ivx,sps)
           right(ivy)   = right(ivy)   + wAdj(0,0,kk+1,ivy,sps)
           right(ivz)   = right(ivz)   + wAdj(0,0,kk+1,ivz,sps)
           right(irhoE) = right(irhoE) + pAdj(0,0,kk+1,sps)

           if( correctForK ) then
             left(itu1)  = left(itu1)  + wAdj(0,0,kk,itu1,sps)
             right(itu1) = right(itu1) + wAdj(0,0,kk+1,itu1,sps)
           endif

           ! Store the normal vector, the porosity and the
           ! mesh velocity if present.

           sx = skAdj(0,0,kk,1,sps); sy = skAdj(0,0,kk,2,sps); sz = skAdj(0,0,kk,3,sps)
           por = porK(i,j,k)
           if( addGridVelocities ) sFace = sFaceKAdj(0,0,kk,sps)

           ! Compute the value of gamma on the face.
           ! Constant gamma for now.

           gammaFace = gammaConstant

           ! Compute the dissipative flux across the interface
           ! and them to dwAdj.

!           call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace)
!           print *,'riemanninputk',left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid
           call riemannFluxAdj(left, right, flux,por,gammaFace,correctForK,sX,sY,sZ,sFace,fineGrid)

           dwAdj(irho,sps)  = dwAdj(irho,sps)  + fact*flux(irho)
           dwAdj(imx,sps)   = dwAdj(imx,sps)   + fact*flux(imx)
           dwAdj(imy,sps)   = dwAdj(imy,sps)   + fact*flux(imy)
           dwAdj(imz,sps)   = dwAdj(imz,sps)   + fact*flux(imz)
           dwAdj(irhoE,sps) = dwAdj(irhoE,sps) + fact*flux(irhoE)

           !print *,'dwupwind',dwadj,'fact',fact,'flux',flux

           ! Update k and set fact to 1 for the second face.

           k    = k + 1
           fact = one

         enddo

       endif orderTest

       end subroutine inviscidUpwindFluxAdj
