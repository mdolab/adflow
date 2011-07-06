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
       subroutine inviscidCentralFlux
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
       use blockPointers
       use cgnsGrid
       use constants
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ind

       real(kind=realType) :: qsp, qsm, rqsp, rqsm, porVel, porFlux
       real(kind=realType) :: pa, fs, sFace, vnp, vnm
       real(kind=realType) :: wx, wy, wz, rvol
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize sFace to zero. This value will be used if the
       ! block is not moving.

       sFace = zero
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the i-direction.                           *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=2,jl
           do i=1,il

             ! Set the dot product of the grid velocity and the
             ! normal in i-direction for a moving face.

             if( addGridVelocities ) sFace = sFaceI(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i+1,j,k,ivx)*sI(i,j,k,1) &
                 + w(i+1,j,k,ivy)*sI(i,j,k,2) &
                 + w(i+1,j,k,ivz)*sI(i,j,k,3)
             vnm = w(i,  j,k,ivx)*sI(i,j,k,1) &
                 + w(i,  j,k,ivy)*sI(i,j,k,2) &
                 + w(i,  j,k,ivz)*sI(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porI(i,j,k) == noFlux)    porFlux = zero
             if(porI(i,j,k) == boundFlux) then
               porVel = zero
               vnp    = sFace
               vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities relative to the grid for
             ! the face as well as the mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i+1,j,k,irho)
             rqsm = qsm*w(i,  j,k,irho)
 
             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i+1,j,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i+1,j,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i+1,j,k,irho) = dw(i+1,j,k,irho) - fs
             dw(i,  j,k,irho) = dw(i,  j,k,irho) + fs

             ind = indFamilyI(i,j,k)
             massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                                                 + factFamilyI(i,j,k)*fs

             fs = rqsp*w(i+1,j,k,ivx) + rqsm*w(i,j,k,ivx) &
                + pa*sI(i,j,k,1)
             dw(i+1,j,k,imx) = dw(i+1,j,k,imx) - fs
             dw(i,  j,k,imx) = dw(i,  j,k,imx) + fs

             fs = rqsp*w(i+1,j,k,ivy) + rqsm*w(i,j,k,ivy) &
                + pa*sI(i,j,k,2)
             dw(i+1,j,k,imy) = dw(i+1,j,k,imy) - fs
             dw(i,  j,k,imy) = dw(i,  j,k,imy) + fs

             fs = rqsp*w(i+1,j,k,ivz) + rqsm*w(i,j,k,ivz) &
                + pa*sI(i,j,k,3)
             dw(i+1,j,k,imz) = dw(i+1,j,k,imz) - fs
             dw(i,  j,k,imz) = dw(i,  j,k,imz) + fs

             fs = qsp*w(i+1,j,k,irhoE) + qsm*w(i,j,k,irhoE) &
                + porFlux*(vnp*p(i+1,j,k) + vnm*p(i,j,k))
             dw(i+1,j,k,irhoE) = dw(i+1,j,k,irhoE) - fs
             dw(i,  j,k,irhoE) = dw(i,  j,k,irhoE) + fs

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the j-direction.                           *
!      *                                                                *
!      ******************************************************************
!
       do k=2,kl
         do j=1,jl
           do i=2,il

             ! Set the dot product of the grid velocity and the
             ! normal in j-direction for a moving face.

             if( addGridVelocities ) sFace = sFaceJ(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i,j+1,k,ivx)*sJ(i,j,k,1) &
                 + w(i,j+1,k,ivy)*sJ(i,j,k,2) &
                 + w(i,j+1,k,ivz)*sJ(i,j,k,3)
             vnm = w(i,j,  k,ivx)*sJ(i,j,k,1) &
                 + w(i,j,  k,ivy)*sJ(i,j,k,2) &
                 + w(i,j,  k,ivz)*sJ(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half
             if(porJ(i,j,k) == noFlux)    porFlux = zero
             if(porJ(i,j,k) == boundFlux) then
               porVel = zero
               vnp    = sFace
               vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities for the face as well as the
             ! mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i,j+1,k,irho)
             rqsm = qsm*w(i,j,  k,irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i,j+1,k) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i,j+1,k. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i,j+1,k,irho) = dw(i,j+1,k,irho) - fs
             dw(i,j,  k,irho) = dw(i,j,  k,irho) + fs

             ind = indFamilyJ(i,j,k)
             massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                                                 + factFamilyJ(i,j,k)*fs

             fs = rqsp*w(i,j+1,k,ivx) + rqsm*w(i,j,k,ivx) &
                + pa*sJ(i,j,k,1)
             dw(i,j+1,k,imx) = dw(i,j+1,k,imx) - fs
             dw(i,j,  k,imx) = dw(i,j,  k,imx) + fs

             fs = rqsp*w(i,j+1,k,ivy) + rqsm*w(i,j,k,ivy) &
                + pa*sJ(i,j,k,2)
             dw(i,j+1,k,imy) = dw(i,j+1,k,imy) - fs
             dw(i,j,  k,imy) = dw(i,j,  k,imy) + fs

             fs = rqsp*w(i,j+1,k,ivz) + rqsm*w(i,j,k,ivz) &
                + pa*sJ(i,j,k,3)
             dw(i,j+1,k,imz) = dw(i,j+1,k,imz) - fs
             dw(i,j,  k,imz) = dw(i,j,  k,imz) + fs

             fs = qsp*w(i,j+1,k,irhoE) + qsm*w(i,j,k,irhoE) &
                + porFlux*(vnp*p(i,j+1,k) + vnm*p(i,j,k))
             dw(i,j+1,k,irhoE) = dw(i,j+1,k,irhoE) - fs
             dw(i,j,  k,irhoE) = dw(i,j,  k,irhoE) + fs

           enddo
         enddo
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Advective fluxes in the k-direction.                           *
!      *                                                                *
!      ******************************************************************
!
       do k=1,kl
         do j=2,jl
           do i=2,il

             ! Set the dot product of the grid velocity and the
             ! normal in k-direction for a moving face.

             if( addGridVelocities ) sFace = sFaceK(i,j,k)

             ! Compute the normal velocities of the left and right state.

             vnp = w(i,j,k+1,ivx)*sK(i,j,k,1) &
                 + w(i,j,k+1,ivy)*sK(i,j,k,2) &
                 + w(i,j,k+1,ivz)*sK(i,j,k,3)
             vnm = w(i,j,k,  ivx)*sK(i,j,k,1) &
                 + w(i,j,k,  ivy)*sK(i,j,k,2) &
                 + w(i,j,k,  ivz)*sK(i,j,k,3)

             ! Set the values of the porosities for this face.
             ! porVel defines the porosity w.r.t. velocity;
             ! porFlux defines the porosity w.r.t. the entire flux.
             ! The latter is only zero for a discontinuous block
             ! block boundary that must be treated conservatively.
             ! The default value of porFlux is 0.5, such that the
             ! correct central flux is scattered to both cells.
             ! In case of a boundFlux the normal velocity is set
             ! to sFace.

             porVel  = one
             porFlux = half

             if(porK(i,j,k) == noFlux)    porFlux = zero
             if(porK(i,j,k) == boundFlux) then
               porVel = zero
               vnp    = sFace
               vnm    = sFace
             endif

             ! Incorporate porFlux in porVel.

             porVel = porVel*porFlux

             ! Compute the normal velocities for the face as well as the
             ! mass fluxes.

             qsp = (vnp - sFace)*porVel
             qsm = (vnm - sFace)*porVel

             rqsp = qsp*w(i,j,k+1,irho)
             rqsm = qsm*w(i,j,k,  irho)

             ! Compute the sum of the pressure multiplied by porFlux.
             ! For the default value of porFlux, 0.5, this leads to
             ! the average pressure.

             pa = porFlux*(p(i,j,k+1) + p(i,j,k))

             ! Compute the fluxes and scatter them to the cells
             ! i,j,k and i,j,k+1. Store the density flux in the
             ! mass flow of the appropriate sliding mesh interface.

             fs = rqsp + rqsm
             dw(i,j,k+1,irho) = dw(i,j,k+1,irho) - fs
             dw(i,j,k,  irho) = dw(i,j,k,  irho) + fs

             ind = indFamilyK(i,j,k)
             massFlowFamilyInv(ind,spectralSol) =       &
                     massFlowFamilyInv(ind,spectralSol) &
                                                 + factFamilyK(i,j,k)*fs

             fs = rqsp*w(i,j,k+1,ivx) + rqsm*w(i,j,k,ivx) &
                + pa*sK(i,j,k,1)
             dw(i,j,k+1,imx) = dw(i,j,k+1,imx) - fs
             dw(i,j,k,  imx) = dw(i,j,k,  imx) + fs

             fs = rqsp*w(i,j,k+1,ivy) + rqsm*w(i,j,k,ivy) &
                + pa*sK(i,j,k,2)
             dw(i,j,k+1,imy) = dw(i,j,k+1,imy) - fs
             dw(i,j,k,  imy) = dw(i,j,k,  imy) + fs

             fs = rqsp*w(i,j,k+1,ivz) + rqsm*w(i,j,k,ivz) &
                + pa*sK(i,j,k,3)
             dw(i,j,k+1,imz) = dw(i,j,k+1,imz) - fs
             dw(i,j,k,  imz) = dw(i,j,k,  imz) + fs

             fs = qsp*w(i,j,k+1,irhoE) + qsm*w(i,j,k,irhoE) &
                + porFlux*(vnp*p(i,j,k+1) + vnm*p(i,j,k))
             dw(i,j,k+1,irhoE) = dw(i,j,k+1,irhoE) - fs
             dw(i,j,k,  irhoE) = dw(i,j,k,  irhoE) + fs

           enddo
         enddo
       enddo

       ! Add the rotational source terms for a moving block in a
       ! steady state computation. These source terms account for the
       ! centrifugal acceleration and the coriolis term. However, as
       ! the the equations are solved in the inertial frame and not
       ! in the moving frame, the form is different than what you
       ! normally find in a text book.

       rotation: if(blockIsMoving .and. equationMode == steady) then

         ! Compute the three nonDimensional angular velocities.

         wx = timeRef*cgnsDoms(nbkGlobal)%rotRate(1)
         wy = timeRef*cgnsDoms(nbkGlobal)%rotRate(2)
         wz = timeRef*cgnsDoms(nbkGlobal)%rotRate(3)

         ! Loop over the internal cells of this block to compute the
         ! rotational terms for the momentum equations.

         do k=2,kl
           do j=2,jl
             do i=2,il
               rvol = w(i,j,k,irho)*vol(i,j,k)

               dw(i,j,k,imx) = dw(i,j,k,imx) &
                             + rvol*(wy*w(i,j,k,ivz) - wz*w(i,j,k,ivy))
               dw(i,j,k,imy) = dw(i,j,k,imy) &
                             + rvol*(wz*w(i,j,k,ivx) - wx*w(i,j,k,ivz))
               dw(i,j,k,imz) = dw(i,j,k,imz) &
                             + rvol*(wx*w(i,j,k,ivy) - wy*w(i,j,k,ivx))
             enddo
           enddo
         enddo

       endif rotation

       end subroutine inviscidCentralFlux
