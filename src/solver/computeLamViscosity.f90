!
!      ******************************************************************
!      *                                                                *
!      * File:          computeLamViscosity.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-10-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeLamViscosity
!
!      ******************************************************************
!      *                                                                *
!      * computeLamViscosity computes the laminar viscosity ratio in    *
!      * the owned cell centers of the given block. Sutherland's law is *
!      * used. It is assumed that the pointes already point to the      *
!      * correct block before entering this subroutine.                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use constants
       use flowVarRefState
       use inputPhysics
       use iteration
       implicit none
!
!      Local parameter.
!
       real(kind=realType), parameter :: twoThird = two*third
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       real(kind=realType)   :: muSuth, TSuth, SSuth, T
       logical               :: correctForK
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if no laminar viscosity needs to be computed.

       if(.not. viscous ) return

       ! Determine whether or not the pressure must be corrected
       ! for the presence of the turbulent kinetic energy.

       if( kPresent ) then
         if((currentLevel <= groundLevel) .or. turbCoupled) then
           correctForK = .true.
         else
           correctForK = .false.
         endif
       else
         correctForK = .false.
       endif

       ! Compute the nonDimensional constants in sutherland's law.

       muSuth = muSuthDim/muRef
       TSuth  = TSuthDim/Tref
       SSuth  = SSuthDim/Tref

       ! Substract 2/3 rho k, which is a part of the normal turbulent
       ! stresses, in case the pressure must be corrected.

       if( correctForK ) then
         do k=2,kl
           do j=2,jl
             do i=2,il
               p(i,j,k) = p(i,j,k) - twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
             enddo
           enddo
         enddo
       endif

       ! Loop over the owned cells of this block and compute the
       ! laminar viscosity ratio.

       do k=2,kl
         do j=2,jl
           do i=2,il

             ! Compute the nonDimensional temperature and the
             ! nonDimensional laminar viscosity.

             T = p(i,j,k)/(RGas*w(i,j,k,irho))
             rlv(i,j,k) = muSuth*((TSuth + SSuth)/(T + SSuth)) &
                        * ((T/TSuth)**1.5_realType)

           enddo
         enddo
       enddo

       ! Add the 2/3 rho k again to the pressure if the pressure was
       ! corrected earlier.

       if( correctForK ) then
         do k=2,kl
           do j=2,jl
             do i=2,il
               p(i,j,k) = p(i,j,k) + twoThird*w(i,j,k,irho)*w(i,j,k,itu1)
             enddo
           enddo
         enddo
       endif

       end subroutine computeLamViscosity
