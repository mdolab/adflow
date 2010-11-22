!
!      ******************************************************************
!      *                                                                *
!      * File:          computeForcesPressureAdj.f90                    *
!      * Author:        Edwin van der Weide,C.A.(Sandy) Mader           *
!      * Starting date: 03-19-2006                                      *
!      * Last modified: 05-06-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computeForceCouplingPressureAdj(wAdj, pAdj)
!
!      ******************************************************************
!      *                                                                *
!      * Simple routine to compute the pressure from the variables w.   *
!      * A calorically perfect gas, i.e. constant gamma, is assumed.    *
!      *                                                                *
!      ******************************************************************
!
       use flowVarRefState
       use inputPhysics
       implicit none
!
!      Subroutine arguments
!
       real(kind=realType), dimension(2,2,2,nw), intent(in) :: wAdj
       real(kind=realType), dimension(2,2,2),intent(out) :: pAdj
!
!      Local variables
!
       integer(kind=intType) :: i, j, k

       real(kind=realType) :: gm1, factK, v2

!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       gm1 = gammaConstant - one

       ! Check the situation.

       if( kPresent ) then

         ! A separate equation for the turbulent kinetic energy is
         ! present. This variable must be taken into account.

         factK = five*third - gammaConstant

         do k=1,2
           do j=1,2
             do i=1,2
               v2 = wAdj(i,j,k,ivx)**2 + wAdj(i,j,k,ivy)**2 &
                  + wAdj(i,j,k,ivz)**2

               pAdj(i,j,k) = gm1*(wAdj(i,j,k,irhoE)    &
                           - half*wAdj(i,j,k,irho)*v2) &
                           + factK*wAdj(i,j,k,irho)*wAdj(i,j,k,itu1)
             enddo
           enddo
         enddo

       else

         ! No separate equation for the turbulent kinetic enery.
         ! Use the standard formula.

          do k=1,2
             do j=1,2
                do i=1,2

               v2 = wAdj(i,j,k,ivx)**2 + wAdj(i,j,k,ivy)**2 &
                  + wAdj(i,j,k,ivz)**2
               
               pAdj(i,j,k) = gm1*(wAdj(i,j,k,irhoE) &
                    - half*wAdj(i,j,k,irho)*v2)
             enddo
           enddo
         enddo

       endif

     end subroutine computeForceCouplingPressureAdj
