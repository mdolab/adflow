!
!      ******************************************************************
!      *                                                                *
!      * File:          computePressureAdj.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-19-2006                                      *
!      * Last modified: 03-20-2006                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine computePressureAdjTS(wAdj, pAdj,nn,level,sps)
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
       use inputTimeSpectral !nIntervalTimespectral
       implicit none
!
!      Subroutine arguments
!

       integer(kind=intType)::nn,level,sps
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nw,nTimeIntervalsSpectral), &
                                                      intent(in) :: wAdj
       real(kind=realType), dimension(-2:2,-2:2,-2:2,nTimeIntervalsSpectral) :: pAdj
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

         do k=-2,2
           do j=-2,2
             do i=-2,2
               v2 = wAdj(i,j,k,ivx,sps)**2 + wAdj(i,j,k,ivy,sps)**2 &
                  + wAdj(i,j,k,ivz,sps)**2

               pAdj(i,j,k,sps) = gm1*(wAdj(i,j,k,irhoE,sps)    &
                           - half*wAdj(i,j,k,irho,sps)*v2) &
                           + factK*wAdj(i,j,k,irho,sps)*wAdj(i,j,k,itu1,sps)
             enddo
           enddo
         enddo

       else

         ! No separate equation for the turbulent kinetic enery.
         ! Use the standard formula.

         do k=-2,2
           do j=-2,2
             do i=-2,2
               v2 = wAdj(i,j,k,ivx,sps)**2 + wAdj(i,j,k,ivy,sps)**2 &
                  + wAdj(i,j,k,ivz,sps)**2

               pAdj(i,j,k,sps) = gm1*(wAdj(i,j,k,irhoE,sps) &
                           - half*wAdj(i,j,k,irho,sps)*v2)
             enddo
           enddo
         enddo

       endif

     end subroutine computePressureAdjTS
