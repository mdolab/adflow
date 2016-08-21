!
!      ******************************************************************
!      *                                                                *
!      * File:          TSBeta.f90                                     *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 11-03-2009                                      *
!      * Last modified: 11-03-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSbeta(degreePolBeta,   coefPolBeta,       &
                        degreeFourBeta,  omegaFourBeta,     &
                        cosCoefFourBeta, sinCoefFourBeta, t)
!
!      ******************************************************************
!      *                                                                *
!      * TSbeta computes the angle of attack for a given Time interval *
!      * in a time spectral solution.                                   *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputPhysics, only : equationMode
       implicit none
!
!      Function type
!
       real(kind=realType) :: TSbeta
!
!      Function arguments.
!
       integer(kind=intType), intent(in) :: degreePolBeta
       integer(kind=intType), intent(in) :: degreeFourBeta

       real(kind=realType), intent(in) :: omegaFourBeta, t

       real(kind=realType), dimension(0:*), intent(in) :: coefPolBeta
       real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourBeta
       real(kind=realType), dimension(*),   intent(in) :: sinCoefFourBeta
!
!      Local variables.
!
       integer(kind=intType) :: nn

       real(kind=realType) :: beta, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSBeta = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       beta = coefPolBeta(0)
       do nn=1,degreePolBeta
         beta = beta + coefPolBeta(nn)*(t**nn)
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       beta = beta + cosCoefFourBeta(0)
       do nn=1,degreeFourBeta
         val = nn*omegaFourBeta*t
         beta = beta + cosCoefFourbeta(nn)*cos(val) &
                   + sinCoefFourbeta(nn)*sin(val)
       enddo

       ! Set TSBeta to phi.

       TSBeta = beta

     end function TSbeta
