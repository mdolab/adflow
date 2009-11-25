!
!      ******************************************************************
!      *                                                                *
!      * File:          TSbetadot.f90                                   *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 11-04-2009                                      *
!      * Last modified: 11-04-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSbetadot(degreePolBeta,   coefPolBeta,       &
                        degreeFourBeta,  omegaFourBeta,     &
                        cosCoefFourBeta, sinCoefFourBeta, t)
!
!      ******************************************************************
!      *                                                                *
!      * TSbeta computes the angle of attack for a given Time interval  *
!      * in a time spectral solution.                                   *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputPhysics
       implicit none
!
!      Function type
!
       real(kind=realType) :: TSbetadot
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

       real(kind=realType) :: betadot, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSBetadot = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       betadot = zero
       do nn=1,degreePolBeta
         betadot = betadot + nn*coefPolBeta(nn)*(t**(nn-1))
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       do nn=1,degreeFourBeta
         val = nn*omegaFourBeta
         betadot = betadot -val* cosCoefFourbeta(nn)*sin(val*t) &
                   +val* sinCoefFourbeta(nn)*cos(val*t)
       enddo

       ! Set TSBeta to phi.

       TSBetadot = betadot

     end function TSbetadot
