!
!      ******************************************************************
!      *                                                                *
!      * File:          TSalphadot.f90                                     *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 11-03-2009                                      *
!      * Last modified: 11-03-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSalphadot(degreePolAlpha,   coefPolAlpha,       &
                        degreeFourAlpha,  omegaFourAlpha,     &
                        cosCoefFourAlpha, sinCoefFourAlpha, t)
!
!      ******************************************************************
!      *                                                                *
!      * TSalpha computes the angle of attack for a given Time interval *
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
       real(kind=realType) :: TSalphadot
!
!      Function arguments.
!
       integer(kind=intType), intent(in) :: degreePolAlpha
       integer(kind=intType), intent(in) :: degreeFourAlpha

       real(kind=realType), intent(in) :: omegaFourAlpha, t

       real(kind=realType), dimension(0:*), intent(in) :: coefPolAlpha
       real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourAlpha
       real(kind=realType), dimension(*),   intent(in) :: sinCoefFourAlpha
!
!      Local variables.
!
       integer(kind=intType) :: nn

       real(kind=realType) :: alphadot, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSAlphadot = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       alphadot = zero
       do nn=1,degreePolAlpha
         alphadot = alphadot + nn*coefPolAlpha(nn)*(t**(nn-1))
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       do nn=1,degreeFourAlpha
         val = nn*omegaFourAlpha
         alphadot = alphadot -val* cosCoefFouralpha(nn)*sin(val*t) &
                   +val* sinCoefFouralpha(nn)*cos(val*t)
       enddo

       ! Set TSAlpha to phi.

       TSAlphadot = alphadot

     end function TSalphadot
