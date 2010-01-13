!
!      ******************************************************************
!      *                                                                *
!      * File:          TSalpha.f90                                     *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 10-30-2009                                      *
!      * Last modified: 10-30-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSalpha(degreePolAlpha,   coefPolAlpha,       &
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
       real(kind=realType) :: TSalpha
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

       real(kind=realType) :: alpha, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSAlpha = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       alpha = coefPolAlpha(0)
       do nn=1,degreePolAlpha
         alpha = alpha + coefPolAlpha(nn)*(t**nn)
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       alpha = alpha + cosCoefFourAlpha(0)
       do nn=1,degreeFourAlpha
         val = nn*omegaFourAlpha*t
         alpha = alpha + cosCoefFouralpha(nn)*cos(val) &
                   + sinCoefFouralpha(nn)*sin(val)
       enddo
       !print *,'inTSalpha',alpha,nn,val,t
       ! Set TSAlpha to phi.

       TSAlpha = alpha

     end function TSalpha
