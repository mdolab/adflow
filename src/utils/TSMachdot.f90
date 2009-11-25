!
!      ******************************************************************
!      *                                                                *
!      * File:          TSMachdot.f90                                     *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 11-04-2009                                      *
!      * Last modified: 11-04-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSMachdot(degreePolMach,   coefPolMach,       &
                        degreeFourMach,  omegaFourMach,     &
                        cosCoefFourMach, sinCoefFourMach, t)
!
!      ******************************************************************
!      *                                                                *
!      * TSmach computes the angle of attack for a given Time interval *
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
       real(kind=realType) :: TSmachdot
!
!      Function arguments.
!
       integer(kind=intType), intent(in) :: degreePolMach
       integer(kind=intType), intent(in) :: degreeFourMach

       real(kind=realType), intent(in) :: omegaFourMach, t

       real(kind=realType), dimension(0:*), intent(in) :: coefPolMach
       real(kind=realType), dimension(0:*), intent(in) :: cosCoefFourMach
       real(kind=realType), dimension(*),   intent(in) :: sinCoefFourMach
!
!      Local variables.
!
       integer(kind=intType) :: nn

       real(kind=realType) :: machdot, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSMachdot = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       machdot = zero
       do nn=1,degreePolMach
         machdot = machdot + nn*coefPolMach(nn)*(t**(nn-1))
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       do nn=1,degreeFourMach
         val = nn*omegaFourMach
         machdot = machdot -val* cosCoefFourmach(nn)*sin(val*t) &
                   +val* sinCoefFourmach(nn)*cos(val*t)
       enddo

       ! Set TSMach to phi.

       TSMachdot = machdot

     end function TSmachdot
