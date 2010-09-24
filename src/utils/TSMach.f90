!
!      ******************************************************************
!      *                                                                *
!      * File:          TSMach.f90                                      *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 10-30-2009                                      *
!      * Last modified: 10-30-2009                                      *
!      *                                                                *
!      ******************************************************************
!
       function TSMach(degreePolMach,   coefPolMach,       &
                        degreeFourMach,  omegaFourMach,     &
                        cosCoefFourMach, sinCoefFourMach, t)
!
!      ******************************************************************
!      *                                                                *
!      * TSMach computes the Mach Number for a given time interval      *
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
       real(kind=realType) :: TSmach
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

       real(kind=realType) :: intervalMach, val
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Return immediately if this is a steady computation.

       if(equationMode == steady) then
         TSMach = zero
         return
       endif

       ! Compute the polynomial contribution. If no polynomial was
       ! specified, the value of index 0 is set to zero automatically.

       intervalMach = coefPolMach(0)
       do nn=1,degreePolMach
         intervalMach = intervalMach + coefPolMach(nn)*(t**nn)
       enddo

       ! Compute the fourier contribution. Again the cosine coefficient
       ! of index 0 is defaulted to zero if not specified.

       intervalMach = intervalMach + cosCoefFourMach(0)
       do nn=1,degreeFourMach
         val = nn*omegaFourMach*t
         intervalMach = intervalMach + cosCoefFourmach(nn)*cos(val) &
                   + sinCoefFourmach(nn)*sin(val)
       enddo
       print *,'inTSMach',intervalMach,nn,val,t
       ! Set TSMach to phi.

       TSMach = intervalMach

     end function TSmach
