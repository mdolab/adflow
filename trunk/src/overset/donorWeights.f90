!
!      ******************************************************************
!      *                                                                *
!      * File:          donorWeights.f90                                *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 05-08-2005                                      *
!      * Last modified: 08-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine getWeights(s, weights)
!
!      ******************************************************************
!      *                                                                *
!      * getWeights converts the parametric coordinates into the actual *
!      * stencil weights used by the solver to interpolate.             *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
!
!      Subroutine arguments
!
       real(kind=realType), dimension(3), intent(in)  :: s
       real(kind=realType), dimension(8), intent(out) :: weights
!
!      Local variables.
!
       real(kind=realType) :: oms1, oms2, oms3
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       oms1 = one - s(1)
       oms2 = one - s(2)
       oms3 = one - s(3)

       weights(1) = oms1*oms2*oms3
       weights(2) = s(1)*oms2*oms3
       weights(3) = oms1*s(2)*oms3
       weights(4) = s(1)*s(2)*oms3
       weights(5) = oms1*oms2*s(3)
       weights(6) = s(1)*oms2*s(3)
       weights(7) = oms1*s(2)*s(3)
       weights(8) = s(1)*s(2)*s(3)

       end subroutine getWeights
