!
!      ******************************************************************
!      *                                                                *
!      * File:          resetssMetric.f90                               *
!      * Author:        Peter Zhoujie Lyu                               *
!      * Starting date: 11-03-2014                                      *
!      * Last modified: 11-03-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetssMetric(nn, ss)
       
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(:,:,:), pointer :: ss


!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! nullify all pointers
       !nullify(ssi, ssj, ssk ,ss)


     end subroutine resetssMetric
