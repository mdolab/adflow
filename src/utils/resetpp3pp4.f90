!
!      ******************************************************************
!      *                                                                *
!      * File:          resetpp3pp4.f90                                 *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetpp3pp4(nn, pp3, pp4)
       
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(:,:),   pointer :: pp3, pp4

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! nullify all pointers
       !nullify(pp3, pp4)


       end subroutine resetpp3pp4
