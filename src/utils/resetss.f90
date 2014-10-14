!
!      ******************************************************************
!      *                                                                *
!      * File:          resetss.f90                                 *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetss(nn, ssi, ssj, ssk, ss)
       
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(:,:,:), pointer :: ssi, ssj, ssk
       real(kind=realType), dimension(:,:,:), pointer :: ss


!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! nullify all pointers
       nullify(ssi, ssj, ssk ,ss)


       end subroutine resetss
