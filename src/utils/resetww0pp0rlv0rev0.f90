!
!      ******************************************************************
!      *                                                                *
!      * File:          resetww0pp0rlv0rev0.f90                           *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetww0pp0rlv0rev0(nn, idim, ddim, ww0, pp0, rlv0, rev0)
       
       use BCTypes
       use blockPointers
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       integer(kind=intType) :: idim, ddim       

       real(kind=realType), dimension(:,:,:), pointer :: ww0
       real(kind=realType), dimension(:,:),   pointer :: pp0
       real(kind=realType), dimension(:,:),   pointer :: rlv0
       real(kind=realType), dimension(:,:),   pointer :: rev0

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and set
       ! the pointers accordinly.
       nullify(ww0, pp0, rlv0, rev0)

       end subroutine resetww0pp0rlv0rev0
