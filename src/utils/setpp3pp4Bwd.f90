!
!      ******************************************************************
!      *                                                                *
!      * File:          setpp3pp4Bwd.f90                                *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setpp3pp4Bwd(nn, pp3, pp4)
       
       use BCTypes
       use blockPointers
!       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
       real(kind=realType), dimension(imaxDim,jmaxDim) :: pp3, pp4

!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and set
       ! the pointers accordinly.

       select case (BCFaceID(nn))
         case (iMin)
           pp3 = p(3,1:,1:)
           pp4 = p(4,1:,1:)
         case (iMax)
           pp3 = p(nx,1:,1:)
           pp4 = p(nx-1,1:,1:)
         case (jMin)
           pp3 = p(1:,3,1:)
           pp4 = p(1:,4,1:)
         case (jMax)
           pp3 = p(1:,ny,1:)
           pp4 = p(1:,ny-1,1:)
         case (kMin)
           pp3 = p(1:,1:,3)
           pp4 = p(1:,1:,4)
         case (kMax)
           pp3 = p(1:,1:,nz)
           pp4 = p(1:,1:,nz-1)
       end select

       end subroutine setpp3pp4Bwd   
       
