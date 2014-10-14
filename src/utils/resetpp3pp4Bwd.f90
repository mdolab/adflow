       !
!      ******************************************************************
!      *                                                                *
!      * File:          resetpp3pp4Bwd.f90                              *
!      * Author:        Eirikur Jonsson                                 *
!      * Starting date: 10-14-2014                                      *
!      * Last modified: 10-14-2014                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resetpp3pp4Bwd(nn, pp3, pp4)
       
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
       print *, "Testing resetpp3pp4bwd"

       select case (BCFaceID(nn))
         case (iMin)
           p(3,1:,1:) = pp3
           p(4,1:,1:) = pp4
           print *, "Testing resetpp3pp4bwd"
         case (iMax)
           p(nx,1:,1:) = pp3
           p(nx-1,1:,1:) = pp4
         case (jMin)
           p(1:,3,1:) = pp3
           p(1:,4,1:) = pp4
         case (jMax)
           p(1:,ny,1:) = pp3
           p(1:,ny-1,1:) = pp4
         case (kMin)
           p(1:,1:,3) = pp3
           p(1:,1:,4) = pp4
         case (kMax)
           p(1:,1:,nz) = pp3
           p(1:,1:,nz-1) = pp4
       end select

       end subroutine resetpp3pp4Bwd


