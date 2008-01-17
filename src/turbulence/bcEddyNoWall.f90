!
!      ******************************************************************
!      *                                                                *
!      * File:          bcEddyNoWall.f90                                *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 06-11-2003                                      *
!      * Last modified: 04-11-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcEddyNoWall(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcEddyNoWall sets the eddy viscosity in the halo cells of      *
!      * subface nn of the block given in blockPointers. The boundary   *
!      * condition on the subface can be anything but a viscous wall.   *
!      * A homogeneous neumann condition is applied, which means that   *
!      * the eddy viscosity is simply copied from the interior cell.    *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      Local variables.
!
       integer(kind=intType) :: i, j

       real(kind=realType), dimension(:,:), pointer :: rev1, rev2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the face id on which the subface is located and
       ! set the pointers rev1 and rev2 accordingly.

       select case (BCFaceID(nn))
         case (iMin)
           rev1 => rev(1, 1:,1:); rev2 => rev(2, 1:,1:)
         case (iMax)
           rev1 => rev(ie,1:,1:); rev2 => rev(il,1:,1:)
         case (jMin)
           rev1 => rev(1:,1 ,1:); rev2 => rev(1:,2 ,1:)
         case (jMax)
           rev1 => rev(1:,je,1:); rev2 => rev(1:,jl,1:)
         case (kMin)
           rev1 => rev(1:,1:,1 ); rev2 => rev(1:,1:,2 )
         case (kMax)
           rev1 => rev(1:,1:,ke); rev2 => rev(1:,1:,kl)
       end select

       ! Loop over the faces of the subface and set the eddy
       ! viscosity in the halo cells.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           rev1(i,j) = rev2(i,j)
         enddo
       enddo

       end subroutine bcEddyNoWall
