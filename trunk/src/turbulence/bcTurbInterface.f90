!
!      ******************************************************************
!      *                                                                *
!      * File:          bcTurbInterface.f90                             *
!      * Author:        Georgi Kalitzin, Edwin van der Weide            *
!      * Starting date: 01-09-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine bcTurbInterface(nn)
!
!      ******************************************************************
!      *                                                                *
!      * bcTurbInterface applies the halo treatment for interface halo  *
!      * cells, sliding mesh interface and domain interface. As these   *
!      * are not really boundary conditions, the variable bvt is simply *
!      * set to keep the current value.                                 *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use BCTypes
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nn
!
!      Local variables.
!
       integer(kind=intType) :: i, j, l

       real(kind=realType), dimension(:,:,:), pointer :: bvt, ww1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointers for bvt and ww, depending on the block face
       ! on which the subface is located.

       select case (BCFaceID(nn))
         case (iMin)
           bvt => bvti1; ww1 => w( 1,1:,1:,:)
         case (iMax)
           bvt => bvti2; ww1 => w(il,1:,1:,:)
         case (jMin)
           bvt => bvtj1; ww1 => w(1:, 1,1:,:)
         case (jMax)
           bvt => bvtj2; ww1 => w(1:,jl,1:,:)
         case (kMin)
           bvt => bvtk1; ww1 => w(1:,1:, 1,:)
         case (kMax)
           bvt => bvtk2; ww1 => w(1:,1:,kl,:)
       end select

       ! Loop over the faces of the subfaces and set the values of
       ! bvt to keep the current value.

       do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
         do i=BCData(nn)%icBeg, BCData(nn)%icEnd
           do l=nt1,nt2
             bvt(i,j,l) = ww1(i,j,l)
           enddo
         enddo
       enddo

       end subroutine bcTurbInterface
