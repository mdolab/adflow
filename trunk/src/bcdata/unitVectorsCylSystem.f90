!
!      ******************************************************************
!      *                                                                *
!      * File:          unitVectorsCylSystem.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-25-2004                                      *
!      * Last modified: 09-27-2004                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine unitVectorsCylSystem(boco)
!
!      ******************************************************************
!      *                                                                *
!      * unitVectorsCylSystem determines the unit vectors of the        *
!      * local coordinate systen of the boundary face defined by the    *
!      * data in BCDataMod. In that local system the axial direction    *
!      * is rotation axis.                                              *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use section
       use BCDataMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: boco
!
!      Local variables.
!
       integer(kind=intType) :: i, j
       real(kind=realType)   :: factInlet, var

       real(kind=realType), dimension(3) :: dir

       real(kind=realType), dimension(:,:,:), pointer :: ss
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the pointers for coordinates and normals of the block
       ! face on which this subface is located. Set factInlet
       ! such that factInlet*normals points into the domain.

       select case (BCFaceID(boco))
         case (iMin)
           xf => x(1,:,:,:);  ss => si(1 ,:,:,:); factInlet =  one
         case (iMax)
           xf => x(il,:,:,:); ss => si(il,:,:,:); factInlet = -one
         case (jMin)
           xf => x(:,1,:,:);  ss => sj(:,1 ,:,:); factInlet =  one
         case (jMax)
           xf => x(:,jl,:,:); ss => sj(:,jl,:,:); factInlet = -one
         case (kMin)
           xf => x(:,:,1,:);  ss => sk(:,:,1 ,:); factInlet =  one
         case (kMax)
           xf => x(:,:,kl,:); ss => sk(:,:,kl,:); factInlet = -one
       end select

       ! Loop over the physical range of the subface to store the sum of
       ! the normals. Note that jBeg, jEnd, iBeg, iEnd cannot be used
       ! here, because they may include the halo faces. Instead the
       ! nodal range is used, which defines the original subface. The
       ! offset of +1 in the start index is there because you need
       ! the face id's.

       dir(1) = zero; dir(2) = zero; dir(3) = zero

       do j=(BCData(boco)%jnBeg+1), BCData(boco)%jnEnd
         do i=(BCData(boco)%inBeg+1), BCData(boco)%inEnd
           dir(1) = dir(1) + ss(i,j,1)
           dir(2) = dir(2) + ss(i,j,2)
           dir(3) = dir(3) + ss(i,j,3)
         enddo
       enddo

       ! Multiply by factInlet to make sure that the normal
       ! is inward pointing.

       dir(1) = dir(1)*factInlet
       dir(2) = dir(2)*factInlet
       dir(3) = dir(3)*factInlet

       ! Determine three unit vectors, which define the local cartesian
       ! coordinate system of the rotation axis. First the axial
       ! direction. If the axis cannot be determined from rotation info,
       ! it is assumed to be the x-axis.

       axis = sections(sectionId)%rotAxis
       var  = axis(1)**2 + axis(2)**2 + axis(3)**2
       if(var < half) then

         ! No rotation axis specified. Assume the x-axis
         ! and set the logical axAssumed to .True.

         axis(1) = one; axis(2) = zero; axis(3) = zero
         axAssumed = .true.
       endif

       ! The axial axis must be such that it points into the
       ! computational domain. If the dot product with dir is
       ! negative the direction of axis should be reversed.

       var = axis(1)*dir(1) + axis(2)*dir(2) + axis(3)*dir(3)
       if(var < zero) then
         axis(1) = -axis(1); axis(2) = -axis(2); axis(3) = -axis(3)
       endif

       ! Two unit vectors define the radial plane. These vectors are
       ! defined up to a constants. Just pick a direction for the second
       ! and create a unit vector normal to axis.

       if(abs(axis(2)) < 0.707107_realType) then
         radVec1(1) = zero; radVec1(2) = one;  radVec1(3) = zero
       else
         radVec1(1) = zero; radVec1(2) = zero; radVec1(3) = one
       endif

       var = radVec1(1)*axis(1) + radVec1(2)*axis(2) &
           + radVec1(3)*axis(3)
       radVec1(1) = radVec1(1) - var*axis(1)
       radVec1(2) = radVec1(2) - var*axis(2)
       radVec1(3) = radVec1(3) - var*axis(3)

       var = one/sqrt(radVec1(1)**2 + radVec1(2)**2 &
           +          radVec1(3)**2)
       radVec1(1) = radVec1(1)*var
       radVec1(2) = radVec1(2)*var
       radVec1(3) = radVec1(3)*var

       ! The second vector of the radial plane is obtained
       ! by taking the cross product of axis and radVec1.

       radVec2(1) = axis(2)*radVec1(3) - axis(3)*radVec1(2)
       radVec2(2) = axis(3)*radVec1(1) - axis(1)*radVec1(3)
       radVec2(3) = axis(1)*radVec1(2) - axis(2)*radVec1(1)

       end subroutine unitVectorsCylSystem
