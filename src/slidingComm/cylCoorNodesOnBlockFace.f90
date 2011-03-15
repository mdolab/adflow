!
!      ******************************************************************
!      *                                                                *
!      * File:          cylCoorNodesOnBlockFace.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-02-2005                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine cylCoorNodesOnBlockFace(xface, iBeg, iEnd,  &
                                          jBeg,  jEnd, color, faceID)
!
!      ******************************************************************
!      *                                                                *
!      * cylCoorNodesOnBlockFace determines the local cylindrical       *
!      * coordinates of the nodes given by the subface (jBeg,jEnd;      *
!      * iBeg,iEnd) on the block face faceid. It is assumed that the    *
!      * pointers in blockPointers, especially the coordinates,         *
!      * already point to the correct block.                            *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use interfaceGroups
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd
       integer(kind=intType), intent(in) :: color, faceID

       real(kind=realType), dimension(iBeg:iEnd,jBeg:jEnd,3), &
                                                  intent(out) :: xface
!
!      Local variables.
!
       integer(kind=intType) :: i, j, ii, jj

       real(kind=realType) :: xx, yy, zz, ax, r1, r2, r

       real(kind=realType), dimension(3) :: rotCenter
       real(kind=realType), dimension(3) :: axis, vec1, vec2

       real(kind=realType), dimension(:,:,:), pointer :: xf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the vectors which define the local cartesian
       ! coordinate system.

       rotCenter = myInterfaces(color)%rotCenter
       axis      = myInterfaces(color)%rotAxis
       vec1      = myInterfaces(color)%radVec1
       vec2      = myInterfaces(color)%radVec2

       ! Set the pointer for xf, depending on the block face to be
       ! looked at. As xf points to a subarray of x, whose indices start
       ! at 0, you must be careful later on when you access the data.
       ! An offset of +1 must be used, because per definition the
       ! subarray starts at 1.

       select case(faceID)
         case (iMin)
           xf => x(1 ,:,:,:)
         case (iMax)
           xf => x(il,:,:,:)
         case (jMin)
           xf => x(:,1, :,:)
         case (jMax)
           xf => x(:,jl,:,:)
         case (kMin)
           xf => x(:,:,1 ,:)
         case (kMax)
           xf => x(:,:,kl,:)
       end select

       ! Loop over the nodes of the subface.

       do j=jBeg,jEnd
         do i=iBeg,iEnd

           ! Determine the cartesian coordinates relative to the center
           ! of rotation. As xf is a pointer to a slice in x, the lower
           ! indices of xf are automatically set to 1. This in contrast
           ! to the original array x, which start at 0.
           ! Therefore i+1 and j+1 must be used.

           xx = xf(i+1,j+1,1) - rotCenter(1)
           yy = xf(i+1,j+1,2) - rotCenter(2)
           zz = xf(i+1,j+1,3) - rotCenter(3)

           ! Determine the axial and two radial components for this point.

           ax = xx*axis(1) + yy*axis(2) + zz*axis(3)
           r1 = xx*vec1(1) + yy*vec1(2) + zz*vec1(3)
           r2 = xx*vec2(1) + yy*vec2(2) + zz*vec2(3)

           ! Compute the radius.

           r = sqrt(r1*r1 + r2*r2)

           ! Compute the polar angle. Be careful when the node is a
           ! polar singularity. In that case the polar angle of the
           ! most appropriate neighbor.

           if(r < eps) then

             ! Polar singularity. Try the i neighbor. Which is either
             ! i-1 or i+1. Again the offset because of the pointer stuff.

             ii = i
             if(ii == iBeg) ii = i + 2

             xx = xf(ii,j+1,1) - rotCenter(1)
             yy = xf(ii,j+1,2) - rotCenter(2)
             zz = xf(ii,j+1,3) - rotCenter(3)

             r1 = xx*vec1(1) + yy*vec1(2) + zz*vec1(3)
             r2 = xx*vec2(1) + yy*vec2(2) + zz*vec2(3)

             ! If r1 and r2 still define a polar singularity check
             ! the j-direction.

             if((abs(r1) < eps) .and. (abs(r2) < eps)) then

               ! Try the j-neighbor.

               jj = j
               if(jj == jBeg) jj = j + 2

               xx = xf(i+1,jj,1) - rotCenter(1)
               yy = xf(i+1,jj,2) - rotCenter(2)
               zz = xf(i+1,jj,3) - rotCenter(3)

               r1 = xx*vec1(1) + yy*vec1(2) + zz*vec1(3)
               r2 = xx*vec2(1) + yy*vec2(2) + zz*vec2(3)

             endif
           endif

           ! Store the polar coordinates of the point.

           xface(i,j,1) = ax
           xface(i,j,2) = r
           xface(i,j,3) = atan2(r2,r1)

         enddo
       enddo

       end subroutine cylCoorNodesOnBlockFace
