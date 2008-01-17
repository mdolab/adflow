!
!      ******************************************************************
!      *                                                                *
!      * File:          determineRotationMatrices.f90                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 11-11-2003                                      *
!      * Last modified: 03-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineRotationMatrices
!
!      ******************************************************************
!      *                                                                *
!      * determineRotationMatrices determines the rotation matrices     *
!      * for the velocities for the periodic transformations of the     *
!      * sliding mesh interfaces.                                       *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use communication
       use commSliding
       use inputPhysics
       use interfaceGroups
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, nn

       real(kind=realType) :: perAngle, rotAngle, cosAngle, sinAngle

       real(kind=realType), dimension(3) :: ax, r1, r2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! The rotation matrices are only needed for unsteady
       ! computations. Threrefore return directly if a steady
       ! computation is performed.

       if(equationMode == steady) return

       ! Determine the number of rotation matrices stored on this
       ! processor. Note that the unity transformation is not stored,
       ! which explains nSlices - 1. As we only consider unsteady
       ! computations here the number of slices on both sides are
       ! identical (this is tested in preprocesssing) and therefore
       ! nSlices1 is taken.

       nRotSliding = 0
       do ii=1,nInterfaceGroups
         nRotSliding = nRotSliding + myInterfaces(ii)%nSlices1 - 1
       enddo

       ! Allocate the memory for the rotation matrices.

       allocate(rotSliding(nRotSliding,3,3), stat=ierr)
       if(ierr /= 0)                                 &
         call terminate("determineRotationMatrices", &
                        "Memory allocation failure for rotSliding.")

       ! Loop over the number of nInterfaceGroups to store the rotation
       ! matrices of the corresponding transformations.

       nn = 0
       nGroups: do ii=1,nInterfaceGroups

         ! Check if rotation matrices must be computed for this
         ! interface.

         testMatrices: if(myInterfaces(ii)%nSlices1 > 1) then

           ! Compute the periodic angle for part 1 and store the three
           ! vectors, which define the axial and the two radial
           ! directions, a bit easier.

           perAngle = two*pi/myInterfaces(ii)%nSlices1
           ax = myInterfaces(ii)%rotAxis
           r1 = myInterfaces(ii)%radVec1
           r2 = myInterfaces(ii)%radVec2

           ! Loop over the possible number of rotation matrices.
           ! Update the counter nn and compute the rotation matrix.

           do jj=1,(myInterfaces(ii)%nSlices1-1)
             nn = nn + 1
             rotAngle = jj*perAngle

             ! Compute the cosine and sine of the rotation angle and
             ! apply the matrix multiplication to obtain the rotation
             ! matrix.

             cosAngle = cos(rotAngle)
             sinAngle = sin(rotAngle)

             rotSliding(nn,1,1) = ax(1)*ax(1)                          &
                                + cosAngle*(r1(1)*r1(1) + r2(1)*r2(1))
             rotSliding(nn,1,2) = ax(1)*ax(2)                          &
                                + cosAngle*(r1(1)*r1(2) + r2(1)*r2(2)) &
                                + sinAngle*(r1(2)*r2(1) - r1(1)*r2(2))
             rotSliding(nn,1,3) = ax(1)*ax(3)                          &
                                + cosAngle*(r1(1)*r1(3) + r2(1)*r2(3)) &
                                + sinAngle*(r1(3)*r2(1) - r1(1)*r2(3))

             rotSliding(nn,2,1) = ax(1)*ax(2)                          &
                                + cosAngle*(r1(1)*r1(2) + r2(1)*r2(2)) &
                                - sinAngle*(r1(2)*r2(1) - r1(1)*r2(2))
             rotSliding(nn,2,2) = ax(2)*ax(2)                          &
                                + cosAngle*(r1(2)*r1(2) + r2(2)*r2(2))
             rotSliding(nn,2,3) = ax(2)*ax(3)                          &
                                + cosAngle*(r1(2)*r1(3) + r2(2)*r2(3)) &
                                + sinAngle*(r1(3)*r2(2) - r1(2)*r2(3))

             rotSliding(nn,3,1) = ax(1)*ax(3)                          &
                                + cosAngle*(r1(1)*r1(3) + r2(1)*r2(3)) &
                                - sinAngle*(r1(3)*r2(1) - r1(1)*r2(3))
             rotSliding(nn,3,2) = ax(2)*ax(3)                          &
                                + cosAngle*(r1(2)*r1(3) + r2(2)*r2(3)) &
                                - sinAngle*(r1(3)*r2(2) - r1(2)*r2(3))
             rotSliding(nn,3,3) = ax(3)*ax(3)                          &
                                + cosAngle*(r1(3)*r1(3) + r2(3)*r2(3))
           enddo

         endif testMatrices

       enddo nGroups

       end subroutine determineRotationMatrices
