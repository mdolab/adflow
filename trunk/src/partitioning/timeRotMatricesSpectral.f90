!
!      ******************************************************************
!      *                                                                *
!      * File:          timeRotMatricesSpectral.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 07-23-2004                                      *
!      * Last modified: 06-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine timeRotMatricesSpectral
!
!      ******************************************************************
!      *                                                                *
!      * timeRotMatricesSpectral determines the rotation matrices       *
!      * used in the time derivatives for the velocity components in    *
!      * the time spectral method. These matrices are the identity      *
!      * matrices for non-rotating sections and something different for *
!      * rotating sections. Therefore the rotation matrices are stored  *
!      * for every section.                                             *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use inputPhysics
       use inputTimeSpectral
       use section
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
       real(kind=realType)   :: tmp, theta, cosTheta, sinTheta

       real(kind=realType), dimension(3) :: xt, yt, zt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! This routine is only used for the spectral solutions. Return
       ! immediately if a different mode is solved.

       if(equationMode /= timeSpectral) return

       ! Allocate the memory for rotMatrixSpectral, which will store
       ! the rotation matrices for all the sections.

       allocate(rotMatrixSpectral(nSections,3,3), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("timeRotMatricesSpectral", &
                        "Memory allocation failure for &
                        &rotMatrixSpectral")

       ! Loop over the number of sections.

       sectionLoop: do nn=1,nSections

         ! Test if the rotation matrix is the unity matrix. This is the
         ! case if this section is not rotating or if the number of
         ! slices is only 1, i.e. the true physical model is computed.

         testUnity: if(.not. sections(nn)%rotating .or. &
                        sections(nn)%nSlices == 1) then

           ! Set the rotation matrix to the unity matrix.

           rotMatrixSpectral(nn,1,1) = one
           rotMatrixSpectral(nn,1,2) = zero
           rotMatrixSpectral(nn,1,3) = zero

           rotMatrixSpectral(nn,2,1) = zero
           rotMatrixSpectral(nn,2,2) = one
           rotMatrixSpectral(nn,2,3) = zero

           rotMatrixSpectral(nn,3,1) = zero
           rotMatrixSpectral(nn,3,2) = zero
           rotMatrixSpectral(nn,3,3) = one

         else testUnity

           ! Section is rotating and only a part of the physical problem
           ! is modelled. Consequently a rotation matrix is present for
           ! the velocity components.

           ! First transform to a frame where the xt-axis points in the
           ! direction of the rotation vector.

           xt(1) = sections(nn)%rotAxis(1)
           xt(2) = sections(nn)%rotAxis(2)
           xt(3) = sections(nn)%rotAxis(3)

           ! Construct the yt axis. It does not matter exactly as long
           ! as it is normal to xt.

           if(abs(xt(2)) < 0.707107_realType) then
             yt(1) = zero
             yt(2) = one
             yt(3) = zero
           else
             yt(1) = zero
             yt(2) = zero
             yt(3) = one
           endif

           ! Make sure that yt is normal to xt.

           tmp   = xt(1)*yt(1) + xt(2)*yt(2) + xt(3)*yt(3)
           yt(1) = yt(1) - tmp*xt(1)
           yt(2) = yt(2) - tmp*xt(2)
           yt(3) = yt(3) - tmp*xt(3)

           ! And create a unit vector.

           tmp   = one/sqrt(yt(1)**2 + yt(2)**2 + yt(3)**2)
           yt(1) = tmp*yt(1)
           yt(2) = tmp*yt(2)
           yt(3) = tmp*yt(3)

           ! Create the vector zt by taking the cross product xt*yt.

           zt(1) = xt(2)*yt(3) - xt(3)*yt(2)
           zt(2) = xt(3)*yt(1) - xt(1)*yt(3)
           zt(3) = xt(1)*yt(2) - xt(2)*yt(1)

           ! Compute the periodic angle theta and its sine and cosine.

           theta     = two*pi/real(sections(nn)%nSlices,realType)
           cosTheta = cos(theta)
           sinTheta = sin(theta)

           ! The rotation matrix in the xt,yt,zt frame is given by
           !
           ! R = | 1      0           0      |
           !     | 0  cos(theta) -sin(theta) |
           !     | 0  sin(theta)  cos(theta) |
           !
           ! The rotation matrix in the standard cartesian frame is then
           ! given by t * r * t^t, where the colums of the transformation
           ! matrix t are the unit vectors xt,yt,zt. One can easily check
           ! this by checking rotation around the y- and z-axis. The
           ! result of this is the expression below.

           rotMatrixSpectral(nn,1,1) = xt(1)*xt(1)                    &
                               + cosTheta*(yt(1)*yt(1) + zt(1)*zt(1))
           rotMatrixSpectral(nn,1,2) = xt(1)*xt(2)                    &
                               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
                               - sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
           rotMatrixSpectral(nn,1,3) = xt(1)*xt(3)                    &
                               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
                               - sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))

           rotMatrixSpectral(nn,2,1) = xt(1)*xt(2)                    &
                               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
                               + sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
           rotMatrixSpectral(nn,2,2) = xt(2)*xt(2)                    &
                               + cosTheta*(yt(2)*yt(2) + zt(2)*zt(2))
           rotMatrixSpectral(nn,2,3) = xt(2)*xt(3)                    &
                               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
                               - sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))

           rotMatrixSpectral(nn,3,1) = xt(1)*xt(3)                    &
                               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
                               + sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))
           rotMatrixSpectral(nn,3,2) = xt(2)*xt(3)                    &
                               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
                               + sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))
           rotMatrixSpectral(nn,3,3) = xt(3)*xt(3)                    &
                               + cosTheta*(yt(3)*yt(3) + zt(3)*zt(3))

         endif testUnity

       enddo sectionLoop

       end subroutine timeRotMatricesSpectral
