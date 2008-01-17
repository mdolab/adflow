!
!      ******************************************************************
!      *                                                                *
!      * File:          copySpectralSolution.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 08-07-2004                                      *
!      * Last modified: 10-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine copySpectralSolution
!
!      ******************************************************************
!      *                                                                *
!      * copySpectralSolution copies the solution of the 1st spectral   *
!      * solution to all spectral solutions. This typically occurs when *
!      * a for the spectral mode a restart is made from a steady or an  *
!      * unsteady solution. Possible rotation effects are taken into    *
!      * account for the velocity components.                           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use flowVarRefState
       use inputTimeSpectral
       use IOModule
       use monitor
       use section
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: sps, spsm1, mm, nn, i, j, k, l

       real(kind=realType) :: dt, tnew, told, tmp
       real(kind=realType) :: theta, cosTheta, sinTheta

       real(kind=realType), dimension(3)   :: rotPoint, uu
       real(kind=realType), dimension(3)   :: xt, yt, zt
       real(kind=realType), dimension(3,3) :: rotMat

       real(kind=realType), dimension(nSections,3,3) :: rotMatSec
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the rotation matrix from one spectral solution to the
       ! other for every section.

       sectionLoop: do nn=1,nSections

         ! Test if the section is rotating.

         testRotating: if( sections(nn)%rotating ) then

           ! Section is rotating. Determine the angle between the
           ! spectral solutions and its sine and cosine.

           i = sections(nn)%nSlices*nTimeIntervalsSpectral
           theta     = two*pi/real(i,realType)
           cosTheta = cos(theta)
           sinTheta = sin(theta)

           ! Transform to a frame where the xt-axis points in the
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

           rotMatSec(nn,1,1) = xt(1)*xt(1)                          &
                              + cosTheta*(yt(1)*yt(1) + zt(1)*zt(1))
           rotMatSec(nn,1,2) = xt(1)*xt(2)                          &
                             + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
                             - sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
           rotMatSec(nn,1,3) = xt(1)*xt(3)                          &
                             + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
                             - sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))

           rotMatSec(nn,2,1) = xt(1)*xt(2)                          &
                             + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
                             + sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
           rotMatSec(nn,2,2) = xt(2)*xt(2)                          &
                             + cosTheta*(yt(2)*yt(2) + zt(2)*zt(2))
           rotMatSec(nn,2,3) = xt(2)*xt(3)                          &
                             + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
                             - sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))

           rotMatSec(nn,3,1) = xt(1)*xt(3)                          &
                             + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
                             + sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))
           rotMatSec(nn,3,2) = xt(2)*xt(3)                          &
                             + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
                             + sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))
           rotMatSec(nn,3,3) = xt(3)*xt(3)                          &
                             + cosTheta*(yt(3)*yt(3) + zt(3)*zt(3))

         else testRotating

           ! Section is not rotating. Set the rotation matrix to the
           ! identity matrix.

           rotMatSec(nn,1,1) = one
           rotMatSec(nn,1,2) = zero
           rotMatSec(nn,1,3) = zero

           rotMatSec(nn,2,1) = zero
           rotMatSec(nn,2,2) = one
           rotMatSec(nn,2,3) = zero

           rotMatSec(nn,3,1) = zero
           rotMatSec(nn,3,2) = zero
           rotMatSec(nn,3,3) = one

         endif testRotating

       enddo sectionLoop

       ! Initialize told to timeUnsteadyRestart. This takes the
       ! possibility into account that the spectral mode is restarted
       ! from an unsteady computation. Although not likely this
       ! possibility is allowed and should therefore be taken into
       ! account. Anyway told corresponds to the time of the 1st
       ! spectral solution. Also determine the time step between the
       ! spectral solutions. Both told and dt are only used to determine
       ! the rigid body motion and if these are specified it is assumed
       ! that there is only one section present in the grid.

       told = timeUnsteadyRestart
       dt   = sections(1)%timePeriod &
            / real(nTimeIntervalsSpectral,realType)

       ! Loop over the number of spectral modes, starting at 2.

       spectralLoop: do sps=2,nTimeIntervalsSpectral

         ! Determine the corresponding time for this spectral solution
         ! and store sps - 1 a bit easier.

         tnew  = told + dt
         spsm1 = sps - 1

         ! Determine the rotation matrix and rotation point between the
         ! told and tnew for the rigid body rotation of the entire mesh.
         ! The rotation point is not needed for the transformation of the
         ! velocities, but rotMatrixRigidBody happens to compute it.

         call rotMatrixRigidBody(tnew, told, rotMat, rotPoint)

         ! Loop over the local number of blocks.

         domains: do nn=1,nDom

           ! Store the section ID of this block a bit easier in mm.

           mm = flowDoms(nn,1,1)%sectionId

           ! Loop over the owned cells of this block. As the number of
           ! cells is identical for all spectral solutions, it does not
           ! matter which mode is taken for the upper dimensions.

           do k=2,flowDoms(nn,1,1)%kl
             do j=2,flowDoms(nn,1,1)%jl
               do i=2,flowDoms(nn,1,1)%il

                 ! Step 1. Copy the solution variables w from
                 ! the previous spectral solution.

                 do l=1,nw
                   IOVar(nn,sps)%w(i,j,k,l) = IOVar(nn,spsm1)%w(i,j,k,l)
                 enddo

                 ! Step 2. Apply the rigid body motion rotation matrix
                 ! to the velocity. Use uu as a temporary storage.

                 uu(1) = rotMat(1,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                       + rotMat(1,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                       + rotMat(1,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                 uu(2) = rotMat(2,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                       + rotMat(2,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                       + rotMat(2,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                 uu(3) = rotMat(3,1)*IOVar(nn,sps)%w(i,j,k,ivx) &
                       + rotMat(3,2)*IOVar(nn,sps)%w(i,j,k,ivy) &
                       + rotMat(3,3)*IOVar(nn,sps)%w(i,j,k,ivz)

                 ! Step 3. Apply the rotation matrix of the section to
                 ! the velocity.

                 IOVar(nn,sps)%w(i,j,k,ivx) = rotMatSec(mm,1,1)*uu(1) &
                                            + rotMatSec(mm,1,2)*uu(2) &
                                            + rotMatSec(mm,1,3)*uu(3)

                 IOVar(nn,sps)%w(i,j,k,ivy) = rotMatSec(mm,2,1)*uu(1) &
                                            + rotMatSec(mm,2,2)*uu(2) &
                                            + rotMatSec(mm,2,3)*uu(3)

                 IOVar(nn,sps)%w(i,j,k,ivz) = rotMatSec(mm,3,1)*uu(1) &
                                            + rotMatSec(mm,3,2)*uu(2) &
                                            + rotMatSec(mm,3,3)*uu(3)
               enddo
             enddo
           enddo

         enddo domains

         ! Set told to tnew for the next spectral solution.

         told = tnew

       enddo spectralLoop

       end subroutine copySpectralSolution
