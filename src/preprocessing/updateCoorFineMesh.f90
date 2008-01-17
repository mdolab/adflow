!
!      ******************************************************************
!      *                                                                *
!      * File:          updateCoorFineMesh.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-02-2004                                      *
!      * Last modified: 06-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine updateCoorFineMesh(dtAdvance, sps)
!
!      ******************************************************************
!      *                                                                *
!      * updateCoorFineMesh updates the coordinates of the              *
!      * moving parts of the current finest mesh by the given amount of *
!      * time, possibly different per section. In unsteady mode all the *
!      * times will be equal, but in time spectral mode they can be     *
!      * different.                                                     *
!      * This routine is called in the full mg cycle to put the fine    *
!      * mesh to the position previously calculated on the coarser      *
!      * grid levels, in the unsteady time loop to advance the          *
!      * coordinates only one time step and in the partitioning part    *
!      * of the spectral mode to compute the coordinates of the given   *
!      * spectral solution sps. As it is used in the full MG cycle,     *
!      * currentLevel points to the correct grid level and not          *
!      * ground level.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use flowVarRefState
       use cgnsGrid
       use inputMotion
       use iteration
       use monitor
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps

       real(kind=realType), dimension(*), intent(in) :: dtAdvance
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, nn

       real(kind=realType) :: displX, displY, displZ
       real(kind=realType) :: tNew, tOld
       real(kind=realType) :: angleX, angleY, angleZ, dx, dy, dz
       real(kind=realType) :: xiX, xiY, xiZ, etaX, etaY, etaZ
       real(kind=realType) :: zetaX, zetaY, zetaZ, xp, yp, zp, t
       real(kind=realType) :: phi, cosPhi, sinPhi, eta, zeta

       real(kind=realType), dimension(3)   :: rotationPoint
       real(kind=realType), dimension(3,3) :: rotationMatrix
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the displacements due to the rigid motion of the mesh.

       displX = zero
       displY = zero
       displZ = zero

       ! Determine the time values of the old and new time level.
       ! It is assumed that the rigid body rotation of the mesh is only
       ! used when only 1 section is present.

       tNew = timeUnsteady + timeUnsteadyRestart
       tOld = tNew - dtAdvance(1)

       ! Compute the rotation matrix of the rigid body rotation as
       ! well as the rotation point; the latter may vary in time due
       ! to rigid body translation.

       call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

       ! Loop over the number of local blocks.

       blockLoop: do nn=1,nDom

         ! Set the pointers for this block on the current level.
         ! Note that currentLevel must be used and not groundLevel,
         ! because groundLevel is 1 level too coarse when this routine
         ! is called in the full mg cycle.

         call setPointers(nn, currentLevel, sps)
!
!        ****************************************************************
!        *                                                              *
!        * The rigid body motion of the entire mesh.                    *
!        *                                                              *
!        ****************************************************************
!
         ! First the rotation.

         do k=1,kl
           do j=1,jl
             do i=1,il

               ! Determine the vector relative to the rotation point.

               xp = x(i,j,k,1) - rotationPoint(1)
               yp = x(i,j,k,2) - rotationPoint(2)
               zp = x(i,j,k,3) - rotationPoint(3)

               ! Apply the transformation matrix to the vector (xp,yp,zp)
               ! and set the new coordinates.

               x(i,j,k,1) = rotationMatrix(1,1)*xp &
                          + rotationMatrix(1,2)*yp &
                          + rotationMatrix(1,3)*zp + rotationPoint(1)
               x(i,j,k,2) = rotationMatrix(2,1)*xp &
                          + rotationMatrix(2,2)*yp &
                          + rotationMatrix(2,3)*zp + rotationPoint(2)
               x(i,j,k,3) = rotationMatrix(3,1)*xp &
                          + rotationMatrix(3,2)*yp &
                          + rotationMatrix(3,3)*zp + rotationPoint(3)
             enddo
           enddo
         enddo

         ! Add the translation.

         do k=1,kl
           do j=1,jl
             do i=1,il
               x(i,j,k,1) = x(i,j,k,1) + displX
               x(i,j,k,2) = x(i,j,k,2) + displY
               x(i,j,k,3) = x(i,j,k,3) + displZ
             enddo
           enddo
         enddo

!
!        ****************************************************************
!        *                                                              *
!        * Determine whether the corresponding cgns block is a rotating *
!        * block. If it is, apply the rotation.                         *
!        * Note that now the section ID of the block is taken into      *
!        * account to allow for different periodic times per section.   *
!        *                                                              *
!        ****************************************************************
!
         if( cgnsDoms(nbkGlobal)%rotatingFrameSpecified ) then

           ! Compute the rotation angles.

           angleX = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(1)
           angleY = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(2)
           angleZ = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(3)

           ! Compute the unit vector in the direction of the rotation
           ! axis, which will be called the xi-direction.

           t    = one/max(eps,sqrt(angleX**2 + angleY**2 + angleZ**2))
           xiX = t*angleX
           xiY = t*angleY
           xiZ = t*angleZ

           ! Determine the rotation angle in xi-direction and its sine
           ! and cosine. Due to the definition of the xi-direction this
           ! angle will always be positive.

           phi    = xiX*angleX + xiY*angleY + xiZ*angleZ
           cosPhi = cos(phi)
           sinPhi = sin(phi)

           ! Loop over the owned coordinates of this block.

           do k=1,kl
             do j=1,jl
               do i=1,il

                 ! Compute the vector relative to center of rotation.

                 dx = x(i,j,k,1) - cgnsDoms(nbkGlobal)%rotCenter(1)
                 dy = x(i,j,k,2) - cgnsDoms(nbkGlobal)%rotCenter(2)
                 dz = x(i,j,k,3) - cgnsDoms(nbkGlobal)%rotCenter(3)

                 ! Compute the coordinates of the point p, which is the
                 ! closest point on the rotation axis.

                 t  = dx*xiX + dy*xiY + dz*xiZ
                 xp = cgnsDoms(nbkGlobal)%rotCenter(1) + t*xiX
                 yp = cgnsDoms(nbkGlobal)%rotCenter(2) + t*xiY
                 zp = cgnsDoms(nbkGlobal)%rotCenter(3) + t*xiZ

                 ! Determine the unit vector in eta direction, which
                 ! is defined from point p to the current point.

                 etaX = x(i,j,k,1) - xp
                 etaY = x(i,j,k,2) - yp
                 etaZ = x(i,j,k,3) - zp

                 eta = sqrt(etaX**2 + etaY**2 + etaZ**2)
                 t   = one/max(eps,eta)

                 etaX = t*etaX
                 etaY = t*etaY
                 etaZ = t*etaZ

                 ! Determine the unit vector in zeta-direction. This is
                 ! the cross product of the unit vectors in xi and in
                 ! eta-direction.

                 zetaX = xiY*etaZ - xiZ*etaY
                 zetaY = xiZ*etaX - xiX*etaZ
                 zetaZ = xiX*etaY - xiY*etaX

                 ! Compute the new eta and zeta coordinates.

                 zeta = eta*sinPhi
                 eta  = eta*cosPhi

                 ! Compute the new cartesian coordinates.

                 x(i,j,k,1) = xp + eta*etaX + zeta*zetaX
                 x(i,j,k,2) = yp + eta*etaY + zeta*zetaY
                 x(i,j,k,3) = zp + eta*etaZ + zeta*zetaZ

               enddo
             enddo
           enddo

         endif

       enddo blockLoop

       end subroutine updateCoorFineMesh
