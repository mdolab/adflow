!
!      ******************************************************************
!      *                                                                *
!      * File:          shiftSolution.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-04-2004                                      *
!      * Last modified: 08-30-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine shiftSolution
!
!      ******************************************************************
!      *                                                                *
!      * shiftSolution shifts the solution of the older time levels,    *
!      * such that a new time step can be started.                      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use flowVarRefState
       use inputTimeSpectral
       use inputUnsteady
       use iteration
       use monitor
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, l, sps, nn, mm, ll

       real(kind=realType) :: tOld, tNew
       real(kind=realType) :: vX, vY, vZ, vXi, vEta, vZeta
       real(kind=realType) :: t, angleX, angleY, angleZ
       real(kind=realType) :: phi, cosPhi, sinPhi
       real(kind=realType) :: xiX, xiY, xiZ, etaX, etaY, etaZ
       real(kind=realType) :: zetaX, zetaY, zetaZ

       real(kind=realType), dimension(3)   :: rotationPoint
       real(kind=realType), dimension(3,3) :: rotationMatrix
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the rotation matrix of the rigid body rotation as well
       ! as the rotation point; the latter is not needed to correct the
       ! velocities, but the routine rotMatrixRigidBody is also used
       ! for coordinates.

       tNew = timeUnsteady + timeUnsteadyRestart
       tOld = tNew         - deltaT

       call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

       ! Loop over the number of spectral solutions and local blocks.
       ! Although this routine is only called in unsteady mode where the
       ! number of spectral solutions is 1, this loop is there just for
       ! consistency.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block on the ground level.

           call setPointers(nn, groundLevel, sps)

           ! Shift the solution already stored in wOld.

           loopOldLevels: do mm=nOldLevels,2,-1

             ! Shift the owned solution variables from level mm-1 to mm.

             ll = mm - 1

             do l=1,nw
               do k=2,kl
                 do j=2,jl
                   do i=2,il
                     wOld(mm,i,j,k,l) = wOld(ll,i,j,k,l)
                   enddo
                 enddo
               enddo
             enddo

           enddo loopOldLevels

           ! Shift the current solution into the 1st level of wOld.
           ! Note that in wOld the conservative flow variables are stored,
           ! while in w the velocity components are stored and not
           ! the momentum. Therefore this must be corrected.
           ! Also the turbulent primitive variables are stored, but this
           ! is okay, because the quasi-linear form of the turbulent
           ! transport equations is solved and not the conservative one.

           do l=1,nw
             do k=2,kl
               do j=2,jl
                 do i=2,il
                   wOld(1,i,j,k,l) = w(i,j,k,l)
                 enddo
               enddo
             enddo
           enddo

           ! Make sure that the momentum variables are stored in wOld.

           do k=2,kl
             do j=2,jl
               do i=2,il
                 wOld(1,i,j,k,ivx) = wOld(1,i,j,k,ivx)*wOld(1,i,j,k,irho)
                 wOld(1,i,j,k,ivy) = wOld(1,i,j,k,ivy)*wOld(1,i,j,k,irho)
                 wOld(1,i,j,k,ivz) = wOld(1,i,j,k,ivz)*wOld(1,i,j,k,irho)
               enddo
             enddo
           enddo

           ! To improve the initial guess of the velocity field the
           ! velocity of rotating parts is rotated. First the rigid
           ! body motion.

           do k=2,kl
             do j=2,jl
               do i=2,il
                 vX = w(i,j,k,ivx)
                 vY = w(i,j,k,ivy)
                 vZ = w(i,j,k,ivz)

                 w(i,j,k,ivx) = rotationMatrix(1,1)*vX &
                              + rotationMatrix(1,2)*vY &
                              + rotationMatrix(1,3)*vZ
                 w(i,j,k,ivy) = rotationMatrix(2,1)*vX &
                              + rotationMatrix(2,2)*vY &
                              + rotationMatrix(2,3)*vZ
                 w(i,j,k,ivz) = rotationMatrix(3,1)*vX &
                              + rotationMatrix(3,2)*vY &
                              + rotationMatrix(3,3)*vZ
               enddo
             enddo
           enddo

           ! Apply an additional correction for the velocity components
           ! if a rotation rate is prescribed for this block.

           rotTest: if( cgnsDoms(nbkGlobal)%rotatingFrameSpecified ) then

             ! Compute the rotation angles.

             angleX = deltaT*cgnsDoms(nbkGlobal)%rotRate(1)
             angleY = deltaT*cgnsDoms(nbkGlobal)%rotRate(2)
             angleZ = deltaT*cgnsDoms(nbkGlobal)%rotRate(3)

             ! Compute the unit vector in the direction of the rotation
             ! axis, which will be called the xi-direction.

             t   = one/max(eps,sqrt(angleX**2 + angleY**2 + angleZ**2))
             xiX = t*angleX
             xiY = t*angleY
             xiZ = t*angleZ

             ! Determine the rotation angle in xi-direction and its sine
             ! and cosine. Due to the definition of the xi-direction this
             ! angle will always be positive.

             phi    = xiX*angleX + xiY*angleY + xiZ*angleZ
             cosPhi = cos(phi)
             sinPhi = sin(phi)

             ! Loop over the cell centers.

             do k=2,kl
               do j=2,jl
                 do i=2,il

                   ! Abbreviate the velocity components a bit easier.

                   vX = w(i,j,k,ivx)
                   vY = w(i,j,k,ivy)
                   vZ = w(i,j,k,ivz)

                   ! Determine the component of the velocity vector
                   ! in xi direction and determine the direction eta,
                   ! the direction of the velocity when the xi component
                   ! is substracted.

                   vXi = vX*xiX + vY*xiY + vZ*xiZ

                   etaX = vX - vXi*xiX
                   etaY = vY - vXi*xiY
                   etaZ = vZ - vXi*xiZ

                   t    = one/max(eps,sqrt(etaX**2 + etaY**2 + etaZ**2))
                   etaX = t*etaX
                   etaY = t*etaY
                   etaZ = t*etaZ

                   ! Determine the velocity component in eta direction.

                   vEta = vX*etaX + vY*etaY + vZ*etaZ

                   ! Determine the unit vector in zeta-direction. This is
                   ! the cross product of the unit vectors in xi and in
                   ! eta-direction.

                   zetaX = xiY*etaZ - xiZ*etaY
                   zetaY = xiZ*etaX - xiX*etaZ
                   zetaZ = xiX*etaY - xiY*etaX

                   ! Determine the velocity components in eta and zeta
                   ! direction after the rotation.

                   vZeta = vEta*sinPhi
                   vEta  = vEta*cosPhi

                   ! Compute the new Cartesian velocity components.

                   w(i,j,k,ivx) = vXi*xiX + vEta*etaX + vZeta*zetaX
                   w(i,j,k,ivy) = vXi*xiY + vEta*etaY + vZeta*zetaY
                   w(i,j,k,ivz) = vXi*xiZ + vEta*etaZ + vZeta*zetaZ

                 enddo
               enddo
             enddo

           endif rotTest

         enddo domains
       enddo spectralLoop

       end subroutine shiftSolution
