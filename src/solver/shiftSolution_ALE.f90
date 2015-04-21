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
       subroutine shiftSolution_ALE
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

       ! Loop over the number of spectral solutions and local blocks.
       ! Although this routine is only called in unsteady mode where the
       ! number of spectral solutions is 1, this loop is there just for
       ! consistency.

       ! *******************************
       ! REMOVED the rigid body rotation part for simplicity
       ! *******************************
       
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

             do l = 1,nw
                do k = 2,kl
                   do j = 2,jl
                      do i = 2,il
                         dwALE(mm,i,j,k,l) = dwALE(ll,i,j,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do l = 1,nwf
                do k = 2,kl
                   do j = 2,jl
                      do i = 2,il
                         fwALE(mm,i,j,k,l) = fwALE(ll,i,j,k,l)
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

           do l = 1,nw
              do k = 2,kl
                 do j = 2,jl
                    do i = 2,il
                       dwALE(1,i,j,k,l) = zero
                    enddo
                 enddo
              enddo
           enddo

           do l = 1,nwf
              do k = 2,kl
                 do j = 2,jl
                    do i = 2,il
                       fwALE(1,i,j,k,l) = zero
                    enddo
                 enddo
              enddo
           enddo

         enddo domains
       enddo spectralLoop

       end subroutine shiftSolution_ALE
