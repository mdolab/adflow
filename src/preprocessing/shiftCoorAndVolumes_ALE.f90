!
!      ******************************************************************
!      *                                                                *
!      * File:          shiftCoorAndVolumes.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-03-2004                                      *
!      * Last modified: 03-23-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine shiftCoorAndVolumes_ALE
!
!      ******************************************************************
!      *                                                                *
!      * shiftCoorAndVolumes shifts the owned coordinates and           *
!      * volumes in case of a deforming mesh for an unsteady            *
!      * computation. In this case the old coordinates are needed to    *
!      * determine the mesh velocities. The loop over the number of     *
!      * spectral solutions is present for consistency, but this number *
!      * will be 1 when this routine is called.                         *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, nn, mm, ll, kk
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions and local blocks.

       spectralLoop: do kk=1,nTimeIntervalsSpectral
         domains: do nn=1,nDom

           ! Set the pointers for this block on the ground level.

           call setPointers(nn, groundLevel,kk)

           ! Shift the coordinates already stored in xOld and the
           ! volumes stored in volOld.

           loopOldLevels: do mm=nOldLevels,2,-1

             ! Shift the coordinates from level mm-1 to mm, including
             ! the halo's.

             ll = mm - 1

             do k=0,ke
               do j=0,je
                 do i=0,ie
                   xOld(mm,i,j,k,1) = xOld(ll,i,j,k,1)
                   xOld(mm,i,j,k,2) = xOld(ll,i,j,k,2)
                   xOld(mm,i,j,k,3) = xOld(ll,i,j,k,3)
                 enddo
               enddo
             enddo

             ! Shift the old volumes from level mm-1 to mm.
             ! Only the owned ones need to be considered.

             do k=2,kl
               do j=2,jl
                 do i=2,il
                   volOld(mm,i,j,k) = volOld(ll,i,j,k)
                 enddo
               enddo
             enddo


           enddo loopOldLevels

           ! Shift the current coordinates into the 1st level of xOld.

           do k=0,ke
             do j=0,je
               do i=0,ie
                 xOld(1,i,j,k,1) = x(i,j,k,1)
                 xOld(1,i,j,k,2) = x(i,j,k,2)
                 xOld(1,i,j,k,3) = x(i,j,k,3)
               enddo
             enddo
           enddo

           ! Shift the current volumes into the 1st level of volOld.

           do k=2,kl
             do j=2,jl
               do i=2,il
                 volOld(1,i,j,k) = vol(i,j,k)
               enddo
             enddo
           enddo




           ! ===========================================================
           !
           ! Added by HDN
           !
           ! Shift intermediate meshes from level 1~(nALEsteps-nALEMeshes)
           ! to (nALEMeshes+1)~nALEsteps
           !
           ! ===========================================================
           loopALELevels: do mm=nALEsteps, nALEMeshes+1, -1

              ll = mm - nALEMeshes

              updateI : do k = 1,ke
                 do j = 1,je
                    do i = 0,ie
                       sFaceIALE(mm,i,j,k) = sFaceIALE(ll,i,j,k)
                       sIALE(mm,i,j,k,1)   = sIALE(ll,i,j,k,1)
                       sIALE(mm,i,j,k,2)   = sIALE(ll,i,j,k,2)
                       sIALE(mm,i,j,k,3)   = sIALE(ll,i,j,k,3)
                    enddo
                 enddo
              enddo updateI

              updateJ : do k = 1,ke
                 do j = 0,je
                    do i = 1,ie
                       sFaceJALE(mm,i,j,k) = sFaceJALE(ll,i,j,k)
                       sJALE(mm,i,j,k,1)   = sJALE(ll,i,j,k,1)
                       sJALE(mm,i,j,k,2)   = sJALE(ll,i,j,k,2)
                       sJALE(mm,i,j,k,3)   = sJALE(ll,i,j,k,3)
                    enddo
                 enddo
              enddo updateJ

              updateK : do k = 0,ke
                 do j = 1,je
                    do i = 1,ie
                       sFaceKALE(mm,i,j,k) = sFaceKALE(ll,i,j,k)
                       sKALE(mm,i,j,k,1)   = sKALE(ll,i,j,k,1)
                       sKALE(mm,i,j,k,2)   = sKALE(ll,i,j,k,2)
                       sKALE(mm,i,j,k,3)   = sKALE(ll,i,j,k,3)
                    enddo
                 enddo
              enddo updateK

           enddo loopALELevels

           ! Shift current mesh configurations into latest positions

           ! updateALELevels: do mm=1, nALEMeshes

           !    updateI2 : do k = 1,ke
           !       do j = 1,je
           !          do i = 0,ie
           !             sFaceIALE(mm,i,j,k) = sFaceI(i,j,k)
           !             sIALE(mm,i,j,k,1)   = sI(i,j,k,1)
           !             sIALE(mm,i,j,k,2)   = sI(i,j,k,2)
           !             sIALE(mm,i,j,k,3)   = sI(i,j,k,3)
           !          enddo
           !       enddo
           !    enddo updateI2

           !    updateJ2 : do k = 1,ke
           !       do j = 0,je
           !          do i = 1,ie
           !             sFaceJALE(mm,i,j,k) = sFaceJ(i,j,k)
           !             sJALE(mm,i,j,k,1)   = sJ(i,j,k,1)
           !             sJALE(mm,i,j,k,2)   = sJ(i,j,k,2)
           !             sJALE(mm,i,j,k,3)   = sJ(i,j,k,3)
           !          enddo
           !       enddo
           !    enddo updateJ2

           !    updateK2 : do k = 0,ke
           !       do j = 1,je
           !          do i = 1,ie
           !             sFaceKALE(mm,i,j,k) = sFaceK(i,j,k)
           !             sKALE(mm,i,j,k,1)   = sK(i,j,k,1)
           !             sKALE(mm,i,j,k,2)   = sK(i,j,k,2)
           !             sKALE(mm,i,j,k,3)   = sK(i,j,k,3)
           !          enddo
           !       enddo
           !    enddo updateK2

           ! enddo updateALELevels
           call setLevelALE(1)

         enddo domains
       enddo spectralLoop

       end subroutine shiftCoorAndVolumes_ALE
