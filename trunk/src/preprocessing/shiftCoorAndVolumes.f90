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
       subroutine shiftCoorAndVolumes
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

         enddo domains
       enddo spectralLoop

       end subroutine shiftCoorAndVolumes
