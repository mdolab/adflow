!
!      ******************************************************************
!      *                                                                *
!      * File:          buildFineBoundaryList.f90                       *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 01-31-2005                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine buildFineBoundaryList(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * buildFineBoundaryList builds the boundary (halo) list for the  *
!      * given fine (or ground) level and spectral mode. This is mostly *
!      * just a straight copy of the information stored in the flowDoms.*
!      * However, if the donors are being treated as guesses then the   *
!      * the interpolants of the list are replaced with the centroid    *
!      * coordinates of the boundary cells.                             *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use boundaryList
       use communication
       use inputOverset
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer                :: ierr
       integer(kind=intType)  :: m, nn, ihalo, i, j, k
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of blocks stored on this processor and
       ! count the total # of overset halos currently stored.

       nHaloOver = 0
       do nn = 1,nDom
         nHaloOver = nHaloOver + flowDoms(nn,level,sps)%nCellsOverset &
                               + flowDoms(nn,level,sps)%nOrphans
       end do

       ! Allocate memory for the halo list for this processor.

       allocate(oversetHalo(nHaloOver), stat=ierr)
       if(ierr /= 0)                            &
        call terminate("buildFineBoundaryList", &
                       "Memory allocation failure for oversetHalo")

       ! Initialize the total halo counter and loop over the blocks.

       ihalo = 0
       domains: do nn = 1,nDom

         ! Set pointers to the current block.

         call setPointers(nn, level, sps)

         ! Change the iblanks on the halos and the boundary to their
         ! appropriate values since only the hole iblanks for the fine
         ! level have been set so far.

         call changeIblanks(.true., 1_intType)
 
         ! Loop over the number of overset cells in the boundary with
         ! donors (not orphans).

         boundaryCopy: do m = 1,nCellsOverset

           ! Update the counter and copy the indices of the boundary
           ! cell and its (possibly approximate) donor cell.
 
           ihalo = ihalo + 1
 
           oversetHalo(ihalo)%myBlock    = nn
           oversetHalo(ihalo)%donorProc  = neighProcOver(m)
           oversetHalo(ihalo)%donorBlock = neighBlockOver(m)
 
           oversetHalo(ihalo)%myI = ibndry(1,m)
           oversetHalo(ihalo)%myJ = ibndry(2,m)
           oversetHalo(ihalo)%myK = ibndry(3,m)
 
           oversetHalo(ihalo)%dI  = idonor(1,m)
           oversetHalo(ihalo)%dJ  = idonor(2,m)
           oversetHalo(ihalo)%dK  = idonor(3,m)

           ! Set levOfInd to 0 so it doesn't corrupt the sort.

           oversetHalo(ihalo)%levOfInd = 0

           ! If the input interpolants for the fine grid are being
           ! ignored then store the centroid, otherwise just copy
           ! the interpolants from the flowDoms.

           if (oversetDonorsAreGuesses) then

             allocate(oversetHalo(ihalo)%interp(3), stat=ierr)
             if(ierr /= 0)                            &
              call terminate("buildFineBoundaryList", &
                             "Memory allocation failure for interp")

             i = oversetHalo(ihalo)%myI
             j = oversetHalo(ihalo)%myJ
             k = oversetHalo(ihalo)%myK

             oversetHalo(ihalo)%interp(1) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,1))
             oversetHalo(ihalo)%interp(2) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,2))
             oversetHalo(ihalo)%interp(3) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,3))

           else

             i = nDonorWeights(oversetInterpType)
             allocate(oversetHalo(ihalo)%interp(i), stat=ierr)
             if(ierr /= 0)                            &
              call terminate("buildFineBoundaryList", &
                             "Memory allocation failure for interp")

             oversetHalo(ihalo)%interp = overint(:,m)

           end if

         end do boundaryCopy

         ! Loop over the number of orphans for this block.

         orphanCopy: do m = 1,nOrphans

           ! Update the counter and copy the indices of the orphan cell.

           ihalo = ihalo + 1

           oversetHalo(ihalo)%myBlock = nn

           i = nCellsOverset + m
           oversetHalo(ihalo)%myI = ibndry(1,i)
           oversetHalo(ihalo)%myJ = ibndry(2,i)
           oversetHalo(ihalo)%myK = ibndry(3,i)

           ! Give the donor info some dummy values. Note that the donor
           ! processor is set to -1 which will maintain the cell as an
           ! orphan. The orphans will be sorted later. 

           oversetHalo(ihalo)%donorProc  = -1
           oversetHalo(ihalo)%donorBlock = 0

           oversetHalo(ihalo)%dI  = 0
           oversetHalo(ihalo)%dJ  = 0
           oversetHalo(ihalo)%dK  = 0

           oversetHalo(ihalo)%levOfInd = 0

           ! If the input interpolants for the fine grid are being
           ! ignored then store the centroid, otherwise just give
           ! the interpolants a value of zero.

           if (oversetDonorsAreGuesses) then

             allocate(oversetHalo(ihalo)%interp(3), stat=ierr)
             if(ierr /= 0)                            &
              call terminate("buildFineBoundaryList", &
                             "Memory allocation failure for interp")

             i = oversetHalo(ihalo)%myI
             j = oversetHalo(ihalo)%myJ
             k = oversetHalo(ihalo)%myK

             oversetHalo(ihalo)%interp(1) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,1))
             oversetHalo(ihalo)%interp(2) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,2))
             oversetHalo(ihalo)%interp(3) = eighth &
                                          * sum(x(i-1:i,j-1:j,k-1:k,3))

           else

             i = nDonorWeights(oversetInterpType)
             allocate(oversetHalo(ihalo)%interp(i), stat=ierr)
             if(ierr /= 0)                            &
              call terminate("buildFineBoundaryList", &
                             "Memory allocation failure for interp")

             oversetHalo(ihalo)%interp = zero

           end if

         end do orphanCopy

         ! If the donors are guesses then the flowDom lists will be 
         ! rebuilt after the donor search process, which may change the
         ! number of boundary cells due to overlap issues; therefore,
         ! deallocate them. Also, add nOrphans to nCellsOverset so that
         ! things are consistent when dealing with orphans after the
         ! donor search process.

         if (oversetDonorsAreGuesses) then
           deallocate(flowDoms(nn,level,sps)%ibndry,         &
                      flowDoms(nn,level,sps)%idonor,         &
                      flowDoms(nn,level,sps)%overint,        &
                      flowDoms(nn,level,sps)%neighBlockOver, &
                      flowDoms(nn,level,sps)%neighProcOver,  &
                      stat=ierr)
           if(ierr /= 0)                            &
            call terminate("buildFineBoundaryList", &
                           "Deallocation failure for overset info")

           flowDoms(nn,level,sps)%nCellsOverset = nCellsOverset &
                                                + nOrphans
         end if
 
       end do domains

       end subroutine buildFineBoundaryList
