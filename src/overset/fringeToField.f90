!
!      ******************************************************************
!      *                                                                *
!      * File:          fringeToField.f90                               *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-17-2005                                      *
!      * Last modified: 08-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine fringeToField(level, sps, blk, fi, fj, fk,       &
                                donorProc, donorBlock, di, dj, dk, &
                                ihalo, n1to1HaloAdd)
!
!      ******************************************************************
!      *                                                                *
!      * fringeToField changes the given fringe cell given by blk and   *
!      * indices (fi,fj,fk) to a field cell because its donor quality   *
!      * may be poor. The residual stencil for the new field cell is    *
!      * then checked to add additional fringe if needed. New fringe    *
!      * is added to the boundary list unless it is a 1-to-1 halo, in   *
!      * which case the iblank is set to an index value to be dealt     *
!      * with later. The input donor processor, block, and indices are  *
!      * used as the starting guesses in the search for any new fringe  *
!      ^ added to the list.                                             *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, blk, fi, fj, fk
       integer(kind=intType), intent(in) :: donorProc, donorBlock
       integer(kind=intType), intent(in) :: di, dj, dk

       integer(kind=intType), intent(inout) :: ihalo, n1to1HaloAdd
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, m, n, del(3)

       logical :: donorInfoNotStored
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize whether the donor info was recorded.

       donorInfoNotStored = .true.

       ! Set the iblank of the new field cell to 1 and update the
       ! number of overset cells for the block.

       flowDoms(blk,level,sps)%iblank(fi,fj,fk) = 1

       flowDoms(blk,level,sps)%nCellsOverset = &
       flowDoms(blk,level,sps)%nCellsOverset - 1

       ! Loop over a layer of cells surrounding the block, which are
       ! needed for bother inviscid and viscous residuals.

       do k = fk-1,fk+1
         do j = fj-1,fj+1
           do i = fi-1,fi+1

             if (flowDoms(blk,level,sps)%iblank(i,j,k) <= 0) then
               call addCellToFringe
             end if

           end do
         end do
       end do

       ! Loop over the additional fringe size in each indical direction
       ! and again look for holes that need to be made fringe.

       do n = 2,fringeSize
         do m = 1,3

           ! Reset the delta index vector for the currnt direction.

           del    = 0
           del(m) = n

           ! Compute and check the indices in the positive direction.

           i = fi + del(1); j = fj + del(2); k = fk + del(3)

           if (flowDoms(blk,level,sps)%iblank(i,j,k) <= 0) then
             call addCellToFringe
           end if

           ! Compute and check the indices in the negative direction.

           i = fi - del(1); j = fj - del(2); k = fk - del(3)

           if (flowDoms(blk,level,sps)%iblank(i,j,k) <= 0) then 
             call addCellToFringe
           end if

         end do
       end do

!      ==================================================================

       contains

         subroutine addCellToFringe
!
!        ****************************************************************
!        *                                                              *
!        * addCellToFringe is a procedure to add the current hole cell  *
!        * given in the loops above to the boundary list, or if it's a  *
!        * a 1-to-1 halo, the iblank for the cell is set to an index    *
!        * signal and the donor info is stored.                         *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local variables.
!
         integer(kind=intType) :: ii

         real(kind=realType), dimension(3) :: xc
!
!        Interfaces.
!
         interface
           subroutine reallocateInteger2(intArray,            &
                                          newSize1, newSize2, &
                                          oldSize1, oldSize2, &
                                          alwaysFreeMem)
             use precision
             implicit none

             integer(kind=intType), dimension(:,:), pointer :: intArray
             integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                  oldSize1, oldSize2
             logical, intent(in) :: alwaysFreeMem
           end subroutine reallocateInteger2
         end interface
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Check whether the cell to be added is a 1-to-1 halo. If yes,
         ! then the donor info is stored and the iblank of the cell is
         ! set to 10 + the index into the donorInfo array. If not, then
         ! simply add it to the halo list for this processor.

         if (flowDoms(blk,level,sps)%iblank(i,j,k) == 0   .and. &
             (i <= 1 .or. i >= flowDoms(blk,level,sps)%ie .or.  &
              j <= 1 .or. j >= flowDoms(blk,level,sps)%je .or.  &
              k <= 1 .or. k >= flowDoms(blk,level,sps)%ke)) then

           ! Check whether the donor info has already been stored.

           if (donorInfoNotStored) then

             ! Update the counter for the number of donor info sets.

             nDonorInfo = nDonorInfo + 1

             ! Check whether the current array size is adequate. If
             ! not, reallocate with more space.

             ii = ubound(donorInfo, 2)

             if (nDonorInfo > ii) then
               call reallocateInteger2(donorInfo, 5_intType, ii+100, &
                                       5_intType, ii, .true.)
             end if

             ! Copy the donor info to the array and reset the logical.

             donorInfo(1,nDonorInfo) = donorProc
             donorInfo(2,nDonorInfo) = donorBlock
             donorInfo(3,nDonorInfo) = di
             donorInfo(4,nDonorInfo) = dj
             donorInfo(5,nDonorInfo) = dk

             donorInfoNotStored = .false.

           end if

           ! Set the iblank to a value that allows easy determination
           ! of its donor estimate info.

           flowDoms(blk,level,sps)%iblank(i,j,k) = 10 + nDonorInfo

           ! Update the total number of 1-to-1 halos that I will need
           ! to communicate to connecting blocks.

           n1to1HaloAdd = n1to1HaloAdd + 1

         else

           ! Check to make sure the size of oversetHalo is sufficient.

           call checkSizeBoundaryList(ihalo, 1_intType)

           ! Update the counter and add the new boundary cell info.

           ihalo = ihalo + 1

           oversetHalo(ihalo)%myBlock = blk
           oversetHalo(ihalo)%myI     = i
           oversetHalo(ihalo)%myJ     = j
           oversetHalo(ihalo)%myK     = k

           ! Compute the cell centroid using the external routine
           ! since the boundary could contain 2nd level halos.

           call cellCentroid(blk, level, sps, i, j, k, &
                             xc(1), xc(2), xc(3))

           oversetHalo(ihalo)%interp = eighth*xc

           ! Add the donor info.

           oversetHalo(ihalo)%donorProc  = donorProc
           oversetHalo(ihalo)%donorBlock = donorBlock
           oversetHalo(ihalo)%dI         = di
           oversetHalo(ihalo)%dJ         = dj
           oversetHalo(ihalo)%dK         = dk

           ! Update the number of overset cells for this block, and
           ! change the iblank to reflect it's now part of the 
           ! boundary.

           flowDoms(blk,level,sps)%nCellsOverset = &
           flowDoms(blk,level,sps)%nCellsOverset + 1

           flowDoms(blk,level,sps)%iblank(i,j,k) = 9

         end if

         end subroutine addCellToFringe

       end subroutine fringeToField
