!
!      ******************************************************************
!      *                                                                *
!      * File:          writeCoorCGNSZone.F90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-20-2004                                      *
!      * Last modified: 10-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeCoorCGNSZone(zone, cgnsZone)
!
!      ******************************************************************
!      *                                                                *
!      * writeCoorCGNSZone writes the coordinates of the given zone     *
!      * to the cgns file(s).                                           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use cgnsNames
       use communication
       use inputIO
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(in) :: cgnsZone

       integer(kind=intType), intent(in) :: zone

#ifdef USE_NO_CGNS
       call terminate("writeCoorCGNSZone", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr, tmp
       integer :: bufSize, realTypeCGNS, cgnsBase, fileInd

       integer, dimension(mpi_status_size) :: status

       integer, dimension(:), allocatable :: proc

       integer(kind=intType) :: i, j, nn, mm, ll, ind
       integer(kind=intType) :: nBlocks, nSubBlocks, offset
       integer(kind=intType) :: sizeCGNSWriteType

       integer(kind=intType), dimension(6)     :: ii
       integer(kind=intType), dimension(nProc) :: nMessages

       integer(kind=intType), dimension(:,:,:), allocatable :: subRanges

       real(kind=realType), dimension(:), allocatable :: buffer

       character, dimension(:), allocatable :: coor

       character(len=maxCGNSNameLen), dimension(3) :: coorNames
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cgns real type depending on the input option.

       select case (precisionGrid)
         case (precisionSingle)
           realTypeCGNS      = RealSingle
           sizeCGNSWriteType = 4
         case (precisionDouble)
           realTypeCGNS      = RealDouble
           sizeCGNSWriteType = 8
       end select

       ! Store the number of local blocks and the offset in
       ! blocksCGNSblock for this zone a bit easier.

       offset  = nBlocksCGNSblock(zone-1)
       nBlocks = nBlocksCGNSblock(zone) - offset

       ! Determine the amount of block parts each processor will send to
       ! processor 0.

       call mpi_gather(nBlocks, 1, sumb_integer, nMessages, 1, &
                       sumb_integer, 0, SUmb_comm_world, ierr)

       ! At the moment the writing of the cgns file is sequential and done
       ! by processor 0. This means that this processor gathers all info
       ! from the other processors and writes it to file.

       rootproc: if(myID == 0) then

         ! I am processor 0 and poor me has to do all the work.

         ! First determine the number of subblocks into the original cgns
         ! block is split.

         nSubBlocks = 0
         do i=1,nProc
           nSubBlocks = nSubBlocks + nMessages(i)
         enddo

         ! Allocate the memory for the ranges and the processor
         ! where the subblock is stored.

         allocate(subRanges(3,2,nSubBlocks), proc(nSubBlocks), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("writeCoorCGNSZone", &
                          "Memory allocation failure for subRanges &
                          &and proc")

         ! Determine the processor ID's where the subRanges are stored.
         ! Note that 1 must be substracted, because the processor
         ! numbering starts at 0.

         nSubBlocks = 0
         do i=1,nProc
           do j=1,nMessages(i)
             nSubBlocks = nSubBlocks + 1
             proc(nSubBlocks) = i - 1
           enddo
         enddo

         ! Determine the subRanges of the subblocks stored on this
         ! processor.  Note that nBlocks can be 0.

         do i=1,nBlocks

           ! Store the local block ID a bit easier in j and copy the
           ! range. This range is identical for spectral solutions and
           ! thus taking the first is okay.

           j = blocksCGNSblock(i+offset)

           subRanges(1,1,i) = flowDoms(j,1,1)%iBegor
           subRanges(1,2,i) = flowDoms(j,1,1)%iEndor

           subRanges(2,1,i) = flowDoms(j,1,1)%jBegor
           subRanges(2,2,i) = flowDoms(j,1,1)%jEndor

           subRanges(3,1,i) = flowDoms(j,1,1)%kBegor
           subRanges(3,2,i) = flowDoms(j,1,1)%kEndor

           ! To avoid duplication of overlap regions, add 1 to the lower
           ! boundaries if the lower boundary is larger than 1.

           if(subRanges(1,1,i) > 1) subRanges(1,1,i) = subRanges(1,1,i) +1
           if(subRanges(2,1,i) > 1) subRanges(2,1,i) = subRanges(2,1,i) +1
           if(subRanges(3,1,i) > 1) subRanges(3,1,i) = subRanges(3,1,i) +1

         enddo

         ! The rest of the block ranges must be obtained by
         ! communication.

         do i=(nBlocks+1),nSubBlocks

           call mpi_recv(ii, 6, sumb_integer, proc(i), proc(i), &
                         SUmb_comm_world, status, ierr)

           subRanges(1,1,i) = ii(1)
           subRanges(1,2,i) = ii(2)
           subRanges(2,1,i) = ii(3)
           subRanges(2,2,i) = ii(4)
           subRanges(3,1,i) = ii(5)
           subRanges(3,2,i) = ii(6)
         enddo

         ! Determine the size of the largest subblock and allocate
         ! the memory for the corresponding buffer.

         bufSize = 0
         do i=1,nSubBlocks
           ll = (subRanges(1,2,i) - subRanges(1,1,i) + 1) &
              * (subRanges(2,2,i) - subRanges(2,1,i) + 1) &
              * (subRanges(3,2,i) - subRanges(3,1,i) + 1)
           bufSize = max(bufSize, ll)
         enddo

         allocate(buffer(bufSize), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("writeCoorCGNSZone", &
                          "Memory allocation failure for buffer")

         ! Allocate the memory for the array used to write the three
         ! coordinates and set the cgns names for them. Note that the
         ! coor array is of type character and therefore the size in
         ! bytes must be allocated.

         ll = cgnsDoms(zone)%il * cgnsDoms(zone)%jl &
            * cgnsDoms(zone)%kl * sizeCGNSWriteType
         allocate(coor(ll), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("writeCoorCGNSZone", &
                          "Memory allocation failure for coor")

         coorNames(1) = cgnsCoorx
         coorNames(2) = cgnsCoory
         coorNames(3) = cgnsCoorz

         ! Loop over the number of grid files to be written.

         gridLoopRoot: do ind=1,nGridsToWrite

           ! Store the file and base ID a bit easier.

           fileInd  = fileIDs(ind)
           cgnsBase = cgnsBases(ind)

           ! Loop over the three coordinates.

           coorLoopRoot: do nn=1,3

             ! Loop over the number of subblocks stored
             ! on this processor.

             do mm=1,nBlocks

               ! Fill buffer with the correct coordinate.

               call storeCoorInBuffer(buffer, zone, ind, nn, &
                                      blocksCGNSblock(mm+offset), tmp)

               ! And store it in coor, depending on the precision used.

               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(coor, buffer,      &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   cgnsDoms(zone)%il, &
                                                   cgnsDoms(zone)%jl, &
                                                   cgnsDoms(zone)%kl, &
                                                   subRanges(1,1,mm))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(coor, buffer,      &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   cgnsDoms(zone)%il, &
                                                   cgnsDoms(zone)%jl, &
                                                   cgnsDoms(zone)%kl, &
                                                   subRanges(1,1,mm))
               end select

             enddo

             ! Loop over the number of subblocks stored on
             ! other processors.

             do mm=(nBlocks+1),nSubBlocks

               ! Receive the range of subblock mm and copy it into coor.

               call mpi_recv(buffer, bufSize, sumb_real, proc(mm), &
                             proc(mm)+1, SUmb_comm_world, status, ierr)

               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(coor, buffer,      &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   cgnsDoms(zone)%il, &
                                                   cgnsDoms(zone)%jl, &
                                                   cgnsDoms(zone)%kl, &
                                                   subRanges(1,1,mm))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(coor, buffer,      &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   1_intType,         &
                                                   cgnsDoms(zone)%il, &
                                                   cgnsDoms(zone)%jl, &
                                                   cgnsDoms(zone)%kl, &
                                                   subRanges(1,1,mm))
               end select

             enddo

             ! Write this coordinate to file; tmp is used to store
             ! the actual number of the coordinate; usually this is
             ! equal to nn.

             call cg_coord_write_f(fileInd, cgnsBase, cgnsZone, &
                                   realTypeCGNS, coorNames(nn), &
                                   coor, tmp, ierr)
             if(ierr /= all_ok)                    &
               call terminate("writeCoorCGNSZone", &
                              "Something wrong when calling &
                              &cg_coord_write_f")

             ! Write the units, if possible.

             if( cgnsDoms(zone)%gridUnitsSpecified ) then

               ! Go to the correct place in the grid file.

               call cg_goto_f(fileInd, cgnsBase, ierr, &
                              "Zone_t", cgnsZone,      &
                              "GridCoordinates_t", 1,  &
                              "DataArray_t", tmp, "end")
               if(ierr /= all_ok)                    &
                 call terminate("writeCoorCGNSZone", &
                                "Something wrong when calling cg_goto_f")

               ! Write the units.

               call cg_units_write_f(cgnsDoms(zone)%mass, &
                                     cgnsDoms(zone)%len,  &
                                     cgnsDoms(zone)%time, &
                                     cgnsDoms(zone)%temp, &
                                     cgnsDoms(zone)%angle, ierr)
               if(ierr /= all_ok)                    &
                 call terminate("writeCoorCGNSZone", &
                                "Something wrong when calling &
                                &cg_units_write_f")
             endif

           enddo coorLoopRoot

         enddo gridLoopRoot

         ! Deallocate the memory which is only allocated on the
         ! root processor.

         deallocate(subRanges, proc, coor, stat=ierr)
         if(ierr /= 0) call terminate("writeCoorCGNSZone", &
                                      "Deallocation error on root proc")

       else rootproc

         ! I am not the root processor and I may have to send data to
         ! the root processor.

         ! Loop over the number of subblocks stored on this processor
         ! to send the size to the root processor. Determine in the
         ! same loop the size of the largest subblock.
         
         bufSize = 0
         do i=1,nBlocks

           ! Store the local block ID a bit easier in j and copy the
           ! range. This range is identical for spectral solutions and
           ! thus taking the first is okay.

           j = blocksCGNSblock(i+offset)

           ii(1) = flowDoms(j,1,1)%iBegor
           ii(2) = flowDoms(j,1,1)%iEndor
           ii(3) = flowDoms(j,1,1)%jBegor
           ii(4) = flowDoms(j,1,1)%jEndor
           ii(5) = flowDoms(j,1,1)%kBegor
           ii(6) = flowDoms(j,1,1)%kEndor

           ! To avoid duplication of overlap regions, add 1 to the lower
           ! boundaries if the lower boundary is larger than 1.

           if(ii(1) > 1) ii(1) = ii(1) +1
           if(ii(3) > 1) ii(3) = ii(3) +1
           if(ii(5) > 1) ii(5) = ii(5) +1

           ! Send the buffer ii to processor 0.

           call mpi_send(ii, 6, sumb_integer, 0, myID, &
                         SUmb_comm_world, ierr)

           ! Check the size of this subblock and update bufSize
           ! if needed.

           ll = (ii(2) - ii(1) + 1) * (ii(4) - ii(3) + 1) &
              * (ii(6) - ii(5) + 1)
           bufSize = max(bufSize, ll)
         enddo

         ! Allocate the memory for buffer.

         allocate(buffer(bufSize), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("writeCoorCGNSZone", &
                          "Memory allocation failure for buffer")

         ! Loop over the number of grids to be written.

         gridLoopOthers: do ind=1,nGridsToWrite

           ! Loop over the three coordinates.

           coorLoopOthers: do nn=1,3

             ! Loop over the number of subblocks stored
             ! on this processor.

             do mm=1,nBlocks

               ! Fill buffer with the correct coordinate and send it to
               ! processor 0.

               call storeCoorInBuffer(buffer, zone, ind, nn, &
                                      blocksCGNSblock(mm+offset), tmp)

               call mpi_send(buffer, tmp, sumb_real, 0, myID+1, &
                             SUmb_comm_world, ierr)
             enddo

           enddo coorLoopOthers
         enddo gridLoopOthers

       endif rootproc

       ! Release the memory of buffer.

       deallocate(buffer, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("writeCoorCGNSZone", &
                        "Deallocation failure for buffer")

#endif

       end subroutine writeCoorCGNSZone

!      ==================================================================

       subroutine storeCoorInBuffer(buffer, zone, ind, coorID, &
                                    blockID, nn)
!
!      ******************************************************************
!      *                                                                *
!      * storeCoorInBuffer stores the given coordinate for the given    *
!      * blockID. The total size of the buffer is returned in nn.       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer, intent(out)              :: nn
       integer(kind=intType), intent(in) :: zone, ind
       integer(kind=intType), intent(in) :: coorID, blockID

       real(kind=realType), dimension(*), intent(out) :: buffer
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: iStart, iEnd, jStart, jEnd, kStart, kEnd

       real(kind=realType) :: LRefInv
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Compute the multiplication factor to obtain the original
       ! coordinates. Note that LRef is corrected to 1.0 when the
       ! coordinates should be written in meters. This happens when
       ! the grid is read.

       LRefInv = one/cgnsDoms(zone)%LRef

       ! Set the range to be stored to the entire block.

       iStart = 1
       iEnd   = flowDoms(blockID,1,1)%il
       jStart = 1
       jEnd   = flowDoms(blockID,1,1)%jl
       kStart = 1
       kEnd   = flowDoms(blockID,1,1)%kl

       ! To avoid duplication of overlap regions, add 1 to the lower
       ! boundaries if the lower boundary in the original block is
       ! larger than 1.

       if(flowDoms(blockID,1,1)%iBegor > 1) iStart = iStart + 1
       if(flowDoms(blockID,1,1)%jBegor > 1) jStart = jStart + 1
       if(flowDoms(blockID,1,1)%kBegor > 1) kStart = kStart + 1

       ! Store the coordinate in the 1D buffer.

       nn = 0
       do k=kStart,kEnd
         do j=jStart,jEnd
           do i=iStart,iEnd
             nn = nn + 1
             buffer(nn) = LRefInv*IOVar(blockID,ind)%w(i,j,k,coorID)
           enddo
         enddo
       enddo

       end subroutine storeCoorInBuffer
