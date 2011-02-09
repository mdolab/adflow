!
!      ******************************************************************
!      *                                                                *
!      * File:          writeSolCGNSZone.F90                            *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 04-12-2003                                      *
!      * Last modified: 10-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeSolCGNSZone(zone, nSolVar, nDiscrVar, solNames)
!
!      ******************************************************************
!      *                                                                *
!      * writeSolCGNSZone writes a volume solution of the given zone    *
!      * to the cgns file(s). In case the solution must be written to a *
!      * separate file, useLinksInCGNS == .true., a link to the zone of *
!      * the grid file is created.                                      *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use cgnsGrid
       use cgnsNames
       use communication
       use flowVarRefState
       use inputIO
       use inputPhysics
       use iteration
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: zone
       integer(kind=intType), intent(in) :: nSolVar, nDiscrVar

       character(len=*), dimension(*), intent(in) :: solNames

#ifdef USE_NO_CGNS

       call terminate("writeSolCGNSZone", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
!
!      Local variables.
!
       integer :: ierr
       integer :: source, bufSize, size, nnVar
       integer :: cgnsInd, cgnsBase, cgnsZone, cgnsSol, realTypeCGNS

       integer, dimension(mpi_status_size) :: status
       integer, dimension(9)               :: sizes

       integer(kind=intType) :: i, j, nn, mm, ll, ind, nVarWritten
       integer(kind=intType) :: nBlocks, nSubblocks, offset
       integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
       integer(kind=intType) :: iBegCGNS, jBegCGNS, kBegCGNS
       integer(kind=intType) :: iEndCGNS, jEndCGNS, kEndCGNS
       integer(kind=intType) :: sizeCGNSWriteType

       integer(kind=intType), dimension(nProc) :: nMessages

       integer(kind=intType), dimension(:), allocatable :: proc

       integer(kind=intType), dimension(:,:,:), allocatable :: subRanges

       real(kind=realType), dimension(:), allocatable :: buffer

       character, dimension(:), allocatable :: sol

       logical :: rindLayerThisSol, unsteadyHigherSol, writeLink

       character(len=maxStringLen) :: linkName, solName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cgns real type depending on the input option.

       select case (precisionSol)
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
!
!        ****************************************************************
!        *                                                              *
!        * I am processor 0 and poor me has to do all the work.         *
!        *                                                              *
!        ****************************************************************
!
         ! Allocate the memory for the array used to write the solution
         ! to file. The size depends whether or not rind layers are to
         ! be written and the precision of the floating point type.

         if( storeRindLayer ) then
           ll = (cgnsDoms(zone)%kl+1) * (cgnsDoms(zone)%jl+1) &
              * (cgnsDoms(zone)%il+1) * sizeCGNSWriteType
         else
           ll = (cgnsDoms(zone)%kl-1) * (cgnsDoms(zone)%jl-1) &
              * (cgnsDoms(zone)%il-1) * sizeCGNSWriteType
         endif

         allocate(sol(ll), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("writeSolCGNSZone", &
                          "Memory allocation failure for sol")

         ! First determine the number of subblocks into the original cgns
         ! block is split.

         nSubblocks = 0
         do i=1,nProc
           nSubblocks = nSubblocks + nMessages(i)
         enddo

         ! Allocate the memory for the ranges and the processor
         ! where the subblock is stored.

         allocate(subRanges(3,2,nSubblocks), proc(nSubblocks), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("writeSolCGNSZone", &
                          "Memory allocation failure for subRanges &
                          &and proc")

         ! Determine the processor ID's where the subRanges are stored.
         ! Note that 1 must be substracted, because the processor numbering
         ! starts at 0.

         nSubblocks = 0
         do i=1,nProc
           do j=1,nMessages(i)
             nSubblocks = nSubblocks + 1
             proc(nSubblocks) = i - 1
           enddo
         enddo

         ! Determine the subranges for the 1st solution.

         rindLayerThisSol = storeRindLayer
         call getSubRangesSol

         ! Allocate the memory for buffer.

         allocate(buffer(bufSize), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("writeSolCGNSZone", &
                          "Memory allocation failure for buffer")

         ! Loop over the number of solutions.

         solLoopRoot: do ind=1,nVolSolToWrite

           ! Determine whether or not we are dealing with an unsteady
           ! higher solution here.

           unsteadyHigherSol = .false.
           if(ind > 1 .and. equationMode == unsteady) &
             unsteadyHigherSol = .true.

           ! For unsteady mode on rigid meshes it is possible that
           ! only the solution must be written; the coordinates are
           ! not written to a file. In case links are used for that
           ! case the link should not be created. The logical
           ! writeLink takes care of that.

           writeLink = useLinksInCGNS
           if(unsteadyHigherSol .and. (.not. deforming_Grid)) &
             writeLink = .false.

           ! Store the file and base ID a bit easier and set
           ! rindLayerThisSol. A rind layer is not written for the
           ! higher solution in unsteady mode.

           cgnsInd          = fileIDs(ind)
           cgnsBase         = cgnsBases(ind)
           rindLayerThisSol = storeRindLayer
           if( unsteadyHigherSol ) rindLayerThisSol = .false.

           ! Check if the subranges must be recomputed. Some dirty stuff
           ! must be done, because Fortran 90/95 does not allow a
           ! comparison between logicals.

           nn = 0; if( storeRindLayer )   nn = 1
           mm = 0; if( rindLayerThisSol ) mm = 1
           if(nn /= mm) call getSubRangesSol

           ! Check whether a zone must be created, which is
           ! true if links are used. Create the zone, if needed.

           createZoneTest: if( useLinksInCGNS ) then

             ! A zone must be created. Use the same name as
             ! in the original grid file.

             sizes(1) = cgnsDoms(zone)%il
             sizes(2) = cgnsDoms(zone)%jl
             sizes(3) = cgnsDoms(zone)%kl
             sizes(4) = cgnsDoms(zone)%nx
             sizes(5) = cgnsDoms(zone)%ny
             sizes(6) = cgnsDoms(zone)%nz
             sizes(7) = 0
             sizes(8) = 0
             sizes(9) = 0

             call cg_zone_write_f(cgnsInd, cgnsBase,              &
                                  cgnsDoms(zone)%zonename, sizes, &
                                  Structured, cgnsZone, ierr)
             if(ierr /= all_ok)                   &
               call terminate("writeSolCGNSZone", &
                              "Something wrong when calling &
                              &cg_zone_write_f")

             ! Check if a link must actually be created.

             writeLinkTest: if( writeLink ) then

               ! Create the link of the coordinates to the zone in the
               ! original grid.  First move to the correct location.

               call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                              cgnsZone, "end")
               if(ierr /= all_ok)                   &
                 call terminate("writeSolCGNSZone", &
                                "Something wrong when calling cg_goto_f")

               ! Determine the link name and write the link to the grid
               ! coordinates.

               linkName = cgnsBasename//"&
                          &/"//cgnsDoms(zone)%zonename//"&
                          &/"//"GridCoordinates"

               call cg_link_write_f("GridCoordinates", &
                                    gridFileNames(ind), linkName, ierr)
               if(ierr /= all_ok)                   &
                 call terminate("writeSolCGNSZone", &
                                "Something wrong when calling &
                                &cg_link_write_f")

               ! If there are overset grids present, then also write a
               ! link to the solution node containing the nodal iblanks.

               if (oversetPresent) then

                 linkName = cgnsBasename//"&
                            &/"//cgnsDoms(zone)%zonename//"&
                            &/"//"Nodal Blanks"

                 call cg_link_write_f("Nodal Blanks", &
                                      gridFileNames(ind), linkName, ierr)
                 if(ierr /= all_ok)                   &
                   call terminate("writeSolCGNSZone", &
                                  "Something wrong when calling &
                                  &cg_link_write_f")

               end if
             endif writeLinkTest

           else createZoneTest

             ! The zone already exists. Simply set cgnsZone to zone.

             cgnsZone = zone

           endif createZoneTest

           ! Create the flow solution node.

           call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone,          &
                               "Flow solution", CellCenter, cgnsSol, &
                               ierr)
           if(ierr /= all_ok)                   &
             call terminate("writeSolCGNSZone", &
                            "Something wrong when calling &
                            &cg_sol_write_f")

           ! Create the rind layers. If rind layers must be stored put
           ! 1 layer on every face of the block; otherwise put 0 layers.
           ! Use sizes as a buffer to store the rind data. The rind data
           ! must be created under the just created solution node.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                          cgnsZone, "FlowSolution_t", cgnsSol, "end")
           if(ierr /= all_ok)                   &
             call terminate("writeSolCGNSZone", &
                            "Something wrong when calling cg_goto_f")

           if( rindLayerThisSol ) then
             sizes(1) = 1; sizes(2) = 1; sizes(3) = 1
             sizes(4) = 1; sizes(5) = 1; sizes(6) = 1
           else
             sizes(1) = 0; sizes(2) = 0; sizes(3) = 0
             sizes(4) = 0; sizes(5) = 0; sizes(6) = 0
           endif

           call cg_rind_write_f(sizes, ierr)
           if(ierr /= all_ok)                   &
             call terminate("writeSolCGNSZone", &
                            "Something wrong when calling &
                            &cg_rind_write_f")

           ! Determine the index range of the solution of the zone
           ! to be written.

           iBegCGNS = 2 - sizes(1)
           jBegCGNS = 2 - sizes(3)
           kBegCGNS = 2 - sizes(5)

           iEndCGNS = cgnsDoms(zone)%il + sizes(2)
           jEndCGNS = cgnsDoms(zone)%jl + sizes(4)
           kEndCGNS = cgnsDoms(zone)%kl + sizes(6)

           ! Determine the number of variables to be written.
           ! For the unsteady higher solutions only the variables
           ! needed for a restart are written.

           nVarWritten = nSolVar+nDiscrVar
           if( unsteadyHigherSol ) nVarWritten = nw

           ! Loop over the number of variables to be written.

           varWriteLoop: do nn=1,nVarWritten

             ! Copy solNames(nn) in solName for later purposes.
             ! Correct this value if the conservative variables
             ! for a consistent unsteady restart must be written.
             ! No need to correct the turbulent variables, because
             ! these names are okay.

             solName = solNames(nn)

             if( unsteadyHigherSol ) then
               nnVar = nn

               select case(nnVar)
                 case (irho)
                   solName = cgnsDensity

                 case (imx)
                   solName = cgnsMomx

                 case (imy)
                   solName = cgnsMomy

                 case (imz)
                   solName = cgnsMomz

                 case (irhoE)
                   solName = cgnsEnergy
               end select
             endif

             ! Loop over the number of subblocks stored on
             ! this processor.

             do mm=1,nBlocks

               ! Set the pointers to the local domain.

               if( unsteadyHigherSol ) then
                 call setPointers(blocksCGNSblock(mm+offset), &
                                  1_intType, 1_intType)
               else
                 call setPointers(blocksCGNSblock(mm+offset), &
                                  1_intType, ind)
               endif

               ! Determine the cell range I have to write. This depends
               ! whether or not halo's must be written.

               iBeg = 2; iEnd = il
               jBeg = 2; jEnd = jl
               kBeg = 2; kEnd = kl

               if(storeRindLayer .and. (.not. unsteadyHigherSol)) then
                 if(iBegor == 1) iBeg = 1
                 if(jBegor == 1) jBeg = 1
                 if(kBegor == 1) kBeg = 1

                 if(iEndor == cgnsDoms(zone)%il) iEnd = ie
                 if(jEndor == cgnsDoms(zone)%jl) jEnd = je
                 if(kEndor == cgnsDoms(zone)%kl) kEnd = ke
               endif

               ! Fill the buffer with the correct solution variable.
               ! The routine called depends on the situation.

               if( unsteadyHigherSol ) then

                 ! Unsteady higher solution. Write the old solution.

                 call storeOldSolInBuffer(buffer, ind, nn, iBeg, iEnd, &
                                          jBeg, jEnd, kBeg, kEnd)
               else

                 ! Standard solution must be written.
 
                 call storeSolInBuffer(buffer, .true., solName, &
                                       iBeg, iEnd, jBeg, jEnd,  &
                                       kBeg, kEnd)
               endif

               ! And store it in sol. The routine called depends on
               ! the desired precision.

               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(sol, buffer,        &
                                                   iBegCGNS, jBegCGNS, &
                                                   kBegCGNS, iEndCGNS, &
                                                   jEndCGNS, kEndCGNS, &
                                                   subRanges(1,1,mm))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(sol, buffer,        &
                                                   iBegCGNS, jBegCGNS, &
                                                   kBegCGNS, iEndCGNS, &
                                                   jEndCGNS, kEndCGNS, &
                                                   subRanges(1,1,mm))
               end select

             enddo

             ! Loop over the number of subblocks stored on
             ! other processors.

             do mm=(nBlocks+1),nSubblocks

               ! Receive the range of subblock mm.

               source = proc(mm)
               call mpi_recv(buffer, bufSize, sumb_real, source, &
                             source+1, SUmb_comm_world, status, ierr)

               ! And store it in sol.

               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(sol, buffer,        &
                                                   iBegCGNS, jBegCGNS, &
                                                   kBegCGNS, iEndCGNS, &
                                                   jEndCGNS, kEndCGNS, &
                                                   subRanges(1,1,mm))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(sol, buffer,        &
                                                   iBegCGNS, jBegCGNS, &
                                                   kBegCGNS, iEndCGNS, &
                                                   jEndCGNS, kEndCGNS, &
                                                   subRanges(1,1,mm))
               end select

             enddo

             ! Write the solution variable to file. Source is just used
             ! as a dummy variable and does not have a meaning.

             call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                                   cgnsSol, realTypeCGNS,       &
                                   solName, sol, source, ierr)

             if(ierr /= all_ok)                   &
               call terminate("writeSolCGNSZone", &
                              "Something wrong when calling &
                              &cg_field_write_f")
           enddo varWriteLoop

         enddo solLoopRoot

         ! Release some memory only allocated on the root processor.

         deallocate(sol, subRanges, proc, stat=ierr)
         if(ierr /= 0) call terminate("writeSolCGNSZone", &
                                      "Deallocation error on root proc")

       else rootproc
!
!        ****************************************************************
!        *                                                              *
!        * I am not the root processor and may have to send some data   *
!        * to the root processor.                                       *
!        *                                                              *
!        ****************************************************************
! 
         ! Determine the subranges for the 1st solution.

         rindLayerThisSol = storeRindLayer
         call getSubRangesSol

         ! Allocate the memory for buffer.

         allocate(buffer(bufSize), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("writeSolCGNSZone", &
                          "Memory allocation failure for buffer")

         ! Loop over the number of solutions.

         solLoopOthers: do ind=1,nVolSolToWrite

           ! Determine whether or not we are dealing with an unsteady
           ! higher solution here.

           unsteadyHigherSol = .false.
           if(ind > 1 .and. equationMode == unsteady) &
             unsteadyHigherSol = .true.

           ! Set rindLayerThisSol. A rind layer is not written for
           ! the higher solutions in unsteady mode.

           rindLayerThisSol = storeRindLayer
           if( unsteadyHigherSol ) rindLayerThisSol = .false.

           ! Check if the subranges must be recomputed. Some dirty stuff
           ! must be done, because Fortran 90/95 does not allow a
           ! comparison between logicals.

           nn = 0; if( storeRindLayer )   nn = 1
           mm = 0; if( rindLayerThisSol ) mm = 1
           if(nn /= mm) call getSubRangesSol

           ! Determine the number of variables to be written.
           ! For the unsteady higher solutions only the variables
           ! needed for a restart are written.

           nVarWritten = nSolVar+nDiscrVar
           if( unsteadyHigherSol ) nVarWritten = nw

           ! Loop over the number of variables to be written.

           do nn=1,nVarWritten

             ! Loop over the number of subblocks stored on
             ! this processor.

             do mm=1,nBlocks

               ! Set the pointers to the local domain.

               if( unsteadyHigherSol ) then
                 call setPointers(blocksCGNSblock(mm+offset), &
                                  1_intType, 1_intType)
               else
                 call setPointers(blocksCGNSblock(mm+offset), &
                                  1_intType, ind)
               endif

               ! Determine the cell range I have to write. This depends
               ! whether or not halo's must be written.

               iBeg = 2; iEnd = il
               jBeg = 2; jEnd = jl
               kBeg = 2; kEnd = kl

               if(storeRindLayer .and. (.not. unsteadyHigherSol)) then
                 if(iBegor == 1) iBeg = 1
                 if(jBegor == 1) jBeg = 1
                 if(kBegor == 1) kBeg = 1

                 if(iEndor == cgnsDoms(zone)%il) iEnd = ie
                 if(jEndor == cgnsDoms(zone)%jl) jEnd = je
                 if(kEndor == cgnsDoms(zone)%kl) kEnd = ke
               endif

               ! Fill the buffer with the correct solution variable.
               ! The routine called depends on the situation.

               if( unsteadyHigherSol ) then

                 ! Unsteady higher solution. Write the old solution.

                 call storeOldSolInBuffer(buffer, ind, nn, iBeg, iEnd, &
                                          jBeg, jEnd, kBeg, kEnd)
               else

                 ! Standard solution must be written.
 
                 call storeSolInBuffer(buffer, .true., solNames(nn), &
                                       iBeg, iEnd, jBeg, jEnd,       &
                                       kBeg, kEnd)
               endif

               ! And send it to processor 0.

               ll   = (iEnd-iBeg+1)*(jEnd-jBeg+1)*(kEnd-kBeg+1)
               size = ll

               call mpi_send(buffer, size, sumb_real, 0, myID+1, &
                             SUmb_comm_world, ierr)
             enddo
           enddo

         enddo solLoopOthers
       endif rootproc

       ! Release some memory.

       deallocate(buffer, stat=ierr)
       if(ierr /= 0)                        &
         call terminate("writeSolCGNSZone", &
                        "Deallocation error for buffer")

       !=================================================================

       contains

         !===============================================================

         subroutine getSubRangesSol
!
!        ****************************************************************
!        *                                                              *
!        * getSubRangesSol determines the subranges of the              *
!        * computational blocks that contribute to the CGNS block which *
!        * is currently written. Also the size of the largest subblock  *
!        * is determined.                                               *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local variables.
!
         integer :: source

         integer(kind=intType) :: i, j, ll
         integer(kind=intType), dimension(6) :: ii
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Initialize bufSize.

         bufSize = 0

         ! Test if I'm the root processor or not.

         testRoot: if(myID == 0) then

           ! I'm the root processor. Determine the subRanges of the 
           ! subblocks stored on locally. Note that nBlocks can be 0.

           do i=1,nBlocks

             ! Store the local block ID a bit easier in j.

             j = blocksCGNSblock(i+offset)

             ! Determine the range; this is the same for all spectral
             ! solutions, so the first one can be used.

             subRanges(1,1,i) = flowDoms(j,1,1)%iBegor + 1
             subRanges(1,2,i) = flowDoms(j,1,1)%iEndor

             subRanges(2,1,i) = flowDoms(j,1,1)%jBegor + 1
             subRanges(2,2,i) = flowDoms(j,1,1)%jEndor

             subRanges(3,1,i) = flowDoms(j,1,1)%kBegor + 1
             subRanges(3,2,i) = flowDoms(j,1,1)%kEndor

             ! Correct in case rind layers must be stored.

             if( rindLayerThisSol ) then

               if(subRanges(1,1,i) == 2) subRanges(1,1,i) = 1
               if(subRanges(2,1,i) == 2) subRanges(2,1,i) = 1
               if(subRanges(3,1,i) == 2) subRanges(3,1,i) = 1

               if(subRanges(1,2,i) == cgnsDoms(zone)%il) &
                  subRanges(1,2,i) =  cgnsDoms(zone)%il + 1
               if(subRanges(2,2,i) == cgnsDoms(zone)%jl) &
                  subRanges(2,2,i) =  cgnsDoms(zone)%jl + 1
               if(subRanges(3,2,i) == cgnsDoms(zone)%kl) &
                  subRanges(3,2,i) =  cgnsDoms(zone)%kl + 1

             endif

           enddo

           ! The rest of the block ranges must be obtained by
           ! communication.

           do i=(nBlocks+1),nSubblocks

             ! Receive the range of subblock i.

             source = proc(i)
             call mpi_recv(ii, 6, sumb_integer, source, source, &
                           SUmb_comm_world, status, ierr)

             subRanges(1,1,i) = ii(1)
             subRanges(1,2,i) = ii(2)
             subRanges(2,1,i) = ii(3)
             subRanges(2,2,i) = ii(4)
             subRanges(3,1,i) = ii(5)
             subRanges(3,2,i) = ii(6)
           enddo

           ! Determine the size of the largest subblock.

           do i=1,nSubBlocks
             ll = (subRanges(1,2,i) - subRanges(1,1,i) + 1) &
                * (subRanges(2,2,i) - subRanges(2,1,i) + 1) &
                * (subRanges(3,2,i) - subRanges(3,1,i) + 1)
             bufSize = max(bufSize, ll)
           enddo

         else testRoot

           ! Loop over the number of subblocks stored on this processor.

           do i=1,nBlocks

             ! Store the local block id a bit easier in j.

             j = blocksCGNSblock(i+offset)

             ! Copy the range of this subblock into the buffer ii.
             ! This is the same for all spectral solutions, so the
             ! first one can be used.

             ii(1) = flowDoms(j,1,1)%iBegor + 1
             ii(2) = flowDoms(j,1,1)%iEndor
             ii(3) = flowDoms(j,1,1)%jBegor + 1
             ii(4) = flowDoms(j,1,1)%jEndor
             ii(5) = flowDoms(j,1,1)%kBegor + 1
             ii(6) = flowDoms(j,1,1)%kEndor

             ! Correct in case rind layers must be stored.

             if( rindLayerThisSol ) then

               if(ii(1) == 2) ii(1) = 1
               if(ii(2) == cgnsDoms(zone)%il) ii(2) = ii(2) + 1
               if(ii(3) == 2) ii(3) = 1
               if(ii(4) == cgnsDoms(zone)%jl) ii(4) = ii(4) + 1
               if(ii(5) == 2) ii(5) = 1
               if(ii(6) == cgnsDoms(zone)%kl) ii(6) = ii(6) + 1

             endif

             ! Send the buffer to processor 0.

             call mpi_send(ii, 6, sumb_integer, 0, myID, &
                           SUmb_comm_world, ierr)

             ! Check the size of this subblock and update bufSize
             ! if needed.

             ll = (ii(2) - ii(1) + 1) * (ii(4) - ii(3) + 1) &
                * (ii(6) - ii(5) + 1)
             bufSize = max(bufSize, ll)

           enddo

         endif testRoot

         end subroutine getSubRangesSol
#endif

       end subroutine writeSolCGNSZone
