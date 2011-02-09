!
!      ******************************************************************
!      *                                                                *
!      * File:          writeSurfsolCGNSZone.F90                        *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      * Starting date: 05-15-2003                                      *
!      * Last modified: 11-04-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writeSurfsolCGNSZone(zone, nBlocks, subface, nSolVar, &
                                       solNames, nZonesWritten, periodic)
!
!      ******************************************************************
!      *                                                                *
!      * writeSurfsolCGNSZone writes a surface solution of the given    *
!      * zone (block) and boundary subface to the cgns surface file(s). *
!      * A distinction must be made between true boundaries and         *
!      * periodic boundaries; the latter are a special kind of internal *
!      * block boundaries. This is indicated by the logical periodic.   *
!      *                                                                *
!      ******************************************************************
!
       use block
       use BCTypes
       use cgnsGrid
       use communication
       use inputIO
       use su_cgns
       use outputMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: zone, nBlocks
       integer(kind=intType), intent(in)    :: subface, nSolVar
       integer(kind=intType), intent(inout) :: nZonesWritten

       character(len=*), dimension(*), intent(in) :: solNames

       logical, intent(in) :: periodic

#ifdef USE_NO_CGNS
       call terminate("writeSurfsolCGNSZone", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       integer :: ierr
       integer :: source, size
       integer :: cgnsBase, cgnsZone, cgnsSol, cgnsInd

       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, offset
       integer(kind=intType) :: mm, mBlocks, faceID, nSubfaces
       integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
       integer(kind=intType) :: il, jl, ind

       integer(kind=intType), dimension(nProc)       :: nMessages
       integer(kind=intType), dimension(3,2,nBlocks) :: nodalRange
       integer(kind=intType), dimension(3,2,nBlocks) :: cellRange

       integer(kind=intType), dimension(:,:,:), allocatable :: rangeNode
       integer(kind=intType), dimension(:,:,:), allocatable :: rangeCell

       real(kind=realType), dimension(:), allocatable :: buffer

       logical :: iOverlap, jOverlap, kOverlap
       logical :: viscousSubface
       logical, dimension(nBlocks) :: contributeToFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the offset in blocksCGNSblock for this zone in offset.

       offset  = nBlocksCGNSblock(zone-1)

       ! Determine the range of this subface and whether or not this
       ! is a viscous subface. For a periodic boundary this info is
       ! retrieved from the internal block boundary; for all others from
       ! the physical boundary.

       if( periodic ) then

         ! Periodic boundary. Set viscousSubface to false.

         viscousSubface = .false.

         ! Store the nodal range of the cgns subface a bit easier.
         ! Make sure that iBeg, jBeg and kBeg contain the lowest values
         ! and iEnd, jEnd and kEnd the highest.

         iBeg = min(cgnsDoms(zone)%conn1to1(subface)%iBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%iEnd)
         jBeg = min(cgnsDoms(zone)%conn1to1(subface)%jBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%jEnd)
         kBeg = min(cgnsDoms(zone)%conn1to1(subface)%kBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%kEnd)
         iEnd = max(cgnsDoms(zone)%conn1to1(subface)%iBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%iEnd)
         jEnd = max(cgnsDoms(zone)%conn1to1(subface)%jBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%jEnd)
         kEnd = max(cgnsDoms(zone)%conn1to1(subface)%kBeg, &
                    cgnsDoms(zone)%conn1to1(subface)%kEnd)

       else

         ! True physical boundary.
         ! If this is an extrapolation boundary (usually singular line),
         ! return. You don't want that info in the solution file.

         if(cgnsDoms(zone)%bocoInfo(subface)%BCType == Extrap) return

         ! Store the nodal range of the cgns subface a bit easier.
         ! Make sure that iBeg, jBeg and kBeg contain the lowest values
         ! and iEnd, jEnd and kEnd the highest.

         iBeg = min(cgnsDoms(zone)%bocoInfo(subface)%iBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%iEnd)
         jBeg = min(cgnsDoms(zone)%bocoInfo(subface)%jBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%jEnd)
         kBeg = min(cgnsDoms(zone)%bocoInfo(subface)%kBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%kEnd)
         iEnd = max(cgnsDoms(zone)%bocoInfo(subface)%iBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%iEnd)
         jEnd = max(cgnsDoms(zone)%bocoInfo(subface)%jBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%jEnd)
         kEnd = max(cgnsDoms(zone)%bocoInfo(subface)%kBeg, &
                    cgnsDoms(zone)%bocoInfo(subface)%kEnd)

         ! Determine whether or not this is a viscous subface.

         viscousSubface = .false.
         if(cgnsDoms(zone)%bocoInfo(subface)%BCType == &
                     NSWallAdiabatic          .or.     &
            cgnsDoms(zone)%bocoInfo(subface)%BCType == &
                     NSWallIsothermal) viscousSubface = .true.

       endif

       ! Update nZonesWritten.

       nZonesWritten = nZonesWritten + 1

       ! Determine the face ID on which the given cgns subface is located.

       if(iBeg == iEnd) then
         faceID = iMax
         if(iBeg == 1) faceID = iMin
       else if(jBeg == jEnd) then
         faceID = jMax
         if(jBeg == 1) faceID = jMin
       else
         faceID = kMax
         if(kBeg == 1) faceID = kMin
       endif

       ! Determine the number of nodes in the two coordinate directions.
       ! These are called il and jl.

       select case (faceID)
         case (iMin,iMax)
           il = jEnd - jBeg + 1
           jl = kEnd - kBeg + 1
         case (jMin,jMax)
           il = iEnd - iBeg + 1
           jl = kEnd - kBeg + 1
         case (kMin,kMax)
           il = iEnd - iBeg + 1
           jl = jEnd - jBeg + 1
       end select

       ! Allocate the memory for buffer, which is used to communicate
       ! the coordinates and solution. Assume that rind layers are
       ! present, such that the solution uses most memory.

       size = (il+1)*(jl+1)
       allocate(buffer(size), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writeSurfsolCGNSZone", &
                        "Memory allocation failure for buffer")

       ! Determine the number of local blocks that actually share the
       ! subface of the original cgns block. Note that nBlocks and
       ! blocksCGNSblock contain information of the entire cgns block,
       ! but not of the subface.

       mBlocks = 0
       do i=1,nBlocks

         ! Store the current local block ID a bit easier.

         mm = blocksCGNSblock(i+offset)

         ! Determine whether or not the cgns subface is (partially)
         ! part of the subblock mm. Initialize the overlaps to .false.

         iOverlap = .false.
         jOverlap = .false.
         kOverlap = .false.

         ! First check the face ID.

         select case (faceID)
           case (iMin)
             if(flowDoms(mm,1,1)%iBegor == iBeg) iOverlap = .true.
           case (iMax)
             if(flowDoms(mm,1,1)%iEndor == iEnd) iOverlap = .true.
           case (jMin)
             if(flowDoms(mm,1,1)%jBegor == jBeg) jOverlap = .true.
           case (jMax)
             if(flowDoms(mm,1,1)%jEndor == jEnd) jOverlap = .true.
           case (kMin)
             if(flowDoms(mm,1,1)%kBegor == kBeg) kOverlap = .true.
           case (kMax)
             if(flowDoms(mm,1,1)%kEndor == kEnd) kOverlap = .true.
         end select

         ! Check the overlap for the other two directions.

         if(iBeg < flowDoms(mm,1,1)%iEndor .and. &
            iEnd > flowDoms(mm,1,1)%iBegor) iOverlap = .true.
         if(jBeg < flowDoms(mm,1,1)%jEndor .and. &
            jEnd > flowDoms(mm,1,1)%jBegor) jOverlap = .true.
         if(kBeg < flowDoms(mm,1,1)%kEndor .and. &
            kEnd > flowDoms(mm,1,1)%kBegor) kOverlap = .true.

         ! If all three directions overlap, this subblock contributes
         ! to the current cgns subface.

         if(iOverlap .and. jOverlap .and. kOverlap) then
           contributeToFace(i) = .true.
           mBlocks = mBlocks +1

           ! Determine the nodal and cell subrange for this subface.

           call determineSubranges

         else
           contributeToFace(i) = .false.
         endif

       enddo

       ! Determine the amount of surface parts each processor will send
       ! to processor 0. The result needs only to be known on
       ! processor 0.

       call mpi_gather(mBlocks, 1, sumb_integer, nMessages, 1, &
                       sumb_integer, 0, SUmb_comm_world, ierr)

       ! At the moment the writing of the cgns file is sequential and done
       ! by processor 0. This means that this processor gathers all info
       ! from the other processors and writes it to file.

       rootproc: if(myID == 0) then

         ! I am processor 0 and poor me has to do all the work.

         ! First determine the number of subfaces into the original cgns
         ! subface is split.

         nSubfaces = 0
         do i=1,nProc
           nSubfaces = nSubfaces + nMessages(i)
         enddo

         ! Allocate the memory for the nodal and cell ranges for each
         ! of the contributing subfaces.

         allocate(rangeNode(3,2,nSubfaces), rangeCell(3,2,nSubfaces), &
                  stat=ierr)
         if(ierr /= 0)                            &
           call terminate("writeSurfsolCGNSZone", &
                          "Memory allocation failure for &
                          &rangeNode, etc")

         ! Store the nodal and cell subranges of all contributions in
         ! rangeNode and rangeCell. Start with my own contributions.
         ! Note that mBlocks could be 0.

         do i=1,mBlocks

           ! Copy the local nodal and cell ranges.

           rangeNode(1,1,i) = nodalRange(1,1,i)
           rangeNode(2,1,i) = nodalRange(2,1,i)
           rangeNode(3,1,i) = nodalRange(3,1,i)

           rangeNode(1,2,i) = nodalRange(1,2,i)
           rangeNode(2,2,i) = nodalRange(2,2,i)
           rangeNode(3,2,i) = nodalRange(3,2,i)

           rangeCell(1,1,i) = cellRange(1,1,i)
           rangeCell(2,1,i) = cellRange(2,1,i)
           rangeCell(3,1,i) = cellRange(3,1,i)

           rangeCell(1,2,i) = cellRange(1,2,i)
           rangeCell(2,2,i) = cellRange(2,2,i)
           rangeCell(3,2,i) = cellRange(3,2,i)

         enddo

         ! The rest of the ranges must be obtained by communication.

         mm = mBlocks + 1
         do i=2,nProc

           ! Check if something must be received from this processor.

           if(nMessages(i) > 0) then

             ! Store the source and size of the messages and receive
             ! the messages. Note that 1 must be substracted from i
             ! to obtain the correct processor id.

             source = i -1
             size   = 6*nMessages(i)

             call mpi_recv(rangeNode(1,1,mm), size, sumb_integer,   &
                           source, source, SUmb_comm_world, status, &
                           ierr)
             call mpi_recv(rangeCell(1,1,mm), size, sumb_integer,     &
                           source, source+1, SUmb_comm_world, status, &
                           ierr)

             ! Update mm.

             mm = mm + nMessages(i)

           endif
         enddo

         ! Loop over the number of solutions to be written.

         solLoopRoot: do ind=1,nSurfSolToWrite

           ! Store the file and base ID a bit easier.

           cgnsInd  = fileIDs(ind)
           cgnsBase = cgnsBases(ind)

           ! Create the surface zone.

           call createSurfaceZone

           ! Write the nodal coordinates.

           call writeSurfaceCoord

           ! Write the cell centered surface solution.

           call writeSurfaceSol

         enddo solLoopRoot

         ! Release the memory of the variables only
         ! processor 0 allocates.

         deallocate(rangeNode, rangeCell, stat=ierr)
         if(ierr /= 0)                               &
           call terminate("writeSurfsolCGNSZone", &
                          "Deallocation error for rangeNode, etc")

       else rootproc

         ! Send the node and cell ranges to processor 0 if a block
         ! contributes to the current cgns subface.

         if(mBlocks > 0) then

           ! Determine the size of the messages and send the nodal
           ! and cell ranges to processor 0.

           size = 6*mBlocks
           call mpi_send(nodalRange, size, sumb_integer, 0, myID, &
                         SUmb_comm_world, ierr)
           call mpi_send(cellRange, size, sumb_integer, 0, myID+1, &
                         SUmb_comm_world, ierr)
         endif

         ! Loop over the number of solutions to be written.

         solLoopOthers: do ind=1,nSurfSolToWrite

           ! Write the nodal coordinates.

           call writeSurfaceCoord

           ! Write the cell centered surface solution.

           call writeSurfaceSol

         enddo solLoopOthers

       endif rootproc

       ! Release the memory of buffer.

       deallocate(buffer, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("writeSurfsolCGNSZone", &
                        "Deallocation error for buffer")

       contains

!        ================================================================

         subroutine determineSubranges
!
!        ****************************************************************
!        *                                                              *
!        * determineSubranges determines the nodal and cell subrange    *
!        * for the given local block ID mm in the current cgns subface. *
!        *                                                              *
!        ****************************************************************
!
         use inputIO
         implicit none
!
!        Local variable
!
         integer(kind=intType) :: ii
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Store mBlocks, the current number of local blocks that
         ! participate to the cgns subface, a bit easier in ii.

         ii = mBlocks

         ! Determine the nodal range of the current subface. Note that
         ! in case multiple blocks contribute to the cgns subface, the
         ! nodes on the interface are stored on both partitions for the
         ! moment. This is corrected later.

         nodalRange(1,1,ii) = max(iBeg,flowDoms(mm,1,1)%iBegor)
         nodalRange(2,1,ii) = max(jBeg,flowDoms(mm,1,1)%jBegor)
         nodalRange(3,1,ii) = max(kBeg,flowDoms(mm,1,1)%kBegor)

         nodalRange(1,2,ii) = min(iEnd,flowDoms(mm,1,1)%iEndor)
         nodalRange(2,2,ii) = min(jEnd,flowDoms(mm,1,1)%jEndor)
         nodalRange(3,2,ii) = min(kEnd,flowDoms(mm,1,1)%kEndor)

         ! The cell range. Step 1, the interior.

         cellRange(1,1,ii) = nodalRange(1,1,ii) +1
         cellRange(2,1,ii) = nodalRange(2,1,ii) +1
         cellRange(3,1,ii) = nodalRange(3,1,ii) +1

         cellRange(1,2,ii) = nodalRange(1,2,ii)
         cellRange(2,2,ii) = nodalRange(2,2,ii)
         cellRange(3,2,ii) = nodalRange(3,2,ii)

         ! Step 2. Correct for possible rind layers.

         if( storeRindLayer ) then

           if(nodalRange(1,1,ii) == iBeg) cellRange(1,1,ii) = iBeg
           if(nodalRange(2,1,ii) == jBeg) cellRange(2,1,ii) = jBeg
           if(nodalRange(3,1,ii) == kBeg) cellRange(3,1,ii) = kBeg

           if(nodalRange(1,2,ii) == iEnd) cellRange(1,2,ii) = iEnd +1
           if(nodalRange(2,2,ii) == jEnd) cellRange(2,2,ii) = jEnd +1
           if(nodalRange(3,2,ii) == kEnd) cellRange(3,2,ii) = kEnd +1

         endif

         ! Step 3. Correct for the face ID.

         select case (faceID)
           case (iMin)
             cellRange(1,1,ii) = 2
             cellRange(1,2,ii) = 2
           case (iMax)
             cellRange(1,1,ii) = iEnd
             cellRange(1,2,ii) = iEnd
           case (jMin)
             cellRange(2,1,ii) = 2
             cellRange(2,2,ii) = 2
           case (jMax)
             cellRange(2,1,ii) = jEnd
             cellRange(2,2,ii) = jEnd
           case (kMin)
             cellRange(3,1,ii) = 2
             cellRange(3,2,ii) = 2
           case (kMax)
             cellRange(3,1,ii) = kEnd
             cellRange(3,2,ii) = kEnd
         end select

         ! Correct the nodal range for possible overlap.

         if(nodalRange(1,1,ii) > iBeg) &
            nodalRange(1,1,ii) = nodalRange(1,1,ii) +1
         if(nodalRange(2,1,ii) > jBeg) &
            nodalRange(2,1,ii) = nodalRange(2,1,ii) +1
         if(nodalRange(3,1,ii) > kBeg) &
            nodalRange(3,1,ii) = nodalRange(3,1,ii) +1

         end subroutine determineSubranges

!        ================================================================

         subroutine createSurfaceZone
!
!        ****************************************************************
!        *                                                              *
!        * createSurfaceZone creates a surface node in the given        *
!        * cgns surface solution file. This routine should only be      *
!        * called by processor 0.                                       *
!        *                                                              *
!        ****************************************************************
!
         use inputIO
         implicit none
!
!        Local variables.
!
         integer, dimension(6) :: sizes

         integer(kind=intType) :: nn

         character(len=maxCGNSNameLen) :: zonename
         character(len=7)              :: integerString
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Determine the sizes of the subface.

         sizes(1) = il
         sizes(2) = jl
         sizes(3) = il -1
         sizes(4) = jl -1
         sizes(5) = 0
         sizes(6) = 0

         ! For all zones a number is added to make to zone name unique.
         ! Create that string here.

         write(integerString,"(i6)") nZonesWritten
         integerString = adjustl(integerString)

         ! Create the zone name. A distinction must be made between
         ! periodic and physical boundaries.

         if( periodic ) then

           zonename = "PeriodicBCZone"//trim(integerString)

         else

           ! True physcical boundary. A distinction is made between zones
           ! that do and don't belong to a family. The basename of the
           ! former boundaries is the family name, such that the entire
           ! family can be easily selected in postprocessing software.

           nn = cgnsDoms(zone)%bocoInfo(subface)%familyID
           if(nn > 0) then

             ! Zone belongs to a family. Add the zone number to the
             ! family name for the zone name.

             zonename = trim(cgnsFamilies(nn)%familyName)//   &
                        trim(integerString)

           else

             ! Zone does not belong to a family. The first part of the
             ! zone name depends on the boundary condition of the cgns
             ! subface.

             select case (cgnsDoms(zone)%bocoInfo(subface)%BCType)
               case (Symm)
                 zonename = "Symmetry"
               case (SymmPolar)
                 zonename = "SymmetryPolar"
               case (NSWallAdiabatic)
                 zonename = "NSWallAdiabatic"
               case (NSWallIsothermal)
                 zonename = "NSWallIsothermal"
               case (EulerWall)
                 zonename = "EulerWall"
               case (FarField)
                 zonename = "FarField"
               case (SupersonicInflow)
                 zonename = "InflowSupersonic"
               case (SubsonicInflow)
                 zonename = "InflowSubsonic"
               case (SupersonicOutflow)
                 zonename = "OutflowSupersonic"
               case (SubsonicOutflow)
                 zonename = "OutflowSubsonic"
               case (MassBleedInflow)
                 zonename = "MassBleedInflow"
               case (MassBleedOutflow)
                 zonename = "MassBleedOutflow"
               case (mDot)
                 zonename = "MDot"
               case (Thrust)
                 zonename = "Thrust"
               case (SlidingInterface)
                 zonename = "Sliding"
               case (OversetOuterBound)
                 zonename = "Overlap"
               case (DomainInterfaceAll)
                 zonename = "DomainAll"
               case (DomainInterfaceRhoUVW)
                 zonename = "DomainRhoUVW"
               case (DomainInterfaceP)
                 zonename = "DomainP"
               case (DomainInterfaceRho)
                 zonename = "DomainRho"
               case (DomainInterfaceTotal)
                 zonename = "DomainTotal"
               case default
                 call terminate("createSurfaceZone", &
                                "Unknown boundary condition")
             end select

             ! Add the number to zone name to make it unique.

             zonename = trim(zoneName)//"BCZone"//trim(integerString)

           endif

         endif

         ! Create the 2D structured zone.

         call cg_zone_write_f(cgnsInd, cgnsBase, zonename, sizes, &
                              Structured, cgnsZone, ierr)
         if(ierr /= all_ok)                    &
           call terminate("createSurfaceZone", &
                          "Something wrong when calling cg_zone_write_f")

         ! Create the flow solution node.

         call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, &
                             "Flow solution", CellCenter, cgnsSol, ierr)
         if(ierr /= all_ok)                    &
           call terminate("createSurfaceZone", &
                          "Something wrong when calling cg_sol_write_f")

         ! Create the rind layers. If rind layers must be stored put
         ! 1 layer on every side of the subface; otherwise put 0 layers.
         ! Use sizes as a buffer to store the rind data. The rind data
         ! must be created under the just created solution node.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
                        cgnsZone, "FlowSolution_t", cgnsSol, "end")
         if(ierr /= all_ok)                    &
           call terminate("createSurfaceZone", &
                          "Something wrong when calling cg_goto_f")

         if( storeRindLayer ) then
           sizes(1) = 1; sizes(2) = 1
           sizes(3) = 1; sizes(4) = 1
         else
           sizes(1) = 0; sizes(2) = 0
           sizes(3) = 0; sizes(4) = 0
         endif

         call cg_rind_write_f(sizes, ierr)
         if(ierr /= all_ok)                    &
           call terminate("createSurfaceZone", &
                          "Something wrong when calling cg_rind_write_f")

         end subroutine createSurfaceZone

!        ================================================================

         subroutine writeSurfaceCoord
!
!        ****************************************************************
!        *                                                              *
!        * WriteSurfaceCoord write the vertex values of the             *
!        * coordinates to the given zone of the cgns surface solution   *
!        * file.                                                        *
!        *                                                              *
!        ****************************************************************
!
         use cgnsNames
         use inputIO
         implicit none
!
!        Local variables.
!
         integer :: realTypeCGNS

         integer(kind=intType) :: i, j, k, kk, ll, mm, ii, jj
         integer(kind=intType) :: lk, lj, li
         integer(kind=intType) :: sizeCGNSWriteType

         real(kind=realType) :: LRefInv

         character, dimension(:), allocatable :: writeBuffer
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
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

         ! Compute the multiplication factor to obtain the original
         ! coordinates. Note that LRef is corrected to 1.0 when the
         ! coordinates should be written in meters. This happens when
         ! the grid is read.

         LRefInv = one/cgnsDoms(zone)%LRef

         ! Processor 0 does the writing and must therefore allocate the
         ! writeBuffer.

         if(myID == 0) then
           mm = (kEnd-kBeg+1) * (jEnd-jBeg+1) * (iEnd-iBeg+1) &
              * sizeCGNSWriteType
           allocate(writeBuffer(mm), stat=ierr)
           if(ierr /= 0)                         &
             call terminate("writeSurfaceCoord", &
                            "Memory allocation failure for writeBuffer")
         endif

         ! Loop over the three coordinates.

         coorLoop: do mm=1,3

           ! Loop over the number of blocks stored on this processor
           ! which may contribute to the subface. Note that
           ! blocksCGNSblock contain all the subblocks part of cgns
           ! block, which is not the same. BlocksCGNSblock cannot be
           ! changed, because it is needed for other subfaces.

           kk = 0
           jj = 0
           do ll=1,nBlocks

             ! Test if the current local block contributes to the
             ! cgns subface to be written.

             if( contributeToFace(ll) ) then

               ! Update the counter kk and store the local block id in ii.

               kk = kk+1
               ii = blocksCGNSblock(ll+offset)

               ! Store the coordinate for this contribution in buffer.
               ! As nodalRange contains the nodal ranges in the
               ! original cgns block, the starting value of the current
               ! subblock must be substracted. These are the actual
               ! indices in the local block and are stored in
               ! lk, lj and li.

               do k=nodalRange(3,1,kk),nodalRange(3,2,kk)
                 lk = k - flowDoms(ii,1,ind)%kBegor + 1
                 do j=nodalRange(2,1,kk),nodalRange(2,2,kk)
                   lj = j - flowDoms(ii,1,ind)%jBegor + 1
                   do i=nodalRange(1,1,kk),nodalRange(1,2,kk)
                     li = i - flowDoms(ii,1,ind)%iBegor + 1

                     ! Update the counter jj and store the coordinate
                     ! in buffer.

                     jj = jj+1
                     buffer(jj) = LRefInv &
                                * flowDoms(ii,1,ind)%x(li,lj,lk,mm)
                   enddo
                 enddo
               enddo

             endif
           enddo

           ! Make a distinction between processor 0 and the other
           ! processors. Processor 0 does all the writing.

           rootproc: if(myID == 0) then

             ! I am processor 0 and must do the writing.

             ! Loop over the other processors and receive possible
             ! messages.

             do ll=2,nProc

               ! Test if processor ll contributes to this subface.

               if(nMessages(ll) > 0) then

                 ! Receive the message. Note that size is an upper limit.
                 ! Furthermore 1 must be substracted from ll to obtain
                 ! the correct processor ID; the processor ID's start
                 ! at 0.

                 size   = il*jl - jj
                 source = ll -1
                 call mpi_recv(buffer(jj+1), size, sumb_real, source, &
                               source, SUmb_comm_world, status, ierr)

                 ! Determine the true size of the message and update
                 ! the counter jj accordingly.

                 call mpi_get_count(status, sumb_real, size, ierr)
                 jj = jj + size
               endif
             enddo

             ! Copy the coordinate from buffer into the correct place
             ! in writeBuffer. The routine called depends on the
             ! desired precision.

             ii = 1
             do kk=1,nSubfaces
               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(writeBuffer,      &
                                                   buffer(ii),       &
                                                   iBeg, jBeg, kBeg, &
                                                   iEnd, jEnd, kEnd, &
                                                   rangeNode(1,1,kk))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(writeBuffer,      &
                                                   buffer(ii),       &
                                                   iBeg, jBeg, kBeg, &
                                                   iEnd, jEnd, kEnd, &
                                                   rangeNode(1,1,kk))
               end select

               ! Update the counter ii for the next subface.

               ii = ii + (rangeNode(1,2,kk) - rangeNode(1,1,kk) + 1) &
                  *      (rangeNode(2,2,kk) - rangeNode(2,1,kk) + 1) &
                  *      (rangeNode(3,2,kk) - rangeNode(3,1,kk) + 1)
             enddo

             ! Write the coordinates, depending on the situation.
             ! In source the actual number is stored; normally this
             ! is equal to mm.

             select case (mm)
               case (1_intType)
                 call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                                       realTypeCGNS, cgnsCoorX,     &
                                       writeBuffer, source, ierr)
               case (2_intType)
                 call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                                       realTypeCGNS, cgnsCoorY,     &
                                       writeBuffer, source, ierr)
               case (3_intType)
                 call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                                       realTypeCGNS, cgnsCoorZ,     &
                                       writeBuffer, source, ierr)
             end select

             if(ierr /= all_ok)                    &
               call terminate("writeSurfaceCoord", &
                              "Something wrong when calling &
                              &cg_coord_write_f")

             ! Write the units, if possible.

             if( cgnsDoms(zone)%gridUnitsSpecified ) then

               ! Go to the correct place in the surface solution file.

               call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                              "Zone_t", cgnsZone,      &
                              "GridCoordinates_t", 1,  &
                              "DataArray_t", source, "end")
               if(ierr /= all_ok)                      &
                 call terminate("writeSurfaceCoord", &
                                "Something wrong when calling cg_goto_f")

               ! Write the units.

               call cg_units_write_f(cgnsDoms(zone)%mass, &
                                     cgnsDoms(zone)%len,  &
                                     cgnsDoms(zone)%time, &
                                     cgnsDoms(zone)%temp, &
                                     cgnsDoms(zone)%angle, ierr)
               if(ierr /= all_ok)                    &
                 call terminate("writeSurfaceCoord", &
                                "Something wrong when calling &
                                &cg_units_write_f")
             endif

           else rootproc

             ! Not the root processor.
             ! Data must be sent to processor 0 if local blocks
             ! contribute to the current cgns subface.

             if( jj > 0 ) then
               size = jj
               call mpi_send(buffer, size, sumb_real, 0, myID, &
                             SUmb_comm_world, ierr)
             endif

           endif rootproc

         enddo coorLoop

         ! Processor 0 must deallocate the writeBuffer.

         if(myID == 0) then
           deallocate(writeBuffer, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("writeSurfaceCoord", &
                            "Deallocation error for writeBuffer")
         endif

         end subroutine writeSurfaceCoord

!        ================================================================

         subroutine writeSurfaceSol
!
!        ****************************************************************
!        *                                                              *
!        * writeSurfaceSol writes the cell centered surface solution    *
!        * to the cgns surface file.                                    *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local variables.
!
         integer :: realTypeCGNS

         integer(kind=intType) :: ii, jj, kk, ll, mm
         integer(kind=intType) :: iiBeg, jjBeg, kkBeg
         integer(kind=intType) :: iiEnd, jjEnd, kkEnd
         integer(kind=intType) :: sizeCGNSWriteType

         character, dimension(:), allocatable :: writeBuffer
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
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

         ! Processor 0 does the writing and must therefore allocate the
         ! writeBuffer.

         if(myID == 0) then

           iiBeg = rangeCell(1,1,1); iiEnd = rangeCell(1,2,1)
           jjBeg = rangeCell(2,1,1); jjEnd = rangeCell(2,2,1)
           kkBeg = rangeCell(3,1,1); kkEnd = rangeCell(3,2,1)
           do ll=2,nSubfaces
             iiBeg = min(iiBeg,rangeCell(1,1,ll))
             jjBeg = min(jjBeg,rangeCell(2,1,ll))
             kkBeg = min(kkBeg,rangeCell(3,1,ll))

             iiEnd = max(iiEnd,rangeCell(1,2,ll))
             jjEnd = max(jjEnd,rangeCell(2,2,ll))
             kkEnd = max(kkEnd,rangeCell(3,2,ll))
           enddo

           mm = (kkEnd-kkBeg+1) * (jjEnd-jjBeg+1) * (iiEnd-iiBeg+1) &
              * sizeCGNSWriteType
           allocate(writeBuffer(mm), stat=ierr)
           if(ierr /= 0)                       &
             call terminate("writeSurfaceSol", &
                            "Memory allocation failure for writeBuffer")
         endif

         ! Loop over the number of solution variables.

         solLoop: do mm=1,nSolVar

           ! Loop over the number of blocks stored on this processor
           ! which may contribute to the subface. Note that
           ! blocksCGNSblock contain all the subblocks part of cgns
           ! block, which is not the same. BlocksCGNSblock cannot be
           ! changed, because it is needed for other subfaces.

           kk = 0
           jj = 0
           do ll=1,nBlocks

             ! Test if the current local block contributes to the
             ! cgns subface to be written.

             if( contributeToFace(ll) ) then

               ! Update the counter kk and store the local block id in ii.

               kk = kk+1
               ii = blocksCGNSblock(ll+offset)

               ! Store the surface solution for this contribution in
               ! buffer. Note that the counter jj is updated in the
               ! routine storeSurfsolInBuffer.

               call storeSurfsolInBuffer(ind, buffer, jj, ii,       &
                                         faceID, cellRange(1,1,kk), &
                                         solNames(mm),              &
                                         viscousSubface)
             endif
           enddo

           ! Make a distinction between processor 0 and the other
           ! processors. Processor 0 does all the writing.

           rootproc: if(myID == 0) then

             ! I am processor 0 and must do the writing.

             ! Loop over the other processors and receive possible
             ! messages.

             do ll=2,nProc

               ! Test if processor ll contributes to this subface.

               if(nMessages(ll) > 0) then

                 ! Receive the message. Note that size is an upper limit.
                 ! Furthermore 1 must be substracted from ll to obtain
                 ! the correct processor id; the processor id's start
                 ! at 0.

                 size   = (il+1)*(jl+1) - jj
                 source = ll -1
                 call mpi_recv(buffer(jj+1), size, sumb_real, source, &
                               source, SUmb_comm_world, status, ierr)

                 ! Determine the true size of the message and update
                 ! the counter jj accordingly.

                 call mpi_get_count(status, sumb_real, size, ierr)
                 jj = jj + size
               endif
             enddo

             ! Copy the variable from buffer into the correct place
             ! in writeBuffer. The routine called depends on the
             ! desired precision.

             ii = 1
             do kk=1,nSubfaces
               select case (precisionGrid)
                 case (precisionSingle)
                   call copyDataBufSinglePrecision(writeBuffer,         &
                                                   buffer(ii),          &
                                                   iiBeg, jjBeg, kkBeg, &
                                                   iiEnd, jjEnd, kkEnd, &
                                                   rangeCell(1,1,kk))
                 case (precisionDouble)
                   call copyDataBufDoublePrecision(writeBuffer,         &
                                                   buffer(ii),          &
                                                   iiBeg, jjBeg, kkBeg, &
                                                   iiEnd, jjEnd, kkEnd, &
                                                   rangeCell(1,1,kk))
               end select

               ! Update the counter ii for the next subface.

               ii = ii + (rangeCell(1,2,kk) - rangeCell(1,1,kk) + 1) &
                  *      (rangeCell(2,2,kk) - rangeCell(2,1,kk) + 1) &
                  *      (rangeCell(3,2,kk) - rangeCell(3,1,kk) + 1)
             enddo

             ! Write the solution variable to file. Source is just used
             ! as a dummy variable and does not have a meaning.

             call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                                   cgnsSol, realTypeCGNS,       &
                                   solNames(mm), writeBuffer,   &
                                   source, ierr)
             if(ierr /= 0)                        &
               call terminate("writeSolCGNSZone", &
                              "Something wrong when &
                              &calling cg_field_write_f")

           else rootproc

             ! Not the root processor.
             ! Data must be sent to processor 0 if local blocks
             ! contribute to the current cgns subface.

             if( jj > 0 ) then
               size = jj
               call mpi_send(buffer, size, sumb_real, 0, myID, &
                             SUmb_comm_world, ierr)
             endif

           endif rootproc

         enddo solLoop

         ! Processor 0 must deallocate the writeBuffer.

         if(myID == 0) then
           deallocate(writeBuffer, stat=ierr)
           if(ierr /= 0)                         &
             call terminate("writeSurfaceSol", &
                            "Deallocation error for writeBuffer")
         endif

         end subroutine writeSurfaceSol

#endif

       end subroutine writeSurfsolCGNSZone
