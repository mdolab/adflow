module writeCGNSSurface

contains

  subroutine writeCGNSSurfaceSol(famList)
    !
    !       writeCGNSSurfaceSol and its subroutines write the surface
    !       solution file(s). The unknowns are stored in the center of the
    !       surface quadrilaterals.
    !
    use cgnsGrid
    use communication
    use su_cgns
    use outputMod
    use inputIteration
    use block
    use blockPointers
    use cgnsNames
    use extraOutput
    use utils, only : terminate, setPointers
    use surfaceFamilies, only : famNames
    use sorting, only : qsortStrings
    implicit none

    ! Input Param
    integer(kind=intType), dimension(:), intent(in) :: famList

    !
    !      Local parameter, the cell dimension.
    !
    integer, parameter :: celldim = 2
    !
    !      Local variables.
    !
    integer :: cgnsInd, ierr

    integer(kind=intType) :: nn, mm, ll, i
    integer(kind=intType) :: nSolVar, nZonesWritten
    character(len=maxStringLen) :: errorMessage
    integer(kind=intType) :: iSurf, nisoSurfVar
    character(len=maxCGNSNameLen), dimension(:), allocatable :: &
         solNames, isoSurfSolNames
    character(len=maxCGNSNameLen) :: contourName
    character(len=maxCGNSNAMELen), dimension(:), allocatable :: famListStr
    ! Determine the number and names of the solution files.

    call surfSolFileNamesWrite

    ! Return immediately if no surface solution files must
    ! be written.

    if(nSurfSolToWrite == 0) return

    ! Write a message that the solution file(s) are being written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "#"
       print "(a,a)", "# Writing surface solution file(s): ", trim(surfSolFileNames(1))
    endif

    ! Allocate the memory for the fileIDs and the bases.

    allocate(fileIDs(nSurfSolToWrite), cgnsBases(nSurfSolToWrite), &
         stat=ierr)
    if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
         "Memory allocation failure for fileIDs &
         &and cgnsBases")

    ! Open the cgns file(s) and write the header. This is only done
    ! by processor 0.

    testRootProc: if(myID == 0) then

       ! Loop over the number of surface solution files to write.

       solLoop: do nn=1,nSurfSolToWrite

          ! Open the cgns file for writing and check if it went okay.
          ! Store the file index for later purposes.
          call cg_open_f(surfSolFileNames(nn), mode_write, cgnsInd, &
               ierr)
          if(ierr /= CG_OK) then
             write(errorMessage,101) trim(surfSolFileNames(nn))
101          format("File",1X,A,1X,"could not be opened by cgns for &
                  &writing")

             call terminate("writeCGNSSurfaceSol", errorMessage)
          endif

          fileIDs(nn) = cgnsInd

          ! Create the base.

          call cg_base_write_f(cgnsInd, "BaseSurfaceSol", celldim, &
               cgnsPhysdim, cgnsBases(nn), ierr)
          if(ierr /= CG_OK)                      &
               call terminate("writeCGNSSurfaceSol", &
               "Something wrong when calling &
               &cg_base_write_f")

          ! Write the header in the cgns file.
          call writeCGNSHeader(cgnsInd, cgnsBases(nn))

       enddo solLoop

    endif testRootProc

    ! Determine the number of variables to be written to the surface
    ! solution file as well as the cgns names.

    call numberOfSurfSolVariables(nSolVar)

    allocate(solNames(nSolVar), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
         "Memory allocation failure for solNames")

    call surfSolNames(solNames)

    ! Based on the famList generate the list of strings corresponding
    ! to the families we want to write.
    allocate(famListStr(size(famList)))
    do i=1,size(famList)
       famListStr(i) = trim(famNames(famList(i)))
    end do

    ! Sort them so we can binary search.
    call qsortStrings(famListStr, size(famListStr))

    ! Loop over the number of cgns blocks and its boundary subfaces
    ! and write the cell centered surface solution of the subface.

    nZonesWritten = 0
    zoneLoop: do nn=1,cgnsNDom

       ! Determine the number of blocks on this processor that belong
       ! to this cgns block.

       mm = nblocksCGNSblock(nn) - nblocksCGNSblock(nn-1)

       ! Loop over the number of boundary subfaces of the original
       ! cgns block and write the cell centered surface solution to
       ! the cgns surface file.

       do ll=1,cgnsDoms(nn)%nBocos

          ! Only write the solution to file if this is a true subface.
          if( cgnsDoms(nn)%bocoInfo(ll)%actualFace )                 &
               call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
               nZonesWritten, .false., famListStr)
       enddo

       ! Loop over the number of internal block boundaries of the
       ! original grid and write the periodic boundaries.

       do ll=1,cgnsDoms(nn)%n1to1

          ! Only periodic boundaries are written; check for this.

          if( cgnsDoms(nn)%conn1to1(ll)%periodic )                   &
               call writeSurfsolCGNSZone(nn, mm, ll, nSolVar, solNames, &
               nZonesWritten, .true., famListStr)
       enddo
    enddo zoneLoop


    ! Check if isosurface will be written. These will be written to
    ! a new base

    testIsoSurafce: if (nIsoSurface > 0)  then

       allocate (cgnsIsoSurfBases(nSurfSolToWrite), stat=ierr)
       testRootProc2: if (myID == 0) then

          ! Loop over the number of surface solution files

          solLoop2: do nn=1,nSurfSolToWrite

             ! Create the new base

             cgnsInd = fileIDs(nn)
             call cg_base_write_f(cgnsInd, "IsoSurfaces", celldim, &
                  cgnsPhysDim, cgnsIsoSurfBases(nn), ierr)
             if (ierr /= CG_OK) &
                  call terminate("WriteCGNSSurfaceSol", &
                  "Something wrong when calling cg_base_write_f for &
                  isoSurface")
          end do solLoop2
       end if testRootProc2

       ! Determine the number of variables to be written to the
       ! isosurface itself well as the cgns names.

       call numberOfIsoSurfVariables(nIsoSurfVar)

       if (nIsoSurfVar > 0) then
          allocate(isoSurfSolNames(nIsoSurfVar), stat=ierr)
          if(ierr /= 0)                           &
               call terminate("writeCGNSSurfaceSol", &
               "Memory allocation failure for isoNames")
          call isoSurfNames(isoSurfSolNames)
       end if

       solLoop3: do ll=1,nSurfSolToWrite ! Numer of spectral instances!
          ! Allocate fn and fc for each domain:
          do nn=1,nDom
             call setPointers(nn, 1, ll)
             allocate(flowDoms(nn, 1, ll)%fn(il, jl, kl))
             allocate(flowDoms(nn, 1, ll)%fc(1:ie, 1:je, 1:ke))
          end do

          ! Finally loop over the required isoSurfaces
          do iSurf=1,nIsoSurface
             call computeIsoVariable(isoSurfaceNames(iSurf), ll, isoValues(iSurf))

11           format(A,A,A,F7.4)
             write(contourName, 11) "Contour ", trim(isoSurfaceNames(iSurf)), "=", isoValues(iSurf)
             call writeIsoSurface(contourName, ll, nIsoSurfVar, isoSurfSolNames)
          end do

          ! deAllocate fn and fc for each domain:
          do nn=1,nDom
             deallocate(flowDoms(nn, 1, ll)%fn, flowDoms(nn, 1, ll)%fc)
          end do
       end do solLoop3

       ! Free memory for bases
       deallocate(cgnsIsoSurfBases, stat=ierr)
       if (nIsoSurfVar > 0) then
          deallocate(isoSurfSolNames)
       end if
    end if testIsoSurafce

    ! Close the cgns file(s). Only processor 0 does this.

    if(myID == 0) then
       do nn=1,nSurfSolToWrite
          call cg_close_f(fileIDs(nn), ierr)
          if(ierr /= CG_OK)                      &
               call terminate("writeCGNSSurfaceSol", &
               "Something wrong when calling cg_close_f")
       enddo
    end if

    ! Deallocate the memory of solNames, fileIDs and cgnsBases.

    deallocate(solNames, fileIDs, cgnsBases, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("writeCGNSSurfaceSol", &
         "Deallocation error for solNames, fileIDs &
         &and cgnsBases")

    ! Wait until all processors (especially processor 0) reach
    ! this point.

    call mpi_barrier(ADflow_comm_world, ierr)

    ! Write a message that the solution file(s) have been written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "# Surface solution file(s) written"
       print "(a)", "#"
    endif
  end subroutine writeCGNSSurfaceSol
  subroutine surfSolFileNamesWrite
    !
    !       surfSolFileNamesWrite determines the names and number of
    !       surface solution files to be written.
    !
    use inputIO
    use inputPhysics
    use inputTimeSpectral
    use monitor
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn

    character(len=7) :: intString

    ! In contrast to the grids and volume solutions, possible states
    ! in the past don't need to be written for the surface. Therefore
    ! the memory allocation can be done independent of the test of
    ! the equation mode we are solving for.

    allocate(surfSolFileNames(nTimeIntervalsSpectral), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("surfSolFileNamesWrite", &
         "Memory allocation failure for surfSolFileNames")

    ! Set the number of surface solution files to be written.

    if( writeSurface ) then
       nSurfSolToWrite = nTimeIntervalsSpectral
    else
       nSurfSolToWrite = 0
    endif

    ! Determine the name(s) of the solution file(s), depending on
    ! the situation.

    select case (equationMode)

    case (steady)

       ! Steady state computation. Possible previous files will
       ! be overwritten.

       surfSolFileNames(1) = surfaceSolFile

       !===============================================================

    case (unsteady)

       ! Unsteady computation. A suffix is added depending on the
       ! time step.

       write(intString,"(i4.4)") timeStepUnsteady + nTimeStepsRestart
       intString = adjustl(intString)

       surfSolFileNames(1) = trim(surfaceSolFile)//"&
            &Timestep"//trim(intString)

       !===============================================================

    case (timeSpectral)

       ! Time spectral computation. A suffix is added depending on
       ! the time instance.

       do nn=1,nTimeIntervalsSpectral
          write(intString,"(i7)") nn
          intString = adjustl(intString)

          surfSolFileNames(nn) = trim(surfaceSolFile)//"&
               &Spectral"//trim(intString)
       enddo

    end select

  end subroutine surfSolFileNamesWrite

  subroutine writeSurfsolCGNSZone(zone, nBlocks, subface, nSolVar, &
       solNames, nZonesWritten, periodic, famListStr)
    !
    !       writeSurfsolCGNSZone writes a surface solution of the given
    !       zone (block) and boundary subface to the cgns surface file(s).
    !       A distinction must be made between true boundaries and
    !       periodic boundaries; the latter are a special kind of internal
    !       block boundaries. This is indicated by the logical periodic.
    !
    use block
    use cgnsGrid
    use communication
    use inputIO
    use su_cgns
    use outputMod
    use utils, only : terminate, convertToLowerCase
    use sorting, only : bsearchStrings
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: zone, nBlocks
    integer(kind=intType), intent(in)    :: subface, nSolVar
    integer(kind=intType), intent(inout) :: nZonesWritten
    character(len=*), dimension(:), intent(in) :: famListStr
    character(len=*), dimension(*), intent(in) :: solNames

    logical, intent(in) :: periodic
    !
    !      Local variables.
    !
    integer :: ierr
    integer :: source, size
    integer :: cgnsBase, cgnsZone, cgnsSol, cgnsInd

    integer, dimension(mpi_status_size) :: mpiStatus

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
    character(len=maxCGNSNameLen) :: tmpStr

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

       ! Possibly do not write this surface if it is not in the famList
       tmpStr = trim(CGNSDoms(zone)%bocoInfo(subface)%wallBCName)
       call convertToLowerCase(tmpStr)
       if (bsearchStrings(tmpStr, famListStr) == 0) then
          return
       end if

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

    call mpi_gather(mBlocks, 1, adflow_integer, nMessages, 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

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

             call mpi_recv(rangeNode(1,1,mm), size, adflow_integer,   &
                  source, source, ADflow_comm_world, mpiStatus, &
                  ierr)
             call mpi_recv(rangeCell(1,1,mm), size, adflow_integer,     &
                  source, source+1, ADflow_comm_world, mpiStatus, &
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
          call mpi_send(nodalRange, size, adflow_integer, 0, myID, &
               ADflow_comm_world, ierr)
          call mpi_send(cellRange, size, adflow_integer, 0, myID+1, &
               ADflow_comm_world, ierr)
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
      !         determineSubranges determines the nodal and cell subrange
      !         for the given local block ID mm in the current cgns subface.
      !
      use inputIO
      implicit none
      !
      !        Local variable
      !
      integer(kind=intType) :: ii

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
      !         createSurfaceZone creates a surface node in the given
      !         cgns surface solution file. This routine should only be
      !         called by processor 0.
      !
      use inputIO
      implicit none
      !
      !        Local variables.
      !
      integer(kind=cgsize_t), dimension(6) :: sizes

      integer(kind=intType) :: nn

      character(len=maxCGNSNameLen) :: zonename
      character(len=7)              :: integerString

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
            case (bcThrust)
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
      if(ierr /= CG_OK)                    &
           call terminate("createSurfaceZone", &
           "Something wrong when calling cg_zone_write_f")

      ! Create the flow solution node.

      call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, &
           "Flow solution", CellCenter, cgnsSol, ierr)
      if(ierr /= CG_OK)                    &
           call terminate("createSurfaceZone", &
           "Something wrong when calling cg_sol_write_f")

      ! Create the rind layers. If rind layers must be stored put
      ! 1 layer on every side of the subface; otherwise put 0 layers.
      ! Use sizes as a buffer to store the rind data. The rind data
      ! must be created under the just created solution node.

      call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
           cgnsZone, "FlowSolution_t", cgnsSol, "end")
      if(ierr /= CG_OK)                    &
           call terminate("createSurfaceZone", &
           "Something wrong when calling cg_goto_f")

      if( storeRindLayer ) then
         sizes(1) = 1; sizes(2) = 1
         sizes(3) = 1; sizes(4) = 1
      else
         sizes(1) = 0; sizes(2) = 0
         sizes(3) = 0; sizes(4) = 0
      endif

      call cg_rind_write_f(int(sizes, intType), ierr)
      if(ierr /= CG_OK)                    &
           call terminate("createSurfaceZone", &
           "Something wrong when calling cg_rind_write_f")

    end subroutine createSurfaceZone

    !        ================================================================

    subroutine writeSurfaceCoord
      !
      !         WriteSurfaceCoord write the vertex values of the
      !         coordinates to the given zone of the cgns surface solution
      !         file.
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

      real(kind=4), dimension(:), allocatable :: writeBuffer4
      real(kind=8), dimension(:), allocatable :: writeBuffer8

      ! Set the cgns real type depending on the input option.

      ! Compute the multiplication factor to obtain the original
      ! coordinates. Note that LRef is corrected to 1.0 when the
      ! coordinates should be written in meters. This happens when
      ! the grid is read.

      LRefInv = one/cgnsDoms(zone)%LRef

      ! Processor 0 does the writing and must therefore allocate the
      ! writeBuffer.

      if(myID == 0) then
         mm = (kEnd-kBeg+1) * (jEnd-jBeg+1) * (iEnd-iBeg+1)
         select case (precisionSurfGrid)
         case (precisionSingle)
            allocate(writeBuffer4(mm), writeBuffer8(0), stat=ierr)
         case (precisionDouble)
            allocate(writeBuffer4(0), writeBuffer8(mm), stat=ierr)
         end select
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
                  call mpi_recv(buffer(jj+1), size, adflow_real, source, &
                       source, ADflow_comm_world, mpiStatus, ierr)

                  ! Determine the true size of the message and update
                  ! the counter jj accordingly.

                  call mpi_get_count(mpiStatus, adflow_real, size, ierr)
                  jj = jj + size
               endif
            enddo

            ! Copy the coordinate from buffer into the correct place
            ! in writeBuffer. The routine called depends on the
            ! desired precision.

            ii = 1
            do kk=1,nSubfaces
               select case (precisionSurfGrid)
               case (precisionSingle)
                  call copyDataBufSinglePrecision(writeBuffer4,      &
                       buffer(ii),       &
                       iBeg, jBeg, kBeg, &
                       iEnd, jEnd, kEnd, &
                       rangeNode(1,1,kk))
               case (precisionDouble)
                  call copyDataBufDoublePrecision(writeBuffer8,      &
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
            select case (precisionSurfGrid)
            case (precisionSingle)
               select case (mm)
               case (1_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realSingle, cgnsCoorX, writeBuffer4, source, ierr)
               case (2_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realSingle, cgnsCoorY, writeBuffer4, source, ierr)
               case (3_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realSingle, cgnsCoorZ, writeBuffer4, source, ierr)
               end select

            case (precisionDouble)
               select case(mm)
               case (1_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realDouble, cgnsCoorX, writeBuffer8, source, ierr)
               case (2_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realDouble, cgnsCoorY, writeBuffer8, source, ierr)
               case (3_intType)
                  call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, &
                       realDouble, cgnsCoorZ, writeBuffer8, source, ierr)
               end select
            end select

            if(ierr /= CG_OK)                    &
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
               if(ierr /= CG_OK)                      &
                    call terminate("writeSurfaceCoord", &
                    "Something wrong when calling cg_goto_f")

               ! Write the units.

               call cg_units_write_f(cgnsDoms(zone)%mass, &
                    cgnsDoms(zone)%len,  &
                    cgnsDoms(zone)%time, &
                    cgnsDoms(zone)%temp, &
                    cgnsDoms(zone)%angle, ierr)
               if(ierr /= CG_OK)                    &
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
               call mpi_send(buffer, size, adflow_real, 0, myID, &
                    ADflow_comm_world, ierr)
            endif

         endif rootproc

      enddo coorLoop

      ! Processor 0 must deallocate the writeBuffer.

      if(myID == 0) then
         deallocate(writeBuffer4, writebuffer8, stat=ierr)
         if(ierr /= 0)                         &
              call terminate("writeSurfaceCoord", &
              "Deallocation error for writeBuffer")
      endif

    end subroutine writeSurfaceCoord

    !        ================================================================

    subroutine writeSurfaceSol
      !
      !         writeSurfaceSol writes the cell centered surface solution
      !         to the cgns surface file.
      !
      implicit none
      !
      !        Local variables.
      !
      integer :: realTypeCGNS

      integer(kind=intType) :: ii, jj, kk, ll, mm
      integer(kind=intType) :: iiBeg, jjBeg, kkBeg
      integer(kind=intType) :: iiEnd, jjEnd, kkEnd

      real(kind=4), dimension(:), allocatable :: writeBuffer4
      real(kind=8), dimension(:), allocatable :: writeBuffer8

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

         mm = (kkEnd-kkBeg+1) * (jjEnd-jjBeg+1) * (iiEnd-iiBeg+1)
         select case (precisionSurfSol)
         case (precisionSingle)
            allocate(writeBuffer4(mm), writeBuffer8(0), stat=ierr)
         case (precisionDouble)
            allocate(writeBuffer4(0), writeBuffer8(mm), stat=ierr)
         end select
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
                    viscousSubface, storeRindLayer)
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
                  call mpi_recv(buffer(jj+1), size, adflow_real, source, &
                       source, ADflow_comm_world, mpiStatus, ierr)

                  ! Determine the true size of the message and update
                  ! the counter jj accordingly.

                  call mpi_get_count(mpiStatus, adflow_real, size, ierr)
                  jj = jj + size
               endif
            enddo

            ! Copy the variable from buffer into the correct place
            ! in writeBuffer. The routine called depends on the
            ! desired precision.

            ii = 1
            do kk=1,nSubfaces
               select case (precisionSurfSol)
               case (precisionSingle)
                  call copyDataBufSinglePrecision(writeBuffer4,        &
                       buffer(ii),          &
                       iiBeg, jjBeg, kkBeg, &
                       iiEnd, jjEnd, kkEnd, &
                       rangeCell(1,1,kk))
               case (precisionDouble)
                  call copyDataBufDoublePrecision(writeBuffer8,        &
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
            select case(precisionSurfSol)
            case (precisionSingle)
               call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                    cgnsSol, realSingle, solNames(mm), writeBuffer4, &
                    source, ierr)
            case (precisionDouble)
               call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, &
                    cgnsSol, realDouble, solNames(mm), writeBuffer8, &
                    source, ierr)
            end select
            if(ierr /= 0)           then
               call terminate("writeSolCGNSZone", &
                    "Something wrong when &
                    &calling cg_field_write_f")
            end if

         else rootproc

            ! Not the root processor.
            ! Data must be sent to processor 0 if local blocks
            ! contribute to the current cgns subface.

            if( jj > 0 ) then
               size = jj
               call mpi_send(buffer, size, adflow_real, 0, myID, &
                    ADflow_comm_world, ierr)
            endif

         endif rootproc

      enddo solLoop

      ! Processor 0 must deallocate the writeBuffer.

      if(myID == 0) then
         deallocate(writeBuffer4, writeBuffer8, stat=ierr)
         if(ierr /= 0)                         &
              call terminate("writeSurfaceSol", &
              "Deallocation error for writeBuffer")
      endif

    end subroutine writeSurfaceSol
  end subroutine writeSurfsolCGNSZone

  subroutine writeIsoSurface(isoName , sps, nIsoSurfVar, isoSurfSolNames)

    ! Implements a marching cubes algrorithm which can be used to
    ! extract iso surfaces or slcies from a solution and store them in a
    ! CGNS surface file.

    use communication
    use blockPointers
    use flowVarRefState
    use inputPhysics
    use su_cgns
    use inputIO
    use outputMod
    use cgnsNames
    use utils, only : reallocateReal2, setPointers, reallocateinteger2, &
         terminate, EChk, pointReduce
    implicit none

    ! Input param
    character(len=*), intent(in)                   :: isoName
    integer(kind=intType) :: sps
    integer(kind=intType) :: nIsoSurfVar
    character(len=*), dimension(*), intent(in) :: isoSurfSolNames
    ! Working param
    integer(kind=intType) :: i, j, k, nn, kk, nMax, icon, iCoor, indexCube, num1, num2
    real(kind=realType), dimension(:, :, :), pointer :: fn
    real(kind=realType), dimension(:, :), pointer :: coords, uniqueCoords, weights
    real(kind=realType), dimension(:), allocatable :: buffer

    integer(kind=intType), dimension(:, :), pointer :: indices
    integer(kind=intType), dimension(:, :), allocatable :: connBuffer, allConn
    integer(kind=intType), dimension(:), allocatable :: link, nPtsProc, nConnProc

    integer(kind=intType) :: ccwOrdering(3, 8), n1(3), n2(3)
    integer(kind=intType) :: matCon1(256, 16), matCon2(12,2), nUnique, ivar
    integer(kind=intType) :: ierr, iProc, tag, cumNodes, cumConn, iCorner, curBlock, idim
    real(kind=realType) :: f(8)

    logical :: logic1
    integer, dimension(mpi_status_size) :: mpiStatus
    integer(kind=intType) :: cgnsInd, cgnsBase, cgnsZOne, coordID, secID, solID, fieldID
    real(kind=realType) :: tol=1e-8 ! Node tol for isosurf pointReduce

    ! Fill up the connecivity matrices
    call getMatCons(matCon1, matCon2, ccwOrdering)

    ! Generate a guess for the size of the iso surface -- sum of the
    ! number of nodes on the faces
    nMax = 0
    do nn=1, nDom
       call setPointers(nn, 1, sps)
       nMax = nMax + 2*il*jl + 2*il*kl + 2*jl*kl
    end do

    ! Allocate size nMax
    allocate(weights(2, nMax), indices(7, nMax))

    iCoor = 0
    do nn=1, nDom
       call setPointers(nn, 1, sps)
       fn => flowDoms(nn, 1, sps)%fn
       ! Now loop over the Cells:
       do k=2, kl
          do j=2, jl
             do i=2, il

                ! Extract value at the corners
                do iCorner=1,8
                   f(iCorner) = fn(i+ccwOrdering(1, iCorner), &
                        j + ccwOrdering(2, iCorner), &
                        k + ccwOrdering(3, iCorner))
                end do

                ! Based on the values at each corner, determine which
                ! type surface we have
                indexcube = 1
                if (f(1) .lt. zero) indexcube = indexcube + 1
                if (f(2) .lt. zero) indexcube = indexcube + 2
                if (f(3) .lt. zero) indexcube = indexcube + 4
                if (f(4) .lt. zero) indexcube = indexcube + 8
                if (f(5) .lt. zero) indexcube = indexcube + 16
                if (f(6) .lt. zero) indexcube = indexcube + 32
                if (f(7) .lt. zero) indexcube = indexcube + 64
                if (f(8) .lt. zero) indexcube = indexcube + 128

                logic1 = .true.

                kk = 1
                do while (logic1)
                   icon = matcon1(indexcube, kk)

                   if (icon == 0) then
                      logic1=.false.
                   else

                      iCoor = iCoor + 1
                      if (iCoor > nMax) then
                         ! Need to realloc the coord array. Make it double the size
                         call reallocateReal2(weights, 2, 2*nMax, 2, nMax, .true.)
                         call reallocateInteger2(indices, 7, 2*nMax, 7, nMax, .true.)
                         nMax = nMax * 2
                      end if

                      num1 = matcon2(icon,1)
                      num2 = matcon2(icon,2)

                      ! Weight factors
                      weights(2, iCoor) = (zero - f(num1))/(f(num2) - f(num1))
                      weights(1, iCoor) = one - weights(2, icoor)

                      ! Indices of nodes
                      n1 = (/i, j, k/) + ccwOrdering(: ,num1)
                      n2 = (/i, j, k/) + ccwOrdering(:, num2)
                      indices(:, iCoor) = (/nn, n1(1), n1(2), n1(3), n2(1), n2(2), n2(3)/)

                      kk = kk + 1
                   end if
                end do

             end do ! I loop
          end do ! J loop
       end do ! K loop
    end do ! Domain loop

    ! We have not actually stored the coordintes; only the positions and
    ! the weights. To compute the coordinates we pass back through and assemble
    allocate(Coords(3, iCoor))

    ! Set pointer to first block
    call setPointers(1, 1, sps)
    curBlock = 1
    do i=1,iCoor

       ! If we've switched blocks, reset points. This stil only calls
       ! setPointer nDom times since there are at most that many
       ! switches
       if (indices(1, i) /= curBlock) then
          call setPointers(indices(1, i), 1, sps)
          curBlock = indices(1, i)
       end if

       ! Computing coordinates is easy; we just juse the weights and the
       ! indices on x
       do idim=1,3
          coords(idim, i) = &
               weights(1, i) * &
               X(indices(2,i), indices(3,i), indices(4,i), idim) + &
               weights(2, i) * &
               X(indices(5,i), indices(6,i), indices(7,i), idim)
       end do
    end do

    ! Now we know the maximum number of coordinates so we can allocate
    ! the unique set and the link array
    allocate(uniqueCoords(3, icoor))
    allocate(link(icoor))

    ! Compute the reduced set of coordinates. The sole purpose of this
    ! is to reduce the filesize. This will typicaly reduce the number of
    ! coordinates by about a factor of 4.

    call pointReduce(coords, iCoor, tol, uniqueCoords, link, nUnique)

    ! Now that we have produced the desired isosurface on each
    ! processor. Communicate the number of number of coordinates and the
    ! number of triangles each proc is going to send to the root:

    allocate(nPtsProc(nProc), nConnProc(nProc))
    nPtsProc(:) = 0_intType
    nConnProc(:) = 0_intType

    call  MPI_Allgather(nUnique, 1, mpi_integer4, nPtsProc, 1, mpi_integer4, &
         adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    call  MPI_Allgather(iCoor/3, 1, mpi_integer4, nConnProc, 1, mpi_integer4, &
         adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)


    if (sum(nPtsProc) > 0) then

       if (myid == 0) then ! Root proc does the writing

          ! Write a new zone:
          cgnsInd = fileIDs(sps)
          cgnsBase = cgnsIsoSurfBases(sps)

          ! Write the unstructured zone
          call cg_zone_write_f(cgnsInd, cgnsBase, isoName, int((/sum(nPtsProc), sum(nConnProc), 0/), cgsize_t), &
               Unstructured, cgnsZone, ierr)
          if (ierr .eq. CG_ERROR) call cg_error_exit_f

          if(ierr /= CG_OK)                    &
               call terminate("writeIsoSurface", &
               "Something wrong when calling cg_zone_write_f")
       end if

    else
       if (myid == 0) then
          ! We don't actually have an isosurface. We will create a zone
          ! that the same structure, but contains only a single triangle
          ! with all the coordinates at zero. This way the zone still
          ! exists and yields a uniform structure which can make
          ! processing easier

          ! Write a new zone:
          cgnsInd = fileIDs(sps)
          cgnsBase = cgnsIsoSurfBases(sps)

          call writeEmptyZone

       end if
       ! Don't forget to deallocate the stuff allocated so far:
       deallocate(nPtsProc, nConnProc, link, uniqueCoords, coords, weights, indices)
       return
    end if

    ! We need to keep track of the cumulative number of nodes since each
    ! proc has done its own ordering
    cumNodes = 0

    ! Communicate and write the coordinates
    do iproc=0, nProc-1

       dataOnProc: if (myid == iproc) then
          allocate(buffer(3*nPtsProc(iProc+1)))

          ! We will swap the order of the coordinates to packed format
          ! since this is what we need for CGNS

          do i=1,nPtsProc(iProc+1)
             buffer(i)                      = uniqueCoords(1, i)
             buffer(1*nPtsProc(iProc+1)+ i) = uniqueCoords(2, i)
             buffer(2*nPtsProc(iProc+1)+ i) = uniqueCoords(3, i)
          end do

       end if dataOnProc

       if (iproc .ne. 0) then
          tag = 13
          if (myid == 0) then
             ! allocate space for the recv
             allocate(buffer(3*nPtsProc(iProc+1)))

             call mpi_recv(buffer, nPtsProc(iProc+1)*3, adflow_real, iProc, tag, &
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if

          if (myid == iProc) then
             call mpi_send(buffer, nPtsProc(iProc+1)*3, adflow_real, 0, tag, &
                  adflow_comm_world, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if
       end if

       if (myid == 0) then
          ! Now do partial writes on the root proc with points we've
          ! received from iProc
          if (nPtsProc(iProc+1) > 0) then

             call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                  'CoordinateX', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                  buffer(1:nPtsProc(iProc+1)), coordID, ierr)

             call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                  'CoordinateY', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                  buffer(nPtsProc(iProc+1)+1:2*nPtsProc(iProc+1)), coordID, ierr)

             call cg_coord_partial_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
                  'CoordinateZ', cumNodes+1, cumNodes+nPtsProc(iProc+1), &
                  buffer(2*nPtsProc(iProc+1)+1:3*nPtsProc(iProc+1)), coordID, ierr)

             if(ierr /= CG_OK)                    &
                  call terminate("writeIsoSurface", &
                  "Something wrong when calling cg_coord_write_f")

             ! Increment by the number of nodes on this proc
             cumNodes = cumNodes + nPtsProc(iProc+1)
          end if
       end if

       ! Buffer was only allocated on root and current iProc
       if (myid == iProc .or. myid == 0) then
          deallocate(buffer)
       end if
    end do

    ! We need to keep track of the cumulative number of nodes since each
    ! proc has done its own ordering
    cumNodes = 0
    cumConn  = 0

    ! The partial write functionality is different between versions 2.5
    ! and 3.1, so we will just gather all the connectivities and do a
    ! final write at the end
    if (myid == 0) then
       allocate(allConn(3, sum(nConnProc)))
    endif

    ! Communicate and write the connectivity
    do iProc=0, nProc-1
       connOnProc: if (myid == iProc) then
          allocate(connBuffer(3,nConnProc(iProc+1)))
          do i=1,nConnProc(iProc+1)
             connBuffer(1, i) = link(3*i-2)
             connBuffer(2, i) = link(3*i-1)
             connBuffer(3, i) = link(3*i  )
          end do
       end if connOnProc

       ! Communication is only necessary if we are not dealing with root
       ! proc:
       if (iproc .ne. 0) then
          tag = 13
          if (myid == 0) then
             ! allocate space for the recv
             allocate(connBuffer(3,nConnProc(iProc+1)))
             call mpi_recv(connBuffer, nConnProc(iProc+1)*3, adflow_integer, iProc, tag, &
                  adflow_comm_world, mpiStatus, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if

          if (myid == iProc) then
             call mpi_send(connBuffer, nConnProc(iProc+1)*3, adflow_integer, 0, tag, &
                  adflow_comm_world, ierr)
             call EChk(ierr, __FILE__, __LINE__)
          end if
       endif

       if (myid == 0) then

          ! Copy into the allCon array and increment the received
          ! local connectivity by cummNodes

          do i=1,nConnProc(iProc+1)
             allConn(:, cumConn + i) = connBuffer(:, i) + cumNodes
          end do

          cumNodes = cumNodes + nPtsProc(iProc+1)
          cumConn  = cumConn  + nConnProc(iProc+1)
       end if

       ! Buffer was only allocated on roto and iProc
       if (myid == iProc .or. myid == 0) then
          deallocate(connBuffer)
       end if
    end do

    ! Finally do the (full) connectivity write
    if (myid == 0) then
       ! Now write on root proc:

       ! Write just the connectively we have in buffer
       call cg_section_write_f(cgnsInd, cgnsBase, cgnsZone, "ELEM", TRI_3, &
            1, sum(nConnProc), 0, allConn, secID, ierr)
       if(ierr /= CG_OK)                    &
            call terminate("writeIsoSurface", &
            "Something wrong when calling cg_section_partial_write_f")

       ! Also free allConn
       deallocate(allConn)
    end if

    ! Finally we have to write solution data for the iso surface
    ! iself. The main reason is that the same code is used for "slices"
    ! as well and in that case, you want to have other data interpolated
    ! on the "isoSurafce" (slice)

    ! Write the solution node:
    if (myid == 0) then
       call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, "isoSurfSolution", Vertex, solID, ierr)
       if(ierr /= CG_OK)                    &
            call terminate("writeIsoSurface", &
            "Something wrong when calling cg_sol_write_f")
    end if

    ! Make the buffer large enough
    allocate(buffer(maxval(nPtsProc)))

    ! Loop over variables to write:
    do iVar=1,nIsoSurfVar

       ! We will reuse the same code as was used for computing the value
       ! onwhich we did the interpolation. However, set 'zero' for the
       ! isovalue such that we get the true value back

       call computeIsoVariable(isoSurfSolNames(iVar), sps, zero)

       ! Set points to first block:
       call setPointers(1, 1, sps)
       curBlock = 1
       fn => flowDoms(1, 1, sps)%fn

       do i=1,iCoor

          ! If we've switched blocks, reset points. This stil only calls
          ! setPointer nDom times since there are at most that many
          ! switches
          if (indices(1, i) /= curBlock) then
             call setPointers(indices(1, i), 1, sps)
             curBlock = indices(1, i)
             fn => flowDoms(curBlock, 1, sps)%fn
          end if

          ! Computing interpolated value is easy using weights:
          ! indices on x
          buffer(link(i)) = weights(1, i) * fn(indices(2,i), indices(3,i), indices(4,i)) + &
               weights(2, i) * fn(indices(5,i), indices(6,i), indices(7,i))
       end do

       cumNodes = 0
       ! Communicate and write the solutions
       do iproc=0, nProc-1

          if (iproc .ne. 0) then
             tag = 13
             if (myid == 0) then
                call mpi_recv(buffer, nPtsProc(iProc+1), adflow_real, iProc, tag, &
                     adflow_comm_world, mpiStatus, ierr)
                call EChk(ierr, __FILE__, __LINE__)
             end if

             if (myid == iProc) then
                call mpi_send(buffer, nPtsProc(iProc+1), adflow_real, 0, tag, &
                     adflow_comm_world, ierr)
                call EChk(ierr, __FILE__, __LINE__)
             end if
          end if

          if (myid == 0) then
             ! Now do partial writes on the root proc with points we've
             ! received from iProc
             if (nPtsProc(iProc+1) > 0) Then
                call cg_field_partial_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
                     isoSurfSolNames(iVar), cumNodes+1, cumNodes + nPtsProc(iProc+1), &
                     buffer, fieldID, ierr)

                if(ierr /= CG_OK)                    &
                     call terminate("writeIsoSurface", &
                     "Something wrong when calling cg_field_partial_write_f")

                ! Increment by the number of nodes on this proc
                cumNodes = cumNodes + nPtsProc(iProc+1)
             end if
          end if
       end do
    end do
    ! Everyone deallocs buffer
    deallocate(buffer)

    ! Clear up temporary allocatable data.
    deallocate(nPtsProc, nConnProc)
    deallocate(coords, uniqueCoords, link)
    deallocate(weights, indices)

  contains

    subroutine writeEmptyZone

      call cg_zone_write_f(cgnsInd, cgnsBase, isoName, int((/3, 1, 0/), cgsize_t), &
           Unstructured, cgnsZone, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
           'CoordinateX', (/zero, zero, zero/), coordID, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
           'CoordinateY', (/zero, zero, zero /), coordID, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
           'CoordinateZ', (/zero, zero, zero/), coordID, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      call cg_section_write_f(cgnsInd, cgnsBase, cgnsZone, "ELEM", TRI_3, &
           1, 1, 0, (/1, 2, 3/), secID, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, "isoSurfSolution", Vertex, solID, ierr)
      if (ierr .eq. CG_ERROR) call cg_error_exit_f

      do iVar = 1, nIsoSurfVar
         call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, isoSurfSolNames(iVar), &
              (/zero, zero, zero/), fieldID, ierr)
         if (ierr .eq. CG_ERROR) call cg_error_exit_f
      end do
    end subroutine writeEmptyZone

  end subroutine writeIsoSurface


  subroutine getMatCons(matcon1, matcon2, ccwOrdering)

    use precision
    implicit none
    integer(kind=intType) :: matCon1(256, 16), matCon2(12, 2), ccwOrdering(3, 8)

    matcon2=RESHAPE((/1,2,3,1,5,6,7,5,1,2,3,4,2,3,4,4,6,7,8,8,5,6,7,8/),(/12,2/))

    ccwOrdering(:, 1) = (/-1, -1, -1/)
    ccwOrdering(:, 2) = (/ 0, -1, -1/)
    ccwOrdering(:, 3) = (/ 0,  0, -1/)
    ccwOrdering(:, 4) = (/-1,  0, -1/)
    ccwOrdering(:, 5) = (/-1, -1,  0/)
    ccwOrdering(:, 6) = (/ 0, -1,  0/)
    ccwOrdering(:, 7) = (/ 0,  0,  0/)
    ccwOrdering(:, 8) = (/-1,  0,  0/)

    matCon1(1,:)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(2,:)=(/1,9,4,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(3,:)=(/1,2,10,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(4,:)=(/2,9,4,10,9,2,0,0,0,0,0,0,0,0,0,0/)
    matCon1(5,:)=(/2,3,11,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(6,:)=(/1,9,4,2,3,11,0,0,0,0,0,0,0,0,0,0/)
    matCon1(7,:)=(/10,3,11,1,3,10,0,0,0,0,0,0,0,0,0,0/)
    matCon1(8,:)=(/3,9,4,3,11,9,11,10,9,0,0,0,0,0,0,0/)
    matCon1(9,:)=(/4,12,3,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(10,:)=(/1,12,3,9,12,1,0,0,0,0,0,0,0,0,0,0/)
    matCon1(11,:)=(/2,10,1,3,4,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(12,:)=(/2,12,3,2,10,12,10,9,12,0,0,0,0,0,0,0/)
    matCon1(13,:)=(/4,11,2,12,11,4,0,0,0,0,0,0,0,0,0,0/)
    matCon1(14,:)=(/1,11,2,1,9,11,9,12,11,0,0,0,0,0,0,0/)
    matCon1(15,:)=(/4,10,1,4,12,10,12,11,10,0,0,0,0,0,0,0/)
    matCon1(16,:)=(/10,9,11,11,9,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(17,:)=(/5,8,9,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(18,:)=(/5,4,1,8,4,5,0,0,0,0,0,0,0,0,0,0/)
    matCon1(19,:)=(/1,2,10,9,5,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(20,:)=(/5,2,10,5,8,2,8,4,2,0,0,0,0,0,0,0/)
    matCon1(21,:)=(/2,3,11,9,5,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(22,:)=(/4,5,8,4,1,5,2,3,11,0,0,0,0,0,0,0/)
    matCon1(23,:)=(/10,3,11,10,1,3,9,5,8,0,0,0,0,0,0,0/)
    matCon1(24,:)=(/3,11,10,3,10,8,3,8,4,8,10,5,0,0,0,0/)
    matCon1(25,:)=(/9,5,8,4,12,3,0,0,0,0,0,0,0,0,0,0/)
    matCon1(26,:)=(/12,5,8,12,3,5,3,1,5,0,0,0,0,0,0,0/)
    matCon1(27,:)=(/10,1,2,9,5,8,3,4,12,0,0,0,0,0,0,0/)
    matCon1(28,:)=(/5,8,12,10,5,12,10,12,3,10,3,2,0,0,0,0/)
    matCon1(29,:)=(/4,11,2,4,12,11,8,9,5,0,0,0,0,0,0,0/)
    matCon1(30,:)=(/2,12,11,2,5,12,2,1,5,8,12,5,0,0,0,0/)
    matCon1(31,:)=(/5,8,9,10,1,12,10,12,11,12,1,4,0,0,0,0/)
    matCon1(32,:)=(/5,8,12,5,12,10,10,12,11,0,0,0,0,0,0,0/)
    matCon1(33,:)=(/10,6,5,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(34,:)=(/10,6,5,1,9,4,0,0,0,0,0,0,0,0,0,0/)
    matCon1(35,:)=(/1,6,5,2,6,1,0,0,0,0,0,0,0,0,0,0/)
    matCon1(36,:)=(/9,6,5,9,4,6,4,2,6,0,0,0,0,0,0,0/)
    matCon1(37,:)=(/2,3,11,10,6,5,0,0,0,0,0,0,0,0,0,0/)
    matCon1(38,:)=(/4,1,9,2,3,11,5,10,6,0,0,0,0,0,0,0/)
    matCon1(39,:)=(/6,3,11,6,5,3,5,1,3,0,0,0,0,0,0,0/)
    matCon1(40,:)=(/3,11,6,4,3,6,4,6,5,4,5,9,0,0,0,0/)
    matCon1(41,:)=(/10,6,5,3,4,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(42,:)=(/1,12,3,1,9,12,5,10,6,0,0,0,0,0,0,0/)
    matCon1(43,:)=(/1,6,5,1,2,6,3,4,12,0,0,0,0,0,0,0/)
    matCon1(44,:)=(/3,2,6,3,6,9,3,9,12,5,9,6,0,0,0,0/)
    matCon1(45,:)=(/11,4,12,11,2,4,10,6,5,0,0,0,0,0,0,0/)
    matCon1(46,:)=(/5,10,6,1,9,2,9,11,2,9,12,11,0,0,0,0/)
    matCon1(47,:)=(/6,5,1,6,1,12,6,12,11,12,1,4,0,0,0,0/)
    matCon1(48,:)=(/6,5,9,6,9,11,11,9,12,0,0,0,0,0,0,0/)
    matCon1(49,:)=(/10,8,9,6,8,10,0,0,0,0,0,0,0,0,0,0/)
    matCon1(50,:)=(/10,4,1,10,6,4,6,8,4,0,0,0,0,0,0,0/)
    matCon1(51,:)=(/1,8,9,1,2,8,2,6,8,0,0,0,0,0,0,0/)
    matCon1(52,:)=(/2,6,4,4,6,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(53,:)=(/10,8,9,10,6,8,11,2,3,0,0,0,0,0,0,0/)
    matCon1(54,:)=(/11,2,3,10,6,1,6,4,1,6,8,4,0,0,0,0/)
    matCon1(55,:)=(/9,1,3,9,3,6,9,6,8,11,6,3,0,0,0,0/)
    matCon1(56,:)=(/3,11,6,3,6,4,4,6,8,0,0,0,0,0,0,0/)
    matCon1(57,:)=(/8,10,6,8,9,10,4,12,3,0,0,0,0,0,0,0/)
    matCon1(58,:)=(/10,6,8,10,8,3,10,3,1,3,8,12,0,0,0,0/)
    matCon1(59,:)=(/3,4,12,1,2,9,2,8,9,2,6,8,0,0,0,0/)
    matCon1(60,:)=(/12,3,2,12,2,8,8,2,6,0,0,0,0,0,0,0/)
    matCon1(61,:)=(/10,6,9,9,6,8,11,2,4,11,4,12,0,0,0,0/)
    matCon1(62,:)=(/6,8,1,6,1,10,8,12,1,2,1,11,12,11,1,0/)
    matCon1(63,:)=(/12,11,1,12,1,4,11,6,1,9,1,8,6,8,1,0/)
    matCon1(64,:)=(/12,11,6,8,12,6,0,0,0,0,0,0,0,0,0,0/)
    matCon1(65,:)=(/11,7,6,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(66,:)=(/1,9,4,6,11,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(67,:)=(/10,1,2,6,11,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(68,:)=(/2,9,4,2,10,9,6,11,7,0,0,0,0,0,0,0/)
    matCon1(69,:)=(/2,7,6,3,7,2,0,0,0,0,0,0,0,0,0,0/)
    matCon1(70,:)=(/2,7,6,2,3,7,4,1,9,0,0,0,0,0,0,0/)
    matCon1(71,:)=(/10,7,6,10,1,7,1,3,7,0,0,0,0,0,0,0/)
    matCon1(72,:)=(/6,10,9,6,9,3,6,3,7,4,3,9,0,0,0,0/)
    matCon1(73,:)=(/3,4,12,11,7,6,0,0,0,0,0,0,0,0,0,0/)
    matCon1(74,:)=(/12,1,9,12,3,1,11,7,6,0,0,0,0,0,0,0/)
    matCon1(75,:)=(/1,2,10,3,4,12,6,11,7,0,0,0,0,0,0,0/)
    matCon1(76,:)=(/6,11,7,2,10,3,10,12,3,10,9,12,0,0,0,0/)
    matCon1(77,:)=(/7,4,12,7,6,4,6,2,4,0,0,0,0,0,0,0/)
    matCon1(78,:)=(/1,9,12,1,12,6,1,6,2,6,12,7,0,0,0,0/)
    matCon1(79,:)=(/4,12,7,1,4,7,1,7,6,1,6,10,0,0,0,0/)
    matCon1(80,:)=(/7,6,10,7,10,12,12,10,9,0,0,0,0,0,0,0/)
    matCon1(81,:)=(/6,11,7,5,8,9,0,0,0,0,0,0,0,0,0,0/)
    matCon1(82,:)=(/5,4,1,5,8,4,7,6,11,0,0,0,0,0,0,0/)
    matCon1(83,:)=(/2,10,1,6,11,7,9,5,8,0,0,0,0,0,0,0/)
    matCon1(84,:)=(/11,7,6,2,10,8,2,8,4,8,10,5,0,0,0,0/)
    matCon1(85,:)=(/7,2,3,7,6,2,5,8,9,0,0,0,0,0,0,0/)
    matCon1(86,:)=(/2,3,6,6,3,7,4,1,5,4,5,8,0,0,0,0/)
    matCon1(87,:)=(/9,5,8,10,1,6,1,7,6,1,3,7,0,0,0,0/)
    matCon1(88,:)=(/8,4,10,8,10,5,4,3,10,6,10,7,3,7,10,0/)
    matCon1(89,:)=(/4,12,3,8,9,5,11,7,6,0,0,0,0,0,0,0/)
    matCon1(90,:)=(/6,11,7,5,8,3,5,3,1,3,8,12,0,0,0,0/)
    matCon1(91,:)=(/1,2,10,5,8,9,3,4,12,6,11,7,0,0,0,0/)
    matCon1(92,:)=(/10,3,2,10,12,3,10,5,12,8,12,5,6,11,7,0/)
    matCon1(93,:)=(/9,5,8,4,12,6,4,6,2,6,12,7,0,0,0,0/)
    matCon1(94,:)=(/6,2,12,6,12,7,2,1,12,8,12,5,1,5,12,0/)
    matCon1(95,:)=(/1,6,10,1,7,6,1,4,7,12,7,4,9,5,8,0/)
    matCon1(96,:)=(/7,6,10,7,10,12,5,8,10,8,12,10,0,0,0,0/)
    matCon1(97,:)=(/11,5,10,7,5,11,0,0,0,0,0,0,0,0,0,0/)
    matCon1(98,:)=(/5,11,7,5,10,11,1,9,4,0,0,0,0,0,0,0/)
    matCon1(99,:)=(/11,1,2,11,7,1,7,5,1,0,0,0,0,0,0,0/)
    matCon1(100,:)=(/9,4,2,9,2,7,9,7,5,7,2,11,0,0,0,0/)
    matCon1(101,:)=(/2,5,10,2,3,5,3,7,5,0,0,0,0,0,0,0/)
    matCon1(102,:)=(/4,1,9,2,3,10,3,5,10,3,7,5,0,0,0,0/)
    matCon1(103,:)=(/1,3,5,5,3,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(104,:)=(/9,4,3,9,3,5,5,3,7,0,0,0,0,0,0,0/)
    matCon1(105,:)=(/11,5,10,11,7,5,12,3,4,0,0,0,0,0,0,0/)
    matCon1(106,:)=(/1,9,3,3,9,12,5,10,11,5,11,7,0,0,0,0/)
    matCon1(107,:)=(/4,12,3,1,2,7,1,7,5,7,2,11,0,0,0,0/)
    matCon1(108,:)=(/7,5,2,7,2,11,5,9,2,3,2,12,9,12,2,0/)
    matCon1(109,:)=(/10,7,5,10,4,7,10,2,4,12,7,4,0,0,0,0/)
    matCon1(110,:)=(/9,12,2,9,2,1,12,7,2,10,2,5,7,5,2,0/)
    matCon1(111,:)=(/4,12,7,4,7,1,1,7,5,0,0,0,0,0,0,0/)
    matCon1(112,:)=(/7,5,9,12,7,9,0,0,0,0,0,0,0,0,0,0/)
    matCon1(113,:)=(/8,11,7,8,9,11,9,10,11,0,0,0,0,0,0,0/)
    matCon1(114,:)=(/1,8,4,1,11,8,1,10,11,7,8,11,0,0,0,0/)
    matCon1(115,:)=(/11,7,8,2,11,8,2,8,9,2,9,1,0,0,0,0/)
    matCon1(116,:)=(/11,7,8,11,8,2,2,8,4,0,0,0,0,0,0,0/)
    matCon1(117,:)=(/2,3,7,2,7,9,2,9,10,9,7,8,0,0,0,0/)
    matCon1(118,:)=(/3,7,10,3,10,2,7,8,10,1,10,4,8,4,10,0/)
    matCon1(119,:)=(/8,9,1,8,1,7,7,1,3,0,0,0,0,0,0,0/)
    matCon1(120,:)=(/8,4,3,7,8,3,0,0,0,0,0,0,0,0,0,0/)
    matCon1(121,:)=(/3,4,12,11,7,9,11,9,10,9,7,8,0,0,0,0/)
    matCon1(122,:)=(/3,1,8,3,8,12,1,10,8,7,8,11,10,11,8,0/)
    matCon1(123,:)=(/2,9,1,2,8,9,2,11,8,7,8,11,3,4,12,0/)
    matCon1(124,:)=(/12,3,2,12,2,8,11,7,2,7,8,2,0,0,0,0/)
    matCon1(125,:)=(/9,10,7,9,7,8,10,2,7,12,7,4,2,4,7,0/)
    matCon1(126,:)=(/1,10,2,12,7,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(127,:)=(/8,9,1,8,1,7,4,12,1,12,7,1,0,0,0,0/)
    matCon1(128,:)=(/8,12,7,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(129,:)=(/8,7,12,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(130,:)=(/4,1,9,12,8,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(131,:)=(/1,2,10,12,8,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(132,:)=(/9,2,10,9,4,2,12,8,7,0,0,0,0,0,0,0/)
    matCon1(133,:)=(/11,2,3,7,12,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(134,:)=(/2,3,11,4,1,9,7,12,8,0,0,0,0,0,0,0/)
    matCon1(135,:)=(/3,10,1,3,11,10,7,12,8,0,0,0,0,0,0,0/)
    matCon1(136,:)=(/7,12,8,3,11,4,11,9,4,11,10,9,0,0,0,0/)
    matCon1(137,:)=(/8,3,4,7,3,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(138,:)=(/8,1,9,8,7,1,7,3,1,0,0,0,0,0,0,0/)
    matCon1(139,:)=(/3,8,7,3,4,8,1,2,10,0,0,0,0,0,0,0/)
    matCon1(140,:)=(/2,7,3,2,9,7,2,10,9,9,8,7,0,0,0,0/)
    matCon1(141,:)=(/11,8,7,11,2,8,2,4,8,0,0,0,0,0,0,0/)
    matCon1(142,:)=(/11,8,7,2,8,11,2,9,8,2,1,9,0,0,0,0/)
    matCon1(143,:)=(/1,4,8,1,8,11,1,11,10,7,11,8,0,0,0,0/)
    matCon1(144,:)=(/8,7,11,8,11,9,9,11,10,0,0,0,0,0,0,0/)
    matCon1(145,:)=(/7,9,5,12,9,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(146,:)=(/4,7,12,4,1,7,1,5,7,0,0,0,0,0,0,0/)
    matCon1(147,:)=(/9,7,12,9,5,7,10,1,2,0,0,0,0,0,0,0/)
    matCon1(148,:)=(/10,5,7,10,7,4,10,4,2,12,4,7,0,0,0,0/)
    matCon1(149,:)=(/7,9,5,7,12,9,3,11,2,0,0,0,0,0,0,0/)
    matCon1(150,:)=(/2,3,11,4,1,12,1,7,12,1,5,7,0,0,0,0/)
    matCon1(151,:)=(/5,12,9,5,7,12,1,3,10,3,11,10,0,0,0,0/)
    matCon1(152,:)=(/11,10,4,11,4,3,10,5,4,12,4,7,5,7,4,0/)
    matCon1(153,:)=(/9,3,4,9,5,3,5,7,3,0,0,0,0,0,0,0/)
    matCon1(154,:)=(/1,5,3,5,7,3,0,0,0,0,0,0,0,0,0,0/)
    matCon1(155,:)=(/2,10,1,3,4,5,3,5,7,5,4,9,0,0,0,0/)
    matCon1(156,:)=(/2,10,5,2,5,3,3,5,7,0,0,0,0,0,0,0/)
    matCon1(157,:)=(/9,2,4,9,7,2,9,5,7,7,11,2,0,0,0,0/)
    matCon1(158,:)=(/11,2,1,11,1,7,7,1,5,0,0,0,0,0,0,0/)
    matCon1(159,:)=(/5,7,4,5,4,9,7,11,4,1,4,10,11,10,4,0/)
    matCon1(160,:)=(/11,10,5,7,11,5,0,0,0,0,0,0,0,0,0,0/)
    matCon1(161,:)=(/5,10,6,8,7,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(162,:)=(/1,9,4,5,10,6,12,8,7,0,0,0,0,0,0,0/)
    matCon1(163,:)=(/6,1,2,6,5,1,8,7,12,0,0,0,0,0,0,0/)
    matCon1(164,:)=(/12,8,7,9,4,5,4,6,5,4,2,6,0,0,0,0/)
    matCon1(165,:)=(/10,6,5,11,2,3,8,7,12,0,0,0,0,0,0,0/)
    matCon1(166,:)=(/7,12,8,2,3,11,1,9,4,5,10,6,0,0,0,0/)
    matCon1(167,:)=(/8,7,12,6,5,11,5,3,11,5,1,3,0,0,0,0/)
    matCon1(168,:)=(/4,5,9,4,6,5,4,3,6,11,6,3,12,8,7,0/)
    matCon1(169,:)=(/8,3,4,8,7,3,6,5,10,0,0,0,0,0,0,0/)
    matCon1(170,:)=(/10,6,5,1,9,7,1,7,3,7,9,8,0,0,0,0/)
    matCon1(171,:)=(/4,7,3,4,8,7,2,6,1,6,5,1,0,0,0,0/)
    matCon1(172,:)=(/7,3,9,7,9,8,3,2,9,5,9,6,2,6,9,0/)
    matCon1(173,:)=(/10,6,5,11,2,7,2,8,7,2,4,8,0,0,0,0/)
    matCon1(174,:)=(/2,7,11,2,8,7,2,1,8,9,8,1,10,6,5,0/)
    matCon1(175,:)=(/5,1,11,5,11,6,1,4,11,7,11,8,4,8,11,0/)
    matCon1(176,:)=(/8,7,11,8,11,9,6,5,11,5,9,11,0,0,0,0/)
    matCon1(177,:)=(/7,10,6,7,12,10,12,9,10,0,0,0,0,0,0,0/)
    matCon1(178,:)=(/4,7,12,1,7,4,1,6,7,1,10,6,0,0,0,0/)
    matCon1(179,:)=(/1,12,9,1,6,12,1,2,6,6,7,12,0,0,0,0/)
    matCon1(180,:)=(/7,12,4,7,4,6,6,4,2,0,0,0,0,0,0,0/)
    matCon1(181,:)=(/2,3,11,10,6,12,10,12,9,12,6,7,0,0,0,0/)
    matCon1(182,:)=(/1,12,4,1,7,12,1,10,7,6,7,10,2,3,11,0/)
    matCon1(183,:)=(/12,9,6,12,6,7,9,1,6,11,6,3,1,3,6,0/)
    matCon1(184,:)=(/7,12,4,7,4,6,3,11,4,11,6,4,0,0,0,0/)
    matCon1(185,:)=(/6,9,10,6,3,9,6,7,3,4,9,3,0,0,0,0/)
    matCon1(186,:)=(/10,6,7,10,7,1,1,7,3,0,0,0,0,0,0,0/)
    matCon1(187,:)=(/2,6,9,2,9,1,6,7,9,4,9,3,7,3,9,0/)
    matCon1(188,:)=(/2,6,7,3,2,7,0,0,0,0,0,0,0,0,0,0/)
    matCon1(189,:)=(/2,4,7,2,7,11,4,9,7,6,7,10,9,10,7,0/)
    matCon1(190,:)=(/11,2,1,11,1,7,10,6,1,6,7,1,0,0,0,0/)
    matCon1(191,:)=(/1,4,9,6,7,11,0,0,0,0,0,0,0,0,0,0/)
    matCon1(192,:)=(/11,6,7,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(193,:)=(/12,6,11,8,6,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(194,:)=(/12,6,11,12,8,6,9,4,1,0,0,0,0,0,0,0/)
    matCon1(195,:)=(/6,12,8,6,11,12,2,10,1,0,0,0,0,0,0,0/)
    matCon1(196,:)=(/11,8,6,11,12,8,10,9,2,9,4,2,0,0,0,0/)
    matCon1(197,:)=(/12,2,3,12,8,2,8,6,2,0,0,0,0,0,0,0/)
    matCon1(198,:)=(/1,9,4,2,3,8,2,8,6,8,3,12,0,0,0,0/)
    matCon1(199,:)=(/10,8,6,10,3,8,10,1,3,3,12,8,0,0,0,0/)
    matCon1(200,:)=(/8,6,3,8,3,12,6,10,3,4,3,9,10,9,3,0/)
    matCon1(201,:)=(/3,6,11,3,4,6,4,8,6,0,0,0,0,0,0,0/)
    matCon1(202,:)=(/9,3,1,9,6,3,9,8,6,11,3,6,0,0,0,0/)
    matCon1(203,:)=(/10,1,2,6,11,4,6,4,8,4,11,3,0,0,0,0/)
    matCon1(204,:)=(/10,9,3,10,3,2,9,8,3,11,3,6,8,6,3,0/)
    matCon1(205,:)=(/2,4,6,4,8,6,0,0,0,0,0,0,0,0,0,0/)
    matCon1(206,:)=(/1,9,8,1,8,2,2,8,6,0,0,0,0,0,0,0/)
    matCon1(207,:)=(/10,1,4,10,4,6,6,4,8,0,0,0,0,0,0,0/)
    matCon1(208,:)=(/10,9,8,6,10,8,0,0,0,0,0,0,0,0,0,0/)
    matCon1(209,:)=(/6,9,5,6,11,9,11,12,9,0,0,0,0,0,0,0/)
    matCon1(210,:)=(/6,1,5,6,12,1,6,11,12,12,4,1,0,0,0,0/)
    matCon1(211,:)=(/1,2,10,9,5,11,9,11,12,11,5,6,0,0,0,0/)
    matCon1(212,:)=(/11,12,5,11,5,6,12,4,5,10,5,2,4,2,5,0/)
    matCon1(213,:)=(/3,6,2,3,9,6,3,12,9,5,6,9,0,0,0,0/)
    matCon1(214,:)=(/1,5,12,1,12,4,5,6,12,3,12,2,6,2,12,0/)
    matCon1(215,:)=(/1,3,6,1,6,10,3,12,6,5,6,9,12,9,6,0/)
    matCon1(216,:)=(/10,5,6,3,12,4,0,0,0,0,0,0,0,0,0,0/)
    matCon1(217,:)=(/3,6,11,4,6,3,4,5,6,4,9,5,0,0,0,0/)
    matCon1(218,:)=(/6,11,3,6,3,5,5,3,1,0,0,0,0,0,0,0/)
    matCon1(219,:)=(/4,11,3,4,6,11,4,9,6,5,6,9,1,2,10,0/)
    matCon1(220,:)=(/6,11,3,6,3,5,2,10,3,10,5,3,0,0,0,0/)
    matCon1(221,:)=(/9,5,6,9,6,4,4,6,2,0,0,0,0,0,0,0/)
    matCon1(222,:)=(/1,5,6,2,1,6,0,0,0,0,0,0,0,0,0,0/)
    matCon1(223,:)=(/9,5,6,9,6,4,10,1,6,1,4,6,0,0,0,0/)
    matCon1(224,:)=(/10,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(225,:)=(/5,12,8,5,10,12,10,11,12,0,0,0,0,0,0,0/)
    matCon1(226,:)=(/1,9,4,5,10,8,10,12,8,10,11,12,0,0,0,0/)
    matCon1(227,:)=(/2,11,12,2,12,5,2,5,1,8,5,12,0,0,0,0/)
    matCon1(228,:)=(/4,2,5,4,5,9,2,11,5,8,5,12,11,12,5,0/)
    matCon1(229,:)=(/5,12,8,10,12,5,10,3,12,10,2,3,0,0,0,0/)
    matCon1(230,:)=(/10,8,5,10,12,8,10,2,12,3,12,2,1,9,4,0/)
    matCon1(231,:)=(/12,8,5,12,5,3,3,5,1,0,0,0,0,0,0,0/)
    matCon1(232,:)=(/12,8,5,12,5,3,9,4,5,4,3,5,0,0,0,0/)
    matCon1(233,:)=(/3,10,11,3,8,10,3,4,8,8,5,10,0,0,0,0/)
    matCon1(234,:)=(/10,11,8,10,8,5,11,3,8,9,8,1,3,1,8,0/)
    matCon1(235,:)=(/4,8,11,4,11,3,8,5,11,2,11,1,5,1,11,0/)
    matCon1(236,:)=(/2,11,3,9,8,5,0,0,0,0,0,0,0,0,0,0/)
    matCon1(237,:)=(/5,10,2,5,2,8,8,2,4,0,0,0,0,0,0,0/)
    matCon1(238,:)=(/5,10,2,5,2,8,1,9,2,9,8,2,0,0,0,0/)
    matCon1(239,:)=(/5,1,4,8,5,4,0,0,0,0,0,0,0,0,0,0/)
    matCon1(240,:)=(/5,9,8,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(241,:)=(/10,11,9,11,12,9,0,0,0,0,0,0,0,0,0,0/)
    matCon1(242,:)=(/4,1,10,4,10,12,12,10,11,0,0,0,0,0,0,0/)
    matCon1(243,:)=(/1,2,11,1,11,9,9,11,12,0,0,0,0,0,0,0/)
    matCon1(244,:)=(/4,2,11,12,4,11,0,0,0,0,0,0,0,0,0,0/)
    matCon1(245,:)=(/2,3,12,2,12,10,10,12,9,0,0,0,0,0,0,0/)
    matCon1(246,:)=(/4,1,10,4,10,12,2,3,10,3,12,10,0,0,0,0/)
    matCon1(247,:)=(/1,3,12,9,1,12,0,0,0,0,0,0,0,0,0,0/)
    matCon1(248,:)=(/4,3,12,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(249,:)=(/3,4,9,3,9,11,11,9,10,0,0,0,0,0,0,0/)
    matCon1(250,:)=(/10,11,3,1,10,3,0,0,0,0,0,0,0,0,0,0/)
    matCon1(251,:)=(/3,4,9,3,9,11,1,2,9,2,11,9,0,0,0,0/)
    matCon1(252,:)=(/2,11,3,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(253,:)=(/2,4,9,10,2,9,0,0,0,0,0,0,0,0,0,0/)
    matCon1(254,:)=(/1,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(255,:)=(/1,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    matCon1(256,:)=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

  end subroutine getMatCons

  subroutine computeIsoVariable(solName, sps, isoVal)
    !
    !       computeIsoVar computes NODE centered values for the given
    !       solName variable. It is essentially equilivent to
    !       sotreSolInBuffer. It is assumed blockPointers are already
    !       set to the correct block.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use cgnsNames
    use flowVarRefState
    use inputPhysics
    use IOModule
    use utils, only : setPointers, terminate
    use flowUtils, only : computePTot
    implicit none
    !
    !      Subroutine arguments.
    character(len=*), intent(in)                   :: solName
    integer(kind=intType), intent(in) :: sps
    real(kind=realType), intent(in) :: isoVal
    !
    !      Local parameters
    !
    real(kind=realType), parameter :: plim   = 0.001_realType
    real(kind=realType), parameter :: rholim = 0.001_realType
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, ii, jj, kk, nn

    real(kind=realType) :: uuy, uuz, vvx, vvz, wwx, wwy, tmp
    real(kind=realType) :: vortx, vorty, vortz, a2, ptotInf, ptot, uova(3), gradP(3), a
    real(kind=realType), dimension(:, :, :), pointer :: fc, fn

    do nn=1,nDom
       call setPointers(nn, 1, sps)
       fc => flowDoms(nn, 1, sps)%fc
       fn => flowDoms(nn, 1, sps)%fn

       select case(solName)

       case (cgnsDensity)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,irho)
                enddo
             enddo
          enddo

       case (cgnsMomx)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivx)
                enddo
             enddo
          enddo

       case (cgnsMomy)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivy)
                enddo
             enddo
          enddo

       case (cgnsMomz)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,irho)*w(i,j,k,ivz)
                enddo
             enddo
          enddo

       case (cgnsEnergy)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,irhoE)
                enddo
             enddo
          enddo

       case (cgnsTurbSaNu,cgnsTurbK)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,itu1)
                enddo
             enddo
          enddo

       case (cgnsTurbOmega,cgnsTurbTau,cgnsTurbEpsilon)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,itu2)
                enddo
             enddo
          enddo

       case (cgnsTurbV2)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,itu3)
                enddo
             enddo
          enddo

       case (cgnsTurbF)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,itu4)
                enddo
             enddo
          enddo

       case (cgnsVelx)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivx)
                enddo
             enddo
          enddo

       case (cgnsVely)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivy)
                enddo
             enddo
          enddo

       case (cgnsVelz)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivz)
                enddo
             enddo
          enddo

       case (cgnsRelVelx)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivx)-s(i,j,k,1)
                enddo
             enddo
          enddo

       case (cgnsRelVely)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivy)-s(i,j,k,2)
                enddo
             enddo
          enddo

       case (cgnsRelVelz)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = w(i,j,k,ivz)-s(i,j,k,3)
                enddo
             enddo
          enddo

       case (cgnsPressure)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = p(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsTemp)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = p(i,j,k)/(RGas*w(i,j,k,irho))
                enddo
             enddo
          enddo

       case (cgnsCp)
          tmp = two/(gammaInf*pInf*MachCoef*MachCoef)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = tmp*(p(i,j,k) - pInf)
                enddo
             enddo
          enddo

       case (cgnsMach)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                        / max(w(i,j,k,irho),rholim)
                   tmp = (w(i,j,k,ivx)**2 + w(i,j,k,ivy)**2 &
                        +  w(i,j,k,ivz)**2)/a2
                   fc(i,j,k) = sqrt(max(zero,tmp))
                enddo
             enddo
          enddo

       case (cgnsRelMach)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   a2  = gamma(i,j,k)*max(p(i,j,k),plim) &
                        / max(w(i,j,k,irho),rholim)
                   tmp = ((w(i,j,k,ivx)-s(i,j,k,1))**2 +&
                        (w(i,j,k,ivy)-s(i,j,k,2))**2 &
                        +(w(i,j,k,ivz)-s(i,j,k,3))**2)/a2
                   fc(i,j,k) = sqrt(max(zero,tmp))
                enddo
             enddo
          enddo


       case (cgnsMachTurb)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   tmp = w(i,j,k,irho)*w(i,j,k,itu1) &
                        / (gamma(i,j,k)*max(p(i,j,k),plim))
                   fc(i,j,k) = sqrt(max(zero,tmp))
                enddo
             enddo
          enddo

       case (cgnsEddy)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = rev(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsEddyRatio)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = rev(i,j,k)/rlv(i,j,k)
                enddo
             enddo
          enddo

       case (cgNSWallDist)
          do k=1,ke
             kk = max(2_intType,k); kk = min(kl,kk)
             do j=1,je
                jj = max(2_intType,j); jj = min(jl,jj)
                do i=1,ie
                   ii = max(2_intType,i); ii = min(il,ii)
                   fc(i,j,k) = d2Wall(ii,jj,kk)
                enddo
             enddo
          enddo

       case (cgnsVortMagn)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   tmp = half/vol(i,j,k)
                   uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                        - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                        + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                        - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                        + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                        - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                   uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                        - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                        + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                        - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                        + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                        - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                   vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                        - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                        + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                        - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                        + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                        - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                   vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                        - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                        + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                        - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                        + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                        - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                   wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                        - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                        + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                        - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                        + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                        - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                   wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                        - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                        + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                        - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                        + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                        - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                   vortx = wwy - vvz; vorty = uuz - wwx; vortz = vvx - uuy

                   fc(i,j,k) = tmp*sqrt(vortx**2 + vorty**2 + vortz**2)
                enddo
             enddo
          enddo

       case (cgnsVortx)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   tmp = half/vol(i,j,k)
                   vvz = si(i,  j,k,3)*w(i+1,j,k,ivy) &
                        - si(i-1,j,k,3)*w(i-1,j,k,ivy) &
                        + sj(i,j,  k,3)*w(i,j+1,k,ivy) &
                        - sj(i,j-1,k,3)*w(i,j-1,k,ivy) &
                        + sk(i,j,k,  3)*w(i,j,k+1,ivy) &
                        - sk(i,j,k-1,3)*w(i,j,k-1,ivy)

                   wwy = si(i,  j,k,2)*w(i+1,j,k,ivz) &
                        - si(i-1,j,k,2)*w(i-1,j,k,ivz) &
                        + sj(i,j,  k,2)*w(i,j+1,k,ivz) &
                        - sj(i,j-1,k,2)*w(i,j-1,k,ivz) &
                        + sk(i,j,k,  2)*w(i,j,k+1,ivz) &
                        - sk(i,j,k-1,2)*w(i,j,k-1,ivz)

                   fc(i,j,k) = tmp*(wwy - vvz)
                enddo
             enddo
          enddo

       case (cgnsVorty)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   tmp = half/vol(i,j,k)
                   uuz = si(i,  j,k,3)*w(i+1,j,k,ivx) &
                        - si(i-1,j,k,3)*w(i-1,j,k,ivx) &
                        + sj(i,j,  k,3)*w(i,j+1,k,ivx) &
                        - sj(i,j-1,k,3)*w(i,j-1,k,ivx) &
                        + sk(i,j,k,  3)*w(i,j,k+1,ivx) &
                        - sk(i,j,k-1,3)*w(i,j,k-1,ivx)

                   wwx = si(i,  j,k,1)*w(i+1,j,k,ivz) &
                        - si(i-1,j,k,1)*w(i-1,j,k,ivz) &
                        + sj(i,j,  k,1)*w(i,j+1,k,ivz) &
                        - sj(i,j-1,k,1)*w(i,j-1,k,ivz) &
                        + sk(i,j,k,  1)*w(i,j,k+1,ivz) &
                        - sk(i,j,k-1,1)*w(i,j,k-1,ivz)

                   fc(i,j,k) = tmp*(uuz - wwx)
                enddo
             enddo
          enddo

       case (cgnsVortz)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   tmp = half/vol(i,j,k)
                   uuy = si(i,  j,k,2)*w(i+1,j,k,ivx) &
                        - si(i-1,j,k,2)*w(i-1,j,k,ivx) &
                        + sj(i,j,  k,2)*w(i,j+1,k,ivx) &
                        - sj(i,j-1,k,2)*w(i,j-1,k,ivx) &
                        + sk(i,j,k,  2)*w(i,j,k+1,ivx) &
                        - sk(i,j,k-1,2)*w(i,j,k-1,ivx)

                   vvx = si(i,  j,k,1)*w(i+1,j,k,ivy) &
                        - si(i-1,j,k,1)*w(i-1,j,k,ivy) &
                        + sj(i,j,  k,1)*w(i,j+1,k,ivy) &
                        - sj(i,j-1,k,1)*w(i,j-1,k,ivy) &
                        + sk(i,j,k,  1)*w(i,j,k+1,ivy) &
                        - sk(i,j,k-1,1)*w(i,j,k-1,ivy)

                   fc(i,j,k) = tmp*(vvx - uuy)
                enddo
             enddo
          enddo

       case (cgnsPtotloss)

          ! Compute the free stream total pressure.

          call computePtot(rhoInf, uInf, zero, zero, &
               pInf, ptotInf)
          ptotInf = one/ptotInf

          ! Loop over the cell centers and compute the
          ! total pressure loss.

          do k=1,ke
             do j=1,je
                do i=1,ie
                   call computePtot(w(i,j,k,irho), w(i,j,k,ivx), &
                        w(i,j,k,ivy),  w(i,j,k,ivz), &
                        p(i,j,k),      ptot)

                   fc(i,j,k) = one - ptot*ptotInf
                enddo
             enddo
          enddo

       case (cgnsResRho)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,irho)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResMomx)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,imx)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResMomy)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,imy)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResMomz)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,imz)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResRhoE)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,irhoE)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResNu,cgnsResK)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,itu1)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResOmega,cgnsResTau,cgnsResEpsilon)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,itu2)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResV2)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,itu3)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsResF)

          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = dw(i,j,k,itu4)/vol(i,j,k)
                enddo
             enddo
          enddo

       case (cgnsShock)

          do k=1,ke
             do j=1,je
                do i=1,ie

                   ! Here we compute U/a <dot> grad P / ||grad P||
                   ! Whre U is the velocity vector, a is the speed of
                   ! sound and P is the pressure.

                   ! U / a
                   a  = sqrt(gamma(i,j,k)*max(p(i,j,k),plim) &
                        / max(w(i,j,k,irho),rholim))

                   if (addGridVelocities) then
                      UovA = (/w(i,j,k,ivx)-s(i,j,k,1), &
                           w(i,j,k,ivy)-s(i,j,k,2), &
                           w(i,j,k,ivz)-s(i,j,k,3)/)/a
                   else
                      UovA = (/w(i,j,k,ivx),w(i,j,k,ivy), w(i,j,k,ivz)/)/a
                   end if
                   ! grad P / ||grad P||

                   gradP(1) = si(i,  j,k,1)*P(i+1,j,k) &
                        - si(i-1,j,k,1)*P(i-1,j,k) &
                        + sj(i,j,  k,1)*P(i,j+1,k) &
                        - sj(i,j-1,k,1)*P(i,j-1,k) &
                        + sk(i,j,k,  1)*P(i,j,k+1) &
                        - sk(i,j,k-1,1)*P(i,j,k-1)

                   gradP(2) = si(i,  j,k,2)*P(i+1,j,k) &
                        - si(i-1,j,k,2)*P(i-1,j,k) &
                        + sj(i,j,  k,2)*P(i,j+1,k) &
                        - sj(i,j-1,k,2)*P(i,j-1,k) &
                        + sk(i,j,k,  2)*P(i,j,k+1) &
                        - sk(i,j,k-1,2)*P(i,j,k-1)

                   gradP(3) = si(i,  j,k,3)*P(i+1,j,k) &
                        - si(i-1,j,k,3)*P(i-1,j,k) &
                        + sj(i,j,  k,3)*P(i,j+1,k) &
                        - sj(i,j-1,k,3)*P(i,j-1,k) &
                        + sk(i,j,k,  3)*P(i,j,k+1) &
                        - sk(i,j,k-1,3)*P(i,j,k-1)

                   ! Protect against divide by zero
                   gradP = gradP / sqrt(gradP(1)**2 + gradP(2)**2 + gradP(3)**2 + 1e-12)

                   ! Dot product
                   fc(i,j,k) = UovA(1)*gradP(1) + UovA(2)*gradP(2) + UovA(3)*gradP(3)
                end do
             end do
          end do

       case (cgnsBlank)
          do k=1,ke
             do j=1,je
                do i=1,ie
                   fc(i,j,k) = real(min(iblank(i,j,k),1_intType),realType)
                enddo
             enddo
          enddo

       case default
          call terminate("computeIsoVariable", &
               "This should not happen")

       end select

       ! We now create nodal values from the cell centered
       ! values. This was the reason for going from 1 to ie
       ! etc.

       do k=1,kl
          do j=1,jl
             do i=1,il
                fn(i,j,k) = eighth*( &
                     fc(i  , j  , k  ) + &
                     fc(i+1, j  , k  ) + &
                     fc(i  , j+1, k  ) + &
                     fc(i+1, j+1, k  ) + &
                     fc(i  , j  , k+1) + &
                     fc(i+1, j  , k+1) + &
                     fc(i  , j+1, k+1) + &
                     fc(i+1, j+1, k+1)) - isoVal
             end do
          end do
       end do
    end do
  end subroutine computeIsoVariable
end module writeCGNSSurface
