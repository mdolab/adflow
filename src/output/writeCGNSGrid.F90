module writeCGNSGrid

contains
  subroutine writeCGNSGridFile
    !
    !       writeCGNSGridFile and its subroutines write the CGNS grid
    !       file(s). Typically this is needed when the coordinates have
    !       changed due to moving parts, deformation or both.
    !
    use cgnsGrid
    use communication
    use IOModule
    use monitor
    use su_cgns
    use outputMod
    use communication
    use inputIteration
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer, dimension(cgnsNDom) :: cgnsZone

    integer(kind=intType) :: nn

    ! Determine the number and names of the grid files.
    ! Also set the pointers for IOVar needed for the general
    ! treatment of the IO.

    call gridFileNamesWrite

    ! Return immediately if no grids have to be written.

    if(nGridsToWrite == 0) return

    ! Allocate the memory for the fileIDs and the bases.

    allocate(fileIDs(nGridsToWrite), cgnsBases(nGridsToWrite), &
         stat=ierr)
    if(ierr /= 0)                         &
         call terminate("writeCGNSGridFile", &
         "Memory allocation failure for fileIDs and &
         &cgnsBases")

    ! Write a message that the grid file(s) are being written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "#"
       print "(a,a)", "# Writing grid file(s): ",trim(gridFileNames(1))
    endif

    ! All grid information is stored on all processors, with the
    ! exception of data which can vary in time or use a large
    ! amount of memory, including the coordinates and
    ! connectivities. Processor 0 writes this information first as
    ! a frame for each file.

    if(myID == 0) then
       do nn=1,nGridsToWrite
          call writeCGNSGridFrame(cgnsZone, nn)
       enddo
    endif

    ! Loop over the number of cgns blocks and write the coordinates
    ! one at a time to conserve memory.

    do nn=1,cgnsNDom
       call writeCoorCGNSZone(nn, cgnsZone(nn))
    enddo

    ! Check if the solution must be written in a different file.
    ! If so the files must be closed the memory for the fileIDs
    ! and bases must be released.

    testGridOnly: if(useLinksInCGNS .or. (.not. writeVolume)) then

       ! Processor 0 closes the files.

       if(myID == 0) then
          do nn=1,nGridsToWrite
             call cg_close_f(fileIDs(nn), ierr)
             if(ierr /= CG_OK)                    &
                  call terminate("writeCGNSGridFile", &
                  "Something wrong when calling cg_close_f")
          enddo
       endif


       ! Release the memory of fileIDs and cgnsBases.

       deallocate(fileIDs, cgnsBases, stat=ierr)
       if(ierr /= 0)                         &
            call terminate("writeCGNSGridFile", &
            "Deallocation failure for fileIDs and &
            &cgnsBases")

    endif testGridOnly

    ! Releases the memory of IOVar.

    deallocate(IOVar, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("writeCGNSGridFile", &
         "Deallocation failure for IOVar")

    ! Wait until all processors (especially processor 0) reach
    ! this point.

    call mpi_barrier(ADflow_comm_world, ierr)

    ! Write a message that the grid file has been written.
    ! Of course only processor 0 does this.

    if(myID == 0 .and. printIterations) then
       print "(a)", "# Grid file(s) written"
       print "(a)", "#"
    endif
  end subroutine writeCGNSGridFile

  subroutine gridFileNamesWrite
    !
    !       gridFileNamesWrite determines the names and number of grid
    !       files to be written. Furthermore, it sets the pointers for
    !       IOVar to make a general treatment of the writing possible.
    !
    use block
    use inputIO
    use inputPhysics
    use inputTimeSpectral
    use IOModule
    use iteration
    use monitor
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, kk, nAvail

    character(len=7) :: intString

    !       Determine the names and number of grid files to be written.
    !
    ! Determine the situation we are having here.

    select case (equationMode)

    case (steady)

       ! Steady state computation. Allocate the memory for the
       ! grid file names. Even if no file needs to be written the
       ! memory is still allocated because the name is always set.

       allocate(gridFileNames(1), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("gridFileNamesWrite", &
            "Memory allocation failure for grid &
            &file names")

       ! Set the number of grid files to be written to either 0 or 1
       ! and set the name accordingly. The name is always set,
       ! because it may be needed for a link in the solution file.

       if( writeGrid ) then
          nGridsToWrite    = 1
          gridFileNames(1) = newGridFile
       else
          nGridsToWrite    = 0
          gridFileNames(1) = gridFile
       endif

       !===============================================================

    case (unsteady)

       ! Unsteady computation. For a consistent restart for a
       ! deforming mesh computation nOldLevels grids must be written.
       ! First determine the number of available solutions.

       nAvail = timeStepUnsteady + nTimeStepsRestart + 1
       nAvail = min(nAvail,nOldLevels)

       ! Allocate the memory for the file names. Note that this is
       ! an upper boundary. It is possible that less files need
       ! to be written.

       allocate(gridFileNames(nAvail), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("gridFileNamesWrite", &
            "Memory allocation failure for &
            &gridFileNames")

       ! Set the names of the files.

       do nn=1,nAvail
          write(intString,"(i7)") timeStepUnsteady + &
               nTimeStepsRestart + 1 - nn
          intString = adjustl(intString)

          gridFileNames(nn) = trim(newGridFile)//"&
               &Timestep"//trim(intString)
       enddo

       ! Determine the number of grid files to be written.
       ! This depends on quite a few things.

       if( writeGrid ) then

          ! Initialize nGridsToWrite to 1. This may change
          ! when the mesh is deforming.

          nGridsToWrite = 1

          if( deforming_Grid ) then

             ! Grids deform during the computation. Check if the
             ! older grids must be written.

             do nn=1,(nAvail-1)
                if(.not. oldSolWritten(nn) ) then
                   nGridsToWrite = nGridsToWrite + 1
                   gridFileNames(nGridsToWrite) = gridFileNames(nn+1)
                endif
             enddo

          endif

       else

          ! No grids need to be written. Correct the grid file name
          ! to the original grid file.

          nGridsToWrite = 0
          gridFileNames(1) = gridFile

       endif

       !===============================================================

    case (timeSpectral)

       ! Time spectral computation. Allocate the file names.
       ! Again this is an upper bound.

       allocate(gridFileNames(nTimeIntervalsSpectral), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("gridFileNamesWrite", &
            "Memory allocation failure for &
            &gridFileNames")

       ! Set the names of the files.

       do nn=1,nTimeIntervalsSpectral
          write(intString,"(i7)") nn
          intString = adjustl(intString)

          gridFileNames(nn) = trim(newGridFile)//"&
               &Spectral"//trim(intString)
       enddo

       ! Set the number of grid files to be written.
       ! This depends on quite a few things.

       ! GKK The logic below is seriously flawed must never had
       ! been tested. If nGridsToWrite is set to 0 which is what
       ! it should be set at if the time spectral grid files are
       ! already written, then in writeCGNSGridFile, fileIDs are
       ! and CGNSbases are NOT allocated at ALL. Then when you try
       ! to write the volume grid, and index into fileIDs and
       ! CGNSbases, you're screwed. nGridsToWrite MUST ALWAYS be
       ! ntimeIntervalsSpectral regardless.


       if( writeGrid ) then

          ! Need some additional checks.

          if( deforming_Grid ) then

             ! Grids deform during the computation and thus
             ! they must be written.

             nGridsToWrite = nTimeIntervalsSpectral

          else if( timeSpectralGridsNotWritten ) then

             ! Grids do not deform, but the time spectral grids have
             ! not been written earlier. So write the grids and set
             ! timeSpectralGridsNotWritten to .false.

             nGridsToWrite = nTimeIntervalsSpectral
             timeSpectralGridsNotWritten = .false.

          else

             ! Although indicated that the grids must be written,
             ! this is not necessary, because they have already been
             ! written earlier and they have not changed.

             nGridsToWrite = nTimeIntervalsSpectral! 0

          endif

       else

          ! It is not needed to write the grid files.

          nGridsToWrite = 0

       endif

    end select
    !
    !       Determine whether or not to use links in CGNS.
    !

    if( writeGrid ) then

       ! Grid file(s) will be written. Compare the (base) names of the
       ! grid and solution files and set useLinksInCGNS accordingly.
       if(newGridFile == solFile) then
          useLinksInCGNS = .false.
       else
          useLinksInCGNS = .true.
       endif

    else

       ! Grid file(s) will not be written. Compare the (base) names of
       ! the original grid and solution files and set useLinksInCGNS
       ! accordingly.

       if(gridFile == solFile) then
          useLinksInCGNS = .false.
       else
          useLinksInCGNS = .true.
       endif

    endif

    !       Set the pointers for IOVar if grid files need to be written.
    !
    testGridsToWrite: if(nGridsToWrite > 0) then

       ! Allocate the memory for IOVar.

       allocate(IOVar(nDom,nGridsToWrite), stat=ierr)
       if(ierr /= 0)                          &
            call terminate("gridFileNamesWrite", &
            "Memory allocation failure for IOVar")

       ! Set the pointer w of IOVar to the correct coordinates.

       select case(equationMode)

       case (steady, timeSpectral)

          ! Steady state or time spectral computation. Simply set the
          ! pointers to the current coordinates.

          do nn=1,nDom
             do mm=1,nGridsToWrite
                IOVar(nn,mm)%pointerOffset = 0
                IOVar(nn,mm)%w => flowDoms(nn,1,mm)%x(1:,1:,1:,:)
             enddo
          enddo

          !=============================================================

       case (unsteady)

          ! Unsteady computation. First coordinates to be written
          ! are the current coordinates.

          do nn=1,nDom
             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)
          enddo

          ! It is possible (for a case with deforming meshes) that
          ! older coordinates need to be written as well.

          if( deforming_Grid ) then
             kk = 1
             do mm=1,(nAvail-1)
                if(.not. oldSolWritten(mm) ) then
                   kk = kk + 1
                   do nn=1,nDom
                      IOVar(nn,kk)%pointerOffset = 0
                      IOVar(nn,kk)%w => flowDoms(nn,1,1)%xOld(mm,1:,1:,1:,:)
                   enddo
                endif
             enddo
          endif

       end select

    endif testGridsToWrite

  end subroutine gridFileNamesWrite

  subroutine writeCGNSGridFrame(cgnsZone, ind)
    !
    !       writeCGNSGridFrame writes the framework for the grid file
    !       gridNames(ind) using the information stored in the module
    !       cgnsGrid. Basically all information but the coordinates is
    !       written by this routine.
    !
    use constants
    use cgnsGrid
    use su_cgns
    use outputMod
    use utils, only : terminate, setCGNSRealType
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)  :: ind
    integer, dimension(*), intent(out) :: cgnsZone

    !
    !      Local variables.
    !
    integer :: ierr, ii, jj, cgnsInd, cgnsBase
    integer(kind=cgsize_t), dimension(9)   :: sizes
    integer(kind=cgsize_t), dimension(3,2) :: zoneRange, donorRange
    integer, dimension(3)   :: transform
    integer(kind=cgsize_t), dimension(:,:), allocatable :: donorData

    integer(kind=intType) :: nn, mm, ll, i, j, k
    integer(kind=intType) :: s1, s2, s3

    real(kind=realType), dimension(3) :: rotCenter, rotRate, translation

    real(kind=realType) :: LRefInv

    character(len=maxStringLen)   :: errorMessage

    type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet

    ! Open the CGNS file for writing and check if it went okay.
    ! Store the file index afterwards.

    call cg_open_f(gridFileNames(ind), mode_write, cgnsInd, ierr)
    if(ierr /= CG_OK) then
       write(errorMessage,*) "File ", trim(gridfileNames(ind)), &
            " could not be opened by cgns &
            &for writing"
       call terminate("writeCGNSGridFrame", errorMessage)
    endif

    fileIDs(ind) = cgnsInd

    ! Create the base. Copy the cell and physical dimensions and
    ! store the base ID for this index.

    call cg_base_write_f(cgnsInd, cgnsBaseName, cgnsCelldim, &
         cgnsPhysdim, cgnsBase, ierr)
    if(ierr /= CG_OK)                     &
         call terminate("writeCGNSGridFrame", &
         "Something wrong when calling cg_base_write_f")

    cgnsBases(ind) = cgnsBase
    !
    !       Write the family info.
    !
    ! Loop over the number of families.

    familyLoop: do nn=1,cgnsNfamilies

       ! Create the family node.

       call cg_family_write_f(cgnsInd, cgnsBase, &
            cgnsFamilies(nn)%familyName, ii, ierr)
       if(ierr /= CG_OK)                     &
            call terminate("writeCGNSGridFrame", &
            "Something wrong when calling &
            &cg_family_write_f")

       ! Write the family BC, if this is present.

       if(cgnsFamilies(nn)%BCTypeCGNS /= Null) then

          call cg_fambc_write_f(cgnsInd, cgnsBase, ii,   &
               cgnsFamilies(nn)%bcName, &
               cgnsFamilies(nn)%BCTypeCGNS, jj, ierr)
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_fambc_write_f")

          ! If the boundary condition is UserDefined add the
          ! description what type of user defined BC it is.

          if(cgnsFamilies(nn)%BCTypeCGNS == UserDefined) then

             ! Ultimately you would like to create the
             ! UserDefinedData_t as a subnode of the family boundary
             ! condition node. However, at the moment CGNS does not
             ! allow this and therefore it is put one level higher.
             ! As only 1 boundary condition per family is allowed,
             ! this is not really problem.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                  "Family_t", ii, "end")
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling cg_goto_f")

             call cg_user_data_write_f(cgnsFamilies(nn)%userDefinedName, &
                  ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_user_data_write_f")
          endif

       endif
    enddo familyLoop
    !
    !       Write all the zone info, except the coordinates.
    !
    ! Loop over the number of zones in the original grid.

    zoneLoop: do nn=1,cgnsNDom

       ! Store the inverse of the scaling factor to meters.

       LRefInv = one/cgnsDoms(nn)%LRef

       ! Store the dimensions of the zone in sizes and create the zone.

       sizes(1) = cgnsDoms(nn)%il
       sizes(2) = cgnsDoms(nn)%jl
       sizes(3) = cgnsDoms(nn)%kl
       sizes(4) = cgnsDoms(nn)%nx
       sizes(5) = cgnsDoms(nn)%ny
       sizes(6) = cgnsDoms(nn)%nz
       sizes(7) = 0
       sizes(8) = 0
       sizes(9) = 0

       call cg_zone_write_f(cgnsInd, cgnsBase,                   &
            cgnsDoms(nn)%zoneName, sizes,        &
            cgnsDoms(nn)%zoneType, cgnsZone(nn), &
            ierr)
       if(ierr /= CG_OK)                     &
            call terminate("writeCGNSGridFrame", &
            "Something wrong when calling &
            &cg_zone_write_f")

       ! Go to the current zone. Needed when family and/or rotating
       ! frame info must be written.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", &
            cgnsZone(nn), "end")
       if(ierr /= CG_OK) &
            call terminate("writeCGNSGridFrame", &
            "Something wrong when calling cg_goto_f")

       ! Check if the zone belongs to a family. If so, write the
       ! family name.

       mm = cgnsDoms(nn)%familyID
       if(mm > 0) then

          call cg_famname_write_f(cgnsFamilies(mm)%familyName, ierr)
          if(ierr /= CG_OK) &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_famname_write_f")
       endif

       ! Write the rotating frame info, if the zone is rotating.

       if( cgnsDoms(nn)%rotatingFrameSpecified ) then

          ! Convert the rotation rate to degrees per second and store
          ! it in a single precision array.

          rotRate(1) = cgnsDoms(nn)%rotRate(1)*180.0_realType/pi
          rotRate(2) = cgnsDoms(nn)%rotRate(2)*180.0_realType/pi
          rotRate(3) = cgnsDoms(nn)%rotRate(3)*180.0_realType/pi

          ! Convert the rotation center to the original units
          ! and also in single precision.

          rotCenter(1) = LRefInv*cgnsDoms(nn)%rotCenter(1)
          rotCenter(2) = LRefInv*cgnsDoms(nn)%rotCenter(2)
          rotCenter(3) = LRefInv*cgnsDoms(nn)%rotCenter(3)

          ! Write the rotation rate and rotation center.

          call cg_rotating_write_f(real(rotRate,4), real(rotCenter,4), ierr)
          if(ierr /= CG_OK) &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_rotating_write_f")

          ! Write the units of the rotation rate of the
          ! rotating frame.

          call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",        &
               cgnsZone(nn), "RotatingCoordinates_t", 1, &
               "DataArray_t", 2, "end")
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling cg_goto_f")

          call cg_dataclass_write_f(Dimensional, ierr)
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_dataclass_write_f")

          call cg_units_write_f(Null, Null, Second, Null, &
               Degree, ierr)
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_units_write_f")
       endif

       ! Loop over all 1 to 1 connectivities of the block and
       ! write the data.

       loop1to1: do mm=1,cgnsDoms(nn)%n1to1

          ! Store the range of the subface in zoneRange and
          ! the range of the donor in donorRange.

          zoneRange(1,1) = cgnsDoms(nn)%conn1to1(mm)%iBeg
          zoneRange(2,1) = cgnsDoms(nn)%conn1to1(mm)%jBeg
          zoneRange(3,1) = cgnsDoms(nn)%conn1to1(mm)%kBeg

          zoneRange(1,2) = cgnsDoms(nn)%conn1to1(mm)%iEnd
          zoneRange(2,2) = cgnsDoms(nn)%conn1to1(mm)%jEnd
          zoneRange(3,2) = cgnsDoms(nn)%conn1to1(mm)%kEnd

          donorRange(1,1) = cgnsDoms(nn)%conn1to1(mm)%diBeg
          donorRange(2,1) = cgnsDoms(nn)%conn1to1(mm)%djBeg
          donorRange(3,1) = cgnsDoms(nn)%conn1to1(mm)%dkBeg

          donorRange(1,2) = cgnsDoms(nn)%conn1to1(mm)%diEnd
          donorRange(2,2) = cgnsDoms(nn)%conn1to1(mm)%djEnd
          donorRange(3,2) = cgnsDoms(nn)%conn1to1(mm)%dkEnd

          ! Check whether the subface is periodic or not.

          periodicTest: if( cgnsDoms(nn)%conn1to1(mm)%periodic ) then

             ! Subface is periodic. Due to the current limitations in
             ! cgns it is not possible to write this info as a 1 to 1
             ! subface and the general connectivity must be used.

             ! First allocate the memory for donorData.

             ll = (abs(donorRange(3,2) - donorRange(3,1)) + 1) &
                  * (abs(donorRange(2,2) - donorRange(2,1)) + 1) &
                  * (abs(donorRange(1,2) - donorRange(1,1)) + 1)

             allocate(donorData(3,ll), stat=ierr)
             if(ierr /= 0) &
                  call terminate("writeCGNSGridFrame", &
                  "Memory allocation failure for &
                  &donorData")

             ! Determine the step for the three directions of
             ! donorData and fill the array.

             s1 = 1; s2 = 1; s3 = 1
             if(donorRange(1,2) < donorRange(1,1)) s1 = -1
             if(donorRange(2,2) < donorRange(2,1)) s2 = -1
             if(donorRange(3,2) < donorRange(3,1)) s3 = -1

             ll = 0
             do k=donorRange(3,1),donorRange(3,2),s3
                do j=donorRange(2,1),donorRange(2,2),s2
                   do i=donorRange(1,1),donorRange(1,2),s1
                      ll = ll+1

                      donorData(1,ll) = i
                      donorData(2,ll) = j
                      donorData(3,ll) = k
                   enddo
                enddo
             enddo

             ! Write the general connectivity.

             ii = ll
             call cg_conn_write_f(cgnsInd, cgnsBase, cgnsZone(nn),       &
                  cgnsDoms(nn)%conn1to1(mm)%connectName, &
                  Vertex, Abutting1to1, PointRange, int(2, cgsize_t),   &
                  zoneRange,                             &
                  cgnsDoms(nn)%conn1to1(mm)%donorName,   &
                  Structured, PointListDonor, Integer,   &
                  int(ii, cgsize_t), donorData, jj, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_conn_write_f")

             ! Deallocate the memory of donorData again.

             deallocate(donorData, stat=ierr)
             if(ierr /= 0)                          &
                  call terminate("writeCGNSGridFrame", &
                  "Deallocation failure for donorData")

             ! Write the periodic info. First transform the rotation
             ! center and translation vector to the original coordinates,
             ! the angles to degrees and store everything in single
             ! precision arrays.

             rotCenter  = cgnsDoms(nn)%conn1to1(mm)%rotationCenter &
                  * LRefInv
             translation = cgnsDoms(nn)%conn1to1(mm)%translation &
                  * LRefInv
             rotRate    = cgnsDoms(nn)%conn1to1(mm)%rotationAngles &
                  * 180.0_realType/pi

             call cg_conn_periodic_write_f(cgnsInd, cgnsBase,           &
                  cgnsZone(nn), jj, real(rotCenter,4), &
                  real(rotRate,4), real(translation,4), ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_conn_periodic_write_f")

             ! Write the units of the periodic rotation.

             call cg_goto_f(cgnsInd, cgnsBase, ierr,         &
                  "Zone_t", cgnsZone(nn),          &
                  "ZoneGridConnectivity_t", 1,     &
                  "GridConnectivity_t", jj,        &
                  "GridConnectivityProperty_t", 1, &
                  "Periodic_t", 1, "DataArray_t", 2, "end")
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling cg_goto_f")

             call cg_dataclass_write_f(Dimensional, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_dataclass_write_f")

             call cg_units_write_f(Null, Null, Null, Null, &
                  Degree, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_units_write_f")

          else periodicTest

             ! Normal 1 to 1 subface. Set the elements for the
             ! abbreviation of the transformation matrix.

             transform(1) = cgnsDoms(nn)%conn1to1(mm)%l1
             transform(2) = cgnsDoms(nn)%conn1to1(mm)%l2
             transform(3) = cgnsDoms(nn)%conn1to1(mm)%l3

             ! Write the connectivity.

             call cg_1to1_write_f(cgnsInd, cgnsBase, cgnsZone(nn),       &
                  cgnsDoms(nn)%conn1to1(mm)%connectName, &
                  cgnsDoms(nn)%conn1to1(mm)%donorName,   &
                  zoneRange, donorRange, transform,      &
                  ii, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_1to1_write_f")

          endif periodicTest

       enddo loop1to1

       ! Loop over the boundary subfaces and write the data.

       loopBocos: do mm=1,cgnsDoms(nn)%nBocos

          ! Check if this is an actual face. If not, continue with
          ! the next face.

          if(.not. cgnsDoms(nn)%bocoInfo(mm)%actualFace) cycle

          ! Store the range of the subface in zoneRange.

          zoneRange(1,1) = cgnsDoms(nn)%bocoInfo(mm)%iBeg
          zoneRange(2,1) = cgnsDoms(nn)%bocoInfo(mm)%jBeg
          zoneRange(3,1) = cgnsDoms(nn)%bocoInfo(mm)%kBeg

          zoneRange(1,2) = cgnsDoms(nn)%bocoInfo(mm)%iEnd
          zoneRange(2,2) = cgnsDoms(nn)%bocoInfo(mm)%jEnd
          zoneRange(3,2) = cgnsDoms(nn)%bocoInfo(mm)%kEnd

          ! Write the boundary condition. As the preprocessing
          ! overwrites the BCType for a family specified BC, the
          ! boundary condition is constructed first and stored in jj.

          jj = cgnsDoms(nn)%bocoInfo(mm)%BCTypeCGNS
          ll = cgnsDoms(nn)%bocoInfo(mm)%familyID
          if(ll > 0) jj = FamilySpecified

          call cg_boco_write_f(cgnsInd, cgnsBase, cgnsZone(nn),    &
               cgnsDoms(nn)%bocoInfo(mm)%bocoName, &
               jj, PointRange, int(2, cgsize_t), zoneRange, ii, ierr)
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling &
               &cg_boco_write_f")

          ! Write the family name for the surface for all boundary
          ! conditions.
          call cg_goto_f(cgnsInd, cgnsBase, ierr, &
               "Zone_t", cgnsZone(nn),  &
               "ZoneBC_t", 1, "BC_t", ii, "end")
          if(ierr /= CG_OK)                     &
               call terminate("writeCGNSGridFrame", &
               "Something wrong when calling cg_goto_f")

          call cg_famname_write_f(cgnsDoms(nn)%bocoINFO(mm)%wallBCName, ierr)

          ! Write the family name if the boundary condition is
          ! specified per family.

          if(ll > 0) then

             ! Go to the current boundary condition and write
             ! the appropriate family name.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                  "Zone_t", cgnsZone(nn),  &
                  "ZoneBC_t", 1, "BC_t", ii, "end")
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling cg_goto_f")

             call cg_famname_write_f(cgnsFamilies(ll)%familyName, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_famname_write_f")
          endif

          ! If the boundary condition is UserDefined, write the
          ! description of what type of user defined BC.

          if(jj == UserDefined) then

             ! Go to the current boundary condition and write
             ! the appropriate data.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                  "Zone_t", cgnsZone(nn),  &
                  "ZoneBC_t", 1, "BC_t", ii, "end")
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling cg_goto_f")

             call cg_user_data_write_f(    &
                  cgnsDoms(nn)%bocoInfo(mm)%userDefinedName, ierr)
             if(ierr /= CG_OK)                     &
                  call terminate("writeCGNSGridFrame", &
                  "Something wrong when calling &
                  &cg_user_data_write_f")
          endif

          ! If this boundary condition has allocated memory for data
          ! sets, write them.

          if( cgnsDoms(nn)%bocoInfo(mm)%dataSetAllocated ) then

             ! Set the pointer for the data sets to make the code
             ! more readable.

             dataSet => cgnsDoms(nn)%bocoInfo(mm)%dataSet

             ! Loop over the number of data sets for this boundary face.

             do ll=1,cgnsDoms(nn)%bocoInfo(mm)%nDataSet

                ! Create the bc dataset node.

                call cg_dataset_write_f(cgnsInd, cgnsBase,       &
                     cgnsZone(nn), ii,        &
                     dataSet(ll)%datasetName, &
                     dataSet(ll)%BCType, jj, ierr)
                if(ierr /= CG_OK)                     &
                     call terminate("writeCGNSGridFrame", &
                     "Something wrong when calling &
                     &cg_dataset_write_f")

                ! Write the Dirichlet and Neumann boundary condition
                ! data sets if present.

                call writeBcdataArrays(dataSet(ll)%ndirichletArrays, &
                     dataSet(ll)%dirichletArrays,  &
                     Dirichlet)

                call writeBcdataArrays(dataSet(ll)%nneumannArrays, &
                     dataSet(ll)%neumannArrays,  &
                     Neumann)
             enddo
          endif

       enddo loopBocos
    enddo zoneLoop

    !=================================================================

  contains

    !===============================================================

    subroutine writeBcdataArrays(narr, arr, DirNeu)
      !
      !         writeBcdataArrays writes the given bc data set arrays,
      !         either of the dirichlet or neumann type, to the correct
      !         position in the CGNS file.
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      integer, intent(in) :: DirNeu

      integer(kind=intType), intent(in) :: narr
      type(cgnsBcdataArray), pointer, dimension(:) :: arr
      !
      !        Local variables.
      !
      integer :: ierr
      integer :: realTypeCGNS

      integer(kind=intType) :: i, j, kk

      real(kind=cgnsRealType), dimension(:), allocatable :: tmp

      ! Return immediately if narr == 0, i.e. if there is nothing
      ! to write.

      if(narr == 0) return

      ! Set the cgns real type.

      realTypeCGNS = setCGNSRealType()

      ! Create the BCData node.

      call cg_bcdata_write_f(cgnsInd, cgnsBase, cgnsZone(nn), &
           ii, jj, DirNeu, ierr)
      if(ierr /= CG_OK)                    &
           call terminate("writeBcdataArrays", &
           "Something wrong when calling &
           &cg_bcdata_write_f")

      ! Loop over the number of data arrays.

      loopDataArrays: do kk=1,narr

         ! Go to the main node of the data arrays.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",       &
              cgnsZone(nn), "ZoneBC_t", 1, "BC_t", ii, &
              "BCDataSet_t", jj, "BCData_t", DirNeu, "end")
         if(ierr /= CG_OK)                    &
              call terminate("writeBcdataArrays", &
              "Something wrong when calling cg_goto_f")

         ! Determine the total size of the prescribed data,
         ! allocate the memory for tmp and copy the data into it.

         j = arr(kk)%dataDim(1)
         do i=2,arr(kk)%nDimensions
            j = j*arr(kk)%dataDim(i)
         enddo

         allocate(tmp(j), stat=ierr)
         if(ierr /= 0)                         &
              call terminate("writeBcdataArrays", &
              "Memory allocation failure for tmp")

         tmp = arr(kk)%dataArr

         ! Write the data array and release the memory of tmp
         ! afterwards.

         call cg_array_write_f(arr(kk)%arrayName,   realTypeCGNS,    &
              arr(kk)%nDimensions, int(arr(kk)%dataDim, cgsize_t), &
              tmp, ierr)
         if(ierr /= CG_OK)                    &
              call terminate("writeBcdataArrays", &
              "Something wrong when calling &
              &cg_array_write_f")

         deallocate(tmp, stat=ierr)
         if(ierr /= 0)                         &
              call terminate("writeBcdataArrays", &
              "Deallocation failure for tmp")

         ! Write the dimensional info for this array.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",       &
              cgnsZone(nn), "ZoneBC_t", 1, "BC_t", ii, &
              "BCDataSet_t", jj, "BCData_t", DirNeu,   &
              "DataArray_t", kk, "end")
         if(ierr /= CG_OK)                    &
              call terminate("writeBcdataArrays", &
              "Something wrong when calling cg_goto_f")

         call cg_dataclass_write_f(Dimensional, ierr)
         if(ierr /= CG_OK)                    &
              call terminate("writeBcdataArrays", &
              "Something wrong when calling &
              &cg_dataclass_write_f")

         call cg_units_write_f(arr(kk)%mass,  arr(kk)%len,  &
              arr(kk)%time,  arr(kk)%temp, &
              arr(kk)%angle, ierr)
         if(ierr /= CG_OK)                    &
              call terminate("writeBcdataArrays", &
              "Something wrong when calling &
              &cg_units_write_f")

      enddo loopDataArrays

    end subroutine writeBcdataArrays
  end subroutine writeCGNSGridFrame

  subroutine writeCoorCGNSZone(zone, cgnsZone)
    !
    !       writeCoorCGNSZone writes the coordinates of the given zone
    !       to the cgns file(s).
    !
    use constants
    use block
    use cgnsGrid
    use cgnsNames
    use communication
    use inputIO
    use su_cgns
    use outputMod
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in) :: cgnsZone

    integer(kind=intType), intent(in) :: zone

    !
    !      Local variables.
    !
    integer :: ierr, tmp
    integer :: bufSize, realTypeCGNS, cgnsBase, fileInd

    integer, dimension(mpi_status_size) :: mpiStatus

    integer, dimension(:), allocatable :: proc

    integer(kind=intType) :: i, j, nn, mm, ll, ind
    integer(kind=intType) :: nBlocks, nSubBlocks, offset
    integer(kind=intType) :: sizeCGNSWriteType

    integer(kind=intType), dimension(6)     :: ii
    integer(kind=intType), dimension(nProc) :: nMessages

    integer(kind=intType), dimension(:,:,:), allocatable :: subRanges

    real(kind=realType), dimension(:), allocatable :: buffer

    real(kind=4), dimension(:) , allocatable :: coor4
    real(kind=8), dimension(:) , allocatable :: coor8

    character(len=maxCGNSNameLen), dimension(3) :: coorNames

    ! Store the number of local blocks and the offset in
    ! blocksCGNSblock for this zone a bit easier.

    offset  = nBlocksCGNSblock(zone-1)
    nBlocks = nBlocksCGNSblock(zone) - offset

    ! Determine the amount of block parts each processor will send to
    ! processor 0.

    call mpi_gather(nBlocks, 1, adflow_integer, nMessages, 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

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

          call mpi_recv(ii, 6, adflow_integer, proc(i), proc(i), &
               ADflow_comm_world, mpiStatus, ierr)

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
       ll = cgnsDoms(zone)%il * cgnsDoms(zone)%jl * cgnsDoms(zone)%kl

       select case (precisionGrid)
       case (precisionSingle)
          allocate(coor4(ll), coor8(0), stat=ierr)

       case (precisionDouble)
          allocate(coor8(ll), coor4(0), stat=ierr)
       end select
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
                   call copyDataBufSinglePrecision(&
                        coor4, buffer,  &
                        1_intType,         &
                        1_intType,         &
                        1_intType,         &
                        cgnsDoms(zone)%il, &
                        cgnsDoms(zone)%jl, &
                        cgnsDoms(zone)%kl, &
                        subRanges(1,1,mm))
                case (precisionDouble)
                   call copyDataBufDoublePrecision(&
                        coor8, buffer, &
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

                call mpi_recv(buffer, bufSize, adflow_real, proc(mm), &
                     proc(mm)+1, ADflow_comm_world, mpiStatus, ierr)

                select case (precisionGrid)
                case (precisionSingle)
                   call copyDataBufSinglePrecision(coor4, buffer,      &
                        1_intType,         &
                        1_intType,         &
                        1_intType,         &
                        cgnsDoms(zone)%il, &
                        cgnsDoms(zone)%jl, &
                        cgnsDoms(zone)%kl, &
                        subRanges(1,1,mm))
                case (precisionDouble)
                   call copyDataBufDoublePrecision(coor8, buffer,      &
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
             select case (precisionGrid)
             case (precisionSingle)
                call cg_coord_write_f(fileInd, cgnsBase, cgnsZone, &
                     realSingle, coorNames(nn), coor4, tmp, ierr)
             case (precisionDouble)
                call cg_coord_write_f(fileInd, cgnsBase, cgnsZone, &
                     realDouble, coorNames(nn), coor8, tmp, ierr)
             end select

             if(ierr /= CG_OK)                    &
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
                if(ierr /= CG_OK)                    &
                     call terminate("writeCoorCGNSZone", &
                     "Something wrong when calling cg_goto_f")

                ! Write the units.

                call cg_units_write_f(cgnsDoms(zone)%mass, &
                     cgnsDoms(zone)%len,  &
                     cgnsDoms(zone)%time, &
                     cgnsDoms(zone)%temp, &
                     cgnsDoms(zone)%angle, ierr)
                if(ierr /= CG_OK)                    &
                     call terminate("writeCoorCGNSZone", &
                     "Something wrong when calling &
                     &cg_units_write_f")
             endif

          enddo coorLoopRoot

       enddo gridLoopRoot

       ! Deallocate the memory which is only allocated on the
       ! root processor.

       deallocate(subRanges, proc, coor4, coor8, stat=ierr)
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

          call mpi_send(ii, 6, adflow_integer, 0, myID, &
               ADflow_comm_world, ierr)

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

                call mpi_send(buffer, tmp, adflow_real, 0, myID+1, &
                     ADflow_comm_world, ierr)
             enddo

          enddo coorLoopOthers
       enddo gridLoopOthers

    endif rootproc

    ! Release the memory of buffer.

    deallocate(buffer, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("writeCoorCGNSZone", &
         "Deallocation failure for buffer")
  end subroutine writeCoorCGNSZone

  !      ==================================================================

  subroutine storeCoorInBuffer(buffer, zone, ind, coorID, &
       blockID, nn)
    !
    !       storeCoorInBuffer stores the given coordinate for the given
    !       blockID. The total size of the buffer is returned in nn.
    !
    use constants
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
end module writeCGNSGrid
