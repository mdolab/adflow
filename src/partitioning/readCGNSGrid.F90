module readCGNSGrid

contains

  subroutine readBlockSizes
    !
    !       readBlockSizes reads the number of blocks and their size
    !       from the given grid file. The data is stored in the module
    !       cgnsGrid.
    !       If multiple grids need to be read for a consistent restart, it
    !       is checked that the number of blocks and the block sizes are
    !       identical.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsFamilies, &
         cgnsCellDim, cgnsPhysDim, cgnsDomsd, cgnsBaseName, &
         cgnsNFamilies
    use communication, only : myid, adflow_comm_world
    use inputPhysics, only: equations, equationMode
    use iteration, only : changing_grid, deforming_grid
    use partitionMod, only: fileIds, gridFiles, nGridsRead
    use utils, only : terminate, nullifyCGNSDomPointers
    use sorting, only : bsearchStrings, qsortStrings
    implicit none
    !
    !      Local variables
    !
    integer :: cgnsInd, cgnsNbases, cgnsBase, cgnsNzones
    integer :: i, nZone
    integer :: ierr
    integer :: nDoubleBoundFaces

    integer(kind=intType) :: ii, nn

    integer(kind=intType), dimension(:), allocatable :: famID

    character(len=2*maxStringLen) :: errorMessage
    character(len=7)              :: integerString

    character(len=maxCGNSNameLen), dimension(:), allocatable :: &
         sortedFamName

    logical :: noUnits

    ! Open the cgns files for reading and check if it went okay.
    ! Later on it is assumed that the 1st file stored in gridFiles
    ! is the "master".

    do nn=1,nGridsRead
       call cg_open_f(gridFiles(nn), mode_read, fileIDs(nn), ierr)
       if(ierr /= CG_OK) then
          write(errorMessage,*) "File ", trim(gridFiles(nn)), &
               " could not be opened for reading"
          call terminate("readBlockSizes", errorMessage)
       endif
    enddo

    cgnsInd = fileIDs(1)

    ! Determine the number of bases in the cgns file.
    ! This must be at least 1.

    call cg_nBases_f(cgnsInd, cgnsNbases, ierr)
    if(ierr /= CG_OK)                 &
         call terminate("readBlockSizes", &
         "Something wrong when calling cg_nBases_f")

    if(cgnsNbases < 1) then
       write(errorMessage,*) "CGNS file ", trim(gridFiles(1)), &
            " does not contain a base"
       call terminate("readBlockSizes", errorMessage)
    endif

    ! Set cgnsBase explicitly to 1, because for the reading
    ! of the connectivity and boundary info only the info
    ! of the first base is needed. Here it is assumed that
    ! base 1 is the primary base.

    cgnsBase = 1

    ! Read the cell and physical dimensions as well as the name for
    ! this base.

    call cg_base_read_f(cgnsInd, cgnsBase, cgnsBasename, &
         cgnsCelldim, cgnsPhysdim, ierr)
    if(ierr /= CG_OK)                 &
         call terminate("readBlockSizes", &
         "Something wrong when calling cg_base_read_f")

    ! Check the cell and physical dimensions. Both must be 3 for
    ! this code to work.

    if(cgnsCelldim /= 3 .or. cgnsPhysdim /= 3) then
       write(errorMessage,100) cgnsCelldim, cgnsPhysdim
100    format("Both the number of cell and physical dimensions should &
            &be 3, not",1X,I1,1X,"and",1X,I1)
       call terminate("readBlockSizes", errorMessage)
    endif

    ! Read the family info, if present.

    call readFamilyInfo(cgnsInd, cgnsBase)

    ! Allocate the memory for famID and sortedFamName.

    allocate(famID(cgnsNFamilies), sortedFamName(cgnsNFamilies), &
         stat=ierr)
    if(ierr /= 0)                      &
         call terminate("readBlockSizes", &
         "Memory allocation failure for sorted family names")

    ! Determine the sorted version of the family names and determine
    ! the entries in the original numbering. FamID is initialized
    ! to -1, which serves as a check later on.

    do i=1,cgnsNFamilies
       sortedFamName(i) = cgnsFamilies(i)%familyName
       famID(i)         = -1
    enddo

    nn = cgnsNFamilies
    call qsortStrings(sortedFamName, nn)

    do i=1,cgnsNFamilies
       ii = bsearchStrings(cgnsFamilies(i)%familyName, sortedFamName)

       if( debug ) then
          if(ii == 0)                        &
               call terminate("readBlockSizes", &
               "Family name not found in sorted &
               &family names.")
       endif

       ! Check if the family name is not already taken. If this is the
       ! case the grid file is not valid.

       if(famID(ii) /= -1)                &
            call terminate("readBlockSizes", &
            "Error occurs when two identical family names &
            &are present")

       famID(ii) = i
    enddo

    ! Initialize whether the grid is changing to whether it is
    ! deforming

    changing_Grid = deforming_Grid


    ! Determine the number of zones/blocks in the grid. Note that the
    ! reading is done using cgnsNzones (an integer variable) which is
    ! then copied into cgnsNDom (an integer(kind=intType) variable).
    ! cgnsNdom cannot be used directly, because a type mismatch may
    ! occur.

    call cg_nZones_f(cgnsInd, cgnsBase, cgnsNzones, ierr)
    if(ierr /= CG_OK)                 &
         call terminate("readBlockSizes", &
         "Something wrong when calling cg_nZones_f")
    cgnsNDom = cgnsNzones

    ! Check if the number of zones for all the grid to be read
    ! are identical.

    do nn=2,nGridsRead
       call cg_nZones_f(fileIDs(nn), cgnsBase, cgnsNzones, ierr)
       if(ierr /= CG_OK)                 &
            call terminate("readBlockSizes", &
            "Something wrong when calling cg_nZones_f")

       if(cgnsNzones /= cgnsNDom) then
          write(errorMessage,*) "File ", trim(gridFiles(nn)), &
               ": Different number of blocks than&
               & in file ", trim(gridFiles(1))
          if(myID == 0) call terminate("readBlockSizes", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)
       endif
    enddo

    ! Allocate the memory for cgnsDoms, the derived data type to
    ! store the info of the cgns grid, and nullify its pointers.

    allocate(cgnsDoms(cgnsNDom), stat=ierr)
    !and a second variable for the forwardAD
    allocate(cgnsDomsd(cgnsNDom), stat=ierr)
    if(ierr /= 0)                      &
         call terminate("readBlockSizes", &
         "Memory allocation failure for cgnsDoms")

    do nn=1,cgnsNDom
       call nullifyCGNSDomPointers(nn)
    enddo

    ! Some initializations before reading the zone info.

    nDoubleBoundFaces = 0
    noUnits           = .false.

    ! Loop over the number of zones.

    zones: do nZone=1, cgnsNDom

       ! Read the zone info

       call readZoneInfo(cgnsBase, nZone, sortedFamName, &
            famID, noUnits)

       ! Count the number of each connectivity for the supported
       ! types (needed for memory allocation and to add up 1-to-1)

       call countConnectivities(cgnsInd, cgnsBase, nZone)

       ! For this zone, read the 1-to-1 block connectivity, the
       ! general connectivities, and the boundary conditions.

       call read1to1Conn(cgnsInd, cgnsBase, nZone)
       call readGeneralConn(cgnsInd, cgnsBase, nZone)
       call readBocos(cgnsInd, cgnsBase, nZone, &
            ndoubleBoundFaces, sortedFamName, famId)

    enddo zones

    ! Release the memory of sortedFamName and famID.

    deallocate(sortedFamName, famID, stat=ierr)
    if(ierr /= 0)                      &
         call terminate("readBlockSizes", &
         "Deallocation error for sortedFamName &
         &and famID")

    ! If there are double boundary faces, print a warning.

    if(myID == 0 .and. nDoubleBoundFaces > 0) then
       write(integerString,"(i6)") nDoubleBoundFaces
       integerString = adjustl(integerString)

       print "(a)", "#"
       print "(a)", "#                      Warning"
       print 110, trim(integerString)
       print "(a)", "# Block connectivity is kept, boundary info &
            &is neglected."
       print "(a)", "#"
110    format("# ",a, " double boundary faces found.")

    endif

    ! If the units could not be determined, print a warning if
    ! the viscous or unsteady equations are solved.

    if(myID == 0 .and. noUnits) then

       if(equations==NSEquations .or. equations==RANSEquations .or. equationMode==unsteady)then

          print "(a)", "#"
          print "(a)", "#                      Warning"
          print "(a)", "# Conversion factor from grid units to &
               &meter not specified and some blocks"
          print "(a)", "# do not have units. It is assumed that &
               &the grid is given in meters."
          print "(a)", "#"

       endif
    endif


  end subroutine readBlockSizes
  subroutine readFamilyInfo(cgnsInd, cgnsBase)
    !
    !       readFamilyInfo determines the number of families in the
    !       given base of the cgns grid and determines their possible
    !       boundary condition, including some user defined ones.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsFamilies, cgnsNFamilies
    use communication, only : myid, adflow_comm_world
    use utils, only: terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase
    !
    !      Local variables.
    !
    integer :: nn, bc, nFamBC, nGeo, nUserData, ierr

    character(len=maxStringLen) :: errorMessage

    ! Determine the number of families in the given base.

    call cg_nfamilies_f(cgnsInd, cgnsBase, nn, ierr)
    if(ierr /= CG_OK)                 &
         call terminate("readFamilyInfo", &
         "Something wrong when calling cg_nfamilies_f")
    cgnsNFamilies = nn

    ! Allocate the memory for cgnsFamilies.

    allocate(cgnsFamilies(nn), stat=ierr)
    if(ierr /= 0)                      &
         call terminate("readFamilyInfo", &
         "Memory allocation failure for cgnsFamilies")

    ! Loop over the number of families and read the info.

    nFam: do nn=1,cgnsNFamilies

       ! Initialize slidingID to 0 to indicate that this family does
       ! not belong to a sliding mesh interface. Idem for the
       ! bleedRegionID.

       cgnsFamilies(nn)%slidingID   = 0
       cgnsFamilies(nn)%bleedRegionID = 0

       ! Initialize the logical to monitor the mass flow to .false.

       cgnsFamilies(nn)%monitorMassFlow = .false.

       ! Nullify the pointer for the prescribed boundary data.

       nullify(cgnsFamilies(nn)%dataSet)

       ! Read the family name and the number of boundary conditions
       ! specified.

       call cg_family_read_f(cgnsInd, cgnsBase, nn,       &
            cgnsFamilies(nn)%familyName, &
            nFamBC, nGeo, ierr)

       if(ierr /= CG_OK)               &
            call terminate("readFamilyInfo", &
            "Something wrong when calling cg_family_read_f")

       ! Determine the boundary condition for this family, if specified.

       select case (nFamBC)

       case (0)
          cgnsFamilies(nn)%BCTypeCGNS = Null
          cgnsFamilies(nn)%BCType     = BCNull
          cgnsFamilies(nn)%bcName     = ""

          !=============================================================

       case (1)
          bc = 1
          call cg_fambc_read_f(cgnsInd, cgnsBase, nn, bc, &
               cgnsFamilies(nn)%bcName,   &
               cgnsFamilies(nn)%BCTypeCGNS, ierr)
          if(ierr /= CG_OK)                 &
               call terminate("readFamilyInfo", &
               "Something wrong when calling &
               &cg_fambc_read_f")

          ! If this is a user defined boundary condition it must
          ! contain more information to determine the internally
          ! used BC.

          testUserDefined: if(cgnsFamilies(nn)%BCTypeCGNS == &
               UserDefined) then

             ! Move to the family and determine the number of
             ! user defined data nodes.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, &
                  "Family_t", nn, "end")
             if(ierr /= CG_OK)                 &
                  call terminate("readFamilyInfo", &
                  "Something wrong when calling cg_goto_f")

             call cg_nuser_data_f(nUserData, ierr)
             if(ierr /= CG_OK)                 &
                  call terminate("readFamilyInfo", &
                  "Something wrong when calling &
                  &cg_nuser_data_f")

             ! nUserData should be 1. Check this.

             if(nUserData /= 1) then
                write(errorMessage,101) trim(cgnsFamilies(nn)%familyName)
                if(myID == 0) &
                     call terminate("readFamilyInfo", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)
             endif

             ! Read the name of the user defined data node.

             call cg_user_data_read_f(nUserData,                        &
                  cgnsFamilies(nn)%userDefinedName, &
                  ierr)
             if(ierr /= CG_OK)                 &
                  call terminate("readFamilyInfo", &
                  "Something wrong when calling &
                  &cg_user_data_read_f")

          else testUserDefined

             ! Set the user defined name to an empty string.

             cgnsFamilies(nn)%userDefinedName = ""

          endif testUserDefined

          ! Determine the internal BC type from the CGNS type and
          ! possibly the user defined name.

          cgnsFamilies(nn)%BCType = &
               internalBC(cgnsFamilies(nn)%BCTypeCGNS, &
               cgnsFamilies(nn)%userDefinedName)

          !=============================================================

       case default
          write(errorMessage,201) trim(cgnsFamilies(nn)%familyName)
          if(myID == 0) &
               call terminate("readFamilyInfo", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)

       end select

    enddo nFam

    ! Format statements.

101 format("Family",1x,a,": Need 1 UserDefinedData_t node for &
         &user defined boundary condition")
102 format("Family",1x,a,": Unknown user defined boundary &
         &condition",1x,a)
201 format("Family",1x,a,": More than 1 boundary condition specified")


  end subroutine readFamilyInfo

  subroutine readZoneInfo(cgnsBase, nZone, sortedFamName, &
       famID, noUnits)
    !
    !       readZoneInfo reads the general information, like zone type
    !       and physical dimensions, for the given zone/block.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsFamilies, cgnsNFamilies
    use communication, only : adflow_comm_world, myid
    use flowVarRefState, only : LRefSpecified, LRef
    use iteration, only : changing_grid
    use partitionMod, only : nGridsRead, fileIDs, gridFiles
    use utils, only : terminate, siAngle, siLen, siAngle
    use sorting, only: bsearchStrings
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsBase, nZone
    character(len=*), dimension(:), intent(in) :: sortedFamName
    integer(kind=intType), dimension(:), intent(in) :: famID

    logical, intent(inout) :: noUnits

    !
    !      Local variables
    !
    integer :: cgnsInd
    integer :: i, ierr, nCoords
    integer :: mass, len, time, temp, angle

    integer(kind=cgsize_t), dimension(9) :: sizesBlock

    integer(kind=intType) :: ii, nn

    real(kind=realType), dimension(3) :: rotCenter, rotRate

    real(kind=realType) :: mult, trans

    character(len=maxCGNSNameLen) :: familyName
    character(len=maxStringLen)   :: errorMessage

    logical :: overwrite


    ! Set the cgns ID for the "master" file and read the size
    ! of the block as well as the zone name.

    cgnsInd = fileIDs(1)

    call cg_zone_read_f(cgnsInd, cgnsBase, nZone, &
         cgnsDoms(nZone)%zoneName, sizesBlock, ierr)
    if(ierr /= CG_OK)               &
         call terminate("readZoneInfo", &
         "Something wrong when calling cg_nZones_f")

    ! Check the zone type.

    call cg_zone_type_f(cgnsInd, cgnsBase, nZone, &
         cgnsDoms(nZone)%zonetype, ierr)
    if(ierr /= CG_OK)               &
         call terminate("readZoneInfo", &
         "Something wrong when calling cg_zone_type_f")

    if(cgnsDoms(nZone)%zonetype /= Structured) then
       write(errorMessage,*) "Zone ",                         &
            trim(cgnsDoms(nZone)%zoneName), &
            " of the grid file is not structured"
       if(myID == 0) call terminate("readZoneInfo", errorMessage)
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Set the values for the number of nodes and cells in i, j and
    ! k-direction.

    cgnsDoms(nZone)%il = sizesBlock(1)
    cgnsDoms(nZone)%jl = sizesBlock(2)
    cgnsDoms(nZone)%kl = sizesBlock(3)

    cgnsDoms(nZone)%nx = sizesBlock(4)
    cgnsDoms(nZone)%ny = sizesBlock(5)
    cgnsDoms(nZone)%nz = sizesBlock(6)

    ! Check the size of this zone for the other grids to be read.
    ! They should be equal. Note that familyName is only used as
    ! a dummy in this call.

    do nn=2,nGridsRead
       call cg_zone_read_f(fileIDs(nn), cgnsBase, nZone, &
            familyName, sizesBlock, ierr)
       if(ierr /= CG_OK)               &
            call terminate("readZoneInfo", &
            "Something wrong when calling cg_nZones_f")

       if(cgnsDoms(nZone)%il /= sizesBlock(1) .or. &
            cgnsDoms(nZone)%jl /= sizesBlock(2) .or. &
            cgnsDoms(nZone)%kl /= sizesBlock(3)) then

          write(errorMessage,*) "File ", trim(gridFiles(nn)), &
               ", zone ", trim(familyName), &
               " Zone dimensions are different&
               & than in file ", trim(gridFiles(1))
          if(myID == 0) call terminate("readBlockSizes", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)
       endif
    enddo

    ! Goto this zone.

    call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, "end")
    if(ierr /= CG_OK)               &
         call terminate("readZoneInfo", &
         "Something wrong when calling cg_goto_f")
    !
    !       Try to read the family name.
    !
    call cg_famname_read_f(familyName, ierr)
    if(ierr == error)                &
         call terminate("readZoneInfo", &
         "Something wrong when calling cg_famname_read_f")

    ! Check if a family name was specified. If so, determine the
    ! corresponding id.

    cgnsDoms(nZone)%familyID = 0
    if(ierr == CG_OK) then

       ! Search the family name in the sorted names. For a valid
       ! grid this name must be found.

       nn = cgnsNFamilies
       ii = bsearchStrings(familyName, sortedFamName)
       if(ii == 0) then

          write(errorMessage,100) trim(familyName)
100       format("Family name",1X,A,1X,"not present in the grid")
          if(myID == 0) call terminate("readZoneInfo", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)

       endif

       ! Set the family number.

       ii = famID(ii)
       cgnsDoms(nZone)%familyID = ii

    endif
    !
    !       Try to determine the units of the coordinates.
    !
    ! Determine the number of coordinates in this zone.

    call cg_ncoords_f(cgnsInd, cgnsBase, nZone, nCoords, ierr)
    if(ierr /= CG_OK)               &
         call terminate("readZoneInfo", &
         "Something wrong when calling cg_ncoords_f")

    ! Check that 3 coordinates are present. If not, terminate.

    if(nCoords /= 3) then
       write(errorMessage,102) trim(cgnsDoms(nZone)%zoneName), nCoords
102    format("The number of coordinates of zone ", a, &
            " of base 1 is", i1, ". This should 3.")

       if(myID == 0) call terminate("readZoneInfo", errorMessage)
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Loop over the three coordinates.

    cgnsDoms(nZone)%gridUnitsSpecified = .false.

    do i=1,3

       ! Go to the correct place in the grid file.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
            "GridCoordinates_t", 1, "DataArray_t", i, &
            "end")
       if(ierr /= CG_OK)               &
            call terminate("readZoneInfo", &
            "Something wrong when calling cg_goto_f")

       call cg_units_read_f(mass, len, time, temp, angle, ierr)
       if(ierr == error)                &
            call terminate("readZoneInfo", &
            "Something wrong when calling cg_units_read_f")

       ! Check if units were specified.

       if(ierr == CG_OK .and. len /= Null) then

          ! Copy the units and set gridUnitsSpecified to .true.

          cgnsDoms(nZone)%mass  = mass
          cgnsDoms(nZone)%len   = len
          cgnsDoms(nZone)%time  = time
          cgnsDoms(nZone)%temp  = temp
          cgnsDoms(nZone)%angle = angle

          cgnsDoms(nZone)%gridUnitsSpecified = .true.

          ! Determine the conversion factor to meters.

          call siLen(len, mult, trans)

          cgnsDoms(nZone)%LRef = mult

       endif

    enddo

    ! Check whether units were specified or not.

    if( cgnsDoms(nZone)%gridUnitsSpecified ) then

       ! Units were specified. Check if a global reference length
       ! was specified as well. If so, compare the conversion factor.
       ! If not identical, processor 0 prints a warning.

       if(myID ==0 .and. LRefSpecified .and. &
            cgnsDoms(nZone)%LRef /= LRef)  then

          print "(a)", "#"
          print "(a)", "#                      Warning"
          print 103, trim(cgnsDoms(nZone)%zoneName)
          print 104
          print "(a)", "#"
103       format("# Zone",1X,A,": Conversion factor to meters &
               &not identical to global conversion factor.")
104       format("# Global conversion factor is ignored.")
       endif

       ! In case no global conversion factor was specified,
       ! set LRef to the LRef of this block.

       if(.not. LRefSpecified) LRef = cgnsDoms(nZone)%LRef

    else

       ! No units specified. Set the reference length of the block
       ! to the global conversion factor.

       cgnsDoms(nZone)%LRef = LRef

       ! If the global reference length was not specified, set
       ! noUnits to .true.

       if(.not. LRefSpecified) noUnits = .true.

    endif
    !
    !       Try to determine the rotation rate and center.
    !
    ! Some initializations.

    cgnsDoms(nZone)%rotatingFrameSpecified = .false.
    cgnsDoms(nZone)%rotRate   = zero
    cgnsDoms(nZone)%rotCenter = zero

    mult = one   ! Assuming radians/s for the rotation rate.

    ! Check if a family overwrite is present.

    overwrite = .false.
    nn = cgnsDoms(nZone)%familyID
    if(nn > 0) then
       if( cgnsFamilies(nn)%rotatingFrameSpecified ) &
            overwrite = .true.
    endif

    testOverwrite: if( overwrite ) then

       ! Rotation rate is specified per family.
       ! Set rotatingFrameSpecified to .true. and copy the data
       ! from the corresponding family. Note that the rotation rate
       ! is already in rad/s and therefore mult does not need to
       ! change.

       cgnsDoms(nZone)%rotatingFrameSpecified = .true.
       cgnsDoms(nZone)%rotRate   = cgnsFamilies(nn)%rotRate
       cgnsDoms(nZone)%rotCenter = cgnsFamilies(nn)%rotCenter

    else testOverwrite

       ! Go to the correct location in the cgns file, where
       ! the rotation rate should be specified.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
            "end")
       if(ierr /= CG_OK)               &
            call terminate("readZoneInfo", &
            "Something wrong when calling cg_goto_f")

       ! No family information specified.
       ! Try to read the rotation rate and center.
       call cg_rotating_read_f(real(rotRate,cgnsPerType), real(rotCenter,cgnsPerType), ierr)
       
       if(ierr == error)                &
            call terminate("readZoneInfo", &
            "Something wrong when calling &
            &cg_rotating_read_f")

       ! Check if a rotating frame is specified.

       if(ierr == CG_OK) then

          ! Set changingGrid to .true.

          changing_grid = .true.

          ! Set rotatingFrameSpecified to .true. and copy
          ! to rotation center and rotation rate.

          cgnsDoms(nZone)%rotatingFrameSpecified = .true.
          cgnsDoms(nZone)%rotRate   = rotRate
          cgnsDoms(nZone)%rotCenter = rotCenter

          ! Determine the conversion factor to rad/s if the
          ! dimensional units are specified. If not, it is assumed
          ! that the dimensions are given in rad/s.

          call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
               "RotatingCoordinates_t", 1, "DataArray_t", 2, &
               "end")
          if(ierr /= CG_OK)               &
               call terminate("readZoneInfo", &
               "Something wrong when calling cg_goto_f")

          call cg_units_read_f(mass, len, time, temp, angle, ierr)
          if(ierr == error)                &
               call terminate("readZoneInfo", &
               "Something wrong when calling &
               &cg_units_read_f")

          ! Check if units were specified. If not assume radians.

          if(ierr == CG_OK .and. angle /= Null) then

             ! Determine the conversion factor to radIans.

             call siAngle(angle, mult, trans)

          else

             ! Angle units not specified. Assume rad/s.
             ! Processor 0 writes a warning to stdout.

             if(myID == 0) then

                print "(a)", "#"
                print "(a)", "#                      Warning"
                print 101, trim(cgnsDoms(nZone)%zoneName)
                print "(a)", "#"
101             format("# Zone",1X,A,": No unit specified for &
                     &the rotation rate, assuming rad/s.")
             endif

          endif
       endif

    endif testOverwrite

    ! Multiply the rotation center with LRef to obtain the correct
    ! coordinates in meters.

    cgnsDoms(nZone)%rotCenter(1) = cgnsDoms(nZone)%LRef &
         * cgnsDoms(nZone)%rotCenter(1)
    cgnsDoms(nZone)%rotCenter(2) = cgnsDoms(nZone)%LRef &
         * cgnsDoms(nZone)%rotCenter(2)
    cgnsDoms(nZone)%rotCenter(3) = cgnsDoms(nZone)%LRef &
         * cgnsDoms(nZone)%rotCenter(3)

    ! Multiply the rotation rate by the conversion factor to obtain
    ! the correct angular velocity.

    cgnsDoms(nZone)%rotRate(1) = cgnsDoms(nZone)%rotRate(1)*mult
    cgnsDoms(nZone)%rotRate(2) = cgnsDoms(nZone)%rotRate(2)*mult
    cgnsDoms(nZone)%rotRate(3) = cgnsDoms(nZone)%rotRate(3)*mult


  end subroutine readZoneInfo

  subroutine countConnectivities(cgnsInd, cgnsBase, nZone)
    !
    !       countConnectivities determines the number of connectivities
    !       for each of the supported types stored in 1to1 and general.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsNonMatchAbuttingConnType
    use communication, only : myid, adflow_comm_world
    use partitionMod, only : subfaceNonMatchType, qsortSubfaceNonMatchType
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in)  :: cgnsInd, cgnsBase, nZone
    !
    !      Local variables.
    !
    integer :: i, ngeneral, ierr
    integer :: n1to1, n1to1General, nNonMatch

    integer :: location, connectType, ptsetType
    integer(kind=cgsize_t) :: npnts, ndataDonor
    integer :: donorZoneType, donorPtsetType, donorDatatype

    integer, dimension(:),     allocatable :: connIDNonMatch
    integer(kind=cgsize_t), dimension(:,:),   allocatable :: donorData
    integer(kind=cgsize_t), dimension(:,:,:), allocatable :: myRangeNonMatch

    integer(kind=intType) :: mm, nn

    integer(kind=intType), dimension(:), allocatable :: multSubfaces

    type(subfaceNonMatchType), dimension(:), allocatable :: &
         subfaceNonMatch

    type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
         connNonMatch

    character(len=maxStringLen)   :: errorMessage
    character(len=maxCGNSNameLen) :: connectName, donorName

    ! Determine the number of 1 to 1 connectivities in this zone.
    ! Note that the reading takes place via an integer type.

    call cg_n1to1_f(cgnsInd, cgnsBase, nZone, i, ierr)
    if(ierr /= CG_OK)                      &
         call terminate("countConnectivities", &
         "Something wrong when calling cg_n1to1_f")

    n1to1 = i

    ! Determine the total number of general connectivities in this
    ! zone.

    call cg_nconns_f(cgnsInd, cgnsBase, nZone, ngeneral, ierr)
    if(ierr /= CG_OK)                      &
         call terminate("countConnectivities", &
         "Something wrong when calling cg_nconns_f")

    ! Allocate the memory for connIDNonMatch and myRangeNonMatch. Note
    ! that this number is an upper bound, because other connectivities
    ! may be present in general connectivities.

    allocate(connIDNonMatch(ngeneral), &
         myRangeNonMatch(3,2,ngeneral),stat=ierr)
    if(ierr /= CG_OK)                      &
         call terminate("countConnectivities", &
         "Memory allocation failure for connIDNonMatch &
         &and myRangeNonMatch")

    ! Loop over ngeneral to find out how many of each supported
    ! types of connectivities are stored here.

    n1to1General = 0
    nNonMatch    = 0

    do i=1,ngeneral

       ! Read the information of this connectivity.

       call cg_conn_info_f(cgnsInd, cgnsBase, nZone, i, connectName, &
            location, connectType, ptsetType, npnts,  &
            donorName, donorZoneType, donorPtsetType, &
            donorDatatype, ndataDonor, ierr)
       if(ierr /= CG_OK)                      &
            call terminate("countConnectivities", &
            "Something wrong when calling cg_conn_info_f")

       ! Check if this is a supported structured connectivity.

       select case (connectType)

       case (Abutting1to1)

          if(location       == Vertex       .and. &
               ptsetType      == PointRange   .and. &
               donorZoneType  == Structured   .and. &
               donorPtsetType == PointListDonor) then

             n1to1General = n1to1General + 1

          else

             ! CGNS format not supported.

             write(errorMessage,101) trim(cgnsDoms(nZone)%zoneName), &
                  trim(connectName)
             if(myID == 0) &
                  call terminate("countConnectivities", errorMessage)
             call mpi_barrier(ADflow_comm_world, ierr)

          endif

          !============================================================

       case (Abutting)

          if(location       == Vertex       .and. &
               ptsetType      == PointRange   .and. &
               donorZoneType  == Structured   .and. &
               donorPtsetType == PointListDonor) then

             nNonMatch = nNonMatch + 1

             connIDNonMatch(nNonMatch) = i

             ! Allocate the memory for donorData and read the info.
             ! Release the memory of donorData after the read, because
             ! only the subface range is needed at the moment.

             allocate(donorData(3,ndataDonor), stat=ierr)
             if(ierr /= 0)                           &
                  call terminate("countConnectivities", &
                  "Memory allocation failure for &
                  &donorData")

             call cg_conn_read_f(cgnsInd, cgnsBase, nZone, i,    &
                  myRangeNonMatch(1,1,nNonMatch), &
                  Integer, donorData, ierr)
             if(ierr /= CG_OK)                      &
                  call terminate("countConnectivities", &
                  "Something wrong when calling &
                  &cg_conn_read_f")

             deallocate(donorData, stat=ierr)
             if(ierr /= 0)                           &
                  call terminate("countConnectivities", &
                  "Deallocation failure for donorData")
          else

             ! CGNS format not supported.

             write(errorMessage,102) trim(cgnsDoms(nZone)%zoneName), &
                  trim(connectName)
             if(myID == 0) &
                  call terminate("countConnectivities", errorMessage)
             call mpi_barrier(ADflow_comm_world, ierr)

          endif

          !============================================================

       case default

          call terminate("countConnectivities", &
               "Unsupportted general connectivity found")

       end select

    enddo

    ! Update the value of n1to1 with the number stored in the
    ! general connectivities.

    n1to1 = n1to1 + n1to1General

    ! Memory allocation for the 1 to 1 connectivities.

    allocate(cgnsDoms(nZone)%conn1to1(n1to1), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("countConnectivities", &
         "Memory allocation failure for conn1to1 and &
         &connOver")
    !
    !       For the non-matching abutting subfaces some more information
    !       needs to be extracted. The reason is that a subface abuts
    !       multiple blocks and in CGNS this info is stored in multiple
    !       connectivities. However it is a lot easier to store that info
    !       together. That's why the non-abbuting subfaces must be sorted
    !       in increasing order to extract this information.
    !
    ! Allocate the memory for subfaceNonMatch and copy the data
    ! from myRangeNonMatch and connIDNonMatch. Release the memory
    ! of these two arrays afterwards. Also allocate the memory
    ! for multSubfaces, which is needed later on to determine the
    ! multiplicity of the subfaces.

    allocate(subfaceNonMatch(nNonMatch), &
         multSubfaces(nNonMatch), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("countConnectivities", &
         "Memory allocation failure for subfaceNonMatch")

    do i=1,nNonMatch
       subfaceNonMatch(i)%iBeg = min(myRangeNonMatch(1,1,i), &
            myRangeNonMatch(1,2,i))
       subfaceNonMatch(i)%jBeg = min(myRangeNonMatch(2,1,i), &
            myRangeNonMatch(2,2,i))
       subfaceNonMatch(i)%kBeg = min(myRangeNonMatch(3,1,i), &
            myRangeNonMatch(3,2,i))

       subfaceNonMatch(i)%iEnd = max(myRangeNonMatch(1,1,i), &
            myRangeNonMatch(1,2,i))
       subfaceNonMatch(i)%jEnd = max(myRangeNonMatch(2,1,i), &
            myRangeNonMatch(2,2,i))
       subfaceNonMatch(i)%kEnd = max(myRangeNonMatch(3,1,i), &
            myRangeNonMatch(3,2,i))

       subfaceNonMatch(i)%connID = connIDNonMatch(i)
    enddo

    deallocate(connIDNonMatch, myRangeNonMatch, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("countConnectivities", &
         "Deallocation failure for connIDNonMatch and &
         &myRangeNonMatch")

    ! Sort subfaceNonMatch in increasing order and determine the
    ! number of different subfaces as well as their multiplicity.

    call qsortSubfaceNonMatchType(subfaceNonMatch, nNonMatch)

    nn           = min(nNonMatch, 1_intType)
    multSubfaces = 1

    do i=2,nNonMatch
       mm = i - 1
       if(subfaceNonMatch(i)%iBeg == subfaceNonMatch(mm)%iBeg .and. &
            subfaceNonMatch(i)%jBeg == subfaceNonMatch(mm)%jBeg .and. &
            subfaceNonMatch(i)%kBeg == subfaceNonMatch(mm)%kBeg .and. &
            subfaceNonMatch(i)%iEnd == subfaceNonMatch(mm)%iEnd .and. &
            subfaceNonMatch(i)%jEnd == subfaceNonMatch(mm)%jEnd .and. &
            subfaceNonMatch(i)%kEnd == subfaceNonMatch(mm)%kEnd) then
          multSubfaces(nn) = multSubfaces(nn) + 1
       else
          nn = nn + 1
       endif
    enddo

    ! Store the number of non-matching connectivities in nNonMatch and
    ! allocate the memory for the non-matching abutting connectivities

    nNonMatch = nn
    allocate(cgnsDoms(nZone)%connNonMatchAbutting(nn), stat=ierr)
    if(ierr /= 0)                           &
         call terminate("countConnectivities", &
         "Memory allocation failure for &
         &connNonMatchAbutting")

    connNonMatch => cgnsDoms(nZone)%connNonMatchAbutting

    ! Loop over the number of non-matching connectivities to copy
    ! the info from subfaceNonMatch.

    nn = 0
    do i=1,nNonMatch

       ! Set the value for the number of donor blocks and allocate
       ! the memory for connectNames, donorNames, donorBlocks and
       ! donorFaceIDs.

       mm = multSubfaces(i)
       connNonMatch(i)%nDonorBlocks = mm

       allocate(connNonMatch(i)%connectNames(mm), &
            connNonMatch(i)%donorNames(mm),   &
            connNonMatch(i)%donorBlocks(mm),  &
            connNonMatch(i)%donorFaceIDs(mm), stat=ierr)
       if(ierr /= 0)                           &
            call terminate("countConnectivities", &
            "Memory allocation failure for connectNames, &
            &donorNames, donorBlocks and donorFaceIDs.")

       ! Loop over the number of donor blocks and copy for the moment
       ! the connectivity ID in donorBlocks.

       do mm=1,connNonMatch(i)%nDonorBlocks
          nn = nn + 1
          connNonMatch(i)%donorBlocks(mm) = subfaceNonMatch(nn)%connID
       enddo

       ! Copy the subface range.

       connNonMatch(i)%iBeg = subfaceNonMatch(nn)%iBeg
       connNonMatch(i)%jBeg = subfaceNonMatch(nn)%jBeg
       connNonMatch(i)%kBeg = subfaceNonMatch(nn)%kBeg

       connNonMatch(i)%iEnd = subfaceNonMatch(nn)%iEnd
       connNonMatch(i)%jEnd = subfaceNonMatch(nn)%jEnd
       connNonMatch(i)%kEnd = subfaceNonMatch(nn)%kEnd

    enddo

    ! Release the memory of subfaceNonMatch again.

    deallocate(subfaceNonMatch, stat=ierr)
    if(ierr /= 0)                           &
         call terminate("countConnectivities", &
         "Deallocation failure for subfaceNonMatch")
    !
    !       Store the number of connectivities in cgnsDoms(nZone).
    !
    cgnsDoms(nZone)%n1to1             = n1to1
    cgnsDoms(nZone)%n1to1General      = n1to1General
    cgnsDoms(nZone)%nNonMatchAbutting = nNonMatch
    !
    !       Format statements.
    !
101 format("Zone",1x,a,", connectivity", 1x,a, ": No support for &
         &this format of an abutting 1 to 1 connectivity")
102 format("Zone",1x,a,", connectivity", 1x,a, ": No support for &
         &this format of a non-matching abutting connectivity")


  end subroutine countConnectivities

  subroutine read1to1Conn(cgnsInd, cgnsBase, nZone)
    !
    !       read1to1Conn reads the 1 to 1 block to block, i.e.
    !       continuous grid lines across block boundaries, connectivities
    !       for the given zone/block.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgns1to1ConnType
    use communication, only : adflow_comm_world, myID
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase, nZone
    !
    !      Local variables
    !
    integer :: cgnsN1to1
    integer :: i, ierr

    integer(kind=cgsize_t), dimension(3,2) :: zoneRange, donorRange
    integer, dimension(3)   :: transform

    character(len=maxCGNSNameLen) :: connectName

    type(cgns1to1ConnType),    pointer, dimension(:) :: conn1to1
    real(kind=realType), dimension(3) :: rotCenter, rotAngles
    real(kind=realType), dimension(3) :: tlation

    ! Determine the number of 1 to 1 connectivities stored in the
    ! CGNS file for this zone.

    cgnsN1to1 = cgnsDoms(nZone)%n1to1 - cgnsDoms(nZone)%n1to1General

    ! Loop over the number of 1 to 1 connectivities stored in the
    ! grid file.

    do i=1,cgnsN1to1

       ! Read the 1 to 1 connectivity info from the cgns file.

       call cg_1to1_read_f(cgnsInd, cgnsBase, nZone, i,             &
            cgnsDoms(nZone)%conn1to1(i)%connectName, &
            cgnsDoms(nZone)%conn1to1(i)%donorName,   &
            zoneRange, donorRange, transform, ierr)
       if(ierr /= CG_OK)                 &
            call terminate("read1to1Conn", &
            "Something wrong when calling cg_1to1_read_f")

       ! Store the zone range and donor range in cgnsDoms.

       cgnsDoms(nZone)%conn1to1(i)%iBeg = zoneRange(1,1)
       cgnsDoms(nZone)%conn1to1(i)%jBeg = zoneRange(2,1)
       cgnsDoms(nZone)%conn1to1(i)%kBeg = zoneRange(3,1)

       cgnsDoms(nZone)%conn1to1(i)%iEnd = zoneRange(1,2)
       cgnsDoms(nZone)%conn1to1(i)%jEnd = zoneRange(2,2)
       cgnsDoms(nZone)%conn1to1(i)%kEnd = zoneRange(3,2)

       cgnsDoms(nZone)%conn1to1(i)%diBeg = donorRange(1,1)
       cgnsDoms(nZone)%conn1to1(i)%djBeg = donorRange(2,1)
       cgnsDoms(nZone)%conn1to1(i)%dkBeg = donorRange(3,1)

       cgnsDoms(nZone)%conn1to1(i)%diEnd = donorRange(1,2)
       cgnsDoms(nZone)%conn1to1(i)%djEnd = donorRange(2,2)
       cgnsDoms(nZone)%conn1to1(i)%dkEnd = donorRange(3,2)

       ! Check the transformation matrices between this zone and the
       ! donor zone and store it in l1, L2 and l3.

       call checkTransform(transform, nZone, i, .true.)

       cgnsDoms(nZone)%conn1to1(i)%l1 = transform(1)
       cgnsDoms(nZone)%conn1to1(i)%l2 = transform(2)
       cgnsDoms(nZone)%conn1to1(i)%l3 = transform(3)

       ! Subface is a normal boundary. Set periodic to .false. and
       ! initialize the periodic data to zero to avoid possible
       ! problems due to uninitialized data.

       cgnsDoms(nZone)%conn1to1(i)%periodic = .false.

       cgnsDoms(nZone)%conn1to1(i)%rotationCenter = zero
       cgnsDoms(nZone)%conn1to1(i)%rotationAngles = zero
       cgnsDoms(nZone)%conn1to1(i)%translation    = zero


       call cg_1to1_periodic_read_f(cgnsInd, cgnsBase, nZone, i, &
            real(rotCenter,cgnsPerType), real(rotAngles,cgnsPerType), real(tlation,cgnsPerType), ierr)
       if(ierr == CG_OK)then
          call readPeriodicSubface1to1(cgnsInd, cgnsBase, nZone, i,    &
               cgnsDoms(nZone)%conn1to1(i)%connectName,                       &
               cgnsDoms(nZone)%conn1to1(i)%periodic,          &
               cgnsDoms(nZone)%conn1to1(i)%rotationCenter,    &
               cgnsDoms(nZone)%conn1to1(i)%rotationAngles,    &
               cgnsDoms(nZone)%conn1to1(i)%translation)
       endif

    enddo



  end subroutine read1to1Conn
  subroutine checkTransform(transform, nZone, n1to1, printWarning)
    !
    !       checkTransform checks the transformation matrix between this
    !       zone and the donor for the given subrange. In case an error is
    !       found it is tried to correct this.
    !
    use constants
    use cgnsGrid, only : cgnsDoms
    use communication, only : myID, adflow_comm_world
    use utils, only : delta, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer, intent(in)                  :: nZone, n1to1
    integer, dimension(3), intent(inout) :: transform
    logical, intent(in)                  :: printWarning
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nDirFace, nDirDonor
    integer(kind=intType) :: sumTransform
    integer(kind=intType) :: l1, L2, l3

    integer(kind=intType), dimension(3) :: haloDir, donorDir
    integer(kind=intType), dimension(3,2) :: zoneRange, donorRange
    integer(kind=intType), dimension(3,3) :: trMat

    character(len=maxCGNSNameLen) :: zoneName, connectName
    character(len=2*maxStringLen)  :: errorMessage
    !
    ! Copy the zoneName and connectName, just for readability later.

    zoneName    = cgnsDoms(nZone)%zoneName
    connectName = cgnsDoms(nZone)%conn1to1(n1to1)%connectName

    ! Copy the zone and donor range into zoneRange and donorRange.

    zoneRange(1,1) = cgnsDoms(nZone)%conn1to1(n1to1)%iBeg
    zoneRange(2,1) = cgnsDoms(nZone)%conn1to1(n1to1)%jBeg
    zoneRange(3,1) = cgnsDoms(nZone)%conn1to1(n1to1)%kBeg

    zoneRange(1,2) = cgnsDoms(nZone)%conn1to1(n1to1)%iEnd
    zoneRange(2,2) = cgnsDoms(nZone)%conn1to1(n1to1)%jEnd
    zoneRange(3,2) = cgnsDoms(nZone)%conn1to1(n1to1)%kEnd

    donorRange(1,1) = cgnsDoms(nZone)%conn1to1(n1to1)%diBeg
    donorRange(2,1) = cgnsDoms(nZone)%conn1to1(n1to1)%djBeg
    donorRange(3,1) = cgnsDoms(nZone)%conn1to1(n1to1)%dkEnd

    donorRange(1,2) = cgnsDoms(nZone)%conn1to1(n1to1)%diEnd
    donorRange(2,2) = cgnsDoms(nZone)%conn1to1(n1to1)%djEnd
    donorRange(3,2) = cgnsDoms(nZone)%conn1to1(n1to1)%dkEnd

    ! Determine the normal direction for the subface and do a trivial
    ! check to see if there is one.

    do nDirFace=1,3
       if(zoneRange(nDirFace,1) == zoneRange(nDirFace,2)) exit
    enddo

    if(nDirFace > 3) then
       if(myID == 0) then
          write(errorMessage,100) trim(connectName), trim(zoneName)
100       format("1 to 1 subface",1X,A,1X,"of zone",1X,A, &
               ": No constant index found")
          call terminate("checkTransform", errorMessage)
       endif

       ! Make sure that other processors wait until they are killed.

       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Idem for the donor subface.

    do nDirDonor=1,3
       if(donorRange(nDirDonor,1) == donorRange(nDirDonor,2)) exit
    enddo

    if(nDirDonor > 3) then
       if(myID == 0) then
          write(errorMessage,110) trim(connectName), trim(zoneName)
110       format("1 to 1 subface",1X,A,1X,"of zone",1X,A, &
               ": No constant index found for donor")
          call terminate("checkTransform", errorMessage)
       endif

       ! Make sure that other processors wait until they are killed.

       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Check if the sum of the absolute values of transform equals 6.
    ! If not, assume that the normal direction is not set correctly.

    sumTransform = abs(transform(1)) + abs(transform(2)) &
         + abs(transform(3))
    if(sumTransform /= 6) then

       ! Change the normal direction of transform and check the
       ! sum again.

       transform(nDirFace) = nDirDonor
       sumTransform = abs(transform(1)) + abs(transform(2)) &
            + abs(transform(3))

       if(sumTransform /= 6) then

          ! Something seriously wrong. I cannot repair this.

          if(myID == 0) then
             write(errorMessage,120) trim(connectName), trim(zoneName)
120          format("1 to 1 subface",1X,A,1X,"of zone",1X,A,  &
                  ": Something seriously wrong with the &
                  &transformation matrix")
             call terminate("checkTransform", errorMessage)
          endif

          ! Make sure that other processors wait until they are killed.

          call mpi_barrier(ADflow_comm_world, ierr)
       else
          ! Repair successful, although the orientation might be wrong.
          ! This will be checked later. Anyway print a warning message
          ! if desired.

          if(myID == 0 .and. printWarning) then
             print "(a)", "#"
             print "(a)", "#                          Warning"
             print 130, trim(connectName), trim(zoneName)
130          format("# 1 to 1 subface",1X,A,1X,"of zone",1X,A, &
                  ": Normal component of the transformation&
                  & matrix successfully corrected.")
             print "(a)", "#"
          endif
       endif
    endif

    ! Create the halo vector for the current zone. This vector is
    ! pointing outwards, i.e. in direction of the donor block.

    haloDir = 0
    haloDir(nDirFace) = 1
    if(zoneRange(nDirFace,1) == 1) haloDir(nDirFace) = -1

    ! Idem for the donor. Also this vector points from the current
    ! block to the donor block, albeit in donor block coordinates.

    donorDir = 0
    donorDir(nDirDonor) = -1
    if(donorRange(nDirDonor,1) == 1) donorDir(nDirDonor) = 1

    ! Determine the full transformation matrix.

    l1 = transform(1)
    L2 = transform(2)
    l3 = transform(3)

    trMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
    trMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
    trMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

    trMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
    trMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
    trMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

    trMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
    trMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
    trMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

    ! Apply the transformation matrix to haloDir.

    l1 = haloDir(1)
    L2 = haloDir(2)
    l3 = haloDir(3)

    haloDir(1) = trMat(1,1)*l1 + trMat(1,2)*l2 + trMat(1,3)*l3
    haloDir(2) = trMat(2,1)*l1 + trMat(2,2)*l2 + trMat(2,3)*l3
    haloDir(3) = trMat(3,1)*l1 + trMat(3,2)*l2 + trMat(3,3)*l3

    ! If the transformation matrix is correct haloDir == donorDir.
    ! If this is not the case, there are two possibilities. Either
    ! the directions are just reversed, which means that the
    ! corresponding element of transform must be reversed, or they
    ! are really different. In the latter case it cannot be
    ! corrected here and the grid file must be adapted.

    if(haloDir(nDirDonor) == 0) then

       ! Something seriously wrong. Exit the program.

       if(myID == 0) then
          write(errorMessage,140) trim(connectName), trim(zoneName)
140       format("1 to 1 subface",1X,A,1X,"of zone",1X,A,  &
               ": Something seriously wrong with the &
               &transformation matrix")
          call terminate("checkTransform", errorMessage)
       endif

       ! Make sure that other processors wait until they are killed.

       call mpi_barrier(ADflow_comm_world, ierr)

    else if(haloDir(nDirDonor)*donorDir(nDirDonor) < 0) then

       ! Simply reverse the sign of the corresponding entry in
       ! transform. Processor 0 prints a warning message if desired.

       transform(nDirFace) = -transform(nDirFace)

       if(myID == 0 .and. printWarning) then
          print "(a)", "#"
          print "(a)", "#                            Warning"
          print 150, trim(connectName), trim(zoneName)
150       format("# 1 to 1 subface",1X,A,1X,"of zone",1X,A, &
               ": Normal component of the transformation&
               & matrix reversed")
          print "(a)", "#"
       endif
   endif

  end subroutine checkTransform
  subroutine readGeneralConn(cgnsInd, cgnsBase, nZone)
    !
    !       readGeneralConn reads and converts the cgns general
    !       connectivities.  Supported connectivites are 1-to-1 and
    !       non-matching abutting.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgns1to1ConnType, &
         cgnsNOnMatchAbuttingConnType
    use communication, only : myid, adflow_comm_world
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in) :: cgnsInd, cgnsBase, nZone

    !
    !      Local variables.
    !
    character(len=maxStringLen) :: errorMessage

    integer :: i, j, nn, nGeneral, n1to1, ierr
    integer :: location, connectType, ptsetType
    integer(kind=cgsize_t) :: npnts
    integer :: donorZoneType, donorPtsetType, donorDatatype
    integer :: id, jj
    integer :: nArrays, dataType, dataDim
    integer(kind=cgsize_t) :: nDataDonor
    integer, dimension(2)   :: dimVector
    integer, dimension(3)   :: ii, transform
    integer(kind=cgsize_t), dimension(3,2) :: myRange

    integer(kind=cgsize_t), dimension(:,:), allocatable :: myData, donorData
    integer, dimension(:,:), allocatable :: map2NonMatch

    real(kind=realType), dimension(3) :: rotationCenter
    real(kind=realType), dimension(3) :: rotationAngles
    real(kind=realType), dimension(3) :: translation

    logical :: periodic, wrongData

    character(len=maxCGNSNameLen) :: connectName, donorName, arrayName

    type(cgns1to1ConnType),    pointer, dimension(:) :: conn1to1

    type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
         connNonMatch
    !
    !      Function definition.
    !
    integer :: setCGNSRealType

    ! Set some pointers for the connectivities to make the code
    ! more readable.

    conn1to1     => cgnsDoms(nZone)%conn1to1
    connNonMatch => cgnsDoms(nZone)%connNonMatchAbutting

    ! Set the counter n1to1 to the currently stored number of 1 to 1
    ! block connectivities, and initialize other counters.

    n1to1 = cgnsDoms(nZone)%n1to1 - cgnsDoms(nZone)%n1to1General

    ! Determine the number of general connectivities.

    call cg_nconns_f(cgnsInd, cgnsBase, nZone, ngeneral, ierr)
    if(ierr /= CG_OK)                  &
         call terminate("readGeneralConn", &
         "Something wrong when calling cg_nconns_f")

    ! Determine the mapping from the general connectivities to the
    ! data structure for the non-matching abutting connectivities.
    ! There are two indices for this mapping, the first is the
    ! index in connNonMatch and the second the neighbor index.
    ! Note that in the data structure a non-matching abutting
    ! subface can have multiple neighbors.
    ! Note that the general connectivity ID is temporarily stored
    ! in connNonMatch(i)%donorBlocks.

    allocate(map2NonMatch(ngeneral,2), stat=ierr)
    if(ierr /= CG_OK)                  &
         call terminate("readGeneralConn", &
         "Memory allocation failure for map2NonMatch")

    do i=1,cgnsDoms(nZone)%nNonMatchAbutting
       do j=1,connNonMatch(i)%nDonorBlocks
          nn = connNonMatch(i)%donorBlocks(j)
          map2NonMatch(nn,1) = i
          map2NonMatch(nn,2) = j
       enddo
    enddo

    ! Loop over the general connectivities.

    nConnLoop: do nn=1,ngeneral

       ! Read the information of this connectivity.

       call cg_conn_info_f(cgnsInd, cgnsBase, nZone, nn, connectName, &
            location, connectType, ptsetType, npnts,   &
            donorName, donorZoneType, donorPtsetType,  &
            donorDatatype, ndataDonor, ierr)
       if(ierr /= CG_OK)                  &
            call terminate("readGeneralConn", &
            "Something wrong when calling cg_conn_info_f")

       ! Read the data based on the type of connectivity.

       connectivityType: select case(connectType)

       case (Abutting1to1)
          !
          !             1-to-1 connectivity stored as a general one. Note that
          !             the check for a valid one has already been done in
          !             countConnectivities.
          !
          ! Update the counter n1to1 and store some info in conn1to1.

          n1to1 = n1to1 + 1

          conn1to1(n1to1)%connectName = connectName
          conn1to1(n1to1)%donorName   = donorName

          ! Allocate the memory for donorData.

          allocate(donorData(3,ndataDonor), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
               "Memory allocation failure for donorData")

          ! Read the ranges of the connectivities.

          call cg_conn_read_f(cgnsInd, cgnsBase, nZone, nn, &
               myRange, Integer, donorData, ierr)
          if(ierr /= CG_OK)                  &
               call terminate("readGeneralConn", &
               "Something wrong when calling &
               &cg_conn_read_f")

          ! Store the range of the current subface and its donor.

          conn1to1(n1to1)%iBeg = myRange(1,1)
          conn1to1(n1to1)%jBeg = myRange(2,1)
          conn1to1(n1to1)%kBeg = myRange(3,1)

          conn1to1(n1to1)%iEnd = myRange(1,2)
          conn1to1(n1to1)%jEnd = myRange(2,2)
          conn1to1(n1to1)%kEnd = myRange(3,2)

          conn1to1(n1to1)%diBeg = donorData(1,1)
          conn1to1(n1to1)%djBeg = donorData(2,1)
          conn1to1(n1to1)%dkBeg = donorData(3,1)

          conn1to1(n1to1)%diEnd = donorData(1,ndataDonor)
          conn1to1(n1to1)%djEnd = donorData(2,ndataDonor)
          conn1to1(n1to1)%dkEnd = donorData(3,ndataDonor)

          ! Determine the transformation matrix between the subface
          ! and the donor subface. Initialize it to 0.

          transform = 0

          ! Determine the fastest changing index in the donor.
          ! Take negative running indices into account.

          ii = donorData(:,2) - donorData(:,1)
          do id=1,3
             if(ii(id) /= 0) exit
          enddo
          if(ii(id) < 0) id = -id

          ! Determine the corresponding index in myRange.

          do jj=1,3
             if(myRange(jj,1) /= myRange(jj,2)) exit
          enddo

          ! Set the corresponding entry of transform; take negative
          ! running indices of myRangle into account.

          if(myRange(jj,1) > myRange(jj,2)) id = -id
          transform(jj) = id

          ! Determine the index in donorData where the second index
          ! of the subface changes for the first time and determine
          ! this second index. Take negative running indices into
          ! account.

          j = abs(myRange(jj,2)-myRange(jj,1)) + 2
          ii = donorData(:,j) - donorData(:,1)
          do id=1,3
             if(ii(id) /= 0) exit
          enddo
          if(ii(id) < 0) id = -id

          ! Determine the corresponding index in myRange and set the
          ! corresponding entry in transform.

          do jj=jj+1,3
             if(myRange(jj,1) /= myRange(jj,2)) exit
          enddo
          if(myRange(jj,1) > myRange(jj,2)) id = -id
          transform(jj) = id

          ! Release the memory of donorData.

          deallocate(donorData, stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
               "Deallocation error for donorData")

          ! Determine the correct third direction of the
          ! transformation matrix. Although not intented,
          ! it also serves as a check.

          call checkTransform(transform, nZone, n1to1, .false.)

          conn1to1(n1to1)%l1 = transform(1)
          conn1to1(n1to1)%l2 = transform(2)
          conn1to1(n1to1)%l3 = transform(3)

          ! Read the periodic info if this is a periodic boundary.

          call readPeriodicSubface(cgnsInd, cgnsBase, nZone, nn,   &
               connectName,                    &
               conn1to1(n1to1)%periodic,       &
               conn1to1(n1to1)%rotationCenter, &
               conn1to1(n1to1)%rotationAngles, &
               conn1to1(n1to1)%translation)

          !=============================================================

       case (Abutting)
          !
          !             Non-matching abutting connectivity. Note that the
          !             check for a valid one has already been done in
          !             countConnectivities.
          !
          ! Determine the indices in connNonMatch where the data of
          ! the current connectivity must be stored.

          i = map2NonMatch(nn,1)
          j = map2NonMatch(nn,2)

          ! Store the names of the connectivity and the donor block.

          connNonMatch(i)%connectNames(j) = connectName
          connNonMatch(i)%donorNames(j)   = donorName

          ! Allocate the memory for donorData.

          allocate(donorData(3,ndataDonor), stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
               "Memory allocation failure for donorData")

          ! Read the ranges of the connectivities.

          call cg_conn_read_f(cgnsInd, cgnsBase, nZone, nn, &
               myRange, Integer, donorData, ierr)
          if(ierr /= CG_OK)                  &
               call terminate("readGeneralConn", &
               "Something wrong when calling &
               &cg_conn_read_f")

          ! Determine the face ID on the abutting donor block.
          ! Note that mu subface range has already been stored in
          ! countConnectivities.

          if(donorData(1,1) == donorData(1,ndataDonor)) then

             if(donorData(1,1) == 1) then
                connNonMatch(i)%donorFaceIDs(j) = iMin
             else
                connNonMatch(i)%donorFaceIDs(j) = iMax
             endif

          else if(donorData(2,1) == donorData(2,ndataDonor)) then

             if(donorData(2,1) == 1) then
                connNonMatch(i)%donorFaceIDs(j) = jMin
             else
                connNonMatch(i)%donorFaceIDs(j) = jMax
             endif

          else if(donorData(3,1) == donorData(3,ndataDonor)) then

             if(donorData(3,1) == 1) then
                connNonMatch(i)%donorFaceIDs(j) = kMin
             else
                connNonMatch(i)%donorFaceIDs(j) = kMax
             endif

          else

             write(errorMessage,100) trim(cgnsDoms(nZone)%zoneName), &
                  trim(connectName)
             if(myID == 0) &
                  call terminate("readGeneralConn", errorMessage)
             call mpi_barrier(ADflow_comm_world, ierr)

          endif

          ! Release the memory of donorData.

          deallocate(donorData, stat=ierr)
          if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
               "Deallocation error for donorData")

          ! Read the periodic info if this is a periodic boundary.

          call readPeriodicSubface(cgnsInd, cgnsBase, nZone, nn,   &
               connectName, periodic,          &
               rotationCenter, rotationAngles, &
               translation)

          ! If j == 1 then the info is simply copied. Otherwise it
          ! is checked if it is okay.

          checkConsistency: if(j == 1) then

             connNonMatch(i)%periodic       = periodic
             connNonMatch(i)%rotationCenter = rotationCenter
             connNonMatch(i)%rotationAngles = rotationAngles
             connNonMatch(i)%translation    = translation

          else checkConsistency

             ! Check for consistency.

             wrongData = .false.

             ! Some dirty stuff to compare the logicals. The FORTRAN
             ! standard does not allow comparisons of logicals.

             ii(1) = 0; if( connNonMatch(i)%periodic ) ii(1) = 1
             ii(2) = 0; if(                 periodic ) ii(2) = 1

             if(ii(1) /= ii(2)) wrongData = .true.

             ! Compare the rotation center.

             if(abs(connNonMatch(i)%rotationCenter(1)              &
                  - rotationCenter(1)) > eps .and. &
                  abs(connNonMatch(i)%rotationCenter(2)              &
                  - rotationCenter(2)) > eps .and. &
                  abs(connNonMatch(i)%rotationCenter(3)              &
                  - rotationCenter(3)) > eps)      &
                  wrongData = .true.

             ! Check the rotation angles.

             if(abs(connNonMatch(i)%rotationAngles(1)              &
                  - rotationAngles(1)) > eps .and. &
                  abs(connNonMatch(i)%rotationAngles(2)              &
                  - rotationAngles(2)) > eps .and. &
                  abs(connNonMatch(i)%rotationAngles(3)              &
                  - rotationAngles(3)) > eps)      &
                  wrongData = .true.

             ! Check the translation vector.

             if(abs(connNonMatch(i)%translation(1)              &
                  - translation(1)) > eps .and. &
                  abs(connNonMatch(i)%translation(2)              &
                  - translation(2)) > eps .and. &
                  abs(connNonMatch(i)%translation(3)              &
                  - translation(3)) > eps)      &
                  wrongData = .true.

             ! Print an error message and exit if inconsistent
             ! data was found.

             write(errorMessage,101) trim(cgnsDoms(nZone)%zoneName), &
                  trim(connectName),              &
                  trim(connNonMatch(i)%connectNames(1))
             if(myID == 0) &
                  call terminate("readGeneralConn", errorMessage)
             call mpi_barrier(ADflow_comm_world, ierr)

          endif checkConsistency

          !=============================================================

       end select connectivityType

    enddo nConnLoop

    ! Release the memory of map2NonMatch again.

    deallocate(map2NonMatch, stat=ierr)
    if(ierr /= CG_OK)                  &
         call terminate("readGeneralConn", &
         "Deallocation failure for map2NonMatch")

    ! Format statements.

100 format("Zone",1x,a,", connectivity", 1x,a, &
         ": Invalid donor subface.")
101 format("Zone",1x,a,", connectivity", 1x,a, &
         ": Inconsistent periodic info compared to connectivity", &
         1x,a,".")
102 format("Zone",1x,a,", connectivity", 1x,a, &
         ": InterpolantsDonor node is missing.")
103 format("Zone",1x,a,", connectivity", 1x,a, &
         ": Wrong number or size of interpolants array.")


  end subroutine readGeneralConn

  subroutine readBocos(cgnsInd, cgnsBase, nZone,                &
       nDoubleBoundFaces, sortedFamName, famID)
    !
    !       ReadBocos reads the boundary condition info for the given
    !       zone/block.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsBcDatasetType, &
         cgnsFamilies, cgnsBcdataArray, cgnsNFamilies
    use communication, only : myID, adflow_comm_world
    use utils, only: terminate, setcgnsRealType
    use sorting, only: bsearchStrings
    implicit none
    !
    !      Subroutine arguments
    !
    integer, intent(in)    :: cgnsInd, cgnsBase, nZone
    integer, intent(inout) :: nDoubleBoundFaces

    character(len=*), dimension(:), intent(in) :: sortedFamName
    integer(kind=intType), dimension(:), intent(in) :: famID
    !
    !      Local variables
    !
    integer :: cgnsNBocos, cgnsNDataSet, nUserData
    integer(kind=cgsize_t) :: cgnsNpnts
    integer :: i, j, match
    integer :: ierr, dummy
    integer :: dirichletFlag, neumannFlag

    integer(kind=cgsize_t), dimension(3,2) :: bcRange

    integer(kind=intType) :: ii, nn

    character(len=maxCGNSNameLen) :: familyName
    character(len=maxStringLen)    :: errorMessage

    logical :: familySpecifiedData

    type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet


    ! Read the number of boundary conditions in this zone/block.
    ! Again the reading takes place via an integer.

    call cg_nbocos_f(cgnsInd, cgnsBase, nZone, cgnsNBocos, ierr)
    if(ierr /= CG_OK)            &
         call terminate("readBocos", &
         "Something wrong when calling cg_nbocos_f")
    cgnsDoms(nZone)%nBocos = cgnsNBocos

    ! Allocate the memory for the boundary condition info for
    ! this zone/block.

    allocate(cgnsDoms(nZone)%bocoInfo(cgnsNBocos), stat=ierr)
    if(ierr /= 0) &
         call terminate("readBocos", &
         "Memory allocation failure for bocoInfo")

    ! Loop over the boundary conditions.

    bocoLoop: do i=1,cgnsNBocos
       !
       !         Read the general info for this boundary condition and set
       !         the dimensions of the subface.
       !
       call cg_boco_info_f(cgnsInd, cgnsBase, nZone, i,                &
            cgnsDoms(nZone)%bocoInfo(i)%bocoName,       &
            cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS,     &
            cgnsDoms(nZone)%bocoInfo(i)%ptsetType,      &
            cgnsNpnts,                                  &
            cgnsDoms(nZone)%bocoInfo(i)%normalIndex,    &
            cgnsDoms(nZone)%bocoInfo(i)%normalListFlag, &
            cgnsDoms(nZone)%bocoInfo(i)%normalDataType, &
            cgnsNDataSet, ierr)

       if(ierr /= CG_OK)           &
            call terminate("readBocos", &
            "Something wrong when calling cg_boco_info_f")

       cgnsDoms(nZone)%bocoInfo(i)%npnts = cgnsNpnts

       ! Nullify the pointer for dataSet.

       nullify(cgnsDoms(nZone)%bocoInfo(i)%dataSet)

       ! Perform some checks.
       if(cgnsDoms(nZone)%bocoInfo(i)%normalListFlag > 0) &
            call terminate("readBocos", &
            "Currently not possible to read &
            &boundary normals")

       ! Check how the boundary conditions are specified. Normally this
       ! is a given point range, but there is limited support for
       ! specifying a point list. The latter is needed for icem grids
       ! when boundary conditions are specified on a point.

       if(cgnsDoms(nZone)%bocoInfo(i)%ptsetType == pointrange) then

          ! Point range specified. The usual situation.
          ! Read the point range for this boundary condition.

          call cg_boco_read_f(cgnsInd, cgnsBase, nZone, i, bcRange, &
               dummy, ierr)
          if(ierr /= CG_OK)             &
               call terminate("readBocos", &
               "Something wrong when calling &
               &cg_boco_read_f")

       else if(cgnsDoms(nZone)%bocoInfo(i)%ptsetType == pointlist) then

          ! List of points specified. This is normally undesired, because
          ! boundary conditions are applied per face. However icem grids
          ! tend to give boundary conditions for corner points and this
          ! must be handled. In this case the number of points specified,
          ! cgnsNpnts, equals 1. In all other situations print an error
          ! message.

          if(cgnsNpnts > 1)             &
               call terminate("readBocos", &
               "Point list with more than 1 point specified")

          ! Read the point index.

          call cg_boco_read_f(cgnsInd, cgnsBase, nZone, i, bcRange, &
               dummy, ierr)
          if(ierr /= CG_OK)             &
               call terminate("readBocos", &
               "Something wrong when calling cg_boco_read_f")

          ! Make sure that bcRange contains a range, i.e. the point
          ! indices are copied to the second column of bcRange.

          bcRange(1,2) = bcRange(1,1)
          bcRange(2,2) = bcRange(2,1)
          bcRange(3,2) = bcRange(3,1)

       else

          ! Unknown ptsetType.

          call terminate("readBocos", "Unknown ptsetType encountered")

       endif

       ! Store the range in the cgns grid type.

       cgnsDoms(nZone)%bocoInfo(i)%iBeg = bcRange(1,1)
       cgnsDoms(nZone)%bocoInfo(i)%jBeg = bcRange(2,1)
       cgnsDoms(nZone)%bocoInfo(i)%kBeg = bcRange(3,1)

       cgnsDoms(nZone)%bocoInfo(i)%iEnd = bcRange(1,2)
       cgnsDoms(nZone)%bocoInfo(i)%jEnd = bcRange(2,2)
       cgnsDoms(nZone)%bocoInfo(i)%kEnd = bcRange(3,2)

       ! Check and see if this is a valid boundary condition or if it
       ! corresponds to either an edge or a point, which must be
       ! ignored by the flow solver.

       match = 0
       if(bcRange(1,1) == bcRange(1,2)) match = match +1
       if(bcRange(2,1) == bcRange(2,2)) match = match +1
       if(bcRange(3,1) == bcRange(3,2)) match = match +1

       if(match == 1) then
          cgnsDoms(nZone)%bocoInfo(i)%actualFace = .true.
       else
          cgnsDoms(nZone)%bocoInfo(i)%actualFace = .false.
       endif

       ! It is possible that the 1 to 1 block connectivities are
       ! repeated in the boundary section. Check for this, but only
       ! if this an actual face.

       if( cgnsDoms(nZone)%bocoInfo(i)%actualFace ) then
          if( checkForDoubleBoundFace(nZone, i) ) then

             ! Face is repeated. Increment nDoubleBoundFaces and set
             ! actualFace to .false.

             nDoubleBoundFaces = nDoubleBoundFaces +1
             cgnsDoms(nZone)%bocoInfo(i)%actualFace = .false.

          endif
       endif
       !
       !         Determine the internally used boundary condition and whether
       !         or not the boundary condition is given on a per family basis.
       !
       cgnsDoms(nZone)%bocoInfo(i)%familyID = 0

       checkActualFace: if( cgnsDoms(nZone)%bocoInfo(i)%actualFace ) then

          ! Determine the type of CGNS boundary condition and act
          ! accordingly.

          select case (cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS)

          case (FamilySpecified)

             ! Boundary condition is specified per family.
             !added to accomodate case where a grid family is specified
             ! but BC's are specified in standard CGNS Format
             cgnsDoms(nZone)%BCFamilies = .True.
             ! Find out the family name to which this boundary
             ! face belongs.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                  "ZoneBC_t", 1, "BC_t", i, "end")
             if(ierr /= CG_OK)             &
                  call terminate("readBocos", &
                  "Something wrong when calling cg_goto_f")

             call cg_famname_read_f(familyName, ierr)

             if(ierr /= CG_OK) then

                write(errorMessage,101) trim(cgnsDoms(nZone)%zoneName), &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
                if(myID == 0) call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)

             endif

             ! Search the family name in the sorted names. For a valid
             ! grid this name must be found.

             nn = cgnsNFamilies
             ii = bsearchStrings(familyName, sortedFamName)
             if(ii == 0) then

                write(errorMessage,102) trim(familyName)
                if(myID == 0) call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)

             endif

             ! Set the family number and the boundary condition types.

             ii = famID(ii)
             cgnsDoms(nZone)%bocoInfo(i)%familyID = ii
             cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS = &
                  cgnsFamilies(ii)%BCTypeCGNS
             cgnsDoms(nZone)%bocoInfo(i)%BCType = &
                  cgnsFamilies(ii)%BCType

             !===========================================================

          case (UserDefined)

             ! A user defined boundary condition is prescribed.
             ! More information should be present. Determine the
             ! number of user defined data nodes.

             call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                  "ZoneBC_t", 1, "BC_t", i, "end")
             if(ierr /= CG_OK)             &
                  call terminate("readBocos", &
                  "Something wrong when calling cg_goto_f")

             call cg_nuser_data_f(nUserData, ierr)
             if(ierr /= CG_OK)            &
                  call terminate("readBocos", &
                  "Something wrong when calling &
                  &cg_nuser_data_f")

             ! nUserData should be 1. Check this.

             if(nUserData /= 1) then
                write(errorMessage,103) trim(cgnsDoms(nZone)%zoneName), &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
                if(myID == 0) call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)
             endif

             ! Read the name of the user defined data node.

             call cg_user_data_read_f(nUserData,                      &
                  cgnsDoms(nZone)%bocoInfo(i)%userDefinedName, &
                  ierr)
             if(ierr /= CG_OK)            &
                  call terminate("readBocos", &
                  "Something wrong when calling &
                  &cg_user_data_read_f")

             ! Determine the corresponding internal boundary
             ! condition from the name just read.

             cgnsDoms(nZone)%bocoInfo(i)%BCType = &
                  internalBC(cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS, &
                  cgnsDoms(nZone)%bocoInfo(i)%userDefinedName)

             ! Print an error message if the BC type was not recognized.

             if(cgnsDoms(nZone)%bocoInfo(i)%BCType == bcNull) then
                write(errorMessage,104) trim(cgnsDoms(nZone)%zoneName), &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName),        &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%userDefinedName)
                if(myID == 0) call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)
             endif

             ! At the moment the domain interfaces as well as the
             ! bleed flows are only possible on a per family basis.

             select case (cgnsDoms(nZone)%bocoInfo(i)%BCType)

             case (MassBleedInflow,    MassBleedOutflow,      &
                  DomainInterfaceAll, DomainInterfaceRhoUVW, &
                  DomainInterfaceP,   DomainInterfaceRho,    &
                  DomainInterfaceTotal)

                write(errorMessage,105) trim(cgnsDoms(nZone)%zoneName), &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName),       &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%userDefinedName)
                if(myID == 0) &
                     call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)

             end select

             ! Try to Read off the family name.
             call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                  "ZoneBC_t", 1, "BC_t", i, "end")

             cgnsDoms(nZone)%bocoInfo(i)%wallBCName = ""
             if (ierr == 0) then ! Node exits
                call cg_famname_read_f(familyName, ierr)
                if (ierr == 0) then
                   cgnsDoms(nZone)%bocoInfo(i)%wallBCName = familyName
                end if
             end if
             !===========================================================

          case default

             ! A standard CGNS boundary condition is used. Determine
             ! the internally used boundary condition.

             cgnsDoms(nZone)%bocoInfo(i)%userDefinedName = ""

             cgnsDoms(nZone)%bocoInfo(i)%BCType = &
                  internalBC(cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS, &
                  cgnsDoms(nZone)%bocoInfo(i)%userDefinedName)

             ! Print an error message if the BC type was not recognized.

             if(cgnsDoms(nZone)%bocoInfo(i)%BCType == bcNull) then
                write(errorMessage,106) trim(cgnsDoms(nZone)%zoneName), &
                     trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
                if(myID == 0) call terminate("readBocos", errorMessage)
                call mpi_barrier(ADflow_comm_world, ierr)
             endif

             !added to accomodate case where a grid family is specified
             ! but BC's are specified in standard CGNS Format
             cgnsDoms(nZone)%BCFamilies = .False.


             ! Try to Read off the Boundary Condition family
             ! name.
             call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                  "ZoneBC_t", 1, "BC_t", i, "end")

             cgnsDoms(nZone)%bocoInfo(i)%wallBCName = ""
             if (ierr == 0) then ! Node exits
                call cg_famname_read_f(familyName, ierr)
                if (ierr == 0) then
                   cgnsDoms(nZone)%bocoInfo(i)%wallBCName = familyName
                end if
             end if

          end select

       endif checkActualFace
       !
       !         Initialize slidingID to 0 to indicate that this boco does
       !         not belong to a sliding mesh interface. If it is, this will
       !         be overwritten after all bocos are read for every zone.
       !
       cgnsDoms(nzone)%bocoInfo(i)%slidingID = 0
       !
       !         Determine the possible rotating rate of the boundary face.
       !
       ! Initialize the rotating center and the rotating rates to
       ! the values of the corresponding block.

       cgnsDoms(nZone)%bocoInfo(i)%rotCenter = &
            cgnsDoms(nZone)%rotCenter
       cgnsDoms(nZone)%bocoInfo(i)%rotRate   = &
            cgnsDoms(nZone)%rotRate

       ! Check if a rotating rate was specified for the family to
       ! which this boundary face belongs.

       ii = cgnsDoms(nZone)%bocoInfo(i)%familyID
       if(cgnsDoms(nZone)%bocoInfo(i)%actualFace .and. ii > 0) then

          ! Check if a rotating rate was specified for the family ii.

          if( cgnsFamilies(ii)%rotatingFrameSpecified ) then

             ! Copy the rotation info from the family to the boundary
             ! face. Note the multiplication with LRef, because the
             ! rotation center for the families has not been scaled
             ! to meters.

             cgnsDoms(nZone)%bocoInfo(i)%rotCenter = &
                  cgnsDoms(nZone)%LRef*cgnsFamilies(ii)%rotCenter
             cgnsDoms(nZone)%bocoInfo(i)%rotRate   = &
                  cgnsFamilies(ii)%rotRate

          endif
       endif
       !
       !         Read and store the prescribed boundary condition data sets.
       !
       ! Initialize dataSetAllocated to .false. and nDataSet to 0.

       cgnsDoms(nZone)%bocoInfo(i)%dataSetAllocated = .false.
       cgnsDoms(nZone)%bocoInfo(i)%nDataSet         = 0

       ! Find out whether data has been specified for the
       ! corresponding family.

       familySpecifiedData = .false.
       ii = cgnsDoms(nZone)%bocoInfo(i)%familyID
       if(ii > 0) then
          if(cgnsFamilies(ii)%nDataSet > 0) &
               familySpecifiedData = .true.
       endif

       ! If family data is specified, set the pointer for the
       ! data sets.

       testFamilySpecified: if( familySpecifiedData ) then

          cgnsDoms(nZone)%bocoInfo(i)%nDataSet = &
               cgnsFamilies(ii)%nDataSet
          cgnsDoms(nZone)%bocoInfo(i)%dataSet => &
               cgnsFamilies(ii)%dataSet

       else testFamilySpecified

          ! No family specified stuff.
          ! Check if data sets are prescribed in the cgns file.

          testDataSets: if(cgnsNDataSet > 0) then

             ! Set dataSetAllocated to true., allocate the memory for
             ! dataSet and set the pointer to make the code more
             ! readable.

             cgnsDoms(nZone)%bocoInfo(i)%dataSetAllocated = .true.
             cgnsDoms(nZone)%bocoInfo(i)%nDataSet = cgnsNDataSet

             allocate(cgnsDoms(nZone)%bocoInfo(i)%dataSet(cgnsNDataSet), &
                  stat=ierr)
             if(ierr /= 0 ) &
                  call terminate("testDataSets", &
                  "Memory allocation failure for dataSet")

             dataSet => cgnsDoms(nZone)%bocoInfo(i)%dataSet

             ! Loop over the number of data sets to extract the data.

             loopDataSet: do j=1,cgnsNDataSet

                ! Find out what kind of data, dirichlet or neumann
                ! (or both), are stored in this data set.

                call cg_dataset_read_f(cgnsInd, cgnsBase, nZone, i, j, &
                     dataSet(j)%datasetName,         &
                     dataSet(j)%BCType,              &
                     dirichletFlag, neumannFlag, ierr)

                ! Nullify the dirichlet and neumann arrays and initialize
                ! the number of presribed data to 0.

                nullify(dataSet(j)%dirichletArrays)
                nullify(dataSet(j)%neumannArrays)

                dataSet(j)%nDirichletArrays = 0
                dataSet(j)%nNeumannArrays   = 0

                ! Read the dirichlet and neumann arrays if data
                ! is present.

                if(dirichletFlag == 1)                               &
                     call readBCDataArrays(dataSet(j)%nDirichletArrays, &
                     dataSet(j)%dirichletArrays,  &
                     Dirichlet)

                if(neumannFlag == 1)                               &
                     call readBCDataArrays(dataSet(j)%nNeumannArrays, &
                     dataSet(j)%neumannArrays,  &
                     Neumann)

             enddo loopDataSet
          endif testDataSets
       endif testFamilySpecified

    enddo bocoLoop

    ! Format statements.

101 format("Zone",1x,a,", boundary face", 1x,a, &
         ": Corresponding family name not given.")
102 format("Family name",1x,a,1x,"not present in the grid")
103 format("Zone",1x,a,", boundary face", 1x,a, &
         ": Need 1 UserDefinedData_t node for user defined &
         &boundary condition")
104 format("Zone",1x,a,", boundary face", 1x,a, &
         ": Unknown user-defined boundary condition", 1x,a)
105 format("Zone",1x,a,", boundary face", 1x,a, &
         ": User-defined boundary condition", 1x,a,1x, &
         "only possible for a family")
106 format("Zone",1x,a,", boundary face", 1x,a, &
         ": boundary condition type missing or not supported")

    !=================================================================

  contains

    !===============================================================

    subroutine readBCDataArrays(nArr, arr, DirNeu)
      !
      !         readBCDataArrays reads the arrays of the given data set
      !         from the cgns file.
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      integer, intent(in) :: DirNeu

      integer(kind=intType), intent(out) :: nArr
      type(cgnsBcdataArray), pointer, dimension(:) :: arr
      !
      !        Local variables.
      !
      integer :: ierr
      integer :: k, l, nArrays, realTypeCGNS
      integer :: mass, len, time, temp, angle, dataType

      integer(kind=intType) :: nn

      real(kind=cgnsRealType), dimension(:), allocatable :: tmp

      logical :: globalUnits

      ! Set the cgns real type.

      realTypeCGNS = setCGNSRealType()

      ! Go to the correct node of the given boundary subface.

      call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
           "ZoneBC_t", 1, "BC_t", i, "BCDataSet_t", j, &
           "BCData_t", DirNeu, "end")
      if(ierr /= CG_OK)                     &
           call terminate("readBCDataArrays", &
           "Something wrong when calling cg_goto_f")

      ! Determine the amount of data arrays present for this node.

      call cg_narrays_f(nArrays, ierr)
      if(ierr /= CG_OK)                   &
           call terminate("readBCDataArrays", &
           "Something wrong when calling cg_narrays_f")
      nArr = nArrays

      ! Allocate the memory for the prescribed data sets.

      allocate(arr(nArr), stat=ierr)
      if(ierr /= 0)                        &
           call terminate("readBCDataArrays", &
           "Memory allocation failure for arr")

      ! Initialize the units to si units.

      do k=1,nArrays
         arr(k)%mass  = Kilogram
         arr(k)%len   = Meter
         arr(k)%time  = Second
         arr(k)%temp  = Kelvin
         arr(k)%angle = Radian
      enddo

      ! Check if this "main" node contains info about the units.
      ! If so set the units of arr to these units.

      globalUnits = .false.
      call cg_units_read_f(mass, len, time, temp, angle, ierr)
      if(ierr == error) &
           call terminate("readBCDataArrays", &
           "Something wrong when calling cg_units_read_f")

      if(ierr == CG_OK) then
         globalUnits = .true.

         do k=1,nArrays
            arr(k)%mass  = mass
            arr(k)%len   = len
            arr(k)%time  = time
            arr(k)%temp  = temp
            arr(k)%angle = angle
         enddo
      endif

      ! Loop over the number of data arrays.

      loopDataArrays: do k=1,nArrays

         ! Go to the main node of the data arrays.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
              "ZoneBC_t", 1, "BC_t", i, "BCDataSet_t", j, &
              "BCData_t", DirNeu, "end")
         if(ierr /= CG_OK)                     &
              call terminate("readBCDataArrays", &
              "Something wrong when calling cg_goto_f")

         ! Determine the name and the dimensions of the array.

         call cg_array_info_f(k, arr(k)%arrayName, dataType,     &
              arr(k)%nDimensions, arr(k)%dataDim, &
              ierr)
         if(ierr /= CG_OK)                     &
              call terminate("readBCDataArrays", &
              "Something wrong when calling &
              &cg_array_info_f")

         ! Determine the total size of the data array and allocate
         ! the memory for dataArr and tmp, which is used to read
         ! the data.

         nn = arr(k)%dataDim(1)
         do l=2,arr(k)%nDimensions
            nn = nn*arr(k)%dataDim(l)
         enddo

         allocate(arr(k)%dataArr(nn), tmp(nn), stat=ierr)
         if(ierr /= 0)                          &
              call terminate("readBCDataArrays", &
              "Memory allocation failure for dataArr &
              &and tmp")

         ! Read the data, copy it from tmp into dataArr and
         ! deallocate tmp again.

         call cg_array_read_as_f(k, realTypeCGNS, tmp, ierr)
         if(ierr /= CG_OK)                     &
              call terminate("readBCDataArrays", &
              "Something wrong when calling &
              &cg_array_read_as_f")

         arr(k)%dataArr = tmp

         deallocate(tmp, stat=ierr)
         if(ierr /= 0) call terminate("loopDataArrays", &
              "Deallocation failure for tmp")

         ! Go to data array node to find out if the dimensions are
         ! specified here.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
              "ZoneBC_t", 1, "BC_t", i, "BCDataSet_t", j, &
              "BCData_t", DirNeu, "DataArray_t", k, "end")
         if(ierr /= CG_OK)                     &
              call terminate("readBCDataArrays", &
              "Something wrong when calling cg_goto_f")

         ! Try to read the units.

         call cg_units_read_f(mass, len, time, temp, angle, ierr)
         if(ierr == error) &
              call terminate("readBCDataArrays", &
              "Something wrong when calling &
              &cg_units_read_f")

         ! If the units are specified overwrite the currently
         ! stored units .

         if(ierr == CG_OK) then

            arr(k)%mass  = mass
            arr(k)%len   = len
            arr(k)%time  = time
            arr(k)%temp  = temp
            arr(k)%angle = angle

         else if(.not. globalUnits) then

            ! No local and global units specified. Processor 0 prints
            ! a warning.

            if(myID == 0) then

               print "(a)", "#"
               print "(a)", "#                      Warning"
               print 100, trim(cgnsDoms(nZone)%zoneName), &
                    trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
               print 101, trim(arr(k)%arrayName)
               print "(a)", "#"
100            format("# Zone ",a,", boundary subface ", a)
101            format("# BC data set ",a, &
                    ": No units specified, assuming SI units")
            endif
         endif

      enddo loopDataArrays

    end subroutine readBCDataArrays


  end subroutine readBocos

  !      ==================================================================

  logical function checkForDoubleBoundFace(nZone, nBound)
    !
    !       CheckForDoubleBoundFace checks whether the given boundary
    !       range for the given zone has already been defined in the 1 to
    !       1 block connectivities. If so .true. is returned, otherwise
    !       .false.
    !
    use cgnsGrid
    implicit none
    !
    !      Function arguments.
    !
    integer, intent(in) :: nZone, nBound
    !
    !      Local variables.
    !
    integer :: i
    integer, dimension(3,2) :: rangeBound, rangeFace

    ! Set the range for the boundary face.

    rangeBound(1,1) = min(cgnsDoms(nZone)%bocoInfo(nBound)%iBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%iEnd)
    rangeBound(1,2) = max(cgnsDoms(nZone)%bocoInfo(nBound)%iBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%iEnd)

    rangeBound(2,1) = min(cgnsDoms(nZone)%bocoInfo(nBound)%jBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%jEnd)
    rangeBound(2,2) = max(cgnsDoms(nZone)%bocoInfo(nBound)%jBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%jEnd)

    rangeBound(3,1) = min(cgnsDoms(nZone)%bocoInfo(nBound)%kBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%kEnd)
    rangeBound(3,2) = max(cgnsDoms(nZone)%bocoInfo(nBound)%kBeg, &
         cgnsDoms(nZone)%bocoInfo(nBound)%kEnd)

    ! Initialize checkForDoubleBoundFace to .false.

    checkForDoubleBoundFace = .false.

    ! Loop over the 1 to 1 block connectivities of this zone.

    do i=1,cgnsDoms(nZone)%n1to1

       ! Set the range for this subface.

       rangeFace(1,1) = min(cgnsDoms(nZone)%conn1to1(i)%iBeg, &
            cgnsDoms(nZone)%conn1to1(i)%iEnd)
       rangeFace(1,2) = max(cgnsDoms(nZone)%conn1to1(i)%iBeg, &
            cgnsDoms(nZone)%conn1to1(i)%iEnd)

       rangeFace(2,1) = min(cgnsDoms(nZone)%conn1to1(i)%jBeg, &
            cgnsDoms(nZone)%conn1to1(i)%jEnd)
       rangeFace(2,2) = max(cgnsDoms(nZone)%conn1to1(i)%jBeg, &
            cgnsDoms(nZone)%conn1to1(i)%jEnd)

       rangeFace(3,1) = min(cgnsDoms(nZone)%conn1to1(i)%kBeg, &
            cgnsDoms(nZone)%conn1to1(i)%kEnd)
       rangeFace(3,2) = max(cgnsDoms(nZone)%conn1to1(i)%kBeg, &
            cgnsDoms(nZone)%conn1to1(i)%kEnd)

       ! And do the check.

       if(rangeBound(1,1) == rangeFace(1,1) .and. &
            rangeBound(1,2) == rangeFace(1,2) .and. &
            rangeBound(2,1) == rangeFace(2,1) .and. &
            rangeBound(2,2) == rangeFace(2,2) .and. &
            rangeBound(3,1) == rangeFace(3,1) .and. &
            rangeBound(3,2) == rangeFace(3,2)) then

          ! Faces are identical. Set checkForDoubleBoundFace to
          ! .true. and exit the loop.

          checkForDoubleBoundFace = .true.
          exit

       endif

    enddo

  end function checkForDoubleBoundFace
  function internalBC(cgnsBocoType, userDefinedName)
    !
    !       internalBC determines the corresponding internally used
    !       boundary condition type for the given CGNS boundary condition.
    !       The flow equations to be solved are taken into account, e.g.
    !       a viscous wall BC for the Euler equations is set to an
    !       inviscid wall.
    !
    use constants
    use su_cgns
    use inputPhysics, only : equations, flowType
    implicit none
    !
    !      Function type.
    !
    integer(kind=intType) :: internalBC
    !
    !      Function argument.
    !
    integer, intent(in) :: cgnsBocoType  ! Note integer and not
    ! integer(intType).
    ! Because of cgns.

    character(len=maxCGNSNameLen), intent(in) :: userDefinedName

    ! Determine the CGNS boundary condition type and set
    ! internalBC accordingly.

    select case (cgnsBocoType)
    case (BCWallInviscid)
       internalBC = EulerWall

    case (BCWall, BCWallViscous, BCWallViscousHeatFlux)
       internalBC = NSWallAdiabatic
       if(equations == EulerEquations) internalBC = EulerWall

    case (BCWallViscousIsothermal)
       internalBC = NSWallIsothermal
       if(equations == EulerEquations) internalBC = EulerWall

    case (BCSymmetryPlane)
       internalBC = Symm

    case (BCSymmetryPolar)
       internalBC = SymmPolar

    case (BCExtrapolate, BCDegenerateLine, BCDegeneratePoint, &
         BCAxisymmetricWedge)
       internalBC = Extrap

    case (BCFarfield, BCInflow, BCOutflow)
       internalBC = FarField

    case (BCInflowSubsonic)
       internalBC = SubsonicInflow

    case (BCInflowSupersonic)
       internalBC = SupersonicInflow

    case (BCOutflowSubsonic)
       internalBC = SubsonicOutflow

    case (BCOutflowSupersonic)
       internalBC = SupersonicOutflow

    case (BCTunnelInflow)
       internalBC = SubsonicInflow

    case (BCTunnelOutflow)
       internalBC = SubsonicOutflow

    case (UserDefined)

       ! Select the internal type base on the user defined name.

       select case (trim(adjustl(userDefinedName)))

       case ("BCMassBleedInflow")
          internalBC = MassBleedInflow

       case ("BCMassBleedOutflow")
          internalBC = MassBleedOutflow

       case ("BCSlidingMesh")
          internalBC = SlidingInterface

       case ("BCOverset")
          internalBC = OversetOuterBound

       case ("BCDomainInterfaceAll")
          internalBC = DomainInterfaceAll

       case ("BCDomainInterfaceRhoUVW")
          internalBC = DomainInterfaceRhoUVW

       case ("BCDomainInterfaceP")
          internalBC = DomainInterfaceP

       case ("BCDomainInterfaceRho")
          internalBC = DomainInterfaceRho

       case ("BCDomainInterfaceTotal")
          internalBC = DomainInterfaceTotal

       case default
          internalBC = bcNull

       end select

    case default
       internalBC = bcNull

    end select

    ! Check whether the boundary conditions is allowed for the
    ! type of flow problem to be solved.

    select case (flowType)

    case (internalFlow)

       ! Internal flow. Not allowed to specify a farfield
       ! boundary condition.

       if(internalBC == FarField) internalBC = bcNotValid

    end select

  end function internalBC
  subroutine readPeriodicSubface(cgnsInd, cgnsBase, zone, conn,  &
       connectName, periodic,          &
       rotationCenter, rotationAngles, &
       translation)
    !
    !       readPeriodicSubface reads the possible periodic info for the
    !       given general subface connectivity.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use communication, only : adflow_comm_world, myid
    use utils, only : siAngle, terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer,          intent(in) :: cgnsInd, cgnsBase, zone, conn
    character(len=*), intent(in) :: connectName

    logical, intent(out) :: periodic

    real(kind=realType), dimension(3), intent(out) :: rotationCenter
    real(kind=realType), dimension(3), intent(out) :: rotationAngles
    real(kind=realType), dimension(3), intent(out) :: translation

    !
    !      Local variables.
    !
    integer :: ierr
    integer :: jj
    integer :: mass, len, time, temp, angle

    real(kind=realType), dimension(3) :: rotCenter, rotAngles
    real(kind=realType), dimension(3) :: tlation

    real(kind=realType) :: mult, trans

    ! Check if this is a periodic boundary.

    call cg_conn_periodic_read_f(cgnsInd, cgnsBase, zone, conn, &
         real(rotCenter,cgnsPerType), real(rotAngles,cgnsPerType), real(tlation,cgnsPerType), ierr)

    testPeriodic: if(ierr == CG_OK) then

       ! Subface is a periodic boundary. Check if the unit for
       ! the rotation angles is specified.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", zone, &
            "ZoneGridConnectivity_t", 1,             &
            "GridConnectivity_t", conn,              &
            "GridConnectivityProperty_t", 1,         &
            "Periodic_t", 1, "DataArray_t", 2, "end")
       if(ierr /= CG_OK)                      &
            call terminate("readPeriodicSubface", &
            "Something wrong when calling cg_goto_f")

       call cg_units_read_f(mass, len, time, temp, angle, ierr)
       if(ierr == error)                   &
            call terminate("readPeriodicSubface", &
            "Something wrong when calling cg_units_read_f")

       ! Check if the angle dimensions were specified.

       if(ierr == CG_OK .and. angle /= Null) then

          ! Determine the conversion factor to radians.

          call siAngle(angle, mult, trans)

       else

          ! Angle units not specified. Assume radians.
          ! Processor 0 writes a warning to stdout.

          if(myID == 0) then

             print "(a)", "#"
             print "(a)", "#                      Warning"
             print 100, trim(cgnsDoms(zone)%zonename), trim(connectName)
             print "(a)", "#"

100          format("# Zone",1X,A,", General connectivity",1X,A, &
                  ": No unit specified for periodic angles, &
                  &assuming radians.")

          endif

          ! Set mult to one.

          mult = one

       endif

       ! Store the info. Convert the rotation center and the
       ! translation vector to meters.

       periodic = .true.

       rotationCenter = rotCenter*cgnsDoms(zone)%LRef
       rotationAngles = rotAngles*mult
       translation    = tlation*cgnsDoms(zone)%LRef

       ! Make sure that the rotation angles are such that it
       ! corresponds to an integer value of the number of
       ! sections per wheel.

       mult = sqrt(rotationAngles(1)**2 &
            +      rotationAngles(2)**2 &
            +      rotationAngles(3)**2)

       if(mult > eps) then

          ! Nonzero angle specified. Determine the number of
          ! sections for the full wheel, which is an integer.

          jj = nint(two*pi/mult)

          ! Store the correction factor for the angles in
          ! mult and correct the periodic angles accordingly.

          mult = two*pi/(jj*mult)

          rotationAngles = mult*rotationAngles

       endif

    else testPeriodic

       ! Subface is a normal boundary. Set periodic to .false. and
       ! initialize the periodic data to zero to avoid possible
       ! problems due to uninitialized data.

       periodic = .false.

       rotationCenter = zero
       rotationAngles = zero
       translation    = zero

    endif testPeriodic

  end subroutine readPeriodicSubface

  subroutine readPeriodicSubface1to1(cgnsInd, cgnsBase, zone, conn,  &
       connectName, periodic,          &
       rotationCenter, rotationAngles, &
       translation)
    !
    !       readPeriodicSubface reads the possible periodic info for the
    !       given general subface connectivity.
    !
    use constants
    use su_cgns
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use communication, only : adflow_comm_world, myid
    use utils, only : siAngle, terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer,          intent(in) :: cgnsInd, cgnsBase, zone, conn
    character(len=*), intent(in) :: connectName

    logical, intent(out) :: periodic

    real(kind=realType), dimension(3), intent(out) :: rotationCenter
    real(kind=realType), dimension(3), intent(out) :: rotationAngles
    real(kind=realType), dimension(3), intent(out) :: translation

    !
    !      Local variables.
    !
    integer :: ierr
    integer :: jj
    integer :: mass, len, time, temp, angle

    real(kind=realType), dimension(3) :: rotCenter, rotAngles
    real(kind=realType), dimension(3) :: tlation

    real(kind=realType) :: mult, trans

    ! Check if this is a periodic boundary.
    call cg_1to1_periodic_read_f(cgnsInd, cgnsBase, zone, conn, &
         real(rotCenter,cgnsPerType), real(rotAngles,cgnsPerType), real(tlation,cgnsPerType), ierr)

    testPeriodic: if(ierr == CG_OK) then

       ! Subface is a periodic boundary. Check if the unit for
       ! the rotation angles is specified.

       call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", zone, &
            "ZoneGridConnectivity_t", 1,             &
            "GridConnectivity1to1_t", conn,              &
            "GridConnectivityProperty_t", 1,         &
            "Periodic_t", 3, "DataArray_t", 3, "end")
       if(ierr /= CG_OK)                      &
            call terminate("readPeriodicSubface1to1", &
            "Something wrong when calling cg_goto_f")

       call cg_units_read_f(mass, len, time, temp, angle, ierr)
       if(ierr == error)                   &
            call terminate("readPeriodicSubface1to1", &
            "Something wrong when calling cg_units_read_f")

       ! Check if the angle dimensions were specified.

       if(ierr == CG_OK .and. angle /= Null) then

          ! Determine the conversion factor to radians.

          call siAngle(angle, mult, trans)

       else

          ! Angle units not specified. Assume radians.
          ! Processor 0 writes a warning to stdout.

          if(myID == 0) then

             print "(a)", "#"
             print "(a)", "#                      Warning"
             print 100, trim(cgnsDoms(zone)%zonename), trim(connectName)
             print "(a)", "#"

100          format("# Zone",1X,A,", 1to1 connectivity",1X,A, &
                  ": No unit specified for periodic angles, &
                  &assuming radians.")

          endif

          ! Set mult to one.

          mult = one

       endif

       ! Store the info. Convert the rotation center and the
       ! translation vector to meters.

       periodic = .true.

       rotationCenter = rotCenter*cgnsDoms(zone)%LRef
       rotationAngles = rotAngles*mult
       translation    = tlation*cgnsDoms(zone)%LRef

       ! Make sure that the rotation angles are such that it
       ! corresponds to an integer value of the number of
       ! sections per wheel.

       mult = sqrt(rotationAngles(1)**2 &
            +      rotationAngles(2)**2 &
            +      rotationAngles(3)**2)

       if(mult > eps) then

          ! Nonzero angle specified. Determine the number of
          ! sections for the full wheel, which is an integer.

          jj = nint(two*pi/mult)

          ! Store the correction factor for the angles in
          ! mult and correct the periodic angles accordingly.

          mult = two*pi/(jj*mult)

          rotationAngles = mult*rotationAngles

       endif

    else testPeriodic

       ! Subface is a normal boundary. Set periodic to .false. and
       ! initialize the periodic data to zero to avoid possible
       ! problems due to uninitialized data.

       periodic = .false.

       rotationCenter = zero
       rotationAngles = zero
       translation    = zero

    endif testPeriodic

  end subroutine readPeriodicSubface1to1

  subroutine readGrid
    !
    !       readGrid reads the coordinates for the blocks or block parts
    !       to be stored on this processor.
    !
    use constants
    use cgnsNames
    use su_cgns
    use block, only : flowDoms, nDom
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use communication, only : adflow_comm_world, myID
    use flowVarRefState, only : LRef
    use inputIO, only : writeCoorMeter
    use IOModule, only : IOVar
    use partitionMod, only : nGridsRead, fileIDs, gridFiles
    use utils, only: setCGNSRealType, terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: cgnsInd, cgnsBase, cgnsZone
    integer :: j, nCoords
    integer :: ierr, realTypeCGNS, datatype

    integer(kind=cgsize_t), dimension(3) :: rangeMin, rangeMax

    integer(kind=intType) :: i, ii, jj, kk, ll, nn
    integer(kind=intType) :: il, jl, kl
    integer(kind=intType) :: typeMismatch

    character(len=7)              :: int1String
    character(len=2*maxStringLen) :: errorMessage
    character(len=maxCGNSNameLen) :: coordname

    real(kind=cgnsRealType), allocatable, dimension(:,:,:) :: buffer

    ! Set the cgns real type and initialize typeMismatch to 0.
    ! Set cgnsBase to 1, because we will always read from base 1;
    ! possible higher bases are ignored.

    realTypeCGNS = setCGNSRealType()
    typeMismatch = 0
    cgnsBase     = 1

    ! Loop over the number of blocks stored on this processor.

    domainLoop: do i=1,nDom

       ! Abbreviate the nodal block dimensions.

       il = flowDoms(i,1,1)%il
       jl = flowDoms(i,1,1)%jl
       kl = flowDoms(i,1,1)%kl

       ! Store the zone number a bit easier and set the range
       ! for reading the coordinates.

       cgnsZone = flowDoms(i,1,1)%cgnsBlockID

       rangeMin(1) = flowDoms(i,1,1)%iBegor
       rangeMin(2) = flowDoms(i,1,1)%jBegor
       rangeMin(3) = flowDoms(i,1,1)%kBegor

       rangeMax(1) = flowDoms(i,1,1)%iEndor
       rangeMax(2) = flowDoms(i,1,1)%jEndor
       rangeMax(3) = flowDoms(i,1,1)%kEndor

       ! Allocate the memory for the read buffer.

       allocate(buffer(il,jl,kl), stat=ierr)
       if(ierr /= 0)                &
            call terminate("readGrid", &
            "Memory allocation error for buffer")

       ! Loop over the number of grids to be read.

       nGridLoop: do nn=1,nGridsRead

          ! Store the file index a bit easier.

          cgnsInd = fileIDs(nn)

          ! Determine the number of coordinates in this zone.

          call cg_ncoords_f(cgnsInd, cgnsBase, cgnsZone, &
               nCoords, ierr)
          if(ierr /= CG_OK)           &
               call terminate("readGrid", &
               "Something wrong when calling cg_ncoords_f")

          ! The coordinates are only read if 3 coordinates are present.

          checkNcoords: if(nCoords == 3) then

             ! Loop over the number of coordinates. Note that the counter j
             ! is an integer. This is for compatibility with cgns.

             coords: do j=1,nCoords

                ! Get the info for this coordinate.

                call cg_coord_info_f(cgnsInd, cgnsBase, cgnsZone, j, &
                     datatype, coordname, ierr)
                if(ierr /= CG_OK)           &
                     call terminate("readGrid", &
                     "Something wrong when calling &
                     &cg_coord_info_f")

                ! Update the value of typeMismatch if the datatype of
                ! the coordinate is not equal to the datatype used in
                ! the solver.

                if(realTypeCGNS /= datatype) &
                     typeMismatch = typeMismatch + 1

                ! Set the value of the counter ll, depending on the name.
                ! Normally the x-coordinate is written first, followed by
                ! the y-coordinate and finally the z-coordinate. But you
                ! never know.

                select case(coordname)
                case (cgnsCoorX)
                   ll = 1
                case (cgnsCoorY)
                   ll = 2
                case (cgnsCoorZ)
                   ll = 3
                case default
                   write(errorMessage,110)                        &
                        trim(cgnsDoms(cgnsZone)%zoneName), &
                        trim(coordname)
110                format("Zone ",a," :Unknown coordinate name, ",a, &
                        ",in grid file")
                   call terminate("readGrid", errorMessage)
                end select

                ! Read the coordinates.

                call cg_coord_read_f(cgnsInd, cgnsBase, cgnsZone, &
                     coordname, realTypeCGNS,     &
                     rangeMin, rangeMax, buffer, ierr)
                if(ierr /= CG_OK)           &
                     call terminate("readGrid", &
                     "Something wrong when calling &
                     &cg_coord_read_f")

                ! Copy the data into IOVar and scale it to meters.
                ! The is effectivly copying into the variable x. w is just temporary through IOVar

                do kk=1,kl
                   do jj=1,jl
                      do ii=1,il
                         IOVar(i,nn)%w(ii,jj,kk,ll) = buffer(ii,jj,kk) &
                              * cgnsDoms(cgnsZone)%LRef
                         !print *,'iovar',IOVar(i,nn)%w(ii,jj,kk,ll),flowdoms(nn,1,1)%x(ii,jj,kk,ll),i,nn
                      enddo
                   enddo
                enddo

             enddo coords

          else checkNcoords

             ! There are not three coordinates present in this base.
             ! An error message is printed and an exit is made.

             write(errorMessage,101) trim(gridFiles(nn)),               &
                  trim(cgnsDoms(cgnsZone)%zoneName), &
                  nCoords
101          format("File ", a, ": The number of coordinates of zone ", &
                  a, " should be 3, not ", i1)

             call terminate("readGrid", errorMessage)

          endif checkNcoords

       enddo nGridLoop

       ! Release the memory of buffer.

       deallocate(buffer, stat=ierr)
       if(ierr /= 0) call terminate("readGrid", &
            "Deallocation error for buffer")
    enddo domainLoop

    ! Close the cgns files.

    do nn=1,nGridsRead
       call cg_close_f(fileIDs(nn), ierr)
       if(ierr /= CG_OK)           &
            call terminate("readGrid", &
            "Something wrong when calling cg_close_f")
    enddo

    ! Determine the global sum of typeMismatch; the result only
    ! needs to be known on processor 0. Use ii as the global buffer
    ! to store the result. If a type mismatch occured,
    ! print a warning.

    call mpi_reduce(typeMismatch, ii, 1, adflow_integer, &
         mpi_sum, 0, ADflow_comm_world, ierr)
    if(myID == 0 .and. ii > 0) then

       write(int1String,"(i6)") ii
       int1String = adjustl(int1String)

       print "(a)", "#"
       print "(a)", "#                      Warning"
       print 120, trim(int1String)
       print "(a)", "#"
120    format("# ",a," type mismatches occured when reading the &
            &coordinates of the blocks")
    endif

    ! If the coordinates in the solution files must be written in
    ! meters, correct this info for all cgns blocks.

    if( writeCoorMeter ) then
       do i=1,cgnsNDom
          cgnsDoms(i)%mass  = Null
          cgnsDoms(i)%len   = Meter
          cgnsDoms(i)%time  = Null
          cgnsDoms(i)%temp  = Null
          cgnsDoms(i)%angle = Null

          cgnsDoms(i)%gridUnitsSpecified = .true.
          cgnsDoms(i)%LRef = one
       enddo

       LRef = one
    endif

  end subroutine readGrid

end module readCGNSGrid
