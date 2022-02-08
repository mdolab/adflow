module partitioning

contains

  subroutine partitionAndReadGrid(partitionOnly)
    !
    !       partitionAndReadGrid determines the partitioning of the
    !       multiblock grid over the processors and reads the grid of the
    !       blocks (or block parts) assigned to this processor. Other
    !       preprocessing activities, such as the proper setup of the halo
    !       communication structure, creation of coarse grids and wall
    !       distance computation, are performed in the preprocessing
    !       library.
    !
    use constants
    use IOModule, only : IOVar
    use partitionMod, only : fileIDs, gridFiles
    use utils, only : terminate
    use loadBalance, only : loadBalanceGrid
    use gridChecking, only : checkFaces, check1to1Subfaces
    use readCGNSGrid, only : readBlockSizes, readGrid
    implicit none

    logical, intent(in) :: partitionOnly
    !
    !      Local variables
    !
    integer :: ierr

    ! Determine the number of grid files that must be read,
    ! as well as the corresponding file names.

    call determineGridFileNames

    ! Read the number of blocks and the block sizes of the grid stored
    ! in the cgns grid file. This info is stored on all processors

    call readBlockSizes
    call determineNeighborIDs

    ! Determine the grid sections.
    call determineSections

    ! If we are just doing a partition test, return
    if (partitionOnly) then
       return
    end if

    ! Determine the number of blocks to be stored on this processor
    ! and the relation to the original grid. Remember that blocks
    ! can be split for load balancing reasons.

    call loadBalanceGrid

    ! Initialize the iblank array for the fine grid domains

    call initFineGridIblank

    ! Allocate the coordinates of the fine grid level and the
    ! derived data type used to read them.

    call allocCoorFineGrid

    ! Read the grid of the blocks (block parts) to be stored
    ! on this processor.

    call readGrid

    ! Determine for the time spectral mode the time of one period,
    ! the rotation matrices for the velocity components and
    ! create the fine grid coordinates of all time spectral locations.

    call timePeriodSpectral
    call timeRotMatricesSpectral
    call fineGridSpectralCoor

    ! Release the memory of fileIDs, gridFiles and IOVar.
    ! They are not needed anymore.

    deallocate(fileIDs, gridFiles, IOVar, stat=ierr)
    if(ierr /= 0)                            &
         call terminate("partitionAndReadGrid", &
         "Deallocation failure for fileIDs, gridFiles &
         &and IOVar")

    ! Check if for all faces on the block boundaries either a
    ! physical boundary condition or a connectivity has been
    ! specified and check if the 1 to 1 subfaces match.

    call checkFaces
    call check1to1Subfaces

  end subroutine partitionAndReadGrid


  subroutine determineGridFileNames
    !
    !       determineGridFileNames determines the number and names of the
    !       files that contain the grids. For steady computations only one
    !       file must be present no matter if a restart is performed or
    !       not. For unsteady the situation is a little more complicated.
    !       If no restart is performed only one file must be present. If a
    !       restart is performed in unsteady or time spectral mode and a
    !       rigid body motion is prescribed again only one grid file is
    !       required; however for a consistent restart with deforming
    !       meshes the grids in the past must be read as well. If this is
    !       not possible only a first order restart can be made in
    !       unsteady mode and some kind of interpolation is used for the
    !       time spectral method.
    !
    use constants
    use communication, only: myID, adflow_comm_world
    use inputIO, only : gridFile
    use inputPhysics, only :equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputUnsteady, only : nOldGridRead
    use iteration, only : deforming_grid, nOldLevels
    use partitionMod, only : gridFiles, interpolSpectral, &
         nGridsRead, fileIDs
    use utils, only : terminate
    implicit none
    !
    !      Local variables
    !
    integer :: ierr

    integer(kind=intType) :: ii, nn, restartID

    character(len=7)            :: integerString
    character(len=maxStringLen) :: tmpName

    ! Initialization of nOldGridRead and interpolSpectral.

    nOldGridRead     = 1
    interpolSpectral = .true.

    ! Determine the desired number of files from which grids at
    ! certain time levels should be read.  This depends on the
    ! equation mode we have to solve for. Also set the corresponding
    ! file names.

    select case(equationMode)

    case (steady)

       ! Steady computation. Only one grid needs to be read.

       nGridsRead = 1
       allocate(fileIDs(nGridsRead), &
            gridFiles(nGridsRead), stat=ierr)
       if(ierr /= 0)                              &
            call terminate("determineGridFileNames", &
            "Memory allocation failure for fileIDs &
            &and gridFiles")

       gridFiles(1) = gridFile

       !===============================================================

    case (unsteady)

       ! Unsteady computation. A further check is required.
       ! EJ: replaced boolean variable restart with .false. for now
       ! Need to refactor this code as well with ALE restart

       testMultipleUnsteady: if(deforming_Grid .and. .false.) then

          ! A restart is made with deforming meshes. For a consistent
          ! restart nOldLevels grids must be read. First determine
          ! the prefix of the grid file and the time step number
          ! from which a restart should be made.

          ii = len_trim(gridFile)
          do
             if(gridFile(ii:ii) < "0" .or. gridFile(ii:ii) > "9") exit
             ii = ii - 1
          enddo

          ! If the last characters of the file name do not contain a
          ! number, the grid file does not come from a previous
          ! unsteady deforming mesh computation and therefore only
          ! one grid will be read.

          if(ii == len_trim(gridFile)) then

             nGridsRead = 1
             allocate(fileIDs(nGridsRead), &
                  gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
                  call terminate("determineGridFileNames", &
                  "Memory allocation failure for fileIDs &
                  &and gridFiles")

             gridFiles(1) = gridFile

          else

             ! Read the integer number from the last characters
             ! of the grid file.

             read(gridFile(ii+1:),*) restartID

             ! Allocate the memory for the file names and set them.

             nGridsRead = nOldLevels
             allocate(fileIDs(nGridsRead), &
                  gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
                  call terminate("determineGridFileNames", &
                  "Memory allocation failure for fileIDs &
                  &and gridFiles")

             do nn=1,nGridsRead
                write(integerString,"(i6)") restartID - nn + 1
                integerString = adjustl(integerString)
                gridFiles(nn) = gridFile(:ii)//trim(integerString)
             enddo

          endif

       else testMultipleUnsteady

          ! The computation either starts from scratch or an unsteady
          ! restart for a rigid grid (possibly moving) is made. In
          ! all cases only one grid file is needed.

          nGridsRead = 1
          allocate(fileIDs(nGridsRead), &
               gridFiles(nGridsRead), stat=ierr)
          if(ierr /= 0)                              &
               call terminate("determineGridFileNames", &
               "Memory allocation failure for fileIDs &
               &and gridFiles")

          gridFiles(1) = gridFile

       endif testMultipleUnsteady

       ! Check if the files can be opened.

       do nn=1,nGridsRead
          open(unit=21,file=gridFiles(nn),status="old",iostat=ierr)
          if(ierr /= 0) exit
          close(unit=21)
       enddo

       ! Possibly correct nGridsRead and set nOldGridRead.
       ! If nOldGridRead == 0, i.e. not a valid grid is found,
       ! print an error message and terminate.

       nGridsRead   = nn - 1
       nOldGridRead = nGridsRead

       if(nOldGridRead == 0) then
          if(myID == 0)                              &
               call terminate("determineGridFileNames", &
               "Grid file(s) could not be opened")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       !===============================================================

    case (timeSpectral)

       ! Time spectral computation. A further check is required.
       ! EJ: replaced boolean variable restart with .false. for now
       ! Need to refactor this code as well with ALE restart

       testMultipleTS: if(deforming_Grid .and. .false.) then

          ! A restart is made with deforming meshes. For a consistent
          ! restart multiple grids must be read. First determine the
          ! the prefix of the grid file from which a restart should
          ! be made.

          ii = len_trim(gridFile)
          do
             if(gridFile(ii:ii) < "0" .or. gridFile(ii:ii) > "9") exit
             ii = ii - 1
          enddo

          ! If the last characters of the file name do not contain a
          ! number, the grid file does not come from a previous
          ! time spectral deforming mesh computation and therefore
          ! only one grid will be read.

          if(ii == len_trim(gridFile)) then

             nGridsRead = 1
             allocate(fileIDs(nGridsRead), &
                  gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
                  call terminate("determineGridFileNames", &
                  "Memory allocation failure for fileIDs &
                  &and gridFiles")

             gridFiles(1) = gridFile

          else

             ! Loop to find out how many time instances were used in
             ! the previous computation from which a restart is made.

             nn = 0
             do
                nn = nn + 1
                write(integerString,"(i6)") nn
                integerString = adjustl(integerString)
                tmpName = gridFile(:ii)//trim(integerString)

                open(unit=21,file=tmpName,status="old",iostat=ierr)
                if(ierr /= 0) exit
                close(unit=21)
             enddo

             nn = nn - 1

             ! Take care of the exceptional situation that nn == 0.
             ! This happens when the restart file ends at with an
             ! integer, but does not correspond to a time spectral
             ! solution. Allocate the memory.

             nGridsRead = max(nn, 1_intType)
             allocate(fileIDs(nGridsRead), &
                  gridFiles(nGridsRead), stat=ierr)
             if(ierr /= 0)                              &
                  call terminate("determineGridFileNames", &
                  "Memory allocation failure for fileIDs &
                  &and gridFiles")

             if(nn == 0) then
                gridFiles(1) = gridFile
             else
                do nn=1,nGridsRead
                   write(integerString,"(i6)") nn
                   integerString = adjustl(integerString)
                   gridFiles(nn) = gridFile(:ii)//trim(integerString)
                enddo
             endif

             ! Check whether or not the coordinates must be interpolated,
             ! i.e. check if nGridsRead == nTimeIntervalsSpectral.

             if(nGridsRead == nTimeIntervalsSpectral) &
                  interpolSpectral = .false.

          endif

       else testMultipleTS

          ! The computation either starts from scratch or a
          ! restart for a rigid grid (possibly moving) is made. In
          ! all cases only one grid file is needed.

          nGridsRead = 1
          allocate(fileIDs(nGridsRead), &
               gridFiles(nGridsRead), stat=ierr)
          if(ierr /= 0)                              &
               call terminate("determineGridFileNames", &
               "Memory allocation failure for fileIDs &
               &and gridFiles")

          gridFiles(1) = gridFile

       endif testMultipleTS

    end select

  end subroutine determineGridFileNames

  subroutine determineNeighborIDs
    !
    !       determineNeighborIDs determines for every internal block
    !       boundary the block ID of the neighbor. In the cgns file only
    !       the zone name is stored, but the ID's are more useful
    !       internally.
    !       Although for this case a quadratic search algorithm is not too
    !       bad (number of blocks are O(1000) maximum), I don't like the
    !       idea of having a quadratic loop in the code. That's why a
    !       O(n log(n)) algorithm is used here.
    !
    use constants
    use cgnsGrid, only : cgnsDoms, cgnsNDOm
    use utils, only : terminate
    use sorting, only: qsortstrings, bsearchstrings
    implicit none
    !
    !      Local variables
    !
    character(len=maxCGNSNameLen), dimension(cgnsNDom) :: zoneNames

    integer(kind=intType), dimension(cgnsNDom) :: zoneNumbers

    integer(kind=intType) :: i, j, k, ii
    !
    ! Copy the zone name from the derived data type into zoneNames.

    do i=1,cgnsNDom
       zoneNames(i) = cgnsDoms(i)%zoneName
    enddo

    ! Sort zoneNames in increasing order.

    call qsortStrings(zoneNames, cgnsNDom)

    ! Initialize zoneNumbers to -1. This serves as a check during
    ! the search.

    zoneNumbers = -1

    ! Find the original zone ids for the sorted zone names.

    do i=1,cgnsNDom
       ii = bsearchStrings(cgnsDoms(i)%zoneName, zoneNames)

       ! Check if the zone number is not already taken. If this is the
       ! case, this means that the grid file contains two identical
       ! zone names.

       if(zoneNumbers(ii) /= -1)                &
            call terminate("determineNeighborIDs", &
            "Error occurs only when two identical zone &
            &names are present")

       ! And set the zone number.

       zoneNumbers(ii) = i
    enddo

    ! Loop over the blocks and its connectivities to find out the
    ! neighbors.

    domains: do i=1,cgnsNDom

       ! The 1-to-1 connectivities.

       do j=1,cgnsDoms(i)%n1to1

          ! Determine the neighbor ID for this internal block boundary.

          ii = bsearchStrings(cgnsDoms(i)%conn1to1(j)%donorName, zoneNames)
          if(ii == 0)                              &
               call terminate("determineNeighborIDs", &
               "donor name not found in sorted zone names")

          cgnsDoms(i)%conn1to1(j)%donorBlock = zoneNumbers(ii)

       enddo

       ! The non-matching abutting connectivities.

       do j=1,cgnsDoms(i)%nNonMatchAbutting

          ! Determine the neighbor ID's for this subface.

          do k=1,cgnsDoms(i)%connNonMatchAbutting(j)%nDonorBlocks

             ii = bsearchStrings(&
                  cgnsDoms(i)%connNonMatchAbutting(j)%donorNames(k), zoneNames)
             if(ii == 0)                              &
                  call terminate("determineNeighborIDs", &
                  "donor name not found in sorted zone names")

             cgnsDoms(i)%connNonMatchAbutting(j)%donorBlocks(k) = &
                  zoneNumbers(ii)
          enddo
       enddo

    enddo domains

  end subroutine determineNeighborIDs


  subroutine determineInterfaceIDs
    !
    !       DetermineInterfaceIDs determines more information for both the
    !       sliding mesh and domain interfaces with other codes, which are
    !       both specified as user defined boundary conditions in CGNS. In
    !       particular the number of sliding mesh interfaces and their
    !       pairings, and number of interfaces with other codes are
    !       determined.
    !       There are some parts in the coupler API routines where
    !       the family-specified domain interfaces are implicitly assumed.
    !       Therefore, it is recommended for the time being that all the
    !       domain interfaces should be family-specified.
    !
    use constants
    use cgnsGrid, only : cgnsDoms, cgnsNDom, cgnsFamilies, &
         famIDsDomainInterfaces, bcIDsDomainInterfaces, famIDsSliding, &
         cgnsNFamilies, cgnsNSliding, cgnsNDomainInterfaces
    use communication, only : adflow_comm_world, myid
    use iteration, only : standAloneMode
    use utils, only : terminate
    use sorting, only: qsortstrings, bsearchstrings

    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: ii, jj, mm, nn
    integer(kind=intType) :: nSlidingFam, nSlidingBC, nSliding
    integer(kind=intType) :: nDomainFam,  nDomainBC

    integer(kind=intType), dimension(:),   allocatable :: famSlidingID
    integer(kind=intType), dimension(:,:), allocatable :: bcSlidingID
    integer(kind=intType), dimension(:),   allocatable :: orID

    character(len=maxStringLen) :: errorMessage

    character(len=maxCGNSNameLen), dimension(:), allocatable :: &
         namesSliding
    character(len=maxCGNSNameLen), dimension(:), allocatable :: &
         namesSorted

    logical :: validInterface
    !
    ! Count the total number of each type of user-defined BC that
    ! needs to be treated here. Note that if a BC is specified by the
    ! the family it is not counted, such that BCs making up a user-
    ! defined family are only counted once.

    nSlidingFam = 0
    nSlidingBC  = 0
    nDomainFam  = 0
    nDomainBC   = 0

    do nn=1,cgnsNFamilies

       select case (cgnsFamilies(nn)%BCType)
       case (SlidingInterface)
          nSlidingFam = nSlidingFam + 1
       case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
            DomainInterfaceP,   DomainInterfaceRho,    &
            DomainInterfaceTotal)
          nDomainFam = nDomainFam + 1
       end select

    enddo

    do nn=1,cgnsNDom
       do mm=1,cgnsDoms(nn)%nBocos
          if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
               cgnsDoms(nn)%bocoInfo(mm)%familyID == 0) then

             select case (cgnsDoms(nn)%bocoInfo(mm)%BCType)
             case (SlidingInterface)
                nSlidingBC = nSlidingBC + 1
             case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                  DomainInterfaceP  , DomainInterfaceRho,    &
                  DomainInterfaceTotal)
                nDomainBC = nDomainBC + 1
             end select

          end if
       end do
    end do

    ! Domain interfaces are only allowed when the code is run in a
    ! multi-disciplinary environment. Check this.

    cgnsNDomainInterfaces = nDomainFam + nDomainBC

    if(standAloneMode .and. cgnsNDomainInterfaces > 0) then
       if(myID == 0) &
            call terminate("determineInterfaceIDs", &
            "Domain interfaces are not allowed in &
            &stand alone mode.")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! The number of sliding mesh interfaces must be even, because each
    ! sliding interface should have two sides. Check this.

    nSliding = nSlidingFam + nSlidingBC
    if(mod(nSliding,2) == 1) then
       if(myID == 0)                              &
            call terminate("determineInterfaceIDs", &
            "Odd number of sliding mesh families found")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Allocate memory to store the names and IDs for the interfaces.

    allocate(famIDsDomainInterfaces(nDomainFam),            &
         bcIDsDomainInterfaces(2,nDomainBC),            &
         namesSliding(nSliding), namesSorted(nSliding), &
         orID(nSliding), famSlidingID(nSlidingFam),     &
         bcSlidingID(2,nSlidingBC), stat=ierr)
    if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
         "Memory allocation failure for names, IDs")

    ! Loop back over the families again and store the names and
    ! IDs this time around. Note the ID is just the index of the
    ! corresponding family.

    ii = 0
    jj = 0

    do nn=1,cgnsNFamilies

       select case (cgnsFamilies(nn)%BCType)
       case (SlidingInterface)
          ii = ii + 1
          namesSliding(ii) = cgnsFamilies(nn)%familyName
          famSlidingID(ii) = nn
       case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
            DomainInterfaceP,   DomainInterfaceRho,    &
            DomainInterfaceTotal)
          jj = jj + 1
          famIDsDomainInterfaces(jj) = nn
       end select

    enddo

    ! Loop back over the boundary conditions again and store the
    ! names and IDs this time around. Note the ID has two parts:
    ! the domain and the index of the BC info in that domain.

    jj = 0

    do nn=1,cgnsNDom
       do mm=1,cgnsDoms(nn)%nBocos
          if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
               cgnsDoms(nn)%bocoInfo(mm)%familyID == 0) then

             select case (cgnsDoms(nn)%bocoInfo(mm)%BCType)
             case (SlidingInterface)
                ii = ii + 1
                namesSliding(ii) = cgnsDoms(nn)%bocoInfo(mm)%bocoName
                bcSlidingID(:,ii-nSlidingFam) = (/ nn, mm /)
             case (DomainInterfaceAll, DomainInterfaceRhoUVW, &
                  DomainInterfaceP,   DomainInterfaceRho,    &
                  DomainInterfaceTotal)
                jj = jj + 1
                bcIDsDomainInterfaces(:,jj) = (/ nn, mm /)
             end select

          end if
       end do
    end do

    ! Initialize orID to -1, which serves as a check later on, and
    ! copy the names of the sliding mesh families in namesSorted.

    do ii=1,nSliding
       orID(ii)        = -1
       namesSorted(ii) = namesSliding(ii)
    enddo

    ! Sort the names of the sliding mesh families in increasing
    ! order and find the corresponding entry in namesSliding.

    call qsortStrings(namesSorted, nSliding)

    do ii=1,nSliding

       ! Search the sorted strings.

       mm = bsearchStrings(namesSliding(ii), namesSorted)

       if(orID(mm) /= -1) then

          ! Family name occurs more than once. This is not allowed.

          write(errorMessage,101) trim(namesSliding(ii))
          if(myID == 0) &
               call terminate("determineInterfaceIDs", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)

       endif

       ! Store the entry in orID.

       orID(mm) = ii

    enddo

    ! Set the number of sliding mesh interfaces and allocate the
    ! memory for famIDsSliding.

    cgnsNSliding = nSliding/2
    allocate(famIDsSliding(cgnsNSliding,2), stat=ierr)
    if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
         "Memory allocation failure for famIDsSliding")

    ! Check if the sorted family names indeed form a
    ! sliding mesh interface.

    do ii=1,cgnsNSliding

       ! Store the entries in namesSorted, which should form a sliding
       ! interface, a bit easier.

       nn = 2*ii
       mm = nn - 1

       ! Store the original values in famIDsSliding.

       jj = orID(nn)/2

       famIDsSliding(jj,1) = famSlidingID(orID(mm))
       famIDsSliding(jj,2) = famSlidingID(orID(nn))

       ! Check if the names form a valid sliding interface.

       validInterface = .true.
       if(len_trim(namesSorted(nn)) /= len_trim(namesSorted(mm))) &
            validInterface = .false.

       jj = len_trim(namesSorted(nn)) - 2
       if(jj <= 0) then
          validInterface = .false.
       else if(namesSorted(nn)(:jj) /= namesSorted(mm)(:jj)) then
          validInterface = .false.
       endif

       ! Print an error message and exit if the two families do not
       ! form a valid sliding interface.

       if(.not. validInterface) then
          write(errorMessage,102) trim(namesSliding(nn)), &
               trim(namesSliding(mm))
          if(myID == 0) &
               call terminate("determineInterfaceIDs", errorMessage)
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! The two names form a valid sliding interface. Store the
       ! sliding interface ID for the original family or BC. Note that
       ! one gets a positive ID and the other one a negative ID.
       ! In this way a distinction is made between the two sides.

       if (orID(mm) > nSlidingFam) then
          jj = bcSlidingID(2,orID(mm))
          mm = bcSlidingID(1,orID(mm))
          cgnsDoms(mm)%bocoInfo(jj)%slidingID = -ii
       else
          mm = famSlidingID(orID(mm))
          cgnsFamilies(mm)%slidingID = -ii
       end if

       if (orID(nn) > nSlidingFam) then
          jj = bcSlidingID(2,orID(nn))
          nn = bcSlidingID(1,orID(nn))
          cgnsDoms(nn)%bocoInfo(jj)%slidingID = ii
       else
          nn = famSlidingID(orID(nn))
          cgnsFamilies(nn)%slidingID = ii
       end if

    enddo

    ! Loop over all the boundary conditions again, this time
    ! copying the sliding interface IDs for those BCs that
    ! were part of a sliding family.

    do nn=1,cgnsNDom
       do mm=1,cgnsDoms(nn)%nBocos
          if (cgnsDoms(nn)%bocoInfo(mm)%actualFace .and.   &
               cgnsDoms(nn)%bocoInfo(mm)%familyID > 0 .and. &
               cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             cgnsDoms(nn)%bocoInfo(mm)%slidingID = &
                  cgnsFamilies(cgnsDoms(nn)%bocoInfo(mm)%familyID)%slidingID

          end if
       end do
    end do

    ! Deallocate the memory used to determine the sliding IDs.

    deallocate(namesSliding, namesSorted, orID, &
         famSlidingID, bcSlidingID, stat=ierr)
    if(ierr /= 0) &
         call terminate("determineInterfaceIDs", &
         "Deallocation failure for names, IDs")

    ! Format statements.

101 format("Family name", 1X,A,1X,"occurs more than once.")
102 format("Family names", 1X,A,1X,"and",1X,A,1X, &
         "do not form a valid sliding mesh interface")

  end subroutine determineInterfaceIDs

  subroutine initFineGridIblank
    !
    !       InitFineGridIblank allocates the fine grid iblank array and
    !       initializes the values for the holes, boundary, and halos. The
    !       holes read into the cgns domains are distributed amongst its
    !       sublocks in the form of iblanks. That is, we do not store a
    !       list of indices for the holes of the flow domains as done in
    !       the CGNS. The number of holes in each domain are also counted.
    !
    use constants
    use block, only : nDom, flowDoms
    use utils, only : terminate
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, l, m, n, cgnsId, sps

    ! Loop over the local blocks.
    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do n = 1,nDom

          ! Allocate memory for the iblank array of this block.

          i = flowDoms(n,1,sps)%ib
          j = flowDoms(n,1,sps)%jb
          k = flowDoms(n,1,sps)%kb
          allocate(flowDoms(n,1,sps)%iblank(0:i,0:j,0:k), &
               flowDoms(n,1,sps)%forcedRecv(0:i,0:j,0:k), &
               flowDoms(n,1,sps)%status(0:i,0:j,0:k), &
               stat=ierr)
          if(ierr /= 0)                             &
               call terminate("initFineGridIblank", &
               "Memory allocation failure for iblank")

          ! Initialize iblank to 1 everywhere, and the number of holes
          ! for this domain to 0.

          flowDoms(n,1,sps)%iblank = 1
          flowDoms(n,1,sps)%forcedRecv = 0
          flowDoms(n,1,sps)%status = 0
       end do domains
    end do spectralLoop

  end subroutine initFineGridIblank

  subroutine timePeriodSpectral
    !
    !       timePeriodSpectral determines the time of one period for the
    !       time spectral method. It is possible that sections have
    !       different periodic times.
    !
    use constants
    use communication, only : myID, adflow_comm_world
    use inputMotion, only : degreeFourXRot, degreeFourYRot, degreeFourZrot, &
         omegaFourAlpha, omegaFourBeta, omegaFourMach, omegaFourXRot, &
         omegaFourYRot, omegaFourZRot, degreeFourMach, degreeFourAlpha, &
         degreeFourBeta
    use inputPhysics, only : equationMOde, flowType
    use inputTimeSpectral, only : omegaFourier
    use section, only : sections, nSections
    use utils, only : terminate
    implicit none
    !
    !      Local parameter.
    !
    real(kind=realType), parameter :: tol  = 1.e-5_realType
    !
    !      Local variables.
    !
    integer :: ierr
    integer(kind=intType) :: nn

    real(kind=realType) :: tt, omega
    real(kind=realType) :: timePeriod

    logical :: timeDetermined

    ! This routine is only used for the spectral solutions. Return
    ! immediately if a different mode is solved.

    if(equationMode /= timeSpectral) return

    ! First check if a rotational frequency has been specified.
    ! Only for external flows.

    timeDetermined = .false.


    externalTest: if(flowType == externalFlow) then

       ! X-rotation.

       if(degreeFourXRot > 0) then
          timePeriod     = two*pi/omegaFourXRot
          timeDetermined = .true.
       endif

       ! Y-rotation.

       if(degreeFourYRot > 0) then
          tt = two*pi/omegaFourYRot

          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       endif

       ! Z-rotation.

       if(degreeFourZRot > 0) then
          tt = two*pi/omegaFourZRot

          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       endif


       ! Alpha
       !print *,'degreeFourAlpha',degreefouralpha,omegafouralpha,sincoeffouralpha
       if(degreeFourAlpha > 0) then
          tt = two*pi/omegaFourAlpha
          !print *,'timePeriod',tt
          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       endif


       ! Beta

       if(degreeFourBeta > 0) then
          tt = two*pi/omegaFourBeta

          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       endif

       ! Mach

       if(degreeFourMach > 0) then
          tt = two*pi/omegaFourMach

          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       endif

       ! aeroelastic case
       if(omegaFourier > 0) then 
         tt = two*pi/omegaFourier

          ! Check if a time period was already determined. If so, try
          ! to determine a common time. Otherwise just copy the data.

          if( timeDetermined ) then
             timePeriod = commonTimeSpectral(timePeriod, tt)
          else
             timePeriod     = tt
             timeDetermined = .true.
          endif
       end if

!!$         ! Altitude.
!!$
!!$         if(degreeFourAltitude > 0) then
!!$           tt = two*pi/omegaFourAltitude
!!$
!!$           ! Check if a time period was already determined. If so, try
!!$           ! to determine a common time. Otherwise just copy the data.
!!$
!!$           if( timeDetermined ) then
!!$             timePeriod = commonTimeSpectral(timePeriod, tt)
!!$           else
!!$             timePeriod     = tt
!!$             timeDetermined = .true.
!!$           endif
!!$         endif

    endif externalTest

    ! If it was possible to determine the time, copy it to the
    ! sections and return.


    if( timeDetermined ) then
       do nn=1,nSections
          sections(nn)%timePeriod = timePeriod/sections(nn)%nSlices
          !print *,'sectionTimePeriod',sections(nn)%timePeriod,nn
       enddo
       return
    endif

    ! Try to determine the periodic time via the rotation rate of the
    ! sections and its number of slices.

    sectionLoop: do nn=1,nSections

       ! Test if the section is rotating, because only for rotating
       ! sections the periodic time can be determined.

       testRotating: if( sections(nn)%rotating ) then

          ! Determine the magnitude of the rotation rate and the
          ! corresponding periodic time period.

          omega = sqrt(sections(nn)%rotRate(1)**2 &
               +      sections(nn)%rotRate(2)**2 &
               +      sections(nn)%rotRate(3)**2)

          tt = two*pi/omega

          ! If a time period was already determined, check if this is
          ! identical to tt. If not print an error message and exit.

          if( timeDetermined ) then

             tt = abs(tt-timePeriod)/timePeriod
             if(tt > tol) then
                if(myID == 0)                          &
                     call terminate("timePeriodSpectral", &
                     "Rotational frequencies of the rotating &
                     &sections are not identical.")
                call mpi_barrier(ADflow_comm_world, ierr)
             endif

          else

             ! Just copy the data.

             timePeriod     = tt
             timeDetermined = .true.

          endif

       endif testRotating
    enddo sectionLoop

    ! Divide the periodic time by the number of slices to get the
    ! characteristic time for every section.

    do nn=1,nSections
       sections(nn)%timePeriod = timePeriod/sections(nn)%nSlices
    enddo

    ! Return if it was possible to determine the time.

    if( timeDetermined ) return

    ! Periodic time could not be determined. Print an error
    ! message and exit.

    if(myID == 0)                          &
         call terminate("timePeriodSpectral", &
         "Not possible to determine the periodic time &
         &for the time spectral method")
    call mpi_barrier(ADflow_comm_world, ierr)

  end subroutine timePeriodSpectral

  function commonTimeSpectral(t1, t2)
    !
    !       The function commonTimeSpectral determines the smallest
    !       possible common time between t1 and t2, such that
    !       tcommon = n1*t1 = n2*t2 and n1, n2 integers.
    !
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Function definition
    !
    real(kind=realType) :: commonTimeSpectral
    !
    !      Function arguments.
    !
    real(kind=realType), intent(in) :: t1, t2
    !
    !      Local parameters.
    !
    integer(kind=intType), parameter :: nMax = 100
    real(kind=realType),   parameter :: tol  = 1.e-5_realType
    !
    !      Local variables.
    !
    integer                :: ierr
    integer(kind=intType) :: n1, n2
    real(kind=realType)   :: tt1, tt2, tt, ratio

    ! Store the largest time in tt1 and the smallest in tt2 and
    ! compute the ratio tt1/tt2, which is >= 1

    tt1   = max(t1, t2)
    tt2   = min(t1, t2)
    ratio = tt1/tt2

    ! Loop to find the smallest integer values of n1 and n2, such
    ! that n1*tt1 = n2*tt2. Note that due to the previous definition
    ! n2 >= n1.

    do n1=1,nMax
       tt = n1*ratio
       n2 = nint(tt)
       tt = abs(tt-n2)
       if(tt <= tol) exit
    enddo

    ! Check if a common time was found

    if(n1 > nMax) then
       if(myID == 0)                          &
            call terminate("commonTimeSpectral", &
            "No common time periodic time found. &
            &Problem may not be periodic")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Set the common time.

    commonTimeSpectral = n1*tt1

  end function commonTimeSpectral
  subroutine timeRotMatricesSpectral
    !
    !       timeRotMatricesSpectral determines the rotation matrices
    !       used in the time derivatives for the velocity components in
    !       the time spectral method. These matrices are the identity
    !       matrices for non-rotating sections and something different for
    !       rotating sections. Therefore the rotation matrices are stored
    !       for every section.
    !
    use constants
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : rotMatrixSpectral
    use section, only : sections, nSections
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn
    real(kind=realType)   :: tmp, theta, cosTheta, sinTheta

    real(kind=realType), dimension(3) :: xt, yt, zt

    ! This routine is only used for the spectral solutions. Return
    ! immediately if a different mode is solved.

    if(equationMode /= timeSpectral) return

    ! Allocate the memory for rotMatrixSpectral, which will store
    ! the rotation matrices for all the sections.

    if( allocated(rotMatrixSpectral)) deallocate(rotMatrixSpectral)
    allocate(rotMatrixSpectral(nSections,3,3), stat=ierr)

    if(ierr /= 0)                               &
         call terminate("timeRotMatricesSpectral", &
         "Memory allocation failure for &
         &rotMatrixSpectral")

    ! Loop over the number of sections.

    sectionLoop: do nn=1,nSections

       ! Test if the rotation matrix is the unity matrix. This is the
       ! case if this section is not rotating or if the number of
       ! slices is only 1, i.e. the true physical model is computed.

       testUnity: if(.not. sections(nn)%rotating .or. &
            sections(nn)%nSlices == 1) then

          ! Set the rotation matrix to the unity matrix.

          rotMatrixSpectral(nn,1,1) = one
          rotMatrixSpectral(nn,1,2) = zero
          rotMatrixSpectral(nn,1,3) = zero

          rotMatrixSpectral(nn,2,1) = zero
          rotMatrixSpectral(nn,2,2) = one
          rotMatrixSpectral(nn,2,3) = zero

          rotMatrixSpectral(nn,3,1) = zero
          rotMatrixSpectral(nn,3,2) = zero
          rotMatrixSpectral(nn,3,3) = one

       else testUnity

          ! Section is rotating and only a part of the physical problem
          ! is modelled. Consequently a rotation matrix is present for
          ! the velocity components.

          ! First transform to a frame where the xt-axis points in the
          ! direction of the rotation vector.

          xt(1) = sections(nn)%rotAxis(1)
          xt(2) = sections(nn)%rotAxis(2)
          xt(3) = sections(nn)%rotAxis(3)

          ! Construct the yt axis. It does not matter exactly as long
          ! as it is normal to xt.

          if(abs(xt(2)) < 0.707107_realType) then
             yt(1) = zero
             yt(2) = one
             yt(3) = zero
          else
             yt(1) = zero
             yt(2) = zero
             yt(3) = one
          endif

          ! Make sure that yt is normal to xt.

          tmp   = xt(1)*yt(1) + xt(2)*yt(2) + xt(3)*yt(3)
          yt(1) = yt(1) - tmp*xt(1)
          yt(2) = yt(2) - tmp*xt(2)
          yt(3) = yt(3) - tmp*xt(3)

          ! And create a unit vector.

          tmp   = one/sqrt(yt(1)**2 + yt(2)**2 + yt(3)**2)
          yt(1) = tmp*yt(1)
          yt(2) = tmp*yt(2)
          yt(3) = tmp*yt(3)

          ! Create the vector zt by taking the cross product xt*yt.

          zt(1) = xt(2)*yt(3) - xt(3)*yt(2)
          zt(2) = xt(3)*yt(1) - xt(1)*yt(3)
          zt(3) = xt(1)*yt(2) - xt(2)*yt(1)

          ! Compute the periodic angle theta and its sine and cosine.

          theta     = two*pi/real(sections(nn)%nSlices,realType)
          cosTheta = cos(theta)
          sinTheta = sin(theta)

          ! The rotation matrix in the xt,yt,zt frame is given by
          !
          ! R = | 1      0           0      |
          !     | 0  cos(theta) -sin(theta) |
          !     | 0  sin(theta)  cos(theta) |
          !
          ! The rotation matrix in the standard cartesian frame is then
          ! given by t * r * t^t, where the colums of the transformation
          ! matrix t are the unit vectors xt,yt,zt. One can easily check
          ! this by checking rotation around the y- and z-axis. The
          ! result of this is the expression below.

          rotMatrixSpectral(nn,1,1) = xt(1)*xt(1)                    &
               + cosTheta*(yt(1)*yt(1) + zt(1)*zt(1))
          rotMatrixSpectral(nn,1,2) = xt(1)*xt(2)                    &
               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
               - sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
          rotMatrixSpectral(nn,1,3) = xt(1)*xt(3)                    &
               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
               - sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))

          rotMatrixSpectral(nn,2,1) = xt(1)*xt(2)                    &
               + cosTheta*(yt(1)*yt(2) + zt(1)*zt(2)) &
               + sinTheta*(yt(1)*zt(2) - yt(2)*zt(1))
          rotMatrixSpectral(nn,2,2) = xt(2)*xt(2)                    &
               + cosTheta*(yt(2)*yt(2) + zt(2)*zt(2))
          rotMatrixSpectral(nn,2,3) = xt(2)*xt(3)                    &
               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
               - sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))

          rotMatrixSpectral(nn,3,1) = xt(1)*xt(3)                    &
               + cosTheta*(yt(1)*yt(3) + zt(1)*zt(3)) &
               + sinTheta*(yt(1)*zt(3) - yt(3)*zt(1))
          rotMatrixSpectral(nn,3,2) = xt(2)*xt(3)                    &
               + cosTheta*(yt(2)*yt(3) + zt(2)*zt(3)) &
               + sinTheta*(yt(2)*zt(3) - yt(3)*zt(2))
          rotMatrixSpectral(nn,3,3) = xt(3)*xt(3)                    &
               + cosTheta*(yt(3)*yt(3) + zt(3)*zt(3))

       endif testUnity

    enddo sectionLoop

  end subroutine timeRotMatricesSpectral


  subroutine fineGridSpectralCoor
    !
    !       fineGridSpectralCoor computes the coordinates of all but
    !       the first spectral solution from the known coordinates of the
    !       first time instance.
    !
    use constants
    use block, only: flowDoms, nDom
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use IOModule, only : IOVar
    use iteration, only : currentLevel
    use monitor, only : timeUnsteady
    use section, only : nSections, sections
    use partitionMod, only : interpolSpectral, nGridsRead
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, ll, i, j, k
    integer(kind=intType) :: il, jl, kl

    real(kind=realType), dimension(nSections) :: dt, t

    ! This routine is only used for the spectral solutions. Return
    ! immediately if a different mode is solved.

    if(equationMode /= timeSpectral) return

    ! Also return immediately if all coordinates were already read
    ! from the grid file.

    if(.not. interpolSpectral) return
    !
    !       Step 1. Perform a rigid body motion of the coordinates of the
    !               1st time instance to the other instances.
    !
    ! Set currentLevel to 1, such that updateCoorFineMesh
    ! updates the coordinates of the correct level.

    currentLevel = 1

    ! Determine the delta t for every section. Remember it is possible
    ! that every section has a different periodic time.

    do nn=1,nSections
       dt(nn) = sections(nn)%timePeriod &
            / real(nTimeIntervalsSpectral,realType)
    enddo

    ! Loop over the number of spectral solutions, starting at 2,
    ! because the first is already known.

    timeUnsteady = zero
    spectralLoop: do ll=2,nTimeIntervalsSpectral

       ! Set the owned coordinates to the coordinates of the first
       ! solution such that updateCoorFineMesh can modify them.

       do nn=1,nDom

          do k=1,flowDoms(nn,1,1)%kl
             do j=1,flowDoms(nn,1,1)%jl
                do i=1,flowDoms(nn,1,1)%il
                   flowDoms(nn,1,ll)%x(i,j,k,1) = flowDoms(nn,1,1)%x(i,j,k,1)
                   flowDoms(nn,1,ll)%x(i,j,k,2) = flowDoms(nn,1,1)%x(i,j,k,2)
                   flowDoms(nn,1,ll)%x(i,j,k,3) = flowDoms(nn,1,1)%x(i,j,k,3)
                enddo
             enddo
          enddo

       enddo

       ! Compute the corresponding times for this spectral solution
       ! and call updateCoorFineMesh to determine the coordinates.

       do nn=1,nSections
          t(nn) = (ll-1)*dt(nn)
       enddo

       call updateCoorFineMesh(t, ll)

    enddo spectralLoop

    ! Return if only one grid has been read.

    if(nGridsRead == 1) return
    !
    !       Step 2. Multiple grids have been read, but the number is not
    !               equal to the number of time instances used in the
    !               computation. As multiple grids have been read this
    !               means that a time spectral computation on a deforming
    !               mesh is performed. Therefore the deformations,
    !               relative to the first grid, must be interpolated.
    !
    ! First allocate the memory of IOVar(..,1)%w.
    ! This will serve as temporary storage for the coordinates of
    ! the 1st spectral instance, because those will be used in the
    ! call to updateCoorFineMesh. In this way this routine can be
    ! used without modification.

    do nn=1,nDom

       kl = flowDoms(nn,1,1)%kl
       jl = flowDoms(nn,1,1)%jl
       il = flowDoms(nn,1,1)%il

       allocate(IOVar(nn,1)%w(il,jl,kl,3), stat=ierr)
       if(ierr /= 0)                            &
            call terminate("fineGridSpectralCoor", &
            "Memory allocation failure for &
            &IOVar(nn,1)%w")

       do k=1,kl
          do j=1,jl
             do i=1,il
                IOVar(nn,1)%w(i,j,k,1) = flowDoms(nn,1,1)%x(i,j,k,1)
                IOVar(nn,1)%w(i,j,k,2) = flowDoms(nn,1,1)%x(i,j,k,2)
                IOVar(nn,1)%w(i,j,k,3) = flowDoms(nn,1,1)%x(i,j,k,3)
             enddo
          enddo
       enddo

    enddo

    ! Determine the delta t for every section. Remember it is possible
    ! that every section has a different periodic time.

    do nn=1,nSections
       dt(nn) = sections(nn)%timePeriod &
            / real(nGridsRead,realType)
    enddo

    ! Loop over the number of spectral solutions read, starting at
    ! 2, and determine the displacements relative to the rigid body
    ! motion of the grid of the 1st time instance.

    timeUnsteady = zero
    spectralLoopRead: do ll=2,nGridsRead

       ! Compute the corresponding times for this spectral solution
       ! and call updateCoorFineMesh to determine the coordinates.

       do nn=1,nSections
          t(nn) = (ll-1)*dt(nn)
       enddo

       call updateCoorFineMesh(t, 1_intType)

       ! Determine the relative displacements for this time instance
       ! and initialize flowDoms(nn,1,1) for the next round.

       do nn=1,nDom

          do k=1,flowDoms(nn,1,1)%kl
             do j=1,flowDoms(nn,1,1)%jl
                do i=1,flowDoms(nn,1,1)%il
                   IOVar(nn,ll)%w(i,j,k,1) = IOVar(nn,ll)%w(i,j,k,1)     &
                        - flowDoms(nn,1,1)%x(i,j,k,1)
                   IOVar(nn,ll)%w(i,j,k,2) = IOVar(nn,ll)%w(i,j,k,2)     &
                        - flowDoms(nn,1,1)%x(i,j,k,2)
                   IOVar(nn,ll)%w(i,j,k,3) = IOVar(nn,ll)%w(i,j,k,3)     &
                        - flowDoms(nn,1,1)%x(i,j,k,3)

                   flowDoms(nn,1,1)%x(i,j,k,1) = IOVar(nn,1)%w(i,j,k,1)
                   flowDoms(nn,1,1)%x(i,j,k,2) = IOVar(nn,1)%w(i,j,k,2)
                   flowDoms(nn,1,1)%x(i,j,k,3) = IOVar(nn,1)%w(i,j,k,3)
                enddo
             enddo
          enddo

       enddo

    enddo spectralLoopRead

    ! The coordinates of IOVar now contain the relative
    ! displacements compared to the rigid body motion of the first
    ! time instance, except for the first time instance.
    ! Set these to zero.

    do nn=1,nDom
       IOVar(nn,1)%w = zero
    enddo

    ! Interpolate the displacements and add them to the currently
    ! stored coordinates.

    call terminate("fineGridSpectralCoor", &
         "Arti should do this interpolation stuff")

    ! Release the memory of the variable w in IOVar.

    do nn=1,nDom
       do ll=1,nGridsRead
          deallocate(IOVar(nn,ll)%w, stat=ierr)
          if(ierr /= 0)                            &
               call terminate("fineGridSpectralCoor", &
               "Deallocation failure for IOVar%w")
       enddo
    enddo

  end subroutine fineGridSpectralCoor
 subroutine updateCoorFineMesh(dtAdvance, sps)
    !
    !       updateCoorFineMesh updates the coordinates of the
    !       moving parts of the current finest mesh by the given amount of
    !       time, possibly different per section. In unsteady mode all the
    !       times will be equal, but in time spectral mode they can be
    !       different.
    !       This routine is called in the full mg cycle to put the fine
    !       mesh to the position previously calculated on the coarser
    !       grid levels, in the unsteady time loop to advance the
    !       coordinates only one time step and in the partitioning part
    !       of the spectral mode to compute the coordinates of the given
    !       spectral solution sps. As it is used in the full MG cycle,
    !       currentLevel points to the correct grid level and not
    !       ground level.
    !
    use constants
    use block
    use blockPointers
    use flowVarRefState
    use cgnsGrid
    use inputMotion
    use iteration
    use monitor
    use utils, only : setPointers, rotMatrixRigidBody
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: sps

    real(kind=realType), dimension(*), intent(in) :: dtAdvance
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, nn

    real(kind=realType) :: displX, displY, displZ
    real(kind=realType) :: tNew, tOld
    real(kind=realType) :: angleX, angleY, angleZ, dx, dy, dz
    real(kind=realType) :: xiX, xiY, xiZ, etaX, etaY, etaZ
    real(kind=realType) :: zetaX, zetaY, zetaZ, xp, yp, zp, t
    real(kind=realType) :: phi, cosPhi, sinPhi, eta, zeta

    real(kind=realType), dimension(3)   :: rotationPoint
    real(kind=realType), dimension(3,3) :: rotationMatrix

    ! Compute the displacements due to the rigid motion of the mesh.

    displX = zero
    displY = zero
    displZ = zero

    ! Determine the time values of the old and new time level.
    ! It is assumed that the rigid body rotation of the mesh is only
    ! used when only 1 section is present.

    tNew = timeUnsteady + timeUnsteadyRestart
    tOld = tNew - dtAdvance(1)

    ! Compute the rotation matrix of the rigid body rotation as
    ! well as the rotation point; the latter may vary in time due
    ! to rigid body translation.

    call rotMatrixRigidBody(tNew, tOld, rotationMatrix, rotationPoint)

    ! Loop over the number of local blocks.

    blockLoop: do nn=1,nDom

       ! Set the pointers for this block on the current level.
       ! Note that currentLevel must be used and not groundLevel,
       ! because groundLevel is 1 level too coarse when this routine
       ! is called in the full mg cycle.

       call setPointers(nn, currentLevel, sps)
       !
       !         The rigid body motion of the entire mesh.
       !
       ! First the rotation.

       do k=1,kl
          do j=1,jl
             do i=1,il

                ! Determine the vector relative to the rotation point.

                xp = x(i,j,k,1) - rotationPoint(1)
                yp = x(i,j,k,2) - rotationPoint(2)
                zp = x(i,j,k,3) - rotationPoint(3)

                ! Apply the transformation matrix to the vector (xp,yp,zp)
                ! and set the new coordinates.

                x(i,j,k,1) = rotationMatrix(1,1)*xp &
                     + rotationMatrix(1,2)*yp &
                     + rotationMatrix(1,3)*zp + rotationPoint(1)
                x(i,j,k,2) = rotationMatrix(2,1)*xp &
                     + rotationMatrix(2,2)*yp &
                     + rotationMatrix(2,3)*zp + rotationPoint(2)
                x(i,j,k,3) = rotationMatrix(3,1)*xp &
                     + rotationMatrix(3,2)*yp &
                     + rotationMatrix(3,3)*zp + rotationPoint(3)
             enddo
          enddo
       enddo

       ! Add the translation.

       do k=1,kl
          do j=1,jl
             do i=1,il
                x(i,j,k,1) = x(i,j,k,1) + displX
                x(i,j,k,2) = x(i,j,k,2) + displY
                x(i,j,k,3) = x(i,j,k,3) + displZ
             enddo
          enddo
       enddo

       !
       !         Determine whether the corresponding cgns block is a rotating
       !         block. If it is, apply the rotation.
       !         Note that now the section ID of the block is taken into
       !         account to allow for different periodic times per section.
       !
       if( cgnsDoms(nbkGlobal)%rotatingFrameSpecified ) then

          ! Compute the rotation angles.

          angleX = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(1)
          angleY = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(2)
          angleZ = dtAdvance(sectionID)*cgnsDoms(nbkGlobal)%rotRate(3)

          ! Compute the unit vector in the direction of the rotation
          ! axis, which will be called the xi-direction.

          t    = one/max(eps,sqrt(angleX**2 + angleY**2 + angleZ**2))
          xiX = t*angleX
          xiY = t*angleY
          xiZ = t*angleZ

          ! Determine the rotation angle in xi-direction and its sine
          ! and cosine. Due to the definition of the xi-direction this
          ! angle will always be positive.

          phi    = xiX*angleX + xiY*angleY + xiZ*angleZ
          cosPhi = cos(phi)
          sinPhi = sin(phi)

          ! Loop over the owned coordinates of this block.

          do k=1,kl
             do j=1,jl
                do i=1,il

                   ! Compute the vector relative to center of rotation.

                   dx = x(i,j,k,1) - cgnsDoms(nbkGlobal)%rotCenter(1)
                   dy = x(i,j,k,2) - cgnsDoms(nbkGlobal)%rotCenter(2)
                   dz = x(i,j,k,3) - cgnsDoms(nbkGlobal)%rotCenter(3)

                   ! Compute the coordinates of the point p, which is the
                   ! closest point on the rotation axis.

                   t  = dx*xiX + dy*xiY + dz*xiZ
                   xp = cgnsDoms(nbkGlobal)%rotCenter(1) + t*xiX
                   yp = cgnsDoms(nbkGlobal)%rotCenter(2) + t*xiY
                   zp = cgnsDoms(nbkGlobal)%rotCenter(3) + t*xiZ

                   ! Determine the unit vector in eta direction, which
                   ! is defined from point p to the current point.

                   etaX = x(i,j,k,1) - xp
                   etaY = x(i,j,k,2) - yp
                   etaZ = x(i,j,k,3) - zp

                   eta = sqrt(etaX**2 + etaY**2 + etaZ**2)
                   t   = one/max(eps,eta)

                   etaX = t*etaX
                   etaY = t*etaY
                   etaZ = t*etaZ

                   ! Determine the unit vector in zeta-direction. This is
                   ! the cross product of the unit vectors in xi and in
                   ! eta-direction.

                   zetaX = xiY*etaZ - xiZ*etaY
                   zetaY = xiZ*etaX - xiX*etaZ
                   zetaZ = xiX*etaY - xiY*etaX

                   ! Compute the new eta and zeta coordinates.

                   zeta = eta*sinPhi
                   eta  = eta*cosPhi

                   ! Compute the new cartesian coordinates.

                   x(i,j,k,1) = xp + eta*etaX + zeta*zetaX
                   x(i,j,k,2) = yp + eta*etaY + zeta*zetaY
                   x(i,j,k,3) = zp + eta*etaZ + zeta*zetaZ

                enddo
             enddo
          enddo

       endif

    enddo blockLoop

  end subroutine updateCoorFineMesh
  subroutine allocCoorFineGrid
    !
    !       allocCoorFineGrid allocates the memory for all the coordinates
    !       of all local blocks. Also the memory for the derived data type
    !       used for the reading is allocated. If an interpolation must be
    !       performed for the time spectral method the variables of this
    !       IO type are allocated as well. For all other cases the pointer
    !       of the variables are set to the appropriate entry in flowDoms.
    !
    use constants
    use block, only : nDom, flowDoms
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use IOModule, only : IOVar
    use iteration, only : nOldLevels, deforming_grid
    use partitionMod, only : interpolSpectral, nGridsRead
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm
    integer(kind=intType) :: il, jl, kl, ie, je, ke

    ! Loop over the local blocks and allocate the memory for the
    ! coordinates.

    blockLoop: do nn=1,nDom

       ! Some abbreviations of the block dimensions.

       il = flowDoms(nn,1,1)%il
       jl = flowDoms(nn,1,1)%jl
       kl = flowDoms(nn,1,1)%kl

       ie = flowDoms(nn,1,1)%ie
       je = flowDoms(nn,1,1)%je
       ke = flowDoms(nn,1,1)%ke

       ! Loop over the number of spectral modes and allocate the
       ! memory of the coordinates, single halo's included. The halo
       ! values will be initialized in the preprocessing part.

       do mm=1,nTimeIntervalsSpectral
          allocate(flowDoms(nn,1,mm)%x(0:ie,0:je,0:ke,3), stat=ierr)
          if(ierr /= 0)                         &
               call terminate("allocCoorFineGrid", &
               "Memory allocation failure for flowDoms%x")

          flowDoms(nn,1,mm)%x = 0.0

          !allocate xInit for all time spectral intervals for meshwarping

          if(ierr /= 0)                         &
               call terminate("allocCoorFineGrid", &
               "Memory allocation failure for flowDoms%xInit")
          !flowDoms(nn,1,mm)%xInit=0.0
          !for the first grid also allocate xPlus and xMinus for the
          !meshwarping verification...
       enddo

       ! For a time accurate computation on deforming meshes, allocate
       ! the memory for the old coordinates. As this is not the time
       ! spectral mode, the third index of flowDoms is always 1.

       if(deforming_Grid .and. equationMode == unsteady) then

          allocate(flowDoms(nn,1,1)%xOld(nOldLevels,0:ie,0:je,0:ke,3), &
               stat=ierr)
          if(ierr /= 0)                         &
               call terminate("allocCoorFineGrid", &
               "Memory allocation failure for xOld")
       endif


       ! Added by HDN
       if(deforming_Grid .and. equationMode == unsteady) then

          allocate(flowDoms(nn,1,1)%xALE(0:ie,0:je,0:ke,3), &
               stat=ierr)
          if(ierr /= 0)                         &
               call terminate("allocCoorFineGrid", &
               "Memory allocation failure for xALE")
       endif

    enddo blockLoop

    ! Allocate the memory for IOVar.

    allocate(IOVar(nDom,nGridsRead), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("allocCoorFineGrid", &
         "Memory allocation failure for IOVar")

    ! Determine the equation mode we are solving and set the pointers
    ! of IOVar accordingly, or even allocate the memory, if needed.

    select case(equationMode)

    case (steady)

       ! Steady computation. Only one grid needs to be read.
       ! Loop over the number of blocks and set the pointer.
       ! No pointer offsets are needed for the coordinates.

       do nn=1,nDom
          IOVar(nn,1)%pointerOffset = 0
          IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)
       enddo

       !===============================================================

    case (unsteady)

       ! Unsteady computation. The first set of coordinates should
       ! be stored in x, other sets (if present) in xOld.
       ! No pointer offsets are needed for the coordinates.

       do nn=1,nDom
          IOVar(nn,1)%pointerOffset = 0
          IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)

          do mm=2,nGridsRead
             IOVar(nn,mm)%pointerOffset = 0
             IOVar(nn,mm)%w => flowDoms(nn,1,1)%xOld(mm-1,1:,1:,1:,:)
          enddo
       enddo

       !===============================================================

    case (timeSpectral)

       ! Time spectral mode. A further check is required.

       testAllocIOVar: if(interpolSpectral .and. &
            nGridsRead > 1) then

          ! A restart is performed for a deforming mesh using a
          ! different number of time instances than the previous
          ! computation. Consequently the coordinates, or better
          ! the deformations, will be interpolated later on. Hence
          ! some additional storage is required for the coordinates
          ! to be read and thus the memory for the variables w of
          ! IOVar is allocated. No halo data is needed here.
          ! Note that for the 1st time instance the pointer is set
          ! to the coordinates of flowDoms, because these are not
          ! interpolated. Only the higher time instances are
          ! interpolated. No pointer offsets are needed for the
          ! coordinates.

          do nn=1,nDom
             il = flowDoms(nn,1,1)%il
             jl = flowDoms(nn,1,1)%jl
             kl = flowDoms(nn,1,1)%kl

             IOVar(nn,1)%pointerOffset = 0
             IOVar(nn,1)%w => flowDoms(nn,1,1)%x(1:,1:,1:,:)

             do mm=2,nGridsRead
                IOVar(nn,mm)%pointerOffset = 0
                allocate(IOVar(nn,mm)%w(il,jl,kl,3), stat=ierr)
                if(ierr /= 0)                         &
                     call terminate("allocCoorFineGrid", &
                     "Memory allocation failure for &
                     &IOVar%w")
             enddo
          enddo

       else testAllocIOVar

          ! One of the following options is true.
          ! - The computation starts from scratch.
          ! - A restart is performed using a rigid grid, possibly
          !   moving. The number of time instances does not have
          !   to be the same.
          ! - A restart with a deforming mesh is performed with the
          !   same number of time instances.
          !
          ! In all these situations the pointers of IOVar are
          ! simply set to the coordinates of flowDoms.

          do nn=1,nDom
             do mm=1,nGridsRead
                IOVar(nn,mm)%pointerOffset = 0
                IOVar(nn,mm)%w => flowDoms(nn,1,mm)%x(1:,1:,1:,:)
             enddo
          enddo

       endif testAllocIOVar

    end select

  end subroutine allocCoorFineGrid
  subroutine checkPartitioning(np,load_inbalance,face_inbalance)

    ! This subroutine runs the load balancing and partitioning algorithm
    ! to determine what the load balancing will be for a given number of
    ! procs np. The output is load_inbalance and face_inbalance.

    use constants
    use communication, only : nProc
    use partitionMod, only : ubvec
    use loadBalance, only : blockDistribution
    implicit none

    integer(kind=intType),intent(in) ::np
    integer(kind=intType) :: nproc_save
    real(kind=realType),intent(out) :: load_inbalance,face_inbalance

    ! Note: This file follows mostly partitionAndReadGrid. See
    ! partitionAndReadGrid.f90 for more infromation

    ! Trick it into thinking we have np processors:
    nproc_save = nproc
    nproc = np

    call blockDistribution

    ! Restore the number of procs
    nproc = nproc_save

    ! Extract the inbalance info:
    load_inbalance = ubvec(1)
    face_inbalance = ubvec(2)

  end subroutine checkPartitioning

  subroutine determineSections
    !
    !       determineSections determines the number of sections, i.e.
    !       grid parts between sliding mesh interfaces, present in the
    !       entire grid.
    !
    use constants
    use block
    use cgnsGrid
    use communication
    use inputTimeSpectral
    use section
    use su_cgns
    use sorting, only : qsortIntegers, bsearchIntegers
    use utils, only : terminate
    implicit none
    !
    !      Local parameter, threshold for allowed angle difference between
    !                       the theoretical and true value, 0.1 degrees.
    !
    real(kind=realType), parameter :: threshold = 0.1_realType
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, ii, jj
    integer(kind=intType) :: nLevel, nSlices, slideID
    integer(kind=intType), dimension(cgnsNDom) :: sectionID, sorted
    integer(kind=intType), dimension(cgnsNsliding,2) :: secSliding

    real(kind=realType) :: cosTheta, cosPhi, cosPsi
    real(kind=realType) :: sinTheta, sinPhi, sinPsi
    real(kind=realType) :: r11, r12, r13, r21, r22, r23
    real(kind=realType) :: r31, r32, r33
    real(kind=realType) :: d1, d2, a1, a0, lamr, lami, angle, dAngle

    logical :: situationChanged

    ! Initialize sectionID to the cgns block number.

    do nn=1,cgnsNDom
       sectionID(nn) = nn
    enddo

    ! Loop to determine the highest cgns block id in every section.

    loopHighestBlock: do

       ! Initialize situationChanged to .false.

       situationChanged = .false.

       ! Loop over the internal block faces of the blocks.

       do nn=1,cgnsNDom

          ! First the 1 to 1 block connectivities.

          do mm=1,cgnsDoms(nn)%n1to1

             ! Store the neighboring block a bit easier and change the
             ! value of sectionID if the neighboring block has a
             ! higher value. In that case situationChanged is .true.

             ii = cgnsDoms(nn)%conn1to1(mm)%donorBlock
             if(sectionID(ii) > sectionID(nn)) then
                sectionID(nn)    = sectionID(ii)
                situationChanged = .true.
             endif

          enddo

          ! No general connectivities yet.

       enddo

       ! Criterion to exit the loop.

       if(.not. situationChanged) exit

    enddo loopHighestBlock

    ! Copy sectionID in sorted and sort sorted in increasing order.

    do nn=1,cgnsNDom
       sorted(nn) = sectionID(nn)
    enddo

    call qsortIntegers(sorted, cgnsNDom)

    ! Determine the number of sections. Note there is at least one
    ! section, because there is at least one block in the grid.

    nSections = 1
    do nn=2,cgnsNDom
       if(sorted(nn) > sorted(nSections)) then
          nSections = nSections +1
          sorted(nSections) = sorted(nn)
       endif
    enddo

    ! Determine the sections to which the owned blocks belong.

    nLevel = ubound(flowDoms,2)
    do nn=1,nDom

       ! Determine the corresponding cgns block id and search its
       ! section id in sorted.

       mm = flowDoms(nn,1,1)%cgnsBlockID
       mm = bsearchIntegers(sectionID(mm), sorted)

       if( debug ) then
          if(mm == 0) call terminate("determineSections", &
               "Entry not found in sorted.")
       endif

       ! Set the section id for all grid levels for all spectral
       ! time intervals.

       do jj=1,nTimeIntervalsSpectral
          do ii=1,nLevel
             flowDoms(nn,ii,jj)%sectionID = mm
          enddo
       enddo

    enddo

    ! Allocate the memory for sections.

    allocate(sections(nSections), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("determineSections", &
         "Memory allocation failure for sections.")

    ! Initialize the number of slices for each of the sections to 1,
    ! periodic and rotating to .false. and rotCenter, rotAxis
    ! and rotRate to zero.

    do nn=1,nSections
       sections(nn)%nSlices   = 1
       sections(nn)%periodic  = .false.
       sections(nn)%rotating  = .false.
       sections(nn)%rotCenter = zero
       sections(nn)%rotAxis   = zero
       sections(nn)%rotRate   = zero
    enddo

    ! Determine the number of slices and the periodic transformation
    ! for the sections.

    loopCGNSDom: do nn=1,cgnsNDom

       ! Search for the corresponding section.

       ii = bsearchIntegers(sectionID(nn), sorted(1:nSections))
       if( debug ) then
          if(ii == 0) call terminate("determineSections", &
               "Entry not found in sorted.")
       endif

       ! It is assumed that periodic info is correct. So if this
       ! section has already been treated, there is no need to do
       ! it again.

       if( sections(ii)%periodic ) cycle

       ! Loop over the 1 to 1 subfaces of the cgns block and try to
       ! find a periodic one.

       do mm=1,cgnsDoms(nn)%n1to1
          if( cgnsDoms(nn)%conn1to1(mm)%periodic ) exit
       enddo

       ! Continue with the next block if this block does not have
       ! periodic subfaces.

       if(mm > cgnsDoms(nn)%n1to1) cycle

       ! Subface mm is a periodic one. Set periodic to .true.

       sections(ii)%periodic = .true.

       ! Set the rotation axis of the section to the rotation
       ! angles of the periodic transformation. This may be
       ! overwritten later on using the rotation rate, but for
       ! some cases this is the only rotation information present.

       sections(ii)%rotAxis = cgnsDoms(nn)%conn1to1(mm)%rotationAngles

       ! Construct the rotation matrix, where it is assumed that the
       ! sequence of rotation is first rotation around the x-axis,
       ! followed by rotation around the y-axis and finally rotation
       ! around the z-axis.

       cosTheta = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(1))
       sinTheta = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(1))

       cosPhi = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(2))
       sinPhi = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(2))

       cosPsi = cos(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(3))
       sinPsi = sin(cgnsDoms(nn)%conn1to1(mm)%rotationAngles(3))

       r11 =  cosPhi*cosPsi
       r21 =  cosPhi*sinPsi
       r31 = -sinPhi

       r12 = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi
       r22 = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi
       r32 = sinTheta*cosPhi

       r13 = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi
       r23 = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi
       r33 = cosTheta*cosPhi

       ! Store the rotation matrix, rotation center and translation
       ! vector for this section.

       sections(ii)%rotCenter  = &
            cgnsDoms(nn)%conn1to1(mm)%rotationCenter
       sections(ii)%translation = &
            cgnsDoms(nn)%conn1to1(mm)%translation

       sections(ii)%rotMatrix(1,1) = r11
       sections(ii)%rotMatrix(2,1) = r21
       sections(ii)%rotMatrix(3,1) = r31

       sections(ii)%rotMatrix(1,2) = r12
       sections(ii)%rotMatrix(2,2) = r22
       sections(ii)%rotMatrix(3,2) = r32

       sections(ii)%rotMatrix(1,3) = r13
       sections(ii)%rotMatrix(2,3) = r23
       sections(ii)%rotMatrix(3,3) = r33

       ! Determine the coefficients of lambda and lambda^2 of the
       ! characteristic polynomial.

       d2 = -r11 - r22 - r33
       d1 =  r11*r22 + r11*r33 + r22*r33 - r12*r21 - r13*r31 - r23*r32

       ! Make use of the fact that one eigenvalue of the transformation
       ! matrix is 1 and determine the coefficients of the quadratic
       ! equation for the other two eigenvalues.

       a1 = d2 + one
       a0 = a1 + d1

       ! Determine the real and imaginary part of the two eigenvalues.
       ! Neglect the factor 1/2 here.

       lamr = -a1
       lami = sqrt(abs(a1*a1 - four*a0))

       ! Determine the angle in the imaginary plane. Due to the
       ! positive definition of lami, this angle will be between
       ! 0 and pi. Take care of the exceptional case that the angle
       ! is zero.

       angle = atan2(lami, lamr)
       if(angle == zero) angle = two*pi

       ! Determine the number of slices.

       nSlices = nint(two*pi/angle)

       ! Determine the angle difference in degrees between
       ! nSlices*angle and a complete rotation. If this is larger than
       ! the threshold processor 0 will print an error message and exit.

       dAngle = abs(180.0_realType*(two*pi - nSlices*angle)/pi)

       if(dAngle >= threshold) then
          if(myID == 0)                         &
               call terminate("determineSections", &
               "Periodic angle not a integer divide of &
               &360 degrees")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! Set the number of slices for this section.

       sections(ii)%nSlices = nSlices

    enddo loopCGNSDom

    ! Again loop over the number of block of the original mesh,
    ! but now determine whether or not the section is rotating.

    do nn=1,cgnsNDom

       ! If the block is rotating, copy that information to
       ! the corresponding section. If the section is not
       ! periodic, also set the rotation center.

       if( cgnsDoms(nn)%rotatingFrameSpecified ) then

          ii = bsearchIntegers(sectionID(nn), sorted(1:nSections))
          sections(ii)%rotating = .true.
          sections(ii)%rotAxis  = cgnsDoms(nn)%rotRate
          sections(ii)%rotRate  = cgnsDoms(nn)%rotRate

          if(.not. sections(ii)%periodic) &
               sections(ii)%rotCenter = cgnsDoms(nn)%rotCenter

       endif
    enddo

    ! Determine the two sections for every sliding mesh
    ! interface.

    secSliding = 0
    do nn=1,cgnsNDom
       do mm=1,cgnsDoms(nn)%nBocos
          if(cgnsDoms(nn)%bocoInfo(mm)%actualFace .and. &
               cgnsDoms(nn)%bocoInfo(mm)%BCType == SlidingInterface) then

             ! Boundary face is part of a sliding mesh interface.
             ! Determine the ID of the interface.

             slideID = abs(cgnsDoms(nn)%bocoInfo(mm)%slidingID)

             ! Determine the section to which this block belongs and
             ! store its id for this sliding interface.

             ii = bsearchIntegers(sectionID(nn), sorted(1:nSections))
             if(secSliding(slideID,1) == 0) then
                secSliding(slideID,1) = ii
             else if(secSliding(slideID,1) /= ii) then
                secSliding(slideID,2) = ii
             endif

          endif
       enddo
    enddo

    ! Loop over the sliding mesh interfaces to set the rotation axis
    ! for non-rotating sections.

    do ii=1,cgnsNsliding

       ! Store the two sections id's a bit easier.

       mm = secSliding(ii,1)
       nn = secSliding(ii,2)

       ! Print an error message if both sections are not rotating.

       if((.not. sections(mm)%rotating) .and. &
            (.not. sections(nn)%rotating) ) then
          if(myID == 0)                         &
               call terminate("determineSections", &
               "Encountered sliding interface between &
               &two non-rotating sections")
          call mpi_barrier(ADflow_comm_world, ierr)
       endif

       ! Set the rotation axis if section mm is not rotating.
       ! If it is not periodic also set the rotation point.

       if(.not. sections(mm)%rotating) then
          sections(mm)%rotAxis = sections(nn)%rotAxis
          if(.not. sections(mm)%periodic) &
               sections(mm)%rotCenter = sections(nn)%rotCenter
       endif

       ! Idem for section nn.

       if(.not. sections(nn)%rotating) then
          sections(nn)%rotAxis = sections(mm)%rotAxis
          if(.not. sections(nn)%periodic) &
               sections(nn)%rotCenter = sections(mm)%rotCenter
       endif

    enddo

    ! Determine the unit rotation axis for the sections.

    do nn=1,nSections
       d1 = one/max(eps,sqrt(sections(nn)%rotAxis(1)**2 &
            +                  sections(nn)%rotAxis(2)**2 &
            +                  sections(nn)%rotAxis(3)**2))

       sections(nn)%rotAxis(1) = d1*sections(nn)%rotAxis(1)
       sections(nn)%rotAxis(2) = d1*sections(nn)%rotAxis(2)
       sections(nn)%rotAxis(3) = d1*sections(nn)%rotAxis(3)
    enddo

  end subroutine determineSections

end module partitioning
