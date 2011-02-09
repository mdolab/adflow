!
!      ******************************************************************
!      *                                                                *
!      * File:          readBlockSizes.F90                              *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-17-2002                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readBlockSizes
!
!      ******************************************************************
!      *                                                                *
!      * readBlockSizes reads the number of blocks and their size       *
!      * from the given grid file. The data is stored in the module     *
!      * cgnsGrid.                                                      *
!      * If multiple grids need to be read for a consistent restart, it *
!      * is checked that the number of blocks and the block sizes are   *
!      * identical.                                                     *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("readBlockSizes", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use cgnsGrid
       use communication
       use constants
       use inputPhysics
       use iteration
       use su_cgns
       use partitionMod
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
!
!      Function definition.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Open the cgns files for reading and check if it went okay.
       ! Later on it is assumed that the 1st file stored in gridFiles
       ! is the "master".

       do nn=1,nGridsRead
         call cg_open_f(gridFiles(nn), mode_read, fileIDs(nn), ierr)
         if(ierr /= all_ok) then
           write(errorMessage,*) "File ", trim(gridFiles(nn)), &
                                " could not be opened for reading"
           call terminate("readBlockSizes", errorMessage)
         endif
       enddo

       cgnsInd = fileIDs(1)

       ! Determine the number of bases in the cgns file.
       ! This must be at least 1.

       call cg_nBases_f(cgnsInd, cgnsNbases, ierr)
       if(ierr /= all_ok)                 &
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
       if(ierr /= all_ok)                 &
         call terminate("readBlockSizes", &
                        "Something wrong when calling cg_base_read_f")

       ! Check the cell and physical dimensions. Both must be 3 for
       ! this code to work.

       if(cgnsCelldim /= 3 .or. cgnsPhysdim /= 3) then
         write(errorMessage,100) cgnsCelldim, cgnsPhysdim
 100     format("Both the number of cell and physical dimensions should &
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
         ii = bsearchStrings(cgnsFamilies(i)%familyName, &
                             sortedFamName, nn)

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
       ! deforming, and possibly overwrite the family specified data.

       changing_Grid = deforming_Grid
       call overwriteFamilyData(sortedFamName, famID)

       ! Determine the number of zones/blocks in the grid. Note that the
       ! reading is done using cgnsNzones (an integer variable) which is
       ! then copied into cgnsNDom (an integer(kind=intType) variable).
       ! cgnsNdom cannot be used directly, because a type mismatch may
       ! occur.

       call cg_nZones_f(cgnsInd, cgnsBase, cgnsNzones, ierr)
       if(ierr /= all_ok)                 &
         call terminate("readBlockSizes", &
                        "Something wrong when calling cg_nZones_f")
       cgnsNDom = cgnsNzones

       ! Check if the number of zones for all the grid to be read
       ! are identical.

       do nn=2,nGridsRead
         call cg_nZones_f(fileIDs(nn), cgnsBase, cgnsNzones, ierr)
         if(ierr /= all_ok)                 &
           call terminate("readBlockSizes", &
                          "Something wrong when calling cg_nZones_f")

         if(cgnsNzones /= cgnsNDom) then
           write(errorMessage,*) "File ", trim(gridFiles(nn)), &
                                 ": Different number of blocks than&
                                 & in file ", trim(gridFiles(1))
           if(myID == 0) call terminate("readBlockSizes", errorMessage)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif
       enddo

       ! Allocate the memory for cgnsDoms, the derived data type to
       ! store the info of the cgns grid, and nullify its pointers.

       allocate(cgnsDoms(cgnsNDom), stat=ierr)
       if(ierr /= 0)                      &
         call terminate("readBlockSizes", &
                        "Memory allocation failure for cgnsDoms")

       do nn=1,cgnsNDom
         call nullifyCGNSDomPointers(nn)
       enddo

       ! Some initializations before reading the zone info.

       nDoubleBoundFaces = 0
       noUnits           = .false.
       oversetPresent    = .false.

       ! Loop over the number of zones.

       zones: do nZone=1, cgnsNDom

         ! Read the zone info

         call readZoneInfo(cgnsBase, nZone, sortedFamName, &
                           famID, noUnits)

         ! Count the number of each connectivity for the supported
         ! types (needed for memory allocation and to add up 1-to-1)

         call countConnectivities(cgnsInd, cgnsBase, nZone,          &
                                  cgnsDoms(nZone)%n1to1,             &
                                  cgnsDoms(nZone)%n1to1General,      &
                                  cgnsDoms(nZone)%nNonMatchAbutting, &
                                  cgnsDoms(nZone)%nOverset)

         ! If there are any overset connectivities, reset
         ! oversetPresent.

         if (cgnsDoms(nZone)%nOverset > 0) oversetPresent = .true.

         ! For this zone, read the 1-to-1 block connectivity, the
         ! general connectivities, the overset holes, and the
         ! boundary conditions.

         call read1to1Conn(cgnsInd, cgnsBase, nZone)
         call readGeneralConn(cgnsInd, cgnsBase, nZone)
         call readOversetHoles(cgnsInd, cgnsBase, nZone)
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
 110     format("# ",a, " double boundary faces found.")

       endif

       ! If the units could not be determined, print a warning if
       ! the viscous or unsteady equations are solved.

       if(myID == 0 .and. noUnits) then

         if(equations == NSEquations .or. equations == RANSEquations &
            .or. equationMode == unsteady) then

           print "(a)", "#"
           print "(a)", "#                      Warning"
           print "(a)", "# Conversion factor from grid units to &
                        &meter not specified and some blocks"
           print "(a)", "# do not have units. It is assumed that &
                        &the grid is given in meters."
           print "(a)", "#"

         endif
       endif

#endif

       end subroutine readBlockSizes
