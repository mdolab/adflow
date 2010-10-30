!
!      ******************************************************************
!      *                                                                *
!      * File:          readConnectivityPlot3D.f90                      *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 02-18-2005                                      *
!      * Last modified: 10-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readConnectivityPlot3D
!
!      ******************************************************************
!      *                                                                *
!      * readConnectivityPlot3D reads the ascii connectivity file,      *
!      * which corresponds to plot3D grid file. The connectivity file   *
!      * contains all the subface information for the block boundaries, *
!      * both for internal and physical boundary subfaces.              *
!      * This file is read by all processors and the information stored *
!      * in this file is also stored on all processors                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use su_cgns
       use flowVarRefState
       use inputIO
       use iteration
       implicit none
!
!      Local parameter.
!
       integer, parameter :: readunit = 35
!
!      Local variables.
!
       integer :: ios, ierr, pos

       integer(kind=intType) :: nn, ii, jj, nSub, nDomConn
       integer(kind=intType) :: nAlloc, nSubfaces, nBocos
       integer(kind=intType) :: n1to1, nNonMatch

       integer(kind=intType), dimension(16)             :: intArray
       integer(kind=intType), dimension(:), allocatable :: familyID

       character(len=maxStringLen) :: value
       character(len=512)          :: string
       character(len=7)            :: int1String, int2String, int3String

       character(len=512), dimension(:), allocatable :: strings, &
                                                        tmpStrings

       character(len=maxCGNSNameLen), dimension(:), allocatable :: names

       type(cgns1to1ConnType), pointer, dimension(:) :: conn1to1
       type(cgnsBocoType),     pointer, dimension(:) :: bocoInfo

       type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
                                                            connNonMatch
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Some initializations.

       cgnsNFamilies = 0
       cgnsCellDim   = 3
       cgnsPhysDim   = 3
       cgnsBaseName  = "Plot3DBase"

       ! Open the plot3D connectivity file for reading. If something
       ! goes wrong processor 0 will print an error message and the
       ! program terminates.

       open(unit=readunit, file=plot3DConnFile, status="old", &
            action="read", iostat=ios)

       if(ios /= 0) then
         write(string,101) trim(plot3DConnFile)
         if(myID == 0) &
           call terminate("readConnectivityPlot3D", string)

         call mpi_barrier(SUmb_comm_world, ierr)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Find out the number of families.                               *
!      *                                                                *
!      ******************************************************************
!
       ! Find the first line, which is not a comment and read the
       ! number of families.

       call findNextLineWithInfo(string)
       read(string,*) cgnsNFamilies
!
!      ******************************************************************
!      *                                                                *
!      * Read the family info.                                          *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the families and some help arrays.

       nn = cgnsNFamilies
       allocate(cgnsFamilies(nn), familyID(nn), names(nn), stat=ierr)
       if(ierr /= 0)                              &
         call terminate("readConnectivityPlot3D", &
                        "Memory allocation failure for cgnsFamilies, &
                        &etc.")

       ! Loop over the number families and read the name and possibly
       ! the boundary condition.

       familyLoop: do nn=1,cgnsNFamilies

         ! Some initializations for this family.

         cgnsFamilies(nn)%slidingID       = 0
         cgnsFamilies(nn)%bleedRegionID   = 0
         cgnsFamilies(nn)%monitorMassFlow = .false.

         nullify(cgnsFamilies(nn)%dataSet)

         ! Find the next line which contains information.

         call findNextLineWithInfo(string)

         ! Determine the first occurance of a space and copy the
         ! family name.

         pos = index(string, " ")
         if(pos == 0) then
           cgnsFamilies(nn)%familyName = string
         else
           cgnsFamilies(nn)%familyName = string(:pos-1)
         endif

         ! Remove the family name from the string and remove the
         ! possible leading spaces.

         string = string(pos+1:)
         string = adjustl(string)
         string = trim(string)

         ! If string does not contain any information anymore, no
         ! boundary condition is specified; otherwise there is.

         bcPresentTest: if(len_trim(string) == 0) then

           cgnsFamilies(nn)%BCTypeCGNS = Null
           cgnsFamilies(nn)%BCType     = BCNull
           cgnsFamilies(nn)%bcName     = ""

         else bcPresentTest

           ! String contains more information. Determine the string
           ! which defines the boundary condition.

           pos = index(string, " ")
           if(pos == 0) then
             cgnsFamilies(nn)%bcName = string
           else
             cgnsFamilies(nn)%bcName = string(:pos-1)
           endif

           ! Determine the corresponding internal and CGNS
           ! boundary condition.

           call cgnsBCForPlot3D(cgnsFamilies(nn)%bcName,     &
                                cgnsFamilies(nn)%BCTypeCGNS, &
                                cgnsFamilies(nn)%BCType)

           ! If this is an unknown boundary condition, terminate.

           if(cgnsFamilies(nn)%BCType == BCNull) then

             write(string,102) trim(cgnsFamilies(nn)%familyName)
             if(myID == 0) &
               call terminate("readConnectivityPlot3D", string)

             call mpi_barrier(SUmb_comm_world, ierr)
           endif

         endif bcPresentTest

       enddo familyLoop

       ! Copy the family names in names and sort them. This is needed
       ! later on in the search to determine the family ID.

       do nn=1,cgnsNFamilies
         names(nn) = cgnsFamilies(nn)%familyName
       enddo

       call qsortStrings(names, cgnsNFamilies)

       ! Determine the family ID's of the sorted names. These are
       ! initialized to -1. This serves as a check to see if all family
       ! names are unique.

       do nn=1,cgnsNFamilies
         familyID(nn) = -1
       enddo

       do nn=1,cgnsNFamilies

         ii = bsearchStrings(cgnsFamilies(nn)%familyName, names, &
                             cgnsNFamilies)

         if(familyID(ii) /= -1) then
           write(string,103) trim(cgnsFamilies(nn)%familyName)
           if(myID == 0) &
             call terminate("readConnectivityPlot3D", string)
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         familyID(ii) = nn

       enddo

       ! Initialize whether the grid is changing to whether it is
       ! deforming, and possibly overwrite the family specified data.

       changing_Grid = deforming_Grid
       call overwriteFamilyData(names, familyID)
!
!      ******************************************************************
!      *                                                                *
!      * Find out and check the number of blocks.                       *
!      *                                                                *
!      ******************************************************************
!
       ! Find the first line, which is not a comment and read the
       ! number of blocks.

       call findNextLineWithInfo(string)
       read(string,*) nDomConn

       ! Check if the number of blocks in the plot3D file and this
       ! connectivity file are equal. If not print an error message
       ! and terminate.

       if(nDomConn /= cgnsNDom) then
         if(myID == 0) &
           call terminate("readConnectivityPlot3D", &
                          "Number of blocks in plot3D grid differs from &
                          &the number in the connectivity file.")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Read the connectivity and boundary condition info for the      *
!      * subfaces.                                                      *
!      *                                                                *
!      ******************************************************************
!
       ! Pre-allocate strings.

       nAlloc = 20
       allocate(strings(nAlloc), stat=ierr)
       if(ierr /= 0) &
         call terminate("readConnectivityPlot3D", &
                        "Memory allocation failure for strings.")

       ! Loop over the number of blocks.

       domainLoop: do nn=1,cgnsNDom

         ! Write the block number to a string and perform some
         ! initializations.

         write(int1String,"(i7)") nn
         int1String = adjustl(int1String)

         cgnsDoms(nn)%zonename = trim(int1String)
         cgnsDoms(nn)%zonetype = Structured

         cgnsDoms(nn)%LRef     = LRef
         cgnsDoms(nn)%familyID = 0

         cgnsDoms(nn)%gridUnitsSpecified     = .false.
         cgnsDoms(nn)%rotatingFrameSpecified = .false.

         cgnsDoms(nn)%rotCenter = zero
         cgnsDoms(nn)%rotRate   = zero

         ! Find the first line, which is not a comment. This line
         ! contains the information about the current block id and
         ! possibly the family it belongs to.

         call findNextLineWithInfo(string)

         ! Create the string value, which should contain the block ID.

         pos = index(string, " ")
         if(pos == 0) then
           value = string
         else
           value = string(:pos-1)
         endif

         ! Read and check the block ID.

         read(value,*) ii
         if(ii /= nn) then
           if(myID == 0) &
             call terminate("readConnectivityPlot3D", &
                            "Block numbers are not numbered &
                            &contiguously in the plot3D &
                            &connectivity file.")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Remove the information just read from the string and remove
         ! the possible leading spaces.

         string = string(pos+1:)
         string = adjustl(string)
         string = trim(string)

         ! Check if a family name was specified.

         if(len_trim(string) > 0) then

           ! Store the family name in value.

           pos = index(string, " ")
           if(pos == 0) then
             value = string
           else
             value = string(:pos-1)
           endif

           ! Search the corresponding entry in the family names.
           ! If it is not found this means that the family name was not
           ! specified; write an error message and terminate.

           ii = bsearchStrings(value, names, cgnsNFamilies)
           if(ii == 0) then
             write(string,106) trim(int1String), trim(value)
             if(myID == 0) &
               call terminate("readConnectivityPlot3D", string)
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           cgnsDoms(nn)%familyID = familyID(ii)

           ! Set the rotation rate of the block if it was specified
           ! for the family. Note the conversion to meters for the
           ! rotation center, because it is assumed that the rotation
           ! center is specified in the original coordinates.

           ii = cgnsDoms(nn)%familyID

           if( cgnsFamilies(ii)%rotatingFrameSpecified ) then
             cgnsDoms(nn)%rotatingFrameSpecified = .true.
             cgnsDoms(nn)%rotRate   = cgnsFamilies(ii)%rotRate
             cgnsDoms(nn)%rotCenter = cgnsFamilies(ii)%rotCenter*LRef
           endif

         endif

         ! Loop over the 6 faces of the block to determine the number
         ! of subfaces of the block.

         nSubfaces = 0
         nBocos    = 0
         n1to1     = 0
         nNonMatch = 0

         faceLoop: do jj=1,6

           ! Find the first line, which is not a comment. This line
           ! contains the information about the current face ID and
           ! the number of subfaces for this face ID. The former can
           ! be ignored.

           call findNextLineWithInfo(string)
           read(string,*) ii, nSub

           ! Check if enough memory has been allocated to store the
           ! subface info in strings. If not, reallocate.

           testReallocate: if((nSub+nSubfaces) > nAlloc) then

             allocate(tmpStrings(nSubfaces), stat=ierr)
             if(ierr /= 0)                              &
               call terminate("readConnectivityPlot3D", &
                              "Memory allocation failure for tmpStrings")

             do ii=1,nSubfaces
               tmpStrings(ii) = strings(ii)
             enddo

             deallocate(strings, stat=ierr)
             if(ierr /= 0)                              &
               call terminate("readConnectivityPlot3D", &
                              "Deallocation failure for strings")

             nAlloc = nAlloc + nSub
             allocate(strings(nAlloc), stat=ierr)
             if(ierr /= 0)                              &
               call terminate("readConnectivityPlot3D", &
                              "Memory allocation failure for strings")

             do ii=1,nSubfaces
               strings(ii) = tmpStrings(ii)
             enddo

             deallocate(tmpStrings, stat=ierr)
             if(ierr /= 0)                              &
               call terminate("readConnectivityPlot3D", &
                              "Deallocation failure for tmpStrings")

           endif testReallocate

           ! Loop over the subfaces of the current face. Store the
           ! string with the information and update the number of
           ! boundary and connectivity subfaces.

           subfaceLoop: do ii=1,nSub

             ! Find the next string which contains information.

             nSubfaces = nSubfaces + 1
             call findNextLineWithInfo(strings(nSubfaces))

             ! Determine the subface type. The first part of the string
             ! is copied into string and converted to lower case, such
             ! that the comparison is easier. It cannot be done directly
             ! in strings, because this may change the family names.

             pos = index(strings(nSubfaces), " ")
             string = strings(nSubfaces)(:pos-1)
             call convertToLowerCase(string)

             select case(string)

               case ("internal1to1", "periodic1to1")
                 n1to1 = n1to1 + 1

               !=========================================================

               case ("internalnonmatch", "periodicnonmatch")
                 nNonMatch = nNonMatch + 1

               !=========================================================

               case default
                 nBocos = nBocos + 1

             end select

           enddo subfaceLoop
         enddo faceLoop

         ! Allocate the memory for the block connectivities and
         ! boundary condition info.

         cgnsDoms(nn)%n1to1General      = 0
         cgnsDoms(nn)%n1to1             = n1to1
         cgnsDoms(nn)%nNonMatchAbutting = nNonMatch
         cgnsDoms(nn)%nBocos            = nBocos

         allocate(cgnsDoms(nn)%conn1to1(n1to1),                 &
                  cgnsDoms(nn)%connNonMatchAbutting(nNonMatch), &
                  cgnsDoms(nn)%bocoInfo(nBocos), stat=ierr)
         if(ierr /= 0)                              &
           call terminate("readConnectivityPlot3D", &
                          "Memory allocation failure for conn1to1, &
                          &connNonMatchAbutting and bocoInfo.")

         ! Set some pointers to make the code more readable.

         conn1to1     => cgnsDoms(nn)%conn1to1
         connNonMatch => cgnsDoms(nn)%connNonMatchAbutting
         bocoInfo     => cgnsDoms(nn)%bocoInfo

         ! Loop over all subfaces and extract the information from
         ! the strings stored.

         nBocos    = 0
         n1to1     = 0
         nNonMatch = 0

         allSubfaceLoop: do ii=1,nSubfaces

           ! Store the subface type in value and remove this part
           ! from the currently active string.

           pos = index(strings(ii), " ")
           value = strings(ii)(:pos-1)

           strings(ii) = strings(ii)(pos+1:)
           strings(ii) = adjustl(strings(ii))
           strings(ii) = trim(strings(ii))

           ! Determine the subface type we have here. To make the
           ! comparison a bit easier value is copied into string
           ! and converted to lower case. The conversion cannot be
           ! done directly on value, because this may change the
           ! family name.

           string = value
           call convertToLowerCase(string)

           subfaceTypeTest: select case(string)

             case ("internal1to1", "periodic1to1")

               ! Internal 1 to 1 subface. Update the counter and
               ! create a name for this connectivity.

               n1to1 = n1to1 + 1
               write(int2String,"(i7)") n1to1
               int2String = adjustl(int2String)

               conn1to1(n1to1)%connectName = "1to1Conn"//trim(int2String)

               ! Read the 16 integers which describe the connectivity
               ! from the string and store it in the appropriate places.

               call readIntegersFromString(strings(ii), intArray, &
                                           16_intType)

               conn1to1(n1to1)%iBeg = intArray(1)
               conn1to1(n1to1)%iEnd = intArray(2)
               conn1to1(n1to1)%jBeg = intArray(3)
               conn1to1(n1to1)%jEnd = intArray(4)
               conn1to1(n1to1)%kBeg = intArray(5)
               conn1to1(n1to1)%kEnd = intArray(6)

               conn1to1(n1to1)%donorBlock = intArray(7)

               conn1to1(n1to1)%diBeg = intArray(8)
               conn1to1(n1to1)%diEnd = intArray(9)
               conn1to1(n1to1)%djBeg = intArray(10)
               conn1to1(n1to1)%djEnd = intArray(11)
               conn1to1(n1to1)%dkBeg = intArray(12)
               conn1to1(n1to1)%dkEnd = intArray(13)

               conn1to1(n1to1)%l1 = intArray(14)
               conn1to1(n1to1)%l2 = intArray(15)
               conn1to1(n1to1)%l3 = intArray(16)

               ! Set the donor name.

               write(int2String,"(i7)") conn1to1(n1to1)%donorBlock
               int2String = adjustl(int2String)

               conn1to1(n1to1)%donorName = trim(int2String)

               ! For a periodic subface some more info is needed.

               if(string == "periodic1to1") then

                 ! Read the periodic info and set the flag periodic
                 ! to .true.

                 call readPeriodicInfo(strings(ii),                    &
                                       conn1to1(n1to1)%rotationCenter, &
                                       conn1to1(n1to1)%rotationAngles, &
                                       conn1to1(n1to1)%translation)

                 ! Set the flag that periodic info was specified
                 ! to .true.

                 conn1to1(n1to1)%periodic = .true.

               else

                 ! No periodic info specified.

                 conn1to1(n1to1)%periodic = .false.

               endif

             !===========================================================

             case ("internalnonmatch", "periodicnonmatch")

               ! Internal non matching block to block connectivity.
               ! Update the counter.

               nNonMatch = nNonMatch + 1

               ! Read the 7 integers which describe the subface and the
               ! number of donor blocks, Store the data in the
               ! appropriate variables.

               call readIntegersFromString(strings(ii), intArray, &
                                           7_intType)

               connNonMatch(nNonMatch)%iBeg = intArray(1)
               connNonMatch(nNonMatch)%iEnd = intArray(2)
               connNonMatch(nNonMatch)%jBeg = intArray(3)
               connNonMatch(nNonMatch)%jEnd = intArray(4)
               connNonMatch(nNonMatch)%kBeg = intArray(5)
               connNonMatch(nNonMatch)%kEnd = intArray(6)

               connNonMatch(nNonMatch)%nDonorBlocks = intArray(7)

               ! Allocate the memory for the connectivity name, the
               ! donor block names, their ID's and the block face
               ! ID which abuts this subface.

               jj = connNonMatch(nNonMatch)%nDonorBlocks
               allocate(connNonMatch(nNonMatch)%connectNames(jj), &
                        connNonMatch(nNonMatch)%donorNames(jj),   &
                        connNonMatch(nNonMatch)%donorBlocks(jj),  &
                        connNonMatch(nNonMatch)%donorFaceIDs(jj), &
                        stat=ierr)
               if(ierr /= 0)                              &
                 call terminate("readConnectivityPlot3D", &
                                "Memory allocation failure for &
                                &connectNames, donorNames, &
                                &donorBlocks and donorFaceIDs")

               ! Store the number of the current connectivity in
               ! int2String to avoid some duplication in the loop.

               write(int2String,"(i7)") nNonMatch
               int2String = adjustl(int2String)

               ! Loop over the donor blocks and read the block ID and
               ! block face ID. Set the name of the donor block as well.

               donorsLoop: do jj=1,connNonMatch(nNonMatch)%nDonorBlocks

                 ! Create a name for this connectivity.

                 write(int3String,"(i7)") jj
                 int3String = adjustl(int3String)

                 connNonMatch(nNonMatch)%connectNames(jj) = &
                  "nonMatchConn"//trim(int2String)//"_"//trim(int3String)

                 ! Read the two integers from the string and set the
                 ! donor block ID and block face ID.

                 call readIntegersFromString(strings(ii), intArray, &
                                             2_intType)

                 connNonMatch(nNonMatch)%donorBlocks(jj) = intArray(1)

                 select case (intArray(2))
                   case (1_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = iMin
                   case (2_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = iMax
                   case (3_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = jMin
                   case (4_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = jMax
                   case (5_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = kMin
                   case (6_intType)
                     connNonMatch(nNonMatch)%donorFaceIDs(jj) = kMax
                   case default
                     write(string,107) trim(int1String), trim(int2String)
                     if(myID == 0) &
                       call terminate("readConnectivityPlot3D", string)

                     call mpi_barrier(SUmb_comm_world, ierr)
                 end select

                 ! Set the donor name.

                 write(int3String,"(i7)") &
                                 connNonMatch(nNonMatch)%donorBlocks(jj)
                 int3String = adjustl(int3String)

                 connNonMatch(nNonMatch)%donorNames(jj) = trim(int3String)

               enddo donorsLoop

               ! For a periodic subface some more info is needed.

               if(string == "periodicnonmatch") then

                 ! Read the periodic info and set the flag periodic
                 ! to .true.

                 call readPeriodicInfo(strings(ii),                    &
                               connNonMatch(nNonMatch)%rotationCenter, &
                               connNonMatch(nNonMatch)%rotationAngles, &
                               connNonMatch(nNonMatch)%translation)

                 ! Set the flag that periodic info was specified
                 ! to .true.

                 connNonMatch(nNonMatch)%periodic = .true.

               else

                 ! No periodic info specified.

                 connNonMatch(nNonMatch)%periodic = .false.

               endif

             !===========================================================

             case default

               ! Boundary subface. Update the counter and create a name
               ! for this boundary.

               nBocos = nBocos + 1

               write(int2String,"(i7)") nBocos
               int2String = adjustl(int2String)

               bocoInfo(nBocos)%bocoName = "Boco"//trim(int2String)

               ! Some initializations, mainly for compatibility
               ! with cgns.

               bocoInfo(nBocos)%ptsetType      = PointRange
               bocoInfo(nBocos)%npnts          = 2
               bocoInfo(nBocos)%normalIndex    = 1  ! just a value
               bocoInfo(nBocos)%normalListFlag = 0
               bocoInfo(nBocos)%normalDataType = RealDouble
               bocoInfo(nBocos)%actualFace     = .true.

               bocoInfo(nBocos)%familyID  = 0
               bocoInfo(nBocos)%slidingID = 0

               bocoInfo(nBocos)%nDataSet         = 0
               bocoInfo(nBocos)%dataSetAllocated = .false.
               nullify(bocoInfo(nBocos)%dataSet)

               ! Initialize the boundary velocity to the velocity
               ! of the corresponding block.

               bocoInfo(nBocos)%rotCenter = cgnsDoms(nn)%rotCenter
               bocoInfo(nBocos)%rotRate   = cgnsDoms(nn)%rotRate

               ! Read the 6 integers which describe to subface range.

               call readIntegersFromString(strings(ii), intArray, &
                                           6_intType)

               bocoInfo(nBocos)%iBeg = intArray(1)
               bocoInfo(nBocos)%iEnd = intArray(2)
               bocoInfo(nBocos)%jBeg = intArray(3)
               bocoInfo(nBocos)%jEnd = intArray(4)
               bocoInfo(nBocos)%kBeg = intArray(5)
               bocoInfo(nBocos)%kEnd = intArray(6)

               ! Determine the boundary condition and/or family
               ! for this subface.

               jj = bsearchStrings(value, names, cgnsNFamilies)

               testFamily: if(jj == 0) then

                 ! Corresponding family not found. Therefore it is
                 ! assumed that a standard boundary condition is
                 ! specified. Convert it to a useable format.

                 call cgnsBCForPlot3D(value,                       &
                                      bocoInfo(nBocos)%BCTypeCGNS, &
                                      bocoInfo(nBocos)%BCType)

                 ! Print an error message if the boundary condition was
                 ! not recognized or if it is a boundary condition that
                 ! should be specified per family when using Plot3D.

                 select case (bocoInfo(nBocos)%BCType)

                   case (bcNull)

                     write(string,108) trim(int1String), trim(int2String)
                     if(myID == 0) &
                       call terminate("readConnectivityPlot3D", string)
                     call mpi_barrier(SUmb_comm_world, ierr)

                   !=====================================================

                   case (MassBleedInflow,       MassBleedOutflow,     &
                         SlidingInterface,      DomainInterfaceAll,   &
                         DomainInterfaceRhoUVW, DomainInterfaceP,     &
                         DomainInterfaceRho,    DomainInterfaceTotal)

                     write(string,109) trim(int1String), trim(int2String)
                     if(myID == 0) &
                       call terminate("readConnectivityPlot3D", string)
                     call mpi_barrier(SUmb_comm_world, ierr)

                 end select

               else testFamily

                 ! Family name is found. Determine the corresponding
                 ! ID and set the boundary conditions to the one of the
                 ! family.

                 jj = familyID(jj)
                 bocoInfo(nBocos)%familyID   = jj
                 bocoInfo(nBocos)%BCTypeCGNS = cgnsFamilies(jj)%BCTypeCGNS
                 bocoInfo(nBocos)%BCType     = cgnsFamilies(jj)%BCType

                 ! Overwrite the rotation center and rate if this was
                 ! specified for the family. Note the multiplication
                 ! with LRef to obtain the coordinates in meter.

                 if( cgnsFamilies(jj)%rotatingFrameSpecified ) then
                   bocoInfo(nBocos)%rotCenter = &
                              LRef*cgnsFamilies(jj)%rotCenter
                   bocoInfo(nBocos)%rotRate   = &
                                   cgnsFamilies(jj)%rotRate
                 endif

                 ! Set the pointer for the data set, if specified.

                 if(cgnsFamilies(jj)%nDataSet > 0) then
                   bocoInfo(nBocos)%nDataSet = cgnsFamilies(jj)%nDataSet
                   bocoInfo(nBocos)%dataSet => cgnsFamilies(jj)%dataSet
                 endif

               endif testFamily

           end select subfaceTypeTest

         enddo allSubfaceLoop

       enddo domainLoop

       ! Release the memory of names and familyID.

       deallocate(names, familyID, stat=ierr)
       if(ierr /= 0)                              &
         call terminate("readConnectivityPlot3D", &
                        "Deallocation failure for names and familyID")

       ! Release the memory of strings.

       deallocate(strings, stat=ierr)
       if(ierr /= 0) &
         call terminate("readConnectivityPlot3D", &
                        "Deallocation failure for strings.")

       ! Close the file.

       close(unit=readunit)

       ! Format statements.

 101   format("Plot3D connectivity file",1X,A,1X,"not found.")
 102   format("Family",1X,A, ": Unknown boundary condition")
 103   format("Family name",1X,A,1X,"occurs more than once.")
 106   format("Block",1X,A,": Unknown family name",1X,A)
 107   format("Block",1X,A,", non-matching abutting subface",1X,A, &
              "Illegal face ID of the donor block.")
 108   format("Block",1X,A,", boundary subface",1X,A, &
              ": Unknown boundary condition type")
 109   format("Block",1X,A,", boundary subface",1X,A, &
              ": Boundary condition is only allowed per family")

       !=================================================================

       contains

         !===============================================================

         subroutine findNextLineWithInfo(infoString)
!
!        ****************************************************************
!        *                                                              *
!        * findNextLineWithInfo finds the next first line in the        *
!        * the readunit, which contains information. As always a        *
!        * comment line is indicated by a "#".                          *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Subroutine arguments.
!
         character(len=*), intent(out) :: infoString
!
!        Local variables.
!
         integer :: ios, ierr, pos
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Loop to find the next line which contains information.

         nextLineLoop: do

           ! Read the next line. If something goes wrong print an
           ! error message and terminate.

           read(unit=readunit, fmt="(a512)", iostat=ios) infoString
           if(ios /= 0) then
             if(myID == 0) &
               call terminate("findNextLineWithInfo", &
                              "Unexpected end of file encountered")
             call mpi_barrier(SUmb_comm_world, ierr)
           endif

           ! Replace all the tab and return characters by spaces and
           ! get rid of the leading and trailing spaces in infoString.

           call replaceTabsAndReturns(infoString)
           infoString = adjustl(infoString)
           infoString = trim(infoString)

           ! In case this is an empty string, continue with the next
           ! line. Idem if this is a comment line

           if(len_trim(infoString) == 0) cycle
           if(infoString(:1) == "#")     cycle

           ! Find a possible comment sign somewhere in the string.
           ! If present the info following the comment sign is ignored.

           pos = index(infoString, "#")
           if(pos > 0) then
             infoString = infoString(:pos-1)
             infoString = trim(infoString)
           endif

           ! Line is found. Exit the loop.

           exit nextLineLoop

         enddo nextLineLoop

         end subroutine findNextLineWithInfo

       end subroutine readConnectivityPlot3D

       !================================================================

       subroutine readPeriodicInfo(string,         rotationCenter, &
                                   rotationAngles, translation)
!
!      ******************************************************************
!      *                                                                *
!      * readPeriodicInfo reads the rotation center, the rotation       *
!      * angles and the translation for the periodic transformation     *
!      * from the given string.                                         *
!      *                                                                *
!      ******************************************************************
!
       use constants
       use flowVarRefState
       implicit none
!
!      Subroutine arguments.
!
       real(kind=realType), dimension(3), intent(out) :: rotationCenter
       real(kind=realType), dimension(3), intent(out) :: rotationAngles
       real(kind=realType), dimension(3), intent(out) :: translation

       character(len=*), intent(inout) :: string
!
!      Local variables.
!
       integer(kind=intType) :: jj

       real(kind=realType) :: mult
       real(kind=realType), dimension(9) :: floatArray
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Read the 9 floating point values from the string and store them
       ! in the appropriate arrays. Note that the angles are assumed to
       ! be specified in degrees and therefore they must be converted to
       ! radians and the lenghts are converted to meters.

       call readFloatsFromString(string, floatArray, 9_intType)

       rotationCenter(1) = floatArray(1)*LRef
       rotationCenter(2) = floatArray(2)*LRef
       rotationCenter(3) = floatArray(3)*LRef

       mult  = pi/180.0_realType
       rotationAngles(1) = floatArray(4)*mult
       rotationAngles(2) = floatArray(5)*mult
       rotationAngles(3) = floatArray(6)*mult

       translation(1) = floatArray(7)*LRef
       translation(2) = floatArray(8)*LRef
       translation(3) = floatArray(9)*LRef

       ! Make sure that the rotation angles are such that it corresponds
       ! to an integer value of the number of sections per wheel.

       mult = sqrt(rotationAngles(1)**2 + rotationAngles(2)**2 &
            +      rotationAngles(3)**2)

       if(mult > eps) then

         ! Nonzero angle specified. Determine the number of sections for
         ! the full wheel, which is an integer.

         jj = nint(two*pi/mult,intType)

         ! Store the correction factor for the angles in mult and correct
         ! the periodic angles accordingly.

         mult = two*pi/(jj*mult)

         rotationAngles = mult*rotationAngles

       endif

       end subroutine readPeriodicInfo

       !=================================================================

       subroutine readIntegersFromString(intString, intArray, nn)
!
!      ******************************************************************
!      *                                                                *
!      * readIntegersFromString reads nn integers from intString and    *
!      * stores them in intArray. The integers are removed from         *
!      * intString afterwards. It is assumed that upon entry the first  *
!      * character of intString is not a space.                         *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType),               intent(in)  :: nn
       integer(kind=intType), dimension(*), intent(out) :: intArray

       character(len=*), intent(inout) :: intString
!
!      Local variables.
!
       integer :: pos, ierr

       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       do i=1,nn

         ! Check if the string still contains some information.
         ! If not terminate.

         if(len_trim(intString) == 0) then
           if(myID == 0) &
             call terminate("readIntegersFromString", &
                            "Unexpected end of string")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Read the integer from intString and remove it from the
         ! string.

         read(intString,*) intArray(i)

         pos = index(intString," ")
         if(pos > 0) intString = intString(pos+1:)

         intString = adjustl(intString)
         intString = trim(intString)

       enddo

       end subroutine readIntegersFromString

       !=================================================================

       subroutine readFloatsFromString(floatString, floatArray, nn)
!
!      ******************************************************************
!      *                                                                *
!      * ReadFloatsFromString reads nn floating point values from       *
!      * floatString and stores them in floatArray. The values are      *
!      * removed from floatString afterwards.                           *
!      * It is assumed that upon entry the first character of           *
!      * floatString is not a space.                                    *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: nn

       real(kind=realType), dimension(*), intent(out) :: floatArray

       character(len=*), intent(inout) :: floatString
!
!      Local variables.
!
       integer :: pos, ierr

       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       do i=1,nn

         ! Check if the string still contains some information.
         ! If not terminate.

         if(len_trim(floatString) == 0) then
           if(myID == 0) &
             call terminate("readFloatsFromString", &
                            "Unexpected end of string")
           call mpi_barrier(SUmb_comm_world, ierr)
         endif

         ! Read the floating point value from floatString and remove
         ! it from the string.

         read(floatString,*) floatArray(i)

         pos = index(floatString," ")
         if(pos > 0) floatString = floatString(pos+1:)

         floatString = adjustl(floatString)
         floatString = trim(floatString)

       enddo

       end subroutine readFloatsFromString
