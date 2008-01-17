!
!      ******************************************************************
!      *                                                                *
!      * File:          readBocos.F90                                   *
!      * Author:        Edwin van der Weide, Steve Repsher,             *
!      *                Seonghyeon Hahn                                 *
!      * Starting date: 12-17-2002                                      *
!      * Last modified: 08-09-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readBocos(cgnsInd, cgnsBase, nZone,                &
                            nDoubleBoundFaces, sortedFamName, famID)
!
!      ******************************************************************
!      *                                                                *
!      * ReadBocos reads the boundary condition info for the given      *
!      * zone/block.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use su_cgns
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in)    :: cgnsInd, cgnsBase, nZone
       integer, intent(inout) :: nDoubleBoundFaces

       character(len=*), dimension(*), intent(in) :: sortedFamName
       integer(kind=intType), dimension(*), intent(in) :: famID
!
!      Local variables
!
       integer :: cgnsNBocos, cgnsNpnts, cgnsNDataSet, nUserData
       integer :: i, j, match
       integer :: ierr, dummy
       integer :: dirichletFlag, neumannFlag

       integer, dimension(3,2) :: bcRange

       integer(kind=intType) :: ii, nn

       character(len=maxCGNSNameLen) :: familyName
       character(len=maxStringLen)    :: errorMessage

       logical :: familySpecifiedData

       type(cgnsBcDatasetType), pointer, dimension(:) :: dataSet
!
!      Function definitions.
!
       integer(kind=intType) :: bsearchStrings
       integer(kind=intType) :: internalBC
       logical               :: checkForDoubleBoundFace
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("readBocos", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Read the number of boundary conditions in this zone/block.
       ! Again the reading takes place via an integer.

       call cg_nbocos_f(cgnsInd, cgnsBase, nZone, cgnsNBocos, ierr)
       if(ierr /= all_ok)            &
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
!        ****************************************************************
!        *                                                              *
!        * Read the general info for this boundary condition and set    *
!        * the dimensions of the subface.                               *
!        *                                                              *
!        ****************************************************************
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
         if(ierr /= all_ok)            &
           call terminate("readBocos", &
                          "Something wrong when calling cg_boco_info_f")

         cgnsDoms(nZone)%bocoInfo(i)%npnts = cgnsNpnts

         ! Nullify the pointer for dataSet.

         nullify(cgnsDoms(nZone)%bocoInfo(i)%dataSet)

         ! Perform some checks.

         if(cgnsDoms(nZone)%bocoInfo(i)%normalListFlag /= 0) &
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
           if(ierr /= all_ok)             &
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
           if(ierr /= all_ok)             &
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
!        ****************************************************************
!        *                                                              *
!        * Determine the internally used boundary condition and whether *
!        * or not the boundary condition is given on a per family basis.*
!        *                                                              *
!        ****************************************************************
!
         cgnsDoms(nZone)%bocoInfo(i)%familyID = 0

         checkActualFace: if( cgnsDoms(nZone)%bocoInfo(i)%actualFace ) then

           ! Determine the type of CGNS boundary condition and act
           ! accordingly.

           select case (cgnsDoms(nZone)%bocoInfo(i)%BCTypeCGNS)

             case (FamilySpecified)

               ! Boundary condition is specified per family.

               ! Find out the family name to which this boundary
               ! face belongs.

               call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                              "ZoneBC_t", 1, "BC_t", i, "end")
               if(ierr /= all_ok)             &
                 call terminate("readBocos", &
                                "Something wrong when calling cg_goto_f")

               call cg_famname_read_f(familyName, ierr)
               if(ierr /= all_ok) then

                 write(errorMessage,101) trim(cgnsDoms(nZone)%zoneName), &
                             trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
                 if(myID == 0) call terminate("readBocos", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)

               endif

               ! Search the family name in the sorted names. For a valid
               ! grid this name must be found.

               nn = cgnsNFamilies
               ii = bsearchStrings(familyName, sortedFamName, nn)
               if(ii == 0) then

                 write(errorMessage,102) trim(familyName)
                 if(myID == 0) call terminate("readBocos", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)

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
               if(ierr /= all_ok)             &
                 call terminate("readBocos", &
                                "Something wrong when calling cg_goto_f")

               call cg_nuser_data_f(nUserData, ierr)
               if(ierr /= all_ok)            &
                 call terminate("readBocos", &
                                "Something wrong when calling &
                                &cg_nuser_data_f")

               ! nUserData should be 1. Check this.

               if(nUserData /= 1) then
                 write(errorMessage,103) trim(cgnsDoms(nZone)%zoneName), &
                             trim(cgnsDoms(nZone)%bocoInfo(i)%bocoName)
                 if(myID == 0) call terminate("readBocos", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

               ! Read the name of the user defined data node.

               call cg_user_data_read_f(nUserData,                      &
                           cgnsDoms(nZone)%bocoInfo(i)%userDefinedName, &
                           ierr)
               if(ierr /= all_ok)            &
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
                 call mpi_barrier(SUmb_comm_world, ierr)
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
                   call mpi_barrier(SUmb_comm_world, ierr)

               end select

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
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

           end select

         endif checkActualFace
!
!        ****************************************************************
!        *                                                              *
!        * Initialize slidingID to 0 to indicate that this boco does    *
!        * not belong to a sliding mesh interface. If it is, this will  *
!        * be overwritten after all bocos are read for every zone.      *
!        *                                                              *
!        ****************************************************************
!
         cgnsDoms(nzone)%bocoInfo(i)%slidingID = 0
!
!        ****************************************************************
!        *                                                              *
!        * Determine the possible rotating rate of the boundary face.   *
!        *                                                              *
!        ****************************************************************
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
!        ****************************************************************
!        *                                                              *
!        * Read and store the prescribed boundary condition data sets.  *
!        *                                                              *
!        ****************************************************************
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

 101   format("Zone",1x,a,", boundary face", 1x,a, &
              ": Corresponding family name not given.")
 102   format("Family name",1x,a,1x,"not present in the grid")
 103   format("Zone",1x,a,", boundary face", 1x,a, &
              ": Need 1 UserDefinedData_t node for user defined &
              &boundary condition")
 104   format("Zone",1x,a,", boundary face", 1x,a, &
              ": Unknown user-defined boundary condition", 1x,a)
 105   format("Zone",1x,a,", boundary face", 1x,a, &
              ": User-defined boundary condition", 1x,a,1x, &
              "only possible for a family")
 106   format("Zone",1x,a,", boundary face", 1x,a, &
              ": boundary condition type missing or not supported")

       !=================================================================

       contains

         !===============================================================

         subroutine readBCDataArrays(nArr, arr, DirNeu)
!
!        ****************************************************************
!        *                                                              *
!        * readBCDataArrays reads the arrays of the given data set      *
!        * from the cgns file.                                          *
!        *                                                              *
!        ****************************************************************
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
!
!        Function definition.
!
         integer :: setCGNSRealType
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
         ! Set the cgns real type.

         realTypeCGNS = setCGNSRealType()

         ! Go to the correct node of the given boundary subface.

         call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", nZone, &
                        "ZoneBC_t", 1, "BC_t", i, "BCDataSet_t", j, &
                        "BCData_t", DirNeu, "end")
         if(ierr /= all_ok)                     &
           call terminate("readBCDataArrays", &
                          "Something wrong when calling cg_goto_f")

         ! Determine the amount of data arrays present for this node.

         call cg_narrays_f(nArrays, ierr)
         if(ierr /= all_ok)                   &
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

         if(ierr == all_ok) then
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
           if(ierr /= all_ok)                     &
             call terminate("readBCDataArrays", &
                            "Something wrong when calling cg_goto_f")

           ! Determine the name and the dimensions of the array.

           call cg_array_info_f(k, arr(k)%arrayName, dataType,     &
                                arr(k)%nDimensions, arr(k)%dataDim, &
                                ierr)
           if(ierr /= all_ok)                     &
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
           if(ierr /= all_ok)                     &
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
           if(ierr /= all_ok)                     &
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

           if(ierr == all_ok) then

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
 100           format("# Zone ",a,", boundary subface ", a)
 101           format("# BC data set ",a, &
                      ": No units specified, assuming SI units")
             endif
           endif

         enddo loopDataArrays

         end subroutine readBCDataArrays

#endif

       end subroutine readBocos

!      ==================================================================

       logical function checkForDoubleBoundFace(nZone, nBound)
!
!      ******************************************************************
!      *                                                                *
!      * CheckForDoubleBoundFace checks whether the given boundary      *
!      * range for the given zone has already been defined in the 1 to  *
!      * 1 block connectivities. If so .true. is returned, otherwise    *
!      * .false.                                                        *
!      *                                                                *
!      ******************************************************************
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
