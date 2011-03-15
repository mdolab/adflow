!
!      ******************************************************************
!      *                                                                *
!      * File:          readGeneralConn.F90                             *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 07-14-2002                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readGeneralConn(cgnsInd, cgnsBase, nZone)
!
!      ******************************************************************
!      *                                                                *
!      * readGeneralConn reads and converts the cgns general            *
!      * connectivities.  Supported connectivites are 1-to-1,           *
!      * non-matching abutting and overset types.                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use su_cgns
       use inputOverset
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in) :: cgnsInd, cgnsBase, nZone

#ifdef USE_NO_CGNS

       call terminate("readGeneralConn", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
!
!      Local variables.
!
       character(len=maxStringLen) :: errorMessage

       integer :: i, j, nn, nGeneral, n1to1, nOver, ierr
       integer :: location, connectType, ptsetType, npnts
       integer :: donorZoneType, donorPtsetType, donorDatatype
       integer :: nDataDonor, id, jj, realTypeCGNS
       integer :: nArrays, dataType, dataDim

       integer, dimension(2)   :: dimVector
       integer, dimension(3)   :: ii, transform
       integer, dimension(3,2) :: myRange

       integer, dimension(:,:), allocatable :: myData, donorData
       integer, dimension(:,:), allocatable :: map2NonMatch

       real(kind=realType), dimension(3) :: rotationCenter
       real(kind=realType), dimension(3) :: rotationAngles
       real(kind=realType), dimension(3) :: translation

       logical :: periodic, wrongData

       character(len=maxCGNSNameLen) :: connectName, donorName, arrayName

       type(cgns1to1ConnType),    pointer, dimension(:) :: conn1to1
       type(cgnsOversetConnType), pointer, dimension(:) :: connOver

       type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
                                                             connNonMatch
!
!      Function definition.
!
       integer :: setCGNSRealType
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the cgns real type for reading overset interpolants.

       realTypeCGNS = setCGNSRealType()

       ! Set some pointers for the connectivities to make the code
       ! more readable.

       conn1to1     => cgnsDoms(nZone)%conn1to1
       connNonMatch => cgnsDoms(nZone)%connNonMatchAbutting
       connOver     => cgnsDoms(nZone)%connOver

       ! Set the counter n1to1 to the currently stored number of 1 to 1
       ! block connectivities, and initialize other counters.

       n1to1 = cgnsDoms(nZone)%n1to1 - cgnsDoms(nZone)%n1to1General
       nOver = 0

       ! Determine the number of general connectivities.

       call cg_nconns_f(cgnsInd, cgnsBase, nZone, ngeneral, ierr)
       if(ierr /= all_ok)                  &
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
       if(ierr /= all_ok)                  &
         call terminate("readGeneralConn", &
                        "Memory allocation failure for map2NonMatch")

       do i=1,cgnsDoms(nZone)%nNonMatchAbutting
         do j=1,connNonMatch(i)%nDonorBlocks
           nn = connNonMatch(i)%donorBlocks(j)
           map2NonMatch(nn,1) = i
           map2NonMatch(nn,2) = j
         enddo
       enddo

       ! Initialize the total number of overset cells to 0.

       cgnsDoms(nZone)%nCellsOverset = 0

       ! Loop over the general connectivities.

       nConnLoop: do nn=1,ngeneral

         ! Read the information of this connectivity.

         call cg_conn_info_f(cgnsInd, cgnsBase, nZone, nn, connectName, &
                             location, connectType, ptsetType, npnts,   &
                             donorName, donorZoneType, donorPtsetType,  &
                             donorDatatype, ndataDonor, ierr)
         if(ierr /= all_ok)                  &
           call terminate("readGeneralConn", &
                          "Something wrong when calling cg_conn_info_f")

         ! Read the data based on the type of connectivity.

         connectivityType: select case(connectType)

           case (Abutting1to1)
!
!            ************************************************************
!            *                                                          *
!            * 1-to-1 connectivity stored as a general one. Note that   *
!            * the check for a valid one has already been done in       *
!            * countConnectivities.                                     *
!            *                                                          *
!            ************************************************************
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
             if(ierr /= all_ok)                  &
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
!            ************************************************************
!            *                                                          *
!            * Non-matching abutting connectivity. Note that the        *
!            * check for a valid one has already been done in           *
!            * countConnectivities.                                     *
!            *                                                          *
!            ************************************************************
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
             if(ierr /= all_ok)                  &
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
               call mpi_barrier(SUmb_comm_world, ierr)

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
               call mpi_barrier(SUmb_comm_world, ierr)

             endif checkConsistency

           !=============================================================

           case (Overset)
!
!            ************************************************************
!            *                                                          *
!            * Overset Connectivity. Note that the check for a valid    *
!            * one has already been done in countConnectivities.        *
!            *                                                          *
!            ************************************************************
!
             ! Update the counter and store some info in connOverset.

             nOver = nOver + 1

             connOver(nover)%connectName = connectName
             connOver(nover)%donorName   = donorName
             connOver(nover)%npnts       = npnts

             ! Update the number of overset cells for this block.

             cgnsDoms(nZone)%ncellsOverset = &
                                cgnsDoms(nZone)%ncellsOverset + npnts

             ! Check to make sure that the number of donor points
             ! equals the number of points for this zone to be
             ! interpolated.

             if(npnts /= ndataDonor)             &
               call terminate("readGeneralConn", &
                              "Number of donor and boundary points &
                              &must equal for overset")

             ! Allocate all the memory for this zone and donorData.

             allocate(myData(3,npnts), donorData(3,npnts), &
                      connOver(nover)%ibndry(3,npnts),     &
                      connOver(nover)%idonor(3,npnts), stat=ierr)
             if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
                              "Memory allocation failure for &
                              &overset data")

             ! Read the indices of the connectivity.

             call cg_conn_read_f(cgnsInd, cgnsBase, nZone, nn, &
                                 myData, Integer, donorData, ierr)
             if(ierr /= all_ok)                  &
               call terminate("readGeneralConn", &
                              "Something wrong when calling &
                              &cg_conn_read_f")

             ! Store the indices of the boundary and donor cells.

             connOver(nover)%ibndry = myData
             connOver(nover)%idonor = donorData

             ! Release the memory of the temporary data arrays.

             deallocate(myData, donorData, stat=ierr)
             if(ierr /= 0)                       &
               call terminate("readGeneralConn", &
                              "Deallocation error for temp &
                              &overset data")

             ! If the input donors are being treated as guesses,
             ! the interpolants are ignored and the array is
             ! allocated with zero size to avoid inconsistencies
             ! in memory release.

             checkReadInterp: if (oversetDonorsAreGuesses) then

               allocate(connOver(nover)%interp(0,npnts), stat=ierr)
               if(ierr /= 0)                       &
                 call terminate("readGeneralConn", &
                                "Memory allocation failure for interp")

             else checkReadInterp

               ! Goto this connectivity's node in the file and read
               ! the number of data arrays.

               call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t",  &
                              nZone, "ZoneGridConnectivity_t", 1, &
                              "GridConnectivity_t", nn, "end")
               if(ierr /= all_ok)                  &
                 call terminate("readGeneralConn", &
                                "Something wrong when calling cg_goto_f")

               call cg_narrays_f(nArrays, ierr)
               if(ierr /= all_ok)                  &
                 call terminate("readGeneralConn", &
                                "Something wrong when calling &
                                &cg_narrays_f")

               ! Loop over the number of data arrays and look for
               ! one with the "InterpolantsDonor" name per the SIDS.

               do j = 1,nArrays
                 call cg_array_info_f(j, arrayName, dataType, &
                                      dataDim, dimVector, ierr)

                 if (trim(arrayName) == "InterpolantsDonor") exit
               enddo

               ! If the interpolants array was not found print
               ! an error.

               if(j > nArrays) then
                 write(errorMessage,102) trim(cgnsDoms(nZone)%zoneName), &
                                         trim(connectName)
                 if(myID == 0) &
                   call terminate("readGeneralConn", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif

               ! Check that the dimensions of the array are compatible
               ! with the interpolation type and the size is correct.
               ! Allocate memory for the interpolants an then read
               ! them if it's okay, otherwise print an error.

               if(dataDim == 2 .and. dimVector(2) == npnts .and.        &
                 (dimVector(1) == nDonorWeights(oversetInterpType) .or. &
                  dimVector(1) == 3)) then

                 jj = dimVector(1)
                 allocate(connOver(nover)%interp(jj,npnts), stat=ierr)
                 if(ierr /= 0)                       &
                   call terminate("readGeneralConn", &
                                  "Memory allocation failure for &
                                  &interpolants")

                 call cg_array_read_as_f(j, realTypeCGNS, &
                                         connOver(nover)%interp, ierr)
                 if(ierr /= all_ok)                  &
                   call terminate("readGeneralConn", &
                                  "Something wrong when calling &
                                  &cg_array_read_as")
               else
                 write(errorMessage,103) trim(cgnsDoms(nZone)%zoneName), &
                                         trim(connectName)
                 if(myID == 0) &
                   call terminate("readGeneralConn", errorMessage)
                 call mpi_barrier(SUmb_comm_world, ierr)
               endif
             endif checkReadInterp

         end select connectivityType

       enddo nConnLoop

       ! Release the memory of map2NonMatch again.

       deallocate(map2NonMatch, stat=ierr)
       if(ierr /= all_ok)                  &
         call terminate("readGeneralConn", &
                        "Deallocation failure for map2NonMatch")

       ! Format statements.

 100   format("Zone",1x,a,", connectivity", 1x,a, &
              ": Invalid donor subface.")
 101   format("Zone",1x,a,", connectivity", 1x,a, &
              ": Inconsistent periodic info compared to connectivity", &
              1x,a,".")
 102   format("Zone",1x,a,", connectivity", 1x,a, &
              ": InterpolantsDonor node is missing.")
 103   format("Zone",1x,a,", connectivity", 1x,a, &
              ": Wrong number or size of interpolants array.")

#endif

       end subroutine readGeneralConn
