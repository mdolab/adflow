!
!      ******************************************************************
!      *                                                                *
!      * File:          countConnectivities.F90                         *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 11-07-2005                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine countConnectivities(cgnsInd, cgnsBase, nZone)
!
!      ******************************************************************
!      *                                                                *
!      * countConnectivities determines the number of connectivities    *
!      * for each of the supported types stored in 1to1 and general.    *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use constants
       use su_cgns
       use partitionMod
       implicit none
!
!      Subroutine arguments
!
       integer, intent(in)  :: cgnsInd, cgnsBase, nZone
!
!      Local variables.
!
       integer :: i, ngeneral, ierr
       integer :: n1to1, n1to1General, nNonMatch, nOverset

       integer :: location, connectType, ptsetType, npnts
       integer :: donorZoneType, donorPtsetType, donorDatatype
       integer :: ndataDonor

       integer, dimension(:),     allocatable :: connIDNonMatch
       integer, dimension(:,:),   allocatable :: donorData
       integer, dimension(:,:,:), allocatable :: myRangeNonMatch

       integer(kind=intType) :: mm, nn

       integer(kind=intType), dimension(:), allocatable :: multSubfaces

       type(subfaceNonMatchType), dimension(:), allocatable :: &
                                                        subfaceNonMatch

       type(cgnsNonMatchAbuttingConnType), pointer, dimension(:) :: &
                                                            connNonMatch

       character(len=maxStringLen)   :: errorMessage
       character(len=maxCGNSNameLen) :: connectName, donorName
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS

       call terminate("countConnectivities", &
                      "Routine should not be called if no cgns support &
                      &is selected.")

#else
       ! Determine the number of 1 to 1 connectivities in this zone.
       ! Note that the reading takes place via an integer type.

       call cg_n1to1_f(cgnsInd, cgnsBase, nZone, i, ierr)
       if(ierr /= all_ok)                      &
         call terminate("countConnectivities", &
                        "Something wrong when calling cg_n1to1_f")

       n1to1 = i

       ! Determine the total number of general connectivities in this
       ! zone.

       call cg_nconns_f(cgnsInd, cgnsBase, nZone, ngeneral, ierr)
       if(ierr /= all_ok)                      &
         call terminate("countConnectivities", &
                        "Something wrong when calling cg_nconns_f")

       ! Allocate the memory for connIDNonMatch and myRangeNonMatch. Note
       ! that this number is an upper bound, because other connectivities
       ! may be present in general connectivities.

       allocate(connIDNonMatch(ngeneral), &
                myRangeNonMatch(3,2,ngeneral),stat=ierr)
       if(ierr /= all_ok)                      &
         call terminate("countConnectivities", &
                        "Memory allocation failure for connIDNonMatch &
                        &and myRangeNonMatch")

       ! Loop over ngeneral to find out how many of each supported
       ! types of connectivities are stored here.

       n1to1General = 0
       nNonMatch    = 0
       nOverset     = 0

       do i=1,ngeneral

         ! Read the information of this connectivity.

         call cg_conn_info_f(cgnsInd, cgnsBase, nZone, i, connectName, &
                             location, connectType, ptsetType, npnts,  &
                             donorName, donorZoneType, donorPtsetType, &
                             donorDatatype, ndataDonor, ierr)
         if(ierr /= all_ok)                      &
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
                call mpi_barrier(SUmb_comm_world, ierr)

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
                if(ierr /= all_ok)                      &
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
                call mpi_barrier(SUmb_comm_world, ierr)

              endif

            !============================================================

            case (Overset)

              if(location       == CellCenter   .and. &
                 ptsetType      == PointList    .and. &
                 donorZoneType  == Structured   .and. &
                 donorPtsetType == CellListDonor) then

                nOverset = nOverset + 1

              else

                ! CGNS format not supported.

                write(errorMessage,103) trim(cgnsDoms(nZone)%zoneName), &
                                        trim(connectName)
                if(myID == 0) &
                  call terminate("countConnectivities", errorMessage)
                call mpi_barrier(SUmb_comm_world, ierr)

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

       ! Memory allocation for the 1 to 1 and overset connectivities.

       allocate(cgnsDoms(nZone)%conn1to1(n1to1), &
                cgnsDoms(nZone)%connOver(nOverset), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("countConnectivities", &
                        "Memory allocation failure for conn1to1 and &
                        &connOver")
!
!      ******************************************************************
!      *                                                                *
!      * For the non-matching abutting subfaces some more information   *
!      * needs to be extracted. The reason is that a subface abuts      *
!      * multiple blocks and in CGNS this info is stored in multiple    *
!      * connectivities. However it is a lot easier to store that info  *
!      * together. That's why the non-abbuting subfaces must be sorted  *
!      * in increasing order to extract this information.               *
!      *                                                                *
!      ******************************************************************
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
!      ******************************************************************
!      *                                                                *
!      * Store the number of connectivities in cgnsDoms(nZone).         *
!      *                                                                *
!      ******************************************************************
!
       cgnsDoms(nZone)%n1to1             = n1to1
       cgnsDoms(nZone)%n1to1General      = n1to1General
       cgnsDoms(nZone)%nNonMatchAbutting = nNonMatch
       cgnsDoms(nZone)%nOverset          = nOverset
!
!      ******************************************************************
!      *                                                                *
!      * Format statements.                                             *
!      *                                                                *
!      ******************************************************************
!
 101   format("Zone",1x,a,", connectivity", 1x,a, ": No support for &
              &this format of an abutting 1 to 1 connectivity")
 102   format("Zone",1x,a,", connectivity", 1x,a, ": No support for &
              &this format of a non-matching abutting connectivity")
 103   format("Zone",1x,a,", connectivity", 1x,a, ": No support for &
              &this format of an overset connectivity")

#endif

       end subroutine countConnectivities
