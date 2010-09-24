!
!      ******************************************************************
!      *                                                                *
!      * File:          getSortedZoneNumbers.F90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-16-2003                                      *
!      * Last modified: 10-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine getSortedZoneNumbers
!
!      ******************************************************************
!      *                                                                *
!      * getSortedZoneNumbers reads the names of the zones of the       *
!      * cgns file given by cgnsInd and cgnsBase. Afterwards the        *
!      * zonenames are sorted in increasing order, such that a binary   *
!      * search algorithm can be employed. The original zone numbers    *
!      * are stored in zoneNumbers.                                     *
!      * If the zone contains a link to a zone containing the           *
!      * coordinates the name of the linked zone is taken.              *
!      *                                                                *
!      ******************************************************************
!
#ifdef USE_NO_CGNS
       call terminate("getSortedZoneNumbers", &
                      "Routine should not be called if no cgns support &
                      &is selected.")
#else
       use cgnsGrid
       use su_cgns
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr
       integer :: zone, zonetype, ncoords, pathLength
       integer :: pos

       integer, dimension(9) :: sizesBlock

       integer(kind=intType) :: nn, ii

       character(len=maxStringLen) :: errorMessage, linkPath
       character(len=maxCGNSNameLen), dimension(cgnsNdom) :: tmpNames

       logical :: nameFound

       character(len=7) :: int1String, int2String
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
       ! Allocate the memory for zoneNames and zoneNumbers.

       allocate(zoneNames(cgnsNdom), zoneNumbers(cgnsNdom), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("getSortedZoneNumbers", &
                        "Memory allocation failure for zoneNames &
                        &and zoneNumbers")

       ! Loop over the number of zones in the file.

       cgnsDomains: do nn=1,cgnsNdom

         ! Initialize nameFound to .false.

         nameFound = .false.

         ! Check if the zone is structured.

         zone = nn
         call cg_zone_type_f(cgnsInd, cgnsBase, zone, zonetype, ierr)
         if(ierr /= all_ok)                       &
           call terminate("getSortedZoneNumbers", &
                          "Something wrong when calling cg_zone_type_f")

         if(zonetype /= structured) then

           write(int1String,"(i7)") cgnsBase
           int1String = adjustl(int1String)
           write(int2String,"(i7)") zone
           int2String = adjustl(int2String)

           write(errorMessage,100) trim(int1String), trim(int2String)
 100       format("Base",1X,A,": Zone",1X,A, " of the cgns restart &
                  &file is not structured")
           call terminate("getSortedZoneNumbers", errorMessage)

         endif

         ! Determine the number of grid coordinates of this zone.

         call cg_ncoords_f(cgnsInd, cgnsBase, zone, ncoords, ierr)
         if(ierr /= all_ok)                       &
           call terminate("getSortedZoneNumbers", &
                          "Something wrong when calling cg_ncoords_f")

         ! If ncoords == 3, there are coordinates present. Then check if
         ! it is a link. If so, take the zone name of the link.

         if(ncoords == 3) then

           ! Go to the coordinates node.

           call cg_goto_f(cgnsInd, cgnsBase, ierr, "Zone_t", zone, &
                          "GridCoordinates_t", 1, "end")
           if(ierr /= all_ok)                       &
             call terminate("getSortedZoneNumbers", &
                            "Something wrong when calling cg_goto_f")

           ! Check if this node is a link.

           call cg_is_link_f(pathLength, ierr)
           if(ierr /= all_ok)                       &
             call terminate("getSortedZoneNumbers", &
                            "Something wrong when calling cg_is_link_f")

           if(pathLength > 0) then

             ! Determine the name of the linkPath.

             call cg_link_read_f(errorMessage, linkPath, ierr)
             if(ierr /= all_ok)                       &
               call terminate("getSortedZoneNumbers", &
                              "Something wrong when calling &
                              &cg_link_read_f")

             ! Find the zone name.
             ! Find, starting from the back, the forward slash.

             pos = index(linkPath, "/", .true.)
             if(pos > 0) then
               linkPath = linkPath(:pos-1)

               ! Find the next forward slash from the back and
               ! remove the leading part from the path name.

               pos = index(linkPath, "/", .true.)
               if(pos > 0) linkPath = linkPath(pos+1:)
             endif

             ! Create the zone name and set nameFound to .true..

             linkPath      = adjustl(linkPath)
             zoneNames(nn) = trim(linkPath)
             nameFound     = .true.

           endif
         endif

         ! If no name was yet found set it to the name of the
         ! current zone.

         if(.not. nameFound) then
           call cg_zone_read_f(cgnsInd, cgnsBase, zone, &
                               zoneNames(nn), sizesBlock, ierr)
           if(ierr /= all_ok)                       &
             call terminate("getSortedZoneNumbers", &
                            "Something wrong when calling &
                            &cg_zone_read_f")
         endif

       enddo cgnsDomains

       ! Set tmpNames to zoneNames and sort the latter
       ! in increasing order.

       do nn=1,cgnsNdom
         tmpNames(nn) = zoneNames(nn)
       enddo

       ! Sort zoneNames in increasing order.

       call qsortStrings(zoneNames, cgnsNdom)

       ! Initialize zoneNumbers to -1. This serves as a check during
       ! the search.

       do nn=1,cgnsNdom
         zoneNumbers(nn) = -1
       enddo

       ! Find the original zone numbers for the sorted zone names.

       do nn=1,cgnsNdom
         ii = bsearchStrings(tmpNames(nn), zoneNames, cgnsNdom)

         ! Check if the zone number is not already taken. If this is the
         ! case, this means that two identical zone names are present.

         if(zoneNumbers(ii) /= -1)                &
           call terminate("getSortedZoneNumbers", &
                          "Error occurs only when two identical zone &
                          &names are present")

         ! And set the zone number.

         zoneNumbers(ii) = nn
       enddo

#endif

       end subroutine getSortedZoneNumbers
