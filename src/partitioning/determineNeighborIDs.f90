!
!      ******************************************************************
!      *                                                                *
!      * File:          determineNeighborIDs.f90                        *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 12-18-2002                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineNeighborIDs
!
!      ******************************************************************
!      *                                                                *
!      * determineNeighborIDs determines for every internal block       *
!      * boundary the block ID of the neighbor. In the cgns file only   *
!      * the zone name is stored, but the ID's are more useful          *
!      * internally.                                                    *
!      *                                                                *
!      * Although for this case a quadratic search algorithm is not too *
!      * bad (number of blocks are O(1000) maximum), I don't like the   *
!      * idea of having a quadratic loop in the code. That's why a      *
!      * O(n log(n)) algorithm is used here.                            *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       implicit none
!
!      Local variables
!
       character(len=maxCGNSNameLen), dimension(cgnsNDom) :: zoneNames

       integer(kind=intType), dimension(cgnsNDom) :: zoneNumbers

       integer(kind=intType) :: i, j, k, ii

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
         ii = bsearchStrings(cgnsDoms(i)%zoneName, zoneNames, cgnsNDom)

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

           ii = bsearchStrings(cgnsDoms(i)%conn1to1(j)%donorName, &
                               zoneNames, cgnsNDom)
           if(ii == 0)                              &
             call terminate("determineNeighborIDs", &
                            "donor name not found in sorted zone names")

           cgnsDoms(i)%conn1to1(j)%donorBlock = zoneNumbers(ii)

         enddo

         ! The non-matching abutting connectivities.

         do j=1,cgnsDoms(i)%nNonMatchAbutting

           ! Determine the neighbor ID's for this subface.

           do k=1,cgnsDoms(i)%connNonMatchAbutting(j)%nDonorBlocks

             ii = bsearchStrings(                                      &
                    cgnsDoms(i)%connNonMatchAbutting(j)%donorNames(k), &
                    zoneNames, cgnsNDom)
             if(ii == 0)                              &
               call terminate("determineNeighborIDs", &
                              "donor name not found in sorted zone names")

             cgnsDoms(i)%connNonMatchAbutting(j)%donorBlocks(k) = &
                                                          zoneNumbers(ii)
           enddo
         enddo

         ! The overset connectivities.

         do j=1,cgnsDoms(i)%nOverset

           ! Determine the neighbor ID for this overset boundary.

           ii = bsearchStrings(cgnsDoms(i)%connOver(j)%donorName, &
                               zoneNames, cgnsNDom)
           if(ii == 0)                              &
             call terminate("determineNeighborIDs", &
                            "donor name not found in sorted zone names")

           cgnsDoms(i)%connOver(j)%donorBlock = zoneNumbers(ii)

         enddo

       enddo domains

       end subroutine determineNeighborIDs
