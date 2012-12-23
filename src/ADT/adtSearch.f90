!
!     ******************************************************************
!     *                                                                *
!     * File:          adtSearch.f90                                   *
!     * Author:        Edwin van der Weide                             *
!     * Starting date: 02-10-2006                                      *
!     * Last modified: 03-17-2006                                      *
!     *                                                                *
!     ******************************************************************
!
      module adtSearch
!
!     ******************************************************************
!     *                                                                *
!     * Module which contains the subroutines for the global search.   *
!     *                                                                *
!     ******************************************************************
!
      use adtLocalSearch
      use adtUtils
      implicit none

      !=================================================================

      contains

        !===============================================================

        subroutine containmentSearch(nCoor,       coor,       &
                                     adtID,       procID,     &
                                     elementType, elementID,  &
                                     uvw,         nInterpol,  &
                                     arrDonor,    arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found the corresponding entry in procID is  *
!       * set to -1 to indicate failure.                               *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor:     Number of coordinates for which the element must  *
!       *            be determined.                                    *
!       * coor:      The coordinates of these points.                  *
!       * adtID:     The ADT to be searched.                           *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point.                       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                             elementType
        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol
!
!       Local variables.
!
        integer :: ierr

        integer(kind=intType) :: jj, nAlloc

        real(kind=realType), dimension(1) :: dummy
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Determine the index in the array ADTs, which stores the given
        ! ID. As the number of trees stored is limited, a linear search
        ! algorithm is okay.

        nAlloc = ubound(ADTs, 1)
        do jj=1,nAlloc
          if(adtID == ADTs(jj)%adtID) exit
        enddo

        ! Check if the ADT to be searched exists. If not stop.
        ! Note that adtTerminate is not called. The reason is that the
        ! processor ID is not known.

        if(jj > nAlloc) stop "ADT to be searched does not exist."

        ! Check if the ADT corresponds to a volume grid. If not terminate.

        if(ADTs(jj)%adtType /= adtVolumeADT) then
          if(ADTs(jj)%myID == 0)                       &
            call adtTerminate(jj, "containmentSearch", &
                              "ADT does not contain a volume mesh.")
          call mpi_barrier(ADTs(jj)%comm, ierr)
        endif

        ! Initialize the search, i.e. determine the number of
        ! coordinates to be searched in each of the local ADT's.

        call initSearch(nCoor, coor, dummy, jj, .true.)

        ! Perform the actual search.

        call search(nCoor,     coor,      procID,   elementType, &
                    elementID, uvw,       dummy,    jj,          &
                    .true.,    nInterpol, arrDonor, arrInterpol)

        end subroutine containmentSearch

        !***************************************************************
        !***************************************************************

        subroutine failSafeSearch(nCoor,       coor,      &
                                  adtID,       procID,    &
                                  elementType, elementID, &
                                  uvw,         dist2,     &
                                  nInterpol,   arrDonor,  &
                                  arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT, which contains that coordinate.    *
!       * If no element is found a minimum distance search is          *
!       * performed, such that always an interpolation can be          *
!       * performed. To indicate that the element does not contain the *
!       * point the element ID is negated.                             *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of coordinates for which the element must be   *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element, i.e. if  *
!       *              a minimum distance search had to be used.       *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights.                   *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                              elementType

        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol

        real(kind=realType), dimension(:), intent(inout) :: dist2
!
!       Local variables.
!
        integer :: ierr
        integer, dimension(:), allocatable :: tmpProcID

        integer(kind=intType) :: i, j, ii, jj, nAlloc, nFail
        integer(kind=intType), dimension(:), allocatable :: tmpElementID
        integer(kind=intType), dimension(:), allocatable :: coorIDs

        integer(kind=adtElementType), dimension(:), allocatable :: &
                                                          tmpElementType

        real(kind=realType), dimension(1) :: dummy

        real(kind=realType), dimension(:),   allocatable :: tmpDist2
        real(kind=realType), dimension(:,:), allocatable :: tmpCoor
        real(kind=realType), dimension(:,:), allocatable :: tmpUVW
        real(kind=realType), dimension(:,:), allocatable :: tmpArrInt
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Determine the index in the array ADTs, which stores the given
        ! ID. As the number of trees stored is limited, a linear search
        ! algorithm is okay.

        nAlloc = ubound(ADTs, 1)
        do jj=1,nAlloc
          if(adtID == ADTs(jj)%adtID) exit
        enddo

        ! Check if the ADT to be searched exists. If not stop.
        ! Note that adtTerminate is not called. The reason is that the
        ! processor ID is not known.

        if(jj > nAlloc) stop "ADT to be searched does not exist."

        ! Check if the ADT corresponds to a volume grid.
        ! If not terminate.

        if(ADTs(jj)%adtType /= adtVolumeADT) then
          if(ADTs(jj)%myID == 0)                    &
            call adtTerminate(jj, "failSafeSearch", &
                              "ADT does not contain a volume mesh.")
          call mpi_barrier(ADTs(jj)%comm, ierr)
        endif

        ! Perform the containment search.

        call initSearch(nCoor, coor, dummy, jj, .true.)

        call search(nCoor,     coor,      procID,   elementType, &
                    elementID, uvw,       dummy,    jj,          &
                    .true.,    nInterpol, arrDonor, arrInterpol)

        ! Determine the number of coordinates for which the containment
        ! search failed. Set for the other coordinates the distance to
        ! zero.

        nFail = 0
        do i=1,nCoor
          if(procID(i) == -1) then
            nFail = nFail + 1
          else
            dist2(i) = adtZero
          endif
        enddo

        ! Determine the global sum of nFail, which is stored in ii.

        call mpi_allreduce(nFail, ii, 1, sumb_integer, mpi_max, &
                           ADTs(jj)%comm, ierr)

        ! Return if ii == 0, because the minimum distance search
        ! is not needed.

        if(ii == 0) return

        ! Allocate the memory for the arrays needed for the minimum
        ! distance search.

        allocate(tmpCoor(3,nFail),      tmpProcID(nFail),           &
                 tmpElementType(nFail), tmpElementID(nFail),        &
                 tmpUVW(3,nFail),       tmpDist2(nFail),            &
                 coorIDs(nFail),        tmpArrInt(nInterpol,nFail), &
                 stat=ierr)
        if(ierr /= 0)                             &
          call adtTerminate(jj, "failSafeSearch", &
                            "Memory allocation failure for the arrays &
                            &for the minimum distance search.")

        ! Store the information needed for minimum distance search in
        ! the arrays tmpCoor, tmpDist2 and coorIDs.

        ii = 0
        do i=1,nCoor
          if(procID(i) == -1) then
            ii = ii + 1

            tmpCoor(1,ii) = coor(1,i)
            tmpCoor(2,ii) = coor(2,i)
            tmpCoor(3,ii) = coor(3,i)
            tmpDist2(ii)  = dist2(i)
            coorIDs(ii)   = i
          endif
        enddo

        ! Perform the minimum distance search for the coordinates for
        ! which the containment search failed.

        call initSearch(nFail, tmpCoor, tmpDist2, jj, .false.)

        call search(nFail,        tmpCoor,   tmpProcID, tmpElementType, &
                    tmpElementID, tmpUVW,    tmpDist2, jj,              &
                    .false.,      nInterpol, arrDonor, tmpArrInt)

        ! Copy for the successful searches the data into the arrays to
        ! be returned by this subroutine. Note that the element IDs are
        ! negated to indicate that the coordinate is outside the element.

        do i=1,nFail
          if(tmpProcID(i) /= -1) then
            ii = coorIDs(i)

            procID(ii)      =  tmpProcID(i)
            elementType(ii) =  tmpElementType(i)
            elementID(ii)   = -tmpElementID(i)
            uvw(1,ii)       =  tmpUVW(1,i)
            uvw(2,ii)       =  tmpUVW(2,i)
            uvw(3,ii)       =  tmpUVW(3,i)
            dist2(ii)       =  tmpDist2(i)

            do j=1,nInterpol
              arrInterpol(j,ii) = tmpArrInt(j,i)
            enddo
          endif
        enddo

        ! Release the memory used in the minimum distance search.

        deallocate(tmpCoor, tmpProcID, tmpElementType, tmpElementID, &
                   tmpUVW,  tmpDist2,  coorIDs,        tmpArrInt,    &
                   stat=ierr)
        if(ierr /= 0)                             &
          call adtTerminate(jj, "failSafeSearch", &
                            "Deallocation failure for the arrays &
                            &for the minimum distance search.")

        end subroutine failSafeSearch

        !***************************************************************
        !***************************************************************

        subroutine initSearch(nCoor, coor, dist2, jj, containmentSearch)
!
!       ****************************************************************
!       *                                                              *
!       * This routine performs the initialization tasks before the    *
!       * actual search takes place. It determines the number and the  *
!       * ID's of the coordinates every local tree may have to         *
!       * interpolate. From this info this routine also determines the *
!       * number of rounds needed in the actual algorithm to avoid a   *
!       * memory bottleneck.                                           *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor:             Number of local coordinates for which the *
!       *                    element must be determined.               *
!       * coor:              The coordinates of these points.          *
!       * jj:                The index in the array ADTs of the tree   *
!       *                    to be searched.                           *
!       * containmentSearch: Whether or not a containment search must  *
!       *                    be performed. If not a minimum distance   *
!       *                    search algorithm is used, which is more   *
!       *                    expensive.                                *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Guaranteed minimum distance squared of the            *
!       *        coordinates to the elements of the ADT. This array    *
!       *        should be initialized by the calling routine,         *
!       *        possibly to a large value. It is only used for a      *
!       *        minimum distance search.                              *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor
        integer(kind=intType), intent(in) :: jj

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:),   intent(inout) :: dist2

        logical, intent(in) :: containmentSearch
!
!       Local variables.
!
        integer :: ierr, nRootLeaves, comm, nProcs, myID
        integer :: myEntryInRootProcs

        integer, dimension(:),   pointer :: rootLeavesProcs

        integer(kind=intType) :: i, j, k, mm, nn

        integer(kind=intType), dimension(:), allocatable :: nCoorPerProc
        integer(kind=intType), dimension(:), allocatable :: nCoorFromProc

        real(kind=realType) :: d1, d2, dx, dy, dz

        real(kind=realType), dimension(:,:,:), pointer :: rootBBoxes
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Set some pointers to make the code more readable.

        nRootLeaves        = ADTs(jj)%nRootLeaves
        myEntryInRootProcs = ADTs(jj)%myEntryInRootProcs
        rootLeavesProcs   => ADTs(jj)%rootLeavesProcs
        rootBBoxes        => ADTs(jj)%rootBBoxes

        comm   = ADTs(jj)%comm
        nProcs = ADTs(jj)%nProcs
        myID   = ADTs(jj)%myID

        ! Determine the global maximum of nCoor. This number will serve
        ! as an upper bound for the number of points to be searched
        ! during a round.

        call mpi_allreduce(nCoor, nCoorMax, 1, sumb_integer, mpi_max, &
                           comm, ierr)
        nCoorMax = max(nCoorMax,nCoorMaxLowerLimit)
!
!       ****************************************************************
!       *                                                              *
!       * Determine for every root leaf the local coordinates, which   *
!       * should be searched in the corresponding ADT. The criterion   *
!       * for a containment search is of course different from a       *
!       * minimum distance search and therefore a distinction must be  *
!       * made between the two methods.                                *
!       *                                                              *
!       ****************************************************************
!
        ! Allocate the memory for nCoorPerRootLeaf and mCoorPerRootLeaf.
        ! Initially the latter is a copy of the former, but its data
        ! will change in the outer loop in the iterative algorithm, see
        ! adtSearch. The numbering starts at 0, because these arrays
        ! will be cumulative storage format arrays.

        allocate(nCoorPerRootLeaf(0:nProcs), &
                 mCoorPerRootLeaf(0:nProcs), stat=ierr)
        if(ierr /= 0)                         &
          call adtTerminate(jj, "initSearch", &
                            "Memory allocation failure for &
                            &nCoorPerRootLeaf and mCoorPerRootLeaf.")

        ! Determine for a minimum distance search the guaranteed minimum
        ! distance squared to each of the root leaves and store the
        ! absolute minimum.

        testMinDistance: if(.not. containmentSearch ) then

          do j=1,nRootLeaves
            do i=1,nCoor
              d1 = abs(coor(1,i) - rootBBoxes(1,1,j))
              d2 = abs(coor(1,i) - rootBBoxes(1,2,j))
              dx = max(d1,d2)

              d1 = abs(coor(2,i) - rootBBoxes(2,1,j))
              d2 = abs(coor(2,i) - rootBBoxes(2,2,j))
              dy = max(d1,d2)

              d1 = abs(coor(3,i) - rootBBoxes(3,1,j))
              d2 = abs(coor(3,i) - rootBBoxes(3,2,j))
              dz = max(d1,d2)

              d2       = dx*dx + dy*dy + dz*dz
              dist2(i) = min(dist2(i), d2)
            enddo
          enddo

        endif testMinDistance

        ! Determine the local number of coordinates to be searched in
        ! each of the local trees, i.e. the array nCoorPerRootLeaf. Note
        ! that the processor storing the ADT stores this number rather
        ! than the root leaf. Furthermore nCoorPerRootLeaf is stored in
        ! cumulative storage format and therefore an outer loop over the
        ! number of root leaves is the most logical thing to do.

        nCoorPerRootLeaf(0) = 0
        mm = 0

        loop1RootLeaves: do j=1,nRootLeaves

          ! Determine the next processor which stores a tree; initialize
          ! in the same loop nCoorPerRootLeaf of the processors in
          ! between. The counter mm contains the entry in
          ! nCoorPerRootLeaf where the number is stored. Remember that
          ! the processor ID's start at 0.

          do k=(mm+1),(rootLeavesProcs(j)+1)
            nCoorPerRootLeaf(k) = nCoorPerRootLeaf(mm)
          enddo

          mm = rootLeavesProcs(j)+1

          ! Make a distinction between a containment and a
          ! minimum distance search.

          test1Containment: if( containmentSearch ) then

            ! Containment search. Loop over the local number of
            ! coordinates and check if the coordinates are within the
            ! bounding box of the current root leaf. If so, update
            ! the counter.

            do i=1,nCoor
              if(coor(1,i) >= rootBBoxes(1,1,j) .and. &
                 coor(1,i) <= rootBBoxes(1,2,j) .and. &
                 coor(2,i) >= rootBBoxes(2,1,j) .and. &
                 coor(2,i) <= rootBBoxes(2,2,j) .and. &
                 coor(3,i) >= rootBBoxes(3,1,j) .and. &
                 coor(3,i) <= rootBBoxes(3,2,j)) &
                 nCoorPerRootLeaf(mm) = nCoorPerRootLeaf(mm) + 1
            enddo

          else test1Containment

            ! Minimum distance search. Loop over the local number of
            ! coordinates and determine the possible minimum distance
            ! squared to the bounding box of the current root leaf.
            ! If less than the currently stored guaranteed minimum
            ! distance squared, update the counter.

            do i=1,nCoor
              if(     coor(1,i) < rootBBoxes(1,1,j)) then
                dx =  coor(1,i) - rootBBoxes(1,1,j)
              else if(coor(1,i) > rootBBoxes(1,2,j)) then
                dx =  coor(1,i) - rootBBoxes(1,2,j)
              else
                dx = adtZero
              endif

              if(     coor(2,i) < rootBBoxes(2,1,j)) then
                dy =  coor(2,i) - rootBBoxes(2,1,j)
              else if(coor(2,i) > rootBBoxes(2,2,j)) then
                dy =  coor(2,i) - rootBBoxes(2,2,j)
              else
                dy = adtZero
              endif

              if(     coor(3,i) < rootBBoxes(3,1,j)) then
                dz =  coor(3,i) - rootBBoxes(3,1,j)
              else if(coor(3,i) > rootBBoxes(3,2,j)) then
                dz =  coor(3,i) - rootBBoxes(3,2,j)
              else
                dz = adtZero
              endif

              d2 = dx*dx + dy*dy + dz*dz
              if(d2 < dist2(i)) &
                nCoorPerRootLeaf(mm) = nCoorPerRootLeaf(mm) + 1

            enddo

          endif test1Containment

        enddo loop1RootLeaves

        ! Fill the rest of the array nCoorPerRootLeaf.

        do k=(mm+1),nProcs
          nCoorPerRootLeaf(k) = nCoorPerRootLeaf(mm)
        enddo

        ! Copy nCoorPerRootLeaf to mCoorPerRootLeaf.

        do j=0,nProcs
          mCoorPerRootLeaf(j) = nCoorPerRootLeaf(j)
        enddo

        ! Allocate the memory for coorPerRootLeaf.

        nn = nCoorPerRootLeaf(nProcs)
        allocate(coorPerRootLeaf(nn), stat=ierr)
        if(ierr /= 0)                                   &
          call adtTerminate(jj, "initSearch",           &
                            "Memory allocation failure for &
                            &coorPerRootLeaf.")

        ! Repeat the loop over the number of root leaves, but now store
        ! the ID's of the local coordinates.

        nn = 0
        loop2RootLeaves: do j=1,nRootLeaves

          ! Make a distinction between a containment and a
          ! minimum distance search.

          test2Containment: if( containmentSearch ) then

            ! Containment search. Store the nodes which are within the
            ! bounding box of the root leaf.

            do i=1,nCoor
              if(coor(1,i) >= rootBBoxes(1,1,j) .and. &
                 coor(1,i) <= rootBBoxes(1,2,j) .and. &
                 coor(2,i) >= rootBBoxes(2,1,j) .and. &
                 coor(2,i) <= rootBBoxes(2,2,j) .and. &
                 coor(3,i) >= rootBBoxes(3,1,j) .and. &
                 coor(3,i) <= rootBBoxes(3,2,j)) then
                 nn = nn + 1
                coorPerRootLeaf(nn) = i
              endif
            enddo

          else test2Containment

            ! Minimum distance search. Store the nodes which have a
            ! smaller possible minimum distance squared than the
            ! currently  stored value.

            do i=1,nCoor
              if(     coor(1,i) < rootBBoxes(1,1,j)) then
                dx =  coor(1,i) - rootBBoxes(1,1,j)
              else if(coor(1,i) > rootBBoxes(1,2,j)) then
                dx =  coor(1,i) - rootBBoxes(1,2,j)
              else
                dx = adtZero
              endif

              if(     coor(2,i) < rootBBoxes(2,1,j)) then
                dy =  coor(2,i) - rootBBoxes(2,1,j)
              else if(coor(2,i) > rootBBoxes(2,2,j)) then
                dy =  coor(2,i) - rootBBoxes(2,2,j)
              else
                dy = adtZero
              endif

              if(     coor(3,i) < rootBBoxes(3,1,j)) then
                dz =  coor(3,i) - rootBBoxes(3,1,j)
              else if(coor(3,i) > rootBBoxes(3,2,j)) then
                dz =  coor(3,i) - rootBBoxes(3,2,j)
              else
                dz = adtZero
              endif

              d2 = dx*dx + dy*dy + dz*dz
              if(d2 < dist2(i)) then
                nn = nn + 1
                coorPerRootLeaf(nn) = i
              endif
            enddo

          endif test2Containment

        enddo loop2RootLeaves
!
!       ****************************************************************
!       *                                                              *
!       * Determine for every tree the number of coordinates from      *
!       * other processors it should search and from that information  *
!       * the number of rounds in the outer loop of the search         *
!       * algorithm in search.                                         *
!       *                                                              *
!       ****************************************************************
!
        ! Allocate the memory for some help arrays.

        allocate(procRecv(nProcs-1),       nCoorProcRecv(nProcs-1),   &
                 nCoorPerProc(0:nProcs-1), nCoorFromProc(0:nProcs-1), &
                 stat=ierr)
        if(ierr /= 0)                         &
          call adtTerminate(jj, "initSearch", &
                            "Memory allocation failure for help arrays.")

        ! Determine the number of coordinates I want to be searched in
        ! the trees of other processors. Store the number I have to
        ! search in my own tree and set it to 0 afterwards; the local
        ! interpolations are handled differently.

        do i=0,(nProcs-1)
          nCoorPerProc(i) = nCoorPerRootLeaf(i+1) - nCoorPerRootLeaf(i)
        enddo

        nLocalInterpol     = nCoorPerProc(myID)
        nCoorPerProc(myID) = 0

        ! Communicate these numbers to the other processors.

        call mpi_alltoall(nCoorPerProc,  1, sumb_integer, &
                          nCoorFromProc, 1, sumb_integer, comm, ierr)

        ! Determine the total number of coordinates I have to interpolate
        ! and the processor ID's I will receive coordinates from.

        nn = 0
        nProcRecv = 0

        do i=0,(nProcs-1)
          if(nCoorFromProc(i) > 0) then
            nProcRecv = nProcRecv + 1
            procRecv(nProcRecv) = i
            nn = nn + nCoorFromProc(i)
          endif
        enddo

        ! Determine the number of rounds in the outer loop of the
        ! interpolation algorithm.

        i = real(nn,realType)/real(nCoorMax,realType)
        if(i*nCoorMax < nn) i = i + 1
        i = max(i, 1_intType)

        call mpi_allreduce(i, nRounds, 1, sumb_integer, mpi_max, &
                           comm, ierr)

        ! Modify the sequence of procRecv, such that number of
        ! coordinates sent from processors during a certain round
        ! is more balanced.

        nn = max(nRootLeaves,1_intType)
        j  = (real(myEntryInRootProcs,realType)/real(nn,realType)) &
           * nProcRecv

        do i=1,nProcRecv
          j = j+1
          if(j > nProcRecv) j = 1

          nn = procRecv(i)
          procRecv(i) = procRecv(j)
          procRecv(j) = nn
        enddo

        ! Store the number of coordinates to be received in
        ! nCoorProcRecv.

        do i=1,nProcRecv
          nCoorProcRecv(i) = nCoorFromProc(procRecv(i))
        enddo

        ! Release the memory of the local help arrays.

        deallocate(nCoorPerProc, nCoorFromProc, stat=ierr)
        if(ierr /= 0)                         &
          call adtTerminate(jj, "initSearch", &
                            "Deallocation failure for nCoorPerProc and &
                            &nCoorFromProc")

        end subroutine initSearch

        !***************************************************************
        !***************************************************************

        subroutine minDistanceSearch(nCoor,       coor,      &
                                     adtID,       procID,    &
                                     elementType, elementID, &
                                     uvw,         dist2,     &
                                     nInterpol,   arrDonor,  &
                                     arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine attempts for every coordinate to find the       *
!       * element in the given ADT which minimizes the distance to     *
!       * this point.                                                  *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor: Number of coordinates for which the element must be   *
!       *        determined.                                           *
!       * coor:  The coordinates of these points.                      *
!       * adtID: The ADT to be searched.                               *
!       * nInterpol: Number of variables to be interpolated.           *
!       * arrDonor:  Array with the donor data; needed to obtain the   *
!       *            interpolated data.                                *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored. If no element is found for a given      *
!       *              point the corresponding entry in procID is set  *
!       *              to -1 to indicate failure. Remember that the    *
!       *              processor ID's start at 0 and not at 1.         *
!       * elementType: The type of element which contains the point.   *
!       * elementID:   The entry in the connectivity of this element   *
!       *              which contains the point. The ID is negative if *
!       *              the coordinate is outside the element.          *
!       * uvw:         The parametric coordinates of the point in the  *
!       *              transformed element; this transformation is     *
!       *              such that every element is transformed into a   *
!       *              standard element in parametric space. The u, v  *
!       *              and w coordinates can be used to determine the  *
!       *              actual interpolation weights. If the tree       *
!       *              corresponds to a surface mesh the third entry   *
!       *              of this array will not be filled.               *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it should be            *
!       *        initialized by the calling program, possibly to a     *
!       *        large value. In this way it is possible to handle     *
!       *        periodic problems as efficiently as possible.         *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor, nInterpol
        character(len=*),         intent(in) :: adtID

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                             elementType

        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol

        real(kind=realType), dimension(:), intent(inout) :: dist2
!
!       Local variables.
!
        integer(kind=intType) :: jj, nAlloc
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Determine the index in the array ADTs, which stores the given
        ! ID. As the number of trees stored is limited, a linear search
        ! algorithm is okay.

        nAlloc = ubound(ADTs, 1)
        do jj=1,nAlloc
          if(adtID == ADTs(jj)%adtID) exit
        enddo

        ! Check if the ADT to be searched exists. If not stop.
        ! Note that adtTerminate is not called. The reason is that the
        ! processor ID is not known.

        if(jj > nAlloc) stop "ADT to be searched does not exist."

        ! Initialize the search, i.e. determine the number of
        ! coordinates to be searched in each of the local ADT's.

        call initSearch(nCoor, coor, dist2, jj, .false.)

        ! Perform the actual search.

        call search(nCoor,     coor,      procID,   elementType, &
                    elementID, uvw,       dist2,    jj,          &
                    .false.,   nInterpol, arrDonor, arrInterpol)

        ! Negate the elementID if the coordinate is outside the element,
        ! i.e. if the distance is larger than zero.

        do jj=1,nCoor
          if(dist2(jj) > adtZero) elementID(jj) = -elementID(jj)
        enddo

        end subroutine minDistanceSearch

        !***************************************************************
        !***************************************************************

        subroutine search(nCoor,       coor,      procID,            &
                          elementType, elementID, uvw,               &
                          dist2,       jj,        containmentSearch, &
                          nInterpol,   arrDonor,  arrInterpol)
!
!       ****************************************************************
!       *                                                              *
!       * This routine implements the parallel part of the search      *
!       * algorithm and calls the appropriate local tree searches.     *
!       *                                                              *
!       * Subroutine intent(in) arguments.                             *
!       * --------------------------------                             *
!       * nCoor:             Number of local coordinates for which the *
!       *                    element must be determined.               *
!       * coor:              The coordinates of these points.          *
!       * jj:                The index in the array ADTs of the tree   *
!       *                    to be searched.                           *
!       * containmentSearch: Whether or not a containment search must  *
!       *                    be performed. If not a minimum distance   *
!       *                    search algorithm is used, which is more   *
!       *                    expensive.                                *
!       * nInterpol:         Number of variables to be interpolated.   *
!       * arrDonor:          Array with the donor data; needed to      *
!       *                    obtain the interpolated data.             *
!       *                                                              *
!       * Subroutine intent(inout) arguments.                          *
!       * -----------------------------------                          *
!       * dist2: Minimum distance squared of the coordinates to the    *
!       *        elements of the ADT. On input it contains the         *
!       *        guarenteed distance squared to one of the root        *
!       *        leaves. On output it contains the distance squared to *
!       *        the nearest element of the global tree. It is only    *
!       *        used for a minimum distance search.                   *
!       *                                                              *
!       * Subroutine intent(out) arguments.                            *
!       * ---------------------------------                            *
!       * procID:      The ID of the processor in the group of the ADT *
!       *              where the element containing the point is       *
!       *              stored.                                         *
!       * elementType: The type of element which contains the point or *
!       *              minimizes the distance to the point.            *
!       * elementID:   The entry in the connectivity of this element.  *
!       * uvw:         The parametric coordinates of (the projection   *
!       *              of) the point in the transformed element; this  *
!       *              transformation is such that every element is    *
!       *              transformed into a standard element in          *
!       *              parametric space. The u, v and w coordinates    *
!       *              can be used to determine the actual             *
!       *              interpolation weights.                          *
!       * arrInterpol: Array with the interpolated data.               *
!       *                                                              *
!       ****************************************************************
!
        implicit none
!
!       Subroutine arguments.
!
        integer(kind=intType), intent(in) :: nCoor
        integer(kind=intType), intent(in) :: jj
        integer(kind=intType), intent(in) :: nInterpol

        real(kind=realType), dimension(:,:), intent(in) :: coor
        real(kind=realType), dimension(:),   intent(inout) :: dist2

        real(kind=realType), dimension(:,:), intent(in) :: arrDonor

        integer,                  dimension(:), intent(out) :: procID
        integer(kind=intType), dimension(:), intent(out) :: elementID

        integer(kind=adtElementType), dimension(:), intent(out) :: &
                                                             elementType

        real(kind=realType), dimension(:,:), intent(out) :: uvw
        real(kind=realType), dimension(:,:), intent(out) :: arrInterpol

        logical, intent(in) :: containmentSearch
!
!       Local variables.
!
        integer :: ierr
        integer :: comm, nProcs, myID
        integer :: nProcRecvCur, nProcSendCur, nVarCoor, nVarUVW
        integer :: startProcRecv, procCur, sizeMessage

        integer, dimension(mpi_status_size) :: status

        integer, dimension(:),   allocatable :: procSendCur
        integer, dimension(:),   allocatable :: sendRequest
        integer, dimension(:,:), allocatable :: sendRecvRequest

        integer(kind=intType) :: i, j, k, k1, l, m, ii, mm, nn
        integer(kind=intType) :: nLocalInterpolRound
        integer(kind=intType) :: iStartLocal, iEndLocal, nCoorRecv

        integer(kind=intType), dimension(:), allocatable :: nCoorPerProc
        integer(kind=intType), dimension(:), allocatable :: nCoorFromProc

        integer(kind=intType), dimension(:,:), allocatable :: intRecv
        integer(kind=intType), dimension(:,:), allocatable :: intBuf

        real(kind=realType), dimension(:,:), allocatable :: coorBuf
        real(kind=realType), dimension(:,:), allocatable :: coorRecv
        real(kind=realType), dimension(:,:), allocatable :: uvwRecv
        real(kind=realType), dimension(:,:), allocatable :: uvwBuf

        logical, dimension(:), allocatable :: coorRequested
!
!       ****************************************************************
!       *                                                              *
!       * Begin execution.                                             *
!       *                                                              *
!       ****************************************************************
!
        ! Some abbreviations to make the code more readable.

        comm   = ADTs(jj)%comm
        nProcs = ADTs(jj)%nProcs
        myID   = ADTs(jj)%myID

        ! Determine the number of variables stored in the coordinate
        ! buffers. For a minimum distance search also the distance
        ! squared is stored, i.e. 4 variables instead of 3.

        if( containmentSearch ) then
          nVarCoor = 3
        else
          nVarCoor = 4
        endif

        ! And the size of the uvw buffers. These contain nVarCoor plus
        ! the number of variables to be interpolated.

        nVarUVW = nVarCoor + max(nInterpol,0_intType)

        ! Initialize procID to -1, which indicates failure.

        do i=1,nCoor
          procID(i) = -1
        enddo

        ! Allocate the memory for some help arrays used in the search
        ! algorithm.

        nn = nCoorPerRootLeaf(nProcs)
        allocate(procSendCur(nProcs-1),       sendRequest(nProcs-1),     &
                 nCoorPerProc(0:nProcs-1),    nCoorFromProc(0:nProcs-1), &
                 sendRecvRequest(2,nProcs-1), coorRequested(nn),         &
                 stat=ierr)
        if(ierr /= 0)                     &
          call adtTerminate(jj, "search", &
                            "Memory allocation failure for help arrays.")

        ! Initialize coorRequested to .false. This indicates that
        ! the corresponding entry in coorPerRootLeaf has not been
        ! requested for interpolation.

        do j=1,nn
          coorRequested(j) = .false.
        enddo

        ! Initialize the starting position in the array procRecv to 1.
        ! This variable indicates the starting position in procRecv for
        ! the current round. Also initializes nCoorFromProc to 0.

        startProcRecv = 1

        do i=0,(nProcs-1)
          nCoorFromProc(i) = 0
        enddo

        ! Determine the number of local interpolations per round and
        ! initialize the iStartLocal and iEndLocal, the start and end
        ! indices for the local interpolation of the current round.

        nn = nLocalInterpol/nRounds
        if(nn*nRounds < nLocalInterpol) nn = nn + 1
        nLocalInterpolRound = nn

        iStartLocal = 0
        iEndLocal   = nLocalInterpolRound
!
!       ****************************************************************
!       *                                                              *
!       * Iterative algorithm to determine the elements containing the *
!       * coordinates or the elements which minimize the distance. The *
!       * algorithm consists of a synchronous outer loop over the      *
!       * number of times the inner loop should be executed. This      *
!       * inner loop is asynchronous and performs the actual ADT       *
!       * search. The outer loop is present to avoid that too much     *
!       * data is communicated to a single processor at once such that *
!       * a memory bottleneck occurs.                                  *
!       *                                                              *
!       ****************************************************************
!
        outerLoop: do mm=1,nRounds

          ! Determine the processors I want data from in this round
          ! as well as the number of coordinates from these nodes.

          nProcRecvCur = 0
          nCoorRecv    = 0

          do i=startProcRecv,nProcRecv

            ! Exit the loop if the maximum number of nodes has been
            ! reached.

            if(nCoorRecv == nCoorMax) exit

            ! Update the number of processors from which I will receive
            ! data during this round.

            nProcRecvCur = nProcRecvCur + 1

            ! Determine the amount I can receive from this processor.
            ! If the upper limit is exceeded cut it off.

            j = nCoorProcRecv(i)
            if((nCoorRecv+j) > nCoorMax) j = nCoorMax - nCoorRecv
            nCoorRecv = nCoorRecv + j

            ! Store the amount j in the appropriate place of
            ! nCoorFromProc and determine the number of nodes to be
            ! sent to this processor in the next round (possibly 0).

            nCoorFromProc(procRecv(i)) = j
            nCoorProcRecv(i) = nCoorProcRecv(i) - j

            ! Exit the loop if still some data should be received from
            ! this processor in a next round. The difference between this
            ! exit and the one in the beginning of the do loop is that
            ! here the counter i is not update yet and therefore
            ! startProcRecv will be set correctly for the next round.

            if(nCoorProcRecv(i) > 0) exit
          enddo

          ! Do an all to all communication such that every processor
          ! knows the amount of data it should send to other processors.

          call mpi_alltoall(nCoorFromProc, 1, sumb_integer, &
                            nCoorPerProc,  1, sumb_integer, comm, ierr)

          ! Set the non-zero entries of nCoorFromProc to zero again
          ! for the next round.

          nn = min(i,nProcRecv)
          do j=startProcRecv, nn
            nCoorFromProc(procRecv(j)) = 0
          enddo

          ! Set the starting index for the next round.

          startProcRecv = i

          ! Determine the number of messages I have to send, the 
          ! corresponding processors and the total number of points to
          ! be sent.

          nProcSendCur = 0
          nn = 0

          do i=0,(nProcs-1)
            if(nCoorPerProc(i) > 0) then
              nProcSendCur = nProcSendCur + 1
              procSendCur(nProcSendCur) = i
              nn = nn + nCoorPerProc(i)
            endif
          enddo

          ! Allocate the memory for the send buffer of the coordinates
          ! and possibly the minimum distance squared.

          allocate(coorBuf(nVarCoor,nn), stat=ierr)
          if(ierr /= 0)                     &
            call adtTerminate(jj, "search", &
                              "Memory allocation failure for coorBuf.")

          ! Send the coordinates to the appropriate processors.
          ! Initialize the counter k in the coordinate buffer to 0.

          k = 0
          sendCoorLoop: do i=1,nProcSendCur

            ! Store the processor ID and the starting entry in
            ! coorPerRootLeaf a bit easier and initialize k1 to k.

            procCur = procSendCur(i)
            nn      = mCoorPerRootLeaf(procCur)
            k1      = k

            ! Loop to fill to buffer to this processor.
            ! A distinction is made between a containment search and a
            ! minimum distance search. In the former case a coordinate is
            ! only sent if it has not been interpolated yet; in the
            ! latter case it is sent if the distance is larger than zero.
            ! For a minimum distance search also the current distance
            ! squared is stored.

            test1Containment: if( containmentSearch ) then

              ! Containment search. Store the coordinates of the points
              ! to be searched in coorBuf. Note that the counter j is
              ! updated inside the loop. This is done to be able to send
              ! the maximum number of coordinates possible.

              j = 0
              do
                ! Exit the loop if the maximum number or if the last
                ! coordinate for this processor has been reached.

                if(j  == nCoorPerProc(procCur) .or. &
                   nn == nCoorPerRootLeaf(procCur+1)) exit

                ! Update the counter nn and check if this coordinate
                ! still needs to be interpolated.

                nn = nn + 1
                l  = coorPerRootLeaf(nn)

                if(procID(l) == -1) then

                  ! Coordinate needs to be interpolated. Update the
                  ! counters j and k and copy the coordinate in the
                  ! send buffer.

                  j = j + 1
                  k = k + 1
                  coorBuf(1,k) = coor(1,l)
                  coorBuf(2,k) = coor(2,l)
                  coorBuf(3,k) = coor(3,l)

                  ! Set the entry in coorRequested to .true. to indicate
                  ! that a request was sent.

                  coorRequested(nn) = .true.

                endif
              enddo

            else test1Containment

              ! Minimum distance search. Even if an earlier minimum
              ! distance was found, it is possible that an even smaller
              ! distance can be found. Therefore only points with a
              ! minimum distance squared of zero will not be sent. Both
              ! the coordinates and the current minimum distance squared
              ! is sent. Again the counter j is updated inside the loop.

              j = 0
              do
                ! Exit the loop if the maximum number or if the last
                ! coordinate for this processor has been reached.

                if(j  == nCoorPerProc(procCur) .or. &
                   nn == nCoorPerRootLeaf(procCur+1)) exit

                ! Update the counter nn and check if this coordinate still
                ! needs to be interpolated.

                nn = nn + 1
                l  = coorPerRootLeaf(nn)

                if(dist2(l) > adtZero) then

                  ! Coordinate needs to be interpolated. Update the
                  ! counters j and k and store the coordinates and the
                  ! distance squared in the send buffer.

                  j = j + 1
                  k = k + 1
                  coorBuf(1,k) = coor(1,l)
                  coorBuf(2,k) = coor(2,l)
                  coorBuf(3,k) = coor(3,l)
                  coorBuf(4,k) = dist2(l)

                  ! Set the entry in coorRequested to .true. to indicate
                  ! that a request was sent.

                  coorRequested(nn) = .true.

                endif
              enddo

            endif test1Containment

            ! Determine the size of the message and send it.
            ! Use nonblocking sends to avoid deadlock.

            sizeMessage = nVarCoor*(k-k1)
            call mpi_isend(coorBuf(1,k1+1), sizeMessage, sumb_real, &
                           procCur,         procCur,     comm,     &
                           sendRequest(i),  ierr)

          enddo sendCoorLoop
!
!         **************************************************************
!         *                                                            *
!         * Perform the local interpolations. This is done here to     *
!         * have an overlap between communication and computation.     *
!         *                                                            *
!         **************************************************************
!
          ! Determine the local number to be interpolated in this round.

          nn = iEndLocal - iStartLocal

          ! Allocate the memory for the coorRecv, the integers for
          ! storing the processor ID, element type and element ID, and
          ! uvwRecv. This memory is also used for the local
          ! interpolation, which explains the max test for intRecv and
          ! uvwRecv. The coordinates are different, because they will be
          ! released after each received message.

          i = max(nn,nCoorRecv)
          allocate(coorRecv(nVarCoor,nn), intRecv(3,i), &
                   uvwRecv(nVarUVW,i),    stat=ierr)
          if(ierr /= 0)                     &
            call adtTerminate(jj, "search", &
                              "Memory allocation failure for &
                              &recv arrays")

          ! Initialize the counters i and j, which are used to fill the
          ! buffer coorRecv.

          j = 0
          i = mCoorPerRootLeaf(myID)

          ! Make a distinction between a containment and a minimum
          ! distance search.

          test2Containment: if( containmentSearch ) then

            ! Containment search. Copy the local coordinates to be
            ! searched in coorRecv. As it is possible that in an earlier
            ! round the coordinate has already been found, only take
            ! coordinates whose element has not been found yet.

            do
              ! Exit the loop if the maximum number or if the last
              ! coordinate of the local interpolation has been reached.

              if(j == nn .or. i == nCoorPerRootLeaf(myID+1)) exit

              ! Update the counter i and check if the corresponding
              ! coordinate still needs to be interpolated. If so, store it
              ! in coorRecv and set the corresponding entry in
              ! coorRequested to .true.

              i = i + 1
              l = coorPerRootLeaf(i)

              if(procID(l) == -1) then
                j = j + 1

                coorRecv(1,j) = coor(1,l)
                coorRecv(2,j) = coor(2,l)
                coorRecv(3,j) = coor(3,l)

                coorRequested(i) = .true.
              endif
            enddo

            ! Perform the local interpolations. Note that j contains the
            ! actual number of coordinates to be searched.

            nn = j
            call containmentTreeSearch(jj,        coorRecv, intRecv, &
                                       uvwRecv,   arrDonor, nn,      &
                                       nInterpol)

            ! Store the interpolation data at the correct location in
            ! the corresponding arrays.

            i = mCoorPerRootLeaf(myID)
            do j=1,nn

              ! Determine the coordinate entry, which corresponds to the
              ! counter j. Remember that in the calls to search routines
              ! only nodes are given which were not interpolated.

              do
                i = i + 1
                if( coorRequested(i) ) exit
              enddo
              l = coorPerRootLeaf(i)

              ! Copy the data if an actual element was found.

              if(intRecv(1,j) >= 0) then
                procID(l)      = intRecv(1,j)
                elementType(l) = intRecv(2,j)
                elementID(l)   = intRecv(3,j)
                uvw(1,l)       = uvwRecv(1,j)
                uvw(2,l)       = uvwRecv(2,j)
                uvw(3,l)       = uvwRecv(3,j)

                do m=1,nInterpol
                  arrInterpol(m,l) = uvwRecv(m+nVarCoor,j)
                enddo
              endif

            enddo

          else test2Containment

            ! Minimum distance search. Store the local coordinates and
            ! distance in coorRecv for the coordinates with non-zero
            ! distance.

            do
              ! Exit the loop if the maximum number or if the last
              ! coordinate of the local interpolation has been reached.

              if(j == nn .or. i == nCoorPerRootLeaf(myID+1)) exit

              ! Update the counter i and check if the corresponding
              ! coordinate still needs to be interpolated. If so, store
              ! its coordinates and distance squared in coorRecv and set
              ! the corresponding entry in coorRequested to .true.

              i = i + 1
              l = coorPerRootLeaf(i)

              if(dist2(l) > adtZero) then
                j = j + 1

                coorRecv(1,j) = coor(1,l)
                coorRecv(2,j) = coor(2,l)
                coorRecv(3,j) = coor(3,l)
                coorRecv(4,j) = dist2(l)

                coorRequested(i) = .true.
              endif
            enddo

            ! Perform the local interpolations. Note that j contains the
            ! actual number of coordinates to be searched.

            nn = j
            call minDistanceTreeSearch(jj,        coorRecv, intRecv, &
                                       uvwRecv,   arrDonor, nn,      &
                                       nInterpol)

            ! Store the interpolation data at the correct location in
            ! the corresponding arrays.

            i = mCoorPerRootLeaf(myID)

            do j=1,nn

              ! Determine the next coordinate entry, whose request was
              ! sent to be interpolated. Remember that it is possible
              ! that some nodes were not sent to be interpolated (if
              ! their distance squared is zero already).

              do
                i = i + 1
                if( coorRequested(i) ) exit
              enddo
              l = coorPerRootLeaf(i)

              ! Copy the data if an actual element was found and if
              ! the corresponding distance squared is less than the
              ! currently stored value.

              if(intRecv(1,j) >= 0 .and. uvwRecv(4,j) < dist2(l)) then

                procID(l)      = intRecv(1,j)
                elementType(l) = intRecv(2,j)
                elementID(l)   = intRecv(3,j)
                uvw(1,l)       = uvwRecv(1,j)
                uvw(2,l)       = uvwRecv(2,j)
                uvw(3,l)       = uvwRecv(3,j)
                dist2(l)       = uvwRecv(4,j)

                do m=1,nInterpol
                  arrInterpol(m,l) = uvwRecv(m+nVarCoor,j)
                enddo
              endif

            enddo

          endif test2Containment

          ! The buffer coorRecv is not needed anymore.
          ! Release the memory.

          deallocate(coorRecv, stat=ierr)
          if(ierr /= 0)                     &
            call adtTerminate(jj, "search", &
                              "Deallocation failure for coorRecv.")

          ! Set iStartLocal and iEndLocal for the next round.
          ! Also update mCoorPerRootLeaf(myID).

          iStartLocal = iEndLocal
          iEndLocal   = iEndLocal + nLocalInterpolRound
          iEndLocal   = max(iEndLocal, nLocalInterpol)

          mCoorPerRootLeaf(myID) = i
!
!         **************************************************************
!         *                                                            *
!         * Perform the interpolations from the other processors.      *
!         * Their coordinates are received in an arbitrary sequence    *
!         * using blocking receives and the interpolated data is sent  *
!         * back using nonblocking sends.                              *
!         *                                                            *
!         **************************************************************
!
          ! Loop over the number of messages I will receive. The counter
          ! ii contains the current starting position for the buffers
          ! to be sent back to the requesting processors.

          ii = 1
          recvSendLoop: do i=1,nProcRecvCur

            ! Block until a message arrives and find the source and size
            ! of the message.

            call mpi_probe(mpi_any_source, myID, comm, status, ierr)

            procCur = status(mpi_source)
            call mpi_get_count(status, sumb_real, sizeMessage, ierr)

            ! Check in debug mode that the message is of correct size.

            if( debug ) then
              if(sizeMessage == mpi_undefined .or. &
                 mod(sizeMessage,nVarCoor) /= 0)   &
                call adtTerminate(jj, "search",    &
                                  "Unexpected size of message")
            endif

            ! Allocate the memory for the coordinates to be received and
            ! receive them using a blocking receive; the message has
            ! already arrived.

            nn = sizeMessage/nVarCoor
            allocate(coorRecv(nVarCoor,nn), stat=ierr)
            if(ierr /= 0)                     &
              call adtTerminate(jj, "search", &
                                "Memory allocation failure for &
                                &coorRecv.")

            call mpi_recv(coorRecv, sizeMessage, sumb_real, procCur, &
                          myID,     comm,        status,   ierr)

            ! Search the corresponding elements in the local tree and
            ! release coorRecv afterwards.

            if( containmentSearch ) then
              call containmentTreeSearch(jj,             coorRecv,       &
                                         intRecv(:,ii:), uvwRecv(:,ii:), &
                                         arrDonor,       nn,             &
                                         nInterpol)
            else
              call minDistanceTreeSearch(jj,             coorRecv,       &
                                         intRecv(:,ii:), uvwRecv(:,ii:), &
                                         arrDonor,       nn,             &
                                         nInterpol)
            endif

            deallocate(coorRecv, stat=ierr)
            if(ierr /= 0)                     &
              call adtTerminate(jj, "search", &
                                "Deallocation failure for coorRecv.")

            ! Send the integer and the floating point information back to
            ! the requesting processor.

            sizeMessage = 3*nn
            call mpi_isend(intRecv(1,ii),        sizeMessage, sumb_integer, &
                           procCur,              procCur+1,   comm,        &
                           sendRecvRequest(1,i), ierr)

            sizeMessage = nVarUVW*nn
            call mpi_isend(uvwRecv(1,ii),        sizeMessage, sumb_real, &
                           procCur,              procCur+2,   comm,     &
                           sendRecvRequest(2,i), ierr)

            ! Update the counter ii for the next message.

            ii = ii + nn

          enddo recvSendLoop

          ! Complete the nonblocking coordinate sends and release the
          ! memory of coorBuf.

          do i=1,nProcSendCur
            call mpi_waitany(nProcSendCur, sendRequest, sizeMessage, &
                             status,       ierr)
          enddo

          deallocate(coorBuf, stat=ierr)
          if(ierr /= 0)                     &
            call adtTerminate(jj, "search", &
                              "Deallocation failure for coorBuf.")

          ! Loop over the number of processors to which I sent requests.
          ! Now it is time to receive the information they interpolated.
          ! The sequence of receiving messages is arbitrary.

          recvLoop: do ii=1,nProcSendCur

            ! Block until an integer message arrives. These messages have
            ! tags of myID+1. Also determine the sending processor and
            ! the size of the message.

            call mpi_probe(mpi_any_source, myID+1, comm, status, ierr)

            procCur = status(mpi_source)
            call mpi_get_count(status, sumb_integer, sizeMessage, ierr)

            ! Check in debug mode that the message is of correct size.

            if( debug ) then
              if(sizeMessage == mpi_undefined .or. &
                 mod(sizeMessage,3) /= 0)          &
                call adtTerminate(jj, "search",    &
                                  "Unexpected size of message")
            endif

            ! Allocate the memory for the integer and uvw buffers, such
            ! that the messages can be received.

            nn = sizeMessage/3
            allocate(intBuf(3,nn), uvwBuf(nVarUVW,nn), stat=ierr)
            if(ierr /= 0)                     &
              call adtTerminate(jj, "search", &
                                "Memory allocation failure for intBuf &
                                &and uvwBuf.")

            call mpi_recv(intBuf, sizeMessage, sumb_integer, procCur, &
                          myID+1, comm,        status,      ierr)

            sizeMessage = nVarUVW*nn
            call mpi_recv(uvwBuf, sizeMessage, sumb_real, procCur, &
                          myID+2, comm,        status,   ierr)

            ! Store the interpolation data at the correct location in the
            ! corresponding arrays. A distinction must be made between
            ! containment search and minimum distance search, because for
            ! the latter it is possible that a better candidate is
            ! already stored.

            i = mCoorPerRootLeaf(procCur)

            test3Containment: if( containmentSearch ) then

              ! Containment search. Loop over the number of points
              ! requested on the other processor.

              do j=1,nn

                ! Determine the next coordinate entry, whose request was
                ! sent to be interpolated. Remember that only nodes are
                ! sent which were not interpolated.

                do
                  i = i + 1
                  if( coorRequested(i) ) exit
                enddo

                ! Copy the data if an actual element was found.

                if(intBuf(1,j) >= 0) then
                  l = coorPerRootLeaf(i)

                  procID(l)      = intBuf(1,j)
                  elementType(l) = intBuf(2,j)
                  elementID(l)   = intBuf(3,j)
                  uvw(1,l)       = uvwBuf(1,j)
                  uvw(2,l)       = uvwBuf(2,j)
                  uvw(3,l)       = uvwBuf(3,j)

                  do m=1,nInterpol
                    arrInterpol(m,l) = uvwBuf(m+nVarCoor,j)
                  enddo
                endif

              enddo

            else test3Containment

              ! Minimum distance search. Loop over the number of points
              ! requested on the other processor.

              do j=1,nn

                ! Determine the next coordinate entry, whose request was
                ! sent to be interpolated. Remember that it is possible
                ! that some nodes were not sent to be interpolated (if
                ! their distance squared is zero already).

                do
                  i = i + 1
                  if( coorRequested(i) ) exit
                enddo

                ! Copy the data if an actual element was found and if
                ! the corresponding distance squared is less than the
                ! currently stored value.

                l = coorPerRootLeaf(i)
                if(intBuf(1,j) >= 0 .and. uvwBuf(4,j) < dist2(l)) then

                  procID(l)      = intBuf(1,j)
                  elementType(l) = intBuf(2,j)
                  elementID(l)   = intBuf(3,j)
                  uvw(1,l)       = uvwBuf(1,j)
                  uvw(2,l)       = uvwBuf(2,j)
                  uvw(3,l)       = uvwBuf(3,j)
                  dist2(l)       = uvwBuf(4,j)

                  do m=1,nInterpol
                    arrInterpol(m,l) = uvwBuf(m+nVarCoor,j)
                  enddo
                endif

              enddo

            endif test3Containment

            ! Update mCoorPerRootLeaf(procCur) for the next round and
            ! release the memory of intBuf and uvwBuf.

            mCoorPerRootLeaf(procCur) = i

            deallocate(intBuf, uvwBuf, stat=ierr)
            if(ierr /= 0)                     &
              call adtTerminate(jj, "search", &
                                "Deallocation failure for intBuf &
                                &and uvwBuf.")
          enddo recvLoop

          ! Complete the nonblocking sends of the interpolated data.

          nProcRecvCur = 2*nProcRecvCur
          do i=1,nProcRecvCur
            call mpi_waitany(nProcRecvCur, sendRecvRequest, sizeMessage, &
                             status,       ierr)
          enddo

          ! Release the memory of the buffers used in the nonblocking
          ! sends of the interpolated data.

          deallocate(intRecv, uvwRecv, stat=ierr)
          if(ierr /= 0)                     &
            call adtTerminate(jj, "search", &
                              "Deallocation failure for intRecv and &
                              &uvwRecv")

          ! Synchronize the processors, because wild cards have been
          ! used in the communication.

          call mpi_barrier(comm, ierr)

        enddo outerLoop

        ! Release the memory of the help arrays allocated in this
        ! routine.

        deallocate(procSendCur,   sendRequest,     nCoorPerProc,  &
                   nCoorFromProc, sendRecvRequest, coorRequested, &
                   stat=ierr)
        if(ierr /= 0)                     &
          call adtTerminate(jj, "search", &
                            "Deallocation failure for help arrays.")

        ! Release the memory of the help arrays stored in the module
        ! adtData.

        deallocate(procRecv,         nCoorProcRecv,   nCoorPerRootLeaf, &
                   mCoorPerRootLeaf, coorPerRootLeaf, stat=ierr)
        if(ierr /= 0)                     &
          call adtTerminate(jj, "search", &
                            "Deallocation failure for help arrays &
                            &stored in the module adtData.")

        end subroutine search

      end module adtSearch
