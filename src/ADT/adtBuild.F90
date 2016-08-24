module adtBuild
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Module which contains all the subroutines for the building of  *
  !     * an ADT, both surface and volume.                               *
  !     *                                                                *
  !     ******************************************************************
  !
  use adtUtils
  implicit none

  !=================================================================

contains

  !===============================================================

  subroutine buildADT(ADT)
    !
    !       ****************************************************************
    !       *                                                              *
    !       * This routine builds the 6 dimensional ADT for the given      *
    !       * ADT. When this routine is called it is assumed that the      *
    !       * bounding boxes of the grid have already been computed; the   *
    !       * ADT for these bounding boxes is built here.                  *
    !       *                                                              *
    !       * Subroutine intent(inout) arguments.                          *
    !       * --------------------------------                             *
    !       * ADT: adt derived type to build                               *
    !       *                                                              *
    !       ****************************************************************
    !
    implicit none
    !
    !       Subroutine arguments.
    !
    type(adtType), intent(inout) :: ADT
    !
    !       Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, ii, kk, ll, mm, nn, nfl, nfr
    integer(kind=intType) :: nLeaves, nBBoxes, splitDir
    integer(kind=intType) :: nLeavesToDivide, nLeavesToDivideNew
    integer(kind=intType) :: nLeavesTot

    integer(kind=intType), dimension(:), pointer :: BB_IDs
    integer(kind=intType), dimension(:), pointer :: BB_IDsNew
    integer(kind=intType), dimension(:), pointer :: nBB_IDs
    integer(kind=intType), dimension(:), pointer :: nBB_IDsNew
    integer(kind=intType), dimension(:), pointer :: curLeaf
    integer(kind=intType), dimension(:), pointer :: curLeafNew
    integer(kind=intType), dimension(:), pointer :: tmpIntPointer

    integer(kind=intType), dimension(0:ADT%nProcs-1) :: tmpArr

    real(kind=realType), dimension(:,:), pointer :: xBBox

    real(kind=realType), dimension(3,2) :: rootLeafBBox
    real(kind=realType), dimension(3,2,0:ADT%nProcs-1) :: &
         rootLeavesBBox

    type(adtLeafType), dimension(:), pointer :: ADTree

    ! Initialize nStack and allocate the corresponding array stack.
    ! These are used in the qsort routine for the bounding boxes.
    ! As this routine is called quite often it is more efficient to
    ! have the stack array available rather than allocate it over
    ! and over again.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                       &
         call adtTerminate(ADT, "buildADT", &
         "Memory allocation failure for stack.")

    ! Determine the number of leaves of the adt. It can be proved
    ! that nLeaves equals nBBoxes - 1 for an optimally balanced
    ! tree. Take the exceptional case of nBBoxes == 0 and
    ! nBBoxes == 1 into account.

    nBBoxes = ADT%nBBoxes
    nLeaves = nBBoxes - 1
    if(nBBoxes <= 1) nLeaves = nLeaves + 1

    ADT%nLeaves = nLeaves

    ! Allocate the memory for the adt.

    allocate(ADT%ADTree(nLeaves), stat=ierr)
    if(ierr /= 0)                       &
         call adtTerminate(ADT, "buildADT", &
         "Memory allocation failure for ADTree.")

    ! Set some pointers to make the code more readable.

    xBBox  => ADT%xBBox
    ADTree => ADT%ADTree

    ! Allocate the memory for the arrays which control the
    ! subdivision of the leaves.

    nn = (nBBoxes+1)/2
    nn = max(nn, 1_intType)

    allocate(BB_IDs(nBBoxes), BB_IDsNew(nBBoxes), &
         nBB_IDs(0:nn),   nBB_IDsNew(0:nn),   &
         curLeaf(nn),     curLeafNew(nn),     stat=ierr)
    if(ierr /= 0)                       &
         call adtTerminate(ADT, "buildADT", &
         "Memory allocation failure for the arrays &
         &used in the subdivision.")

    ! Initialize the arrays BB_IDs, nBB_IDs and curLeaf, such that
    ! all bounding boxes belong to the root leaf. Also set the
    ! counters nLeavesToDivide and nLeavesTot, depending on the
    ! situation

    nBB_IDs(0) = 0; nBB_IDs(1) = nBBoxes
    curLeaf(1) = 1

    do i=1,nBBoxes
       BB_IDs(i) = i
    enddo

    nLeavesToDivide = min(nLeaves, 1_intType)
    nLeavesTot      = nLeavesToDivide

    ! Initialize splitDir to 0, such that the first time it will
    ! split in direction 1.

    splitDir = 0

    ! Loop to subdivide the leaves. The division is such that the
    ! adt is optimally balanced.

    leafDivision: do

       ! Criterion to exit the loop.

       if(nLeavesToDivide == 0) exit

       ! Initializations for the next round of subdivisions and
       ! increment splitDir.

       nLeavesToDivideNew = 0
       nBB_IDsNew(0) = 0

       splitdir = splitDir + 1
       if(splitDir > 6) splitDir = 1

       ! Loop over the current number of leaves to be divided.

       currentLeavesLoop: do i=1,nLeavesToDivide

          ! Store the number of bounding boxes present in the leaf
          ! in nn, the current leaf number in mm and i-1 in ii.

          ii = i-1
          nn = nBB_IDs(i) - nBB_IDs(ii)
          mm = curLeaf(i)

          ! Determine the bounding box coordinates of this leaf.

          ll = BB_IDs(nBB_IDs(ii)+1)
          ADTree(mm)%xMin(1) = xBBox(1,ll)
          ADTree(mm)%xMin(2) = xBBox(2,ll)
          ADTree(mm)%xMin(3) = xBBox(3,ll)
          ADTree(mm)%xMin(4) = xBBox(4,ll)
          ADTree(mm)%xMin(5) = xBBox(5,ll)
          ADTree(mm)%xMin(6) = xBBox(6,ll)

          ADTree(mm)%xMax(1) = xBBox(1,ll)
          ADTree(mm)%xMax(2) = xBBox(2,ll)
          ADTree(mm)%xMax(3) = xBBox(3,ll)
          ADTree(mm)%xMax(4) = xBBox(4,ll)
          ADTree(mm)%xMax(5) = xBBox(5,ll)
          ADTree(mm)%xMax(6) = xBBox(6,ll)

          do j=(nBB_IDs(ii)+2),nBB_IDs(i)
             ll = BB_IDs(j)

             ADTree(mm)%xMin(1) = min(ADTree(mm)%xMin(1), xBBox(1,ll))
             ADTree(mm)%xMin(2) = min(ADTree(mm)%xMin(2), xBBox(2,ll))
             ADTree(mm)%xMin(3) = min(ADTree(mm)%xMin(3), xBBox(3,ll))
             ADTree(mm)%xMin(4) = min(ADTree(mm)%xMin(4), xBBox(4,ll))
             ADTree(mm)%xMin(5) = min(ADTree(mm)%xMin(5), xBBox(5,ll))
             ADTree(mm)%xMin(6) = min(ADTree(mm)%xMin(6), xBBox(6,ll))

             ADTree(mm)%xMax(1) = max(ADTree(mm)%xMax(1), xBBox(1,ll))
             ADTree(mm)%xMax(2) = max(ADTree(mm)%xMax(2), xBBox(2,ll))
             ADTree(mm)%xMax(3) = max(ADTree(mm)%xMax(3), xBBox(3,ll))
             ADTree(mm)%xMax(4) = max(ADTree(mm)%xMax(4), xBBox(4,ll))
             ADTree(mm)%xMax(5) = max(ADTree(mm)%xMax(5), xBBox(5,ll))
             ADTree(mm)%xMax(6) = max(ADTree(mm)%xMax(6), xBBox(6,ll))
          enddo

          ! Determine the situation of the leaf. It is either a
          ! terminal leaf or a leaf that must be subdivided.

          terminalTest: if(nn <= 2) then

             ! Terminal leaf. Store the ID's of the bounding boxes with
             ! negative numbers in children.

             ADTree(mm)%children(1) = -BB_IDs(nBB_IDs(ii)+1)
             ADTree(mm)%children(2) = -BB_IDs(nBB_IDs(i))

          else terminalTest

             ! Leaf must be divided. Sort the bounding boxes of the
             ! current leaf in increasing order; the sorting is based
             ! on the coordinate in the split direction.

             call qsortBBoxes(BB_IDs(nBB_IDs(ii)+1:), nn, ADT, splitDir)

             ! Determine the number of bounding boxes in the left leaf.
             ! This number is at least 2. The actual number stored in
             ! kk is this number plus an offset. Also initialize the
             ! counter nfl, which is used to store the bounding boxes
             ! in the arrays for the new round.

             kk  = (nn+1)/2 + nBB_IDs(ii)
             nfl = nBB_IDsNew(nLeavesToDivideNew)

             ! Copy the ID's of the left bounding boxes into BB_IDsNew.
             ! Also update nLeavesToDivideNew and the corresponding
             ! entry in nBB_IDsNew.

             do k=(nBB_IDs(ii)+1),kk
                nfl = nfl + 1
                BB_IDsNew(nfl) = BB_IDs(k)
             enddo

             nLeavesToDivideNew = nLeavesToDivideNew + 1
             nBB_IDsNew(nLeavesToDivideNew) = nfl

             ! Update the total number of leaves and store this number
             ! in child 1 of the current leaf and in the current leaves
             ! for the next round.

             nLeavesTot = nLeavesTot + 1
             ADTree(mm)%children(1) = nLeavesTot
             curLeafNew(nLeavesToDivideNew) = nLeavesTot

             ! The right leaf will only be created if it has more than
             ! one bounding box in it, i.e. if the original leaf has
             ! more than three bounding boxes. If the new leaf only has
             ! one bounding box in it, it is not created; instead the
             ! bounding box is stored in the current leaf.

             if(nn == 3) then

                ! Only three bounding boxes present in the current leaf.
                ! The right leaf is not created and the last bounding
                ! box is stored as the second child of the current leaf.

                ADTree(mm)%children(2) = -BB_IDs(nBB_IDs(i))

             else

                ! More than 3 bounding boxes are present and thus the
                ! right leaf is created. Copy the ID's from BB_IDs into
                ! BB_IDsNew and update the counters for the new round.

                nfr = nBB_IDsNew(nLeavesToDivideNew)
                do k=(kk+1),nBB_IDs(i)
                   nfr = nfr + 1
                   BB_IDsNew(nfr) = BB_IDs(k)
                enddo

                nLeavesToDivideNew = nLeavesToDivideNew + 1
                nBB_IDsNew(nLeavesToDivideNew) = nfr

                ! Update the total number of leaves and store this number
                ! in child 2 of the current leaf and in the current
                ! leaves for the next round.

                nLeavesTot = nLeavesTot + 1
                ADTree(mm)%children(2) = nLeavesTot
                curLeafNew(nLeavesToDivideNew) = nLeavesTot

             endif

          endif terminalTest

       enddo currentLeavesLoop

       ! Swap the pointers for the next round.

       nLeavesToDivide = nLeavesToDivideNew

       tmpIntPointer => BB_IDs
       BB_IDs        => BB_IDsNew
       BB_IDsNew     => tmpIntPointer

       tmpIntPointer => nBB_IDs
       nBB_IDs       => nBB_IDsNew
       nBB_IDsNew    => tmpIntPointer

       tmpIntPointer => curLeaf
       curLeaf       => curLeafNew
       curLeafNew    => tmpIntPointer

    enddo leafDivision

    ! Deallocate the arrays used to build the local tree.

    deallocate(stack, BB_IDs, BB_IDsNew, nBB_IDs, nBB_IDsNew, &
         curLeaf, curLeafNew, stat=ierr)
    if(ierr /= 0)                       &
         call adtTerminate(ADT, "buildADT", &
         "Deallocation failure for the local arrays.")
    !
    !       ****************************************************************
    !       *                                                              *
    !       * Local tree has been built. Now determine the global          *
    !       * information for this tree.                                   *
    !       *                                                              *
    !       ****************************************************************
    !
    ! Determine the number and processor ID's of the non-empty local
    ! trees by gathering the data from the participating processors.

    ii = 1
    if(nBBoxes == 0) ii = 0

    call mpi_allgather(ii, 1, sumb_integer, tmpArr, 1, sumb_integer, &
         ADT%comm, ierr)

    ii = 0
    do i=0,(ADT%nProcs-1)
       ii = ii + tmpArr(i)
    enddo

    ADT%nRootLeaves = ii

    ! Allocate the memory for both the processor ID's and the
    ! 3D bounding box of the root leaf.

    allocate(ADT%rootLeavesProcs(ii), &
         ADT%rootBBoxes(3,2,ii), stat=ierr)
    if(ierr /= 0)                       &
         call adtTerminate(ADT, "buildADT", &
         "Memory allocation failure for &
         &rootLeavesProcs and rootBBoxes.")

    ! Determine the processor ID's of the non-empty trees.

    ii = 0
    ADT%myEntryInRootProcs = 0

    do i=0,(ADT%nProcs-1)
       if(tmpArr(i) > 0) then
          ii = ii + 1
          ADT%rootLeavesProcs(ii) = i
          if(ADT%myID == i) ADT%myEntryInRootProcs = ii
       endif
    enddo

    ! Determine the local 3D bounding box of the root leaf.
    ! If no local tree is present, just set it to zero to avoid
    ! problems. The values do not matter.

    if(nBBoxes == 0) then
       rootLeafBBox = zero
    else
       rootLeafBBox(1,1) = ADTree(1)%xMin(1)
       rootLeafBBox(2,1) = ADTree(1)%xMin(2)
       rootLeafBBox(3,1) = ADTree(1)%xMin(3)

       rootLeafBBox(1,2) = ADTree(1)%xMax(4)
       rootLeafBBox(2,2) = ADTree(1)%xMax(5)
       rootLeafBBox(3,2) = ADTree(1)%xMax(6)
    endif

    ! Gather the data of the root leaves.

    call mpi_allgather(rootLeafBBox, 6, sumb_real, rootLeavesBBox, &
         6, sumb_real, ADT%comm, ierr)

    ! Store the 3D root bounding boxes of the non-empty trees in
    ! the data structure for the current ADT.

    ii = 0
    do i=0,(ADT%nProcs-1)
       if(tmpArr(i) > 0) then
          ii = ii + 1

          ADT%rootBBoxes(1,1,ii) = rootLeavesBBox(1,1,i)
          ADT%rootBBoxes(2,1,ii) = rootLeavesBBox(2,1,i)
          ADT%rootBBoxes(3,1,ii) = rootLeavesBBox(3,1,i)

          ADT%rootBBoxes(1,2,ii) = rootLeavesBBox(1,2,i)
          ADT%rootBBoxes(2,2,ii) = rootLeavesBBox(2,2,i)
          ADT%rootBBoxes(3,2,ii) = rootLeavesBBox(3,2,i)
       endif
    enddo

  end subroutine buildADT

  !***************************************************************
  !***************************************************************

  subroutine buildSurfaceADT(nTria,  nQuads,   nNodes,    &
       coor,   triaConn, quadsConn, &
       BBox,   useBBox,  comm,      &
       adtID)
    !
    !       ****************************************************************
    !       *                                                              *
    !       * This routine builds the 6 dimensional ADT, which stores the  *
    !       * given surface grid. The memory intensive part of these       *
    !       * arguments, the arrays with the coordinates and               *
    !       * connectivities, are not copied. Instead pointers are set to  *
    !       * these arrays. It is therefore the responsibility of the user *
    !       * not to deallocate this memory before all the searches have   *
    !       * been performed.                                              *
    !       *                                                              *
    !       * Subroutine intent(in) arguments.                             *
    !       * --------------------------------                             *
    !       * nNodes:    Number of local nodes in the given grid.          *
    !       * nTria:     Idem for the triangles.                           *
    !       * nQuads:    Idem for the quadrilaterals.                      *
    !       * BBox(3,2): The possible bounding box. Only elements within   *
    !       *            this box will be stored in the ADT.               *
    !       * useBBox:   Whether or not to use the bounding box.           *
    !       * comm:      MPI-communicator for the global ADT.              *
    !       * adtID:     The ID of the ADT.                                *
    !       *                                                              *
    !       * Subroutine intent(in), target arguments.                     *
    !       * ----------------------------------------                     *
    !       * coor(3,nNodes):      Nodal coordinates of the local grid.    *
    !       * triaConn(3,nTria):   Local connectivity of the triangles.    *
    !       * quadsConn(4,nQuads): Idem for the quadrilaterals.            *
    !       *                                                              *
    !       ****************************************************************
    !
    implicit none
    !
    !       Subroutine arguments.
    !
    integer, intent(in)          :: comm
    character(len=*), intent(in) :: adtID

    integer(kind=intType), intent(in) :: nTria
    integer(kind=intType), intent(in) :: nQuads
    integer(kind=intType), intent(in) :: nNodes

    logical, intent(in) :: useBBox

    integer(kind=intType), dimension(:,:), intent(in), &
         target :: triaConn
    integer(kind=intType), dimension(:,:), intent(in), &
         target :: quadsConn

    real(kind=realType), dimension(3,2), intent(in) :: BBox

    real(kind=realType), dimension(:,:), intent(in), &
         target :: coor
    !
    !       Local variables.
    !
    integer :: ierr, ll, nNPE, jj

    integer(kind=adtElementType) :: elType

    integer(kind=intType) :: i, j, ii, mm, nn, nElem

    integer(kind=intType), dimension(:,:), pointer :: conn

    real(kind=realType), dimension(3) :: xMin, xMax

    logical, dimension(:), allocatable :: elementWithinBBox

    type(adtType), pointer :: ADT

    ! ! Allocate or reallocate the memory for ADTs. This depends
    ! ! whether or not this is the first ADT to be built.

    if( allocated(ADTs) ) then
       call reallocateADTs(adtID, jj)
    else
       call allocateADTs
       jj = 1
    endif

    ! The ADT we are currently working with
    ADT => ADTs(jj)

    ! Make sure the ADT is active and store the ID of this ADT.

    ADT%isActive = .true.
    ADT%adtID    = adtID

    ! Copy the communicator and determine the number of processors
    ! and my processor ID in this group.

    ADT%comm = comm
    call mpi_comm_rank(comm, ADT%myID,   ierr)
    call mpi_comm_size(comm, ADT%nProcs, ierr)

    ! Set the ADT type, which is a surface ADT.

    ADT%adtType = adtSurfaceADT

    ! Copy the number of nodes and surface elements and set the
    ! number of volume elements to 0; only a surface grid has been
    ! given.

    ADT%nNodes = nNodes
    ADT%nTria  = nTria
    ADT%nQuads = nQuads

    ADT%nTetra  = 0
    ADT%nPyra   = 0
    ADT%nPrisms = 0
    ADT%nHexa   = 0

    ! Set the pointers for the coordinates and the
    ! surface connectivities.

    ADT%coor      => coor
    ADT%triaConn  => triaConn
    ADT%quadsConn => quadsConn

    ! Determine the number of elements to be stored in the ADT.
    ! This depends whether or not the global bounding box should be
    ! used when building the ADT.

    testBBox: if( useBBox ) then

       ! Global bounding box is used. Allocate the memory for the
       ! logical elementWithinBBox.

       nn = nTria + nQuads
       allocate(elementWithinBBox(nn), stat=ierr)
       if(ierr /= 0)                              &
            call adtTerminate(ADT, "buildSurfaceADT", &
            "Memory allocation failure for &
            &elementWithinBBox.")

       ! Loop over the number of element types.

       ii = 0
       elementLoop1: do ll=1,2

          ! Set the correct pointers for this element.

          call setSurfacePointers(ll)

          ! Loop over the elements and determine the bounding box of
          ! each element.

          do i=1,nElem
             ii = ii + 1

             mm = conn(1,i)
             xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
             xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
             xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

             do j=2,nNPE
                mm = conn(j,i)

                xMin(1) = min(xMin(1),coor(1,mm))
                xMin(2) = min(xMin(2),coor(2,mm))
                xMin(3) = min(xMin(3),coor(3,mm))

                xMax(1) = max(xMax(1),coor(1,mm))
                xMax(2) = max(xMax(2),coor(2,mm))
                xMax(3) = max(xMax(3),coor(3,mm))
             enddo

             ! Check if the bounding box is (partially) inside the
             ! global bounding box. If so, set elementWithinBBox
             ! to .true.; otherwise set it to .false.

             if(xMax(1) >= BBox(1,1) .and. xMin(1) <= BBox(1,2) .and. &
                  xMax(2) >= BBox(2,1) .and. xMin(2) <= BBox(2,2) .and. &
                  xMax(3) >= BBox(3,1) .and. xMin(3) <= BBox(3,2)) then
                elementWithinBBox(ii) = .true.
             else
                elementWithinBBox(ii) = .false.
             endif

          enddo
       enddo elementLoop1

       ! Determine the local number of elements within the global
       ! bounding box.

       ii = 0
       do i=1,nn
          if( elementWithinBBox(i) ) ii = ii + 1
       enddo

       ADT%nBBoxes = ii

       ! Allocate the memory for the bounding box coordinates, the
       ! corresponding element type and the index in the connectivity.

       allocate(ADT%xBBox(6,ii), ADT%elementType(ii), &
            ADT%elementID(ii), stat=ierr)
       if(ierr /= 0)                              &
            call adtTerminate(ADT, "buildSurfaceADT", &
            "Memory allocation failure for bounding &
            &box data.")

       ! Repeat the loop over the all the elements, but now store
       ! the bounding boxes in the ADT.

       ii = 0
       nn = 0
       elementLoop2: do ll=1,2

          ! Set the correct pointers for this element.

          call setSurfacePointers(ll)

          ! Loop over the elements and store the bounding box info,
          ! if needed.

          do i=1,nElem
             ii = ii + 1
             testWithin: if( elementWithinBBox(ii) ) then

                nn = nn + 1

                ADT%elementType(nn) = elType
                ADT%elementID(nn)   = i

                mm = conn(1,i)
                xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
                xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
                xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

                do j=2,nNPE
                   mm = conn(j,i)

                   xMin(1) = min(xMin(1),coor(1,mm))
                   xMin(2) = min(xMin(2),coor(2,mm))
                   xMin(3) = min(xMin(3),coor(3,mm))

                   xMax(1) = max(xMax(1),coor(1,mm))
                   xMax(2) = max(xMax(2),coor(2,mm))
                   xMax(3) = max(xMax(3),coor(3,mm))
                enddo

                ADT%xBBox(1,nn) = xMin(1)
                ADT%xBBox(2,nn) = xMin(2)
                ADT%xBBox(3,nn) = xMin(3)

                ADT%xBBox(4,nn) = xMax(1)
                ADT%xBBox(5,nn) = xMax(2)
                ADT%xBBox(6,nn) = xMax(3)

             endif testWithin
          enddo
       enddo elementLoop2

       ! Deallocate the memory for elementWithinBBox.

       deallocate(elementWithinBBox, stat=ierr)
       if(ierr /= 0)                              &
            call adtTerminate(ADT, "buildSurfaceADT", &
            "Deallocation failure for &
            &elementWithinBBox.")

    else testBBox

       ! No global bounding box. The number of local bounding boxes
       ! to be stored is the total number of local surface elements.

       ii = nTria + nQuads
       ADT%nBBoxes = ii

       ! Allocate the memory for the bounding box coordinates, the
       ! corresponding element type and the index in the connectivity.

       allocate(ADT%xBBox(6,ii), ADT%elementType(ii), &
            ADT%elementID(ii), stat=ierr)
       if(ierr /= 0)                              &
            call adtTerminate(ADT, "buildSurfaceADT", &
            "Memory allocation failure for bounding &
            &box data.")

       ! Loop over the number of element types present, i.e. 2,
       ! to store the bounding boxes; nn is the counter.

       nn = 0
       elementLoop3: do ll=1,2

          ! Set the correct pointers for this element.

          call setSurfacePointers(ll)

          ! Loop over the number of elements and store the bounding
          ! box info.

          do i=1,nElem
             nn = nn + 1

             ADT%elementType(nn) = elType
             ADT%elementID(nn)   = i

             mm = conn(1,i)
             xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
             xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
             xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

             do j=2,nNPE
                mm = conn(j,i)

                xMin(1) = min(xMin(1),coor(1,mm))
                xMin(2) = min(xMin(2),coor(2,mm))
                xMin(3) = min(xMin(3),coor(3,mm))

                xMax(1) = max(xMax(1),coor(1,mm))
                xMax(2) = max(xMax(2),coor(2,mm))
                xMax(3) = max(xMax(3),coor(3,mm))
             enddo

             ADT%xBBox(1,nn) = xMin(1)
             ADT%xBBox(2,nn) = xMin(2)
             ADT%xBBox(3,nn) = xMin(3)

             ADT%xBBox(4,nn) = xMax(1)
             ADT%xBBox(5,nn) = xMax(2)
             ADT%xBBox(6,nn) = xMax(3)
          enddo

       enddo elementLoop3

    endif testBBox

    ! Build the ADT from the now known boundary boxes.

    call buildADT(ADT)

    !===============================================================

  contains

    !=============================================================

    subroutine setSurfacePointers(ll)
      !
      !         **************************************************************
      !         *                                                            *
      !         * This internal subroutine sets the pointers to the correct  *
      !         * surface element, such that a loop over the element types   *
      !         * can be used.                                               *
      !         *                                                            *
      !         * Subroutine intent(in) arguments.                           *
      !         * --------------------------------                           *
      !         * ll: Element type for which the pointers must be used.      *
      !         *                                                            *
      !         **************************************************************
      !
      implicit none
      !
      !         Subroutine arguments.
      !
      integer, intent(in) :: ll

      select case (ll)
      case (1)
         elType = adtTriangle;      nElem = nTria;  nNPE = 3
         conn => triaConn
      case (2)
         elType = adtQuadrilateral; nElem = nQuads; nNPE = 4
         conn => quadsConn
      case (3)
      end select

    end subroutine setSurfacePointers

  end subroutine buildSurfaceADT

  !***************************************************************
  !***************************************************************

  subroutine buildVolumeADT(nTetra,    nPyra,    nPrisms,    &
       nHexa,     nNodes,   coor,       &
       tetraConn, pyraConn, prismsConn, &
       hexaConn,  BBox,     useBBox,    &
       comm,      adtID)
    !
    !       ****************************************************************
    !       *                                                              *
    !       * This routine builds the 6 dimensional ADT, which stores the  *
    !       * given volume grid. The memory intensive part of these        *
    !       * arguments, the arrays with the coordinates and               *
    !       * connectivities, are not copied. Instead pointers are set to  *
    !       * these arrays. It is therefore the responsibility of the user *
    !       * not to deallocate this memory before all the searches have   *
    !       * been performed.                                              *
    !       *                                                              *
    !       * Subroutine intent(in) arguments.                             *
    !       * --------------------------------                             *
    !       * nNodes:    Number of local nodes in the given grid.          *
    !       * nTetra:    Idem for the tetrahedra.                          *
    !       * nPyra:     Idem for the pyramids.                            *
    !       * nPrisms:   Idem for the prisms.                              *
    !       * nHexa:     Idem for the hexahedra.                           *
    !       * BBox(3,2): The possible bounding box. Only elements within   *
    !       *            this box will be stored in the ADT.               *
    !       * useBBox:   Whether or not to use the bounding box.           *
    !       * comm:      MPI-communicator for the global ADT.              *
    !       * adtID:     The ID of the ADT.                                *
    !       *                                                              *
    !       * Subroutine intent(in), target arguments.                     *
    !       * ----------------------------------------                     *
    !       * coor(3,nNodes):        Nodal coordinates of the local grid.  *
    !       * tetraConn(4,nTetra):   Local connectivity of the tetrahedra. *
    !       * pyraConn(5,nPyra):     Idem for the pyramids.                *
    !       * prismsConn(6,nPrisms): Idem for the prisms.                  *
    !       * hexaConn(8,nHexa):     Idem for the hexahedra.               *
    !       *                                                              *
    !       ****************************************************************
    !
    implicit none
    !
    !       Subroutine arguments.
    !
    integer, intent(in)          :: comm
    character(len=*), intent(in) :: adtID

    integer(kind=intType), intent(in) :: nTetra
    integer(kind=intType), intent(in) :: nPyra
    integer(kind=intType), intent(in) :: nPrisms
    integer(kind=intType), intent(in) :: nHexa
    integer(kind=intType), intent(in) :: nNodes

    logical, intent(in) :: useBBox

    integer(kind=intType), dimension(:,:), intent(in), &
         target :: tetraConn
    integer(kind=intType), dimension(:,:), intent(in), &
         target :: pyraConn
    integer(kind=intType), dimension(:,:), intent(in), &
         target :: prismsConn
    integer(kind=intType), dimension(:,:), intent(in), &
         target :: hexaConn

    real(kind=realType), dimension(3,2), intent(in) :: BBox

    real(kind=realType), dimension(:,:), intent(in), &
         target :: coor
    !
    !       Local variables.
    !
    integer :: ierr, ll, nNPE

    integer(kind=adtElementType) :: elType

    integer(kind=intType) :: i, j, ii, jj, mm, nn, nElem

    integer(kind=intType), dimension(:,:), pointer :: conn

    real(kind=realType), dimension(3) :: xMin, xMax

    logical, dimension(:), allocatable :: elementWithinBBox

    type(adtType), pointer :: ADT

    ! Allocate or reallocate the memory for ADTs. This depends
    ! whether or not this is the first ADT to be built.

    if( allocated(ADTs) ) then
       call reallocateADTs(adtID, jj)
    else
       call allocateADTs
       jj = 1
    endif

    ADT => ADTs(jj)

    ! Make sure the ADT is active and store the ID of this ADT.

    ADT%isActive = .true.
    ADT%adtID    = adtID

    ! Copy the communicator and determine the number of processors
    ! and my processor ID in this group.

    ADT%comm = comm
    call mpi_comm_rank(comm, ADT%myID,   ierr)
    call mpi_comm_size(comm, ADT%nProcs, ierr)

    ! Set the ADT type, which is a volume ADT.

    ADT%adtType = adtVolumeADT

    ! Copy the number of nodes and volume elements and set the number
    ! of surface elements to 0; only a volume grid has been given.

    ADT%nNodes  = nNodes
    ADT%nTetra  = nTetra
    ADT%nPyra   = nPyra
    ADT%nPrisms = nPrisms
    ADT%nHexa   = nHexa

    ADT%nTria  = 0
    ADT%nQuads = 0

    ! Set the pointers for the coordinates and the
    ! volume connectivities.

    ADT%coor       => coor
    ADT%tetraConn  => tetraConn
    ADT%pyraConn   => pyraConn
    ADT%prismsConn => prismsConn
    ADT%hexaConn   => hexaConn

    ! Determine the number of elements to be stored in the ADT.
    ! This depends whether or not the global bounding box should be
    ! used when building the ADT.

    testBBox: if( useBBox ) then

       ! Global bounding box is used. Allocate the memory for the
       ! logical elementWithinBBox.

       nn = nTetra + nPyra + nPrisms + nHexa
       allocate(elementWithinBBox(nn), stat=ierr)
       if(ierr /= 0)                             &
            call adtTerminate(ADT, "buildVolumeADT", &
            "Memory allocation failure for &
            &elementWithinBBox.")

       ! Loop over the number of element types.

       ii = 0
       elementLoop1: do ll=1,4

          ! Set the correct pointers for this element.

          call setVolumePointers(ll)

          ! Loop over the elements and determine the bounding box of
          ! each element.

          do i=1,nElem
             ii = ii + 1

             mm = conn(1,i)
             xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
             xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
             xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

             do j=2,nNPE
                mm = conn(j,i)

                xMin(1) = min(xMin(1),coor(1,mm))
                xMin(2) = min(xMin(2),coor(2,mm))
                xMin(3) = min(xMin(3),coor(3,mm))

                xMax(1) = max(xMax(1),coor(1,mm))
                xMax(2) = max(xMax(2),coor(2,mm))
                xMax(3) = max(xMax(3),coor(3,mm))
             enddo

             ! Check if the bounding box is (partially) inside the
             ! global bounding box. If so, set elementWithinBBox
             ! to .true.; otherwise set it to .false.

             if(xMax(1) >= BBox(1,1) .and. xMin(1) <= BBox(1,2) .and. &
                  xMax(2) >= BBox(2,1) .and. xMin(2) <= BBox(2,2) .and. &
                  xMax(3) >= BBox(3,1) .and. xMin(3) <= BBox(3,2)) then
                elementWithinBBox(ii) = .true.
             else
                elementWithinBBox(ii) = .false.
             endif

          enddo
       enddo elementLoop1

       ! Determine the local number of elements within the global
       ! bounding box.

       ii = 0
       do i=1,nn
          if( elementWithinBBox(i) ) ii = ii + 1
       enddo

       ADT%nBBoxes = ii

       ! Allocate the memory for the bounding box coordinates, the
       ! corresponding element type and the index in the connectivity.

       allocate(ADT%xBBox(6,ii), ADT%elementType(ii), &
            ADT%elementID(ii), stat=ierr)
       if(ierr /= 0)                             &
            call adtTerminate(ADT, "buildVolumeADT", &
            "Memory allocation failure for bounding &
            &box data.")

       ! Repeat the loop over the all the elements, but now store
       ! the bounding boxes in the ADT.

       ii = 0
       nn = 0
       elementLoop2: do ll=1,4

          ! Set the correct pointers for this element.

          call setVolumePointers(ll)

          ! Loop over the elements and store the bounding box info,
          ! if needed.

          do i=1,nElem
             ii = ii + 1
             testWithin: if( elementWithinBBox(ii) ) then

                nn = nn + 1

                ADT%elementType(nn) = elType
                ADT%elementID(nn)   = i

                mm = conn(1,i)
                xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
                xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
                xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

                do j=2,nNPE
                   mm = conn(j,i)

                   xMin(1) = min(xMin(1),coor(1,mm))
                   xMin(2) = min(xMin(2),coor(2,mm))
                   xMin(3) = min(xMin(3),coor(3,mm))

                   xMax(1) = max(xMax(1),coor(1,mm))
                   xMax(2) = max(xMax(2),coor(2,mm))
                   xMax(3) = max(xMax(3),coor(3,mm))
                enddo

                ADT%xBBox(1,nn) = xMin(1)
                ADT%xBBox(2,nn) = xMin(2)
                ADT%xBBox(3,nn) = xMin(3)

                ADT%xBBox(4,nn) = xMax(1)
                ADT%xBBox(5,nn) = xMax(2)
                ADT%xBBox(6,nn) = xMax(3)

             endif testWithin
          enddo
       enddo elementLoop2

       ! Deallocate the memory for elementWithinBBox.

       deallocate(elementWithinBBox, stat=ierr)
       if(ierr /= 0)                             &
            call adtTerminate(ADT, "buildVolumeADT", &
            "Deallocation failure for &
            &elementWithinBBox.")

    else testBBox

       ! No global bounding box. The number of local bounding boxes
       ! to be stored is the total number of local volume elements.

       ii = nTetra + nPyra + nPrisms + nHexa
       ADT%nBBoxes = ii

       ! Allocate the memory for the bounding box coordinates, the
       ! corresponding element type and the index in the connectivity.

       allocate(ADT%xBBox(6,ii), ADT%elementType(ii), &
            ADT%elementID(ii), stat=ierr)
       if(ierr /= 0)                             &
            call adtTerminate(ADT, "buildVolumeADT", &
            "Memory allocation failure for bounding &
            &box data.")

       ! Loop over the number of element types present, i.e. 4,
       ! to store the bounding boxes; nn is the counter.

       nn = 0
       elementLoop3: do ll=1,4

          ! Set the correct pointers for this element.

          call setVolumePointers(ll)

          ! Loop over the number of elements and store the bounding
          ! box info.

          do i=1,nElem
             nn = nn + 1

             ADT%elementType(nn) = elType
             ADT%elementID(nn)   = i

             mm = conn(1,i)
             xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
             xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
             xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

             do j=2,nNPE
                mm = conn(j,i)

                xMin(1) = min(xMin(1),coor(1,mm))
                xMin(2) = min(xMin(2),coor(2,mm))
                xMin(3) = min(xMin(3),coor(3,mm))

                xMax(1) = max(xMax(1),coor(1,mm))
                xMax(2) = max(xMax(2),coor(2,mm))
                xMax(3) = max(xMax(3),coor(3,mm))
             enddo

             ADT%xBBox(1,nn) = xMin(1)
             ADT%xBBox(2,nn) = xMin(2)
             ADT%xBBox(3,nn) = xMin(3)

             ADT%xBBox(4,nn) = xMax(1)
             ADT%xBBox(5,nn) = xMax(2)
             ADT%xBBox(6,nn) = xMax(3)
          enddo

       enddo elementLoop3

    endif testBBox

    ! Build the ADT from the now known boundary boxes.

    call buildADT(ADT)

    !===============================================================

  contains

    !=============================================================

    subroutine setVolumePointers(ll)
      !
      !         **************************************************************
      !         *                                                            *
      !         * This internal subroutine sets the pointers to the correct  *
      !         * volume element, such that a loop over the element types    *
      !         * can be used.                                               *
      !         *                                                            *
      !         * Subroutine intent(in) arguments.                           *
      !         * --------------------------------                           *
      !         * ll: Element type for which the pointers must be used.      *
      !         *                                                            *
      !         **************************************************************
      !
      implicit none
      !
      !         Subroutine arguments.
      !
      integer, intent(in) :: ll

      select case (ll)
      case (1)
         elType = adtTetrahedron; nElem = nTetra;  nNPE = 4
         conn => tetraConn
      case (2)
         elType = adtPyramid;     nElem = nPyra;   nNPE = 5
         conn => pyraConn
      case (3)
         elType = adtPrism;       nElem = nPrisms; nNPE = 6
         conn => prismsConn
      case (4)
         elType = adtHexahedron;  nElem = nHexa;   nNPE = 8
         conn => hexaConn
      end select

    end subroutine setVolumePointers

  end subroutine buildVolumeADT

  subroutine buildSerialHex(nHexa, nNodes, coor, hexaConn, ADT)
    !
    !       ****************************************************************
    !       *                                                              *
    !       * This a specialized routine that builds and ADT tree for      *
    !       * hex volumes only and only in serial. Also, this routine does *
    !       * use adtDats's ADTs() array list...the user must supply the   *
    !       * adtType to use and is responsible for all data management of *
    !       * this type.                                                   *
    !       * The memory intensive part of these                           *
    !       * arguments, the arrays with the coordinates and               *
    !       * connectivities, are not copied. Instead pointers are set to  *
    !       * these arrays. It is therefore the responsibility of the user *
    !       * not to deallocate this memory before all the searches have   *
    !       * been performed.                                              *
    !       *                                                              *
    !       * Subroutine intent(in) arguments.                             *
    !       * --------------------------------                             *
    !       * nHexa :    Number of hexa cells
    !       * nNodes:    Number of nodes in the given grid.                *
    !       *                                                              *
    !       * Subroutine intent(in), target arguments.                     *
    !       * ----------------------------------------                     *
    !       * coor(3,nNodes):        Nodal coordinates of the local grid.  *
    !       * hexaConn(8,nHexa):     Idem for the hexahedra.               *
    !       *                                                              *
    !       * Subroutine intent(out), arguments.                           *
    !       * ----------------------------------------                     *
    !       * ADT : The newly completed ADT                                *
    !       *                                                              *
    !       ****************************************************************
    !
    implicit none
    !
    !       Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nHexa
    integer(kind=intType), intent(in) :: nNodes

    integer(kind=intType), dimension(:,:), intent(in), &
         target :: hexaConn

    real(kind=realType), dimension(:,:), intent(in), &
         target :: coor
    type(adtType), intent(out) :: ADT
    !
    !       Local variables.
    !
    integer :: ierr, ll, nNPE
    integer(kind=intType) :: i, j, mm
    real(kind=realType), dimension(3) :: xMin, xMax

    ! We need to set comm...explictly mpi_comm_self
    ADT%comm = MPI_COMM_SELF
    ADT%nProcs = 1
    ADT%myID = 0
    ! Set the ADT type, which is a volume ADT.

    ADT%adtType = adtVolumeADT

    ! Copy the number of nodes and volume elements and set the number
    ! of surface elements to 0; only a volume grid has been given.

    ADT%nNodes  = nNodes
    ADT%nHexa   = nHexa
    ADT%nTetra  = 0
    ADT%nPyra   = 0
    ADT%nPrisms = 0
    ADT%nTria  = 0
    ADT%nQuads = 0

    ! Set the pointers for the coordinates and the
    ! volume connectivities.

    ADT%coor       => coor
    ADT%hexaConn   => hexaConn
    nullify(ADT%tetraConn, ADT%pyraConn, ADT%prismsConn)

    ADT%nBBoxes = nHexa

    ! Allocate the memory for the bounding box coordinates, the
    ! corresponding element type and the index in the connectivity.

    allocate(ADT%xBBox(6, nHexa))
    allocate(ADT%elementType(nHexa))
    allocate(ADT%elementID(nHexa))

    ! All hexas
    ADT%elementType = adtHexahedron


    ! Loop over the number of elements and store the bounding
    ! box info.
    nNPE = 8

    do i=1,nHexa

       mm = hexaConn(1,i)
       xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
       xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
       xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

       do j=2,nNPE
          mm = hexaConn(j,i)

          xMin(1) = min(xMin(1),coor(1,mm))
          xMin(2) = min(xMin(2),coor(2,mm))
          xMin(3) = min(xMin(3),coor(3,mm))

          xMax(1) = max(xMax(1),coor(1,mm))
          xMax(2) = max(xMax(2),coor(2,mm))
          xMax(3) = max(xMax(3),coor(3,mm))
       enddo

       ADT%xBBox(1,i) = xMin(1)
       ADT%xBBox(2,i) = xMin(2)
       ADT%xBBox(3,i) = xMin(3)

       ADT%xBBox(4,i) = xMax(1)
       ADT%xBBox(5,i) = xMax(2)
       ADT%xBBox(6,i) = xMax(3)

       ! elementID is just sequential since we only have 1 element type
       ADT%elementID(i) = i

    enddo

    ! Build the ADT from the now known boundary boxes.

    call buildADT(ADT)

  end subroutine buildSerialHex

  subroutine destroySerialHex(ADT)
    ! Deallocate the data allocated from the ADT

    implicit none
    type(adtType), intent(inout) :: ADT

    deallocate(ADT%xBBox)
    deallocate(ADT%elementType)
    deallocate(ADT%elementID)
    deallocate(ADT%ADTree)
  end subroutine destroySerialHex

  subroutine buildSerialQuad(nQuad, nNodes, coor, quadsConn, ADT)
    !
    !       ****************************************************************
    !       *                                                              *
    !       * This a specialized routine that builds and ADT tree for      *
    !       * quad surface meshes and only in serial. Also, this routine   *
    !       * does not                                                     *
    !       * use adtDats's ADTs() array list...the user must supply the   *
    !       * adtType to use and is responsible for all data management of *
    !       * this type.                                                   *
    !       * The memory intensive part of these                           *
    !       * arguments, the arrays with the coordinates and               *
    !       * connectivities, are not copied. Instead pointers are set to  *
    !       * these arrays. It is therefore the responsibility of the user *
    !       * not to deallocate this memory before all the searches have   *
    !       * been performed.                                              *
    !       *                                                              *
    !       * Subroutine intent(in) arguments.                             *
    !       * --------------------------------                             *
    !       * nQuad :    Number of quad cells
    !       * nNodes:    Number of nodes in the given grid.                *
    !       *                                                              *
    !       * Subroutine intent(in), target arguments.                     *
    !       * ----------------------------------------                     *
    !       * coor(3,nNodes):        Nodal coordinates of the local grid.  *
    !       * quadsConn(8,nHexa):     Connectivity for quad cells           *
    !       *                                                              *
    !       * Subroutine intent(out), arguments.                           *
    !       * ----------------------------------------                     *
    !       * ADT : The newly completed ADT                                *
    !       *                                                              *
    !       ****************************************************************
    !
    implicit none
    !
    !       Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nQuad
    integer(kind=intType), intent(in) :: nNodes

    integer(kind=intType), dimension(:,:), intent(in), &
         target :: quadsConn

    real(kind=realType), dimension(:,:), intent(in), &
         target :: coor
    type(adtType), intent(out) :: ADT
    !
    !       Local variables.
    !
    integer :: ierr, ll, nNPE
    integer(kind=intType) :: i, j, mm
    real(kind=realType), dimension(3) :: xMin, xMax

    ! We need to set comm...explictly mpi_comm_self
    ADT%comm = MPI_COMM_SELF
    ADT%nProcs = 1
    ADT%myID = 0
    ! Set the ADT type, which is a surface ADT.

    ADT%adtType = adtSurfaceADT

    ! Copy the number of nodes and volume elements and set the number
    ! of surface elements to 0; only a volume grid has been given.

    ADT%nNodes  = nNodes
    ADT%nHexa   = 0
    ADT%nTetra  = 0
    ADT%nPyra   = 0
    ADT%nPrisms = 0
    ADT%nTria  = 0
    ADT%nQuads = nQuad

    ! Set the pointers for the coordinates and the
    ! volume connectivities.

    ADT%coor       => coor
    ADT%quadsConn   => quadsConn
    nullify(ADT%triaConn)

    ADT%nBBoxes = nQuad

    ! Allocate the memory for the bounding box coordinates, the
    ! corresponding element type and the index in the connectivity.

    allocate(ADT%xBBox(6, nQuad))
    allocate(ADT%elementType(nQuad))
    allocate(ADT%elementID(nQuad))

    ! All hexas
    ADT%elementType = adtQuadrilateral

    ! Loop over the number of elements and store the bounding
    ! box info.
    nNPE = 4

    do i=1, nQuad

       mm = quadsConn(1,i)
       xMin(1) = coor(1,mm); xMax(1) = coor(1,mm)
       xMin(2) = coor(2,mm); xMax(2) = coor(2,mm)
       xMin(3) = coor(3,mm); xMax(3) = coor(3,mm)

       do j=2,nNPE
          mm = quadsConn(j,i)

          xMin(1) = min(xMin(1),coor(1,mm))
          xMin(2) = min(xMin(2),coor(2,mm))
          xMin(3) = min(xMin(3),coor(3,mm))

          xMax(1) = max(xMax(1),coor(1,mm))
          xMax(2) = max(xMax(2),coor(2,mm))
          xMax(3) = max(xMax(3),coor(3,mm))
       enddo

       ADT%xBBox(1,i) = xMin(1)
       ADT%xBBox(2,i) = xMin(2)
       ADT%xBBox(3,i) = xMin(3)

       ADT%xBBox(4,i) = xMax(1)
       ADT%xBBox(5,i) = xMax(2)
       ADT%xBBox(6,i) = xMax(3)

       ! elementID is just sequential since we only have 1 element type
       ADT%elementID(i) = i

    enddo

    ! Build the ADT from the now known boundary boxes.

    call buildADT(ADT)

  end subroutine buildSerialQuad

  subroutine destroySerialQuad(ADT)
    ! Deallocate the data allocated from the ADT

    implicit none
    type(adtType), intent(inout) :: ADT

    deallocate(ADT%xBBox)
    deallocate(ADT%elementType)
    deallocate(ADT%elementID)
    deallocate(ADT%ADTree)
  end subroutine destroySerialQuad
end module adtBuild
