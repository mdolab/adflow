!
!      ******************************************************************
!      *                                                                *
!      * File:          determineCommPattern.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-20-2003                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineCommPattern(level)
!
!      ******************************************************************
!      *                                                                *
!      * determineCommPattern determines the communication pattern      *
!      * for the indicated grid level from the given block distribution *
!      * and corresponding halo info. Both the first and second level   *
!      * cell halo communication pattern as well as the first level     *
!      * nodal halo communication pattern is determined.                *
!      * A recursive algorithm is used. First the face halo's are       *
!      * determined and from those the indirect halo's can be obtained  *
!      * by looping over the level of indirectness.                     *
!      *                                                                *
!      * This routine controls the creation of the communication        *
!      * pattern and basically contains the function calls to the       *
!      * subtasks.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use haloList
       use periodicInfo
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer               :: ierr
       integer(kind=intType) :: i
!
!      Interfaces
!
       interface
       subroutine determineIndirectHalos(nHalo, iihalo, entityHalo,  &
                                         transform, entityIndex,     &
                                         start, nLevel, offset,      &
                                         gridLevel)
       use haloList
       use indirectHalo
       use communication
       implicit none
       integer(kind=intType), intent(in) :: nHalo, start, nLevel, offset
       integer(kind=intType), intent(in) :: gridLevel
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType), dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
       end subroutine determineIndirectHalos

       !=================================================================

       subroutine determinePeriodicData(entityHalo,   nHalo, &
                                        externalComm, internalComm)
       use communication
       use haloList
       use periodicInfo
       implicit none
       integer(kind=intType), intent(in)            :: nHalo
       type(haloListType), dimension(:), intent(in) :: entityHalo

       type(commType),         intent(inout) :: externalComm
       type(internalCommType), intent(inout) :: internalComm

       end subroutine determinePeriodicData
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of periodic faces.

       call determinePeriodicFaces

       ! Determine the amount of 1st and 2nd level cell halo's and 1st
       ! level node halo's.

       call determineNumberOfHalos(level)

       ! Allocate the memory for the help variables needed to
       ! determine the communication pattern.

       call allocMemHaloList(level)

       ! Determine the face halo's for both the 1st level cell and the
       ! 1st level node halo's.

       call determineFaceHalos(level)

       ! Determine the indirect 1st level node halo's.

       call determineIndirectHalos(nNodeHalo1st, iinode1st,    &
                                   nodeHalo1st, transformNode, &
                                   nodeIndex,   1_intType,     &
                                   1_intType,   0_intType, level)

       ! Release the memory of transformNode and nodeIndex.

       do i=1,nDom
         deallocate(nodeIndex(i)%entryList, stat=ierr)
         if(ierr /= 0)                           &
          call terminate("determineCommPattern", &
                         "Deallocation error for &
                         &nodeIndex(i)%entryList")
       enddo

       deallocate(nodeIndex, transformNode, stat=ierr)
       if(ierr /= 0)                           &
        call terminate("determineCommPattern", &
                       "Deallocation error for nodeIndex &
                       &and transformNode")

       ! Determine the indirect 1st level cell halo's.

       call determineIndirectHalos(nCellHalo1st, iicell1st,    &
                                   cellHalo1st, transformCell, &
                                   cellIndex,   2_intType,     &
                                   1_intType,    0_intType, level)

       ! Initialize the 2nd level cell halo list. Basically the 1st
       ! level cell halo list is copied.

       call init2ndLevelCellHalos

       ! Determine the indirect 2nd level cell halo's. As the indirect
       ! 1st level halo's are already treated an offset of 1 is passed.

       call determineIndirectHalos(nCellHalo2nd, iicell2nd,    &
                                   cellHalo2nd, transformCell, &
                                   cellIndex,   2_intType,     &
                                   2_intType,   1_intType, level)

       ! Release the memory of transformCell and cellIndex.

       do i=1,nDom
         deallocate(cellIndex(i)%entryList, stat=ierr)
         if(ierr /= 0)                             &
          call terminate("determineCommPattern", &
                         "Deallocation error for &
                         &cellIndex(i)%entryList")
       enddo

       deallocate(cellIndex, transformCell, stat=ierr)
       if(ierr /= 0)                           &
        call terminate("determineCommPattern", &
                       "Deallocation error for cellIndex &
                       &and transformCell")

       ! Sort the three lists in increasing order.

       call qsortHaloListType(nodeHalo1st, nNodeHalo1st)
       call qsortHaloListType(cellHalo1st, nCellHalo1st)
       call qsortHaloListType(cellHalo2nd, nCellHalo2nd)

       ! Determine the final communication data structures to store the
       ! halo info.

       call finalCommStructures(nodeHalo1st, nNodeHalo1st,  &
                                commPatternNode_1st(level), &
                                internalNode_1st(level), 0_intType)

       call finalCommStructures(cellHalo1st, nCellHalo1st,  &
                                commPatternCell_1st(level), &
                                internalCell_1st(level), 0_intType)

       call finalCommStructures(cellHalo2nd, nCellHalo2nd,  &
                                commPatternCell_2nd(level), &
                                internalCell_2nd(level), 0_intType)

       ! Determine the transformation for periodic halo's.

       call determinePeriodicData(nodeHalo1st, nNodeHalo1st,  &
                                  commPatternNode_1st(level), &
                                  internalNode_1st(level))

       call determinePeriodicData(cellHalo1st, nCellHalo1st,  &
                                  commPatternCell_1st(level), &
                                  internalCell_1st(level))

       call determinePeriodicData(cellHalo2nd, nCellHalo2nd,  &
                                  commPatternCell_2nd(level), &
                                  internalCell_2nd(level))

       ! Deallocate the memory for the 3 halo lists.

       call deallocatePointersHaloList(nodeHalo1st, nNodeHalo1st)
       call deallocatePointersHaloList(cellHalo1st, nCellHalo1st)
       call deallocatePointersHaloList(cellHalo2nd, nCellHalo2nd)

       deallocate(nodeHalo1st, cellHalo1st, cellHalo2nd, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("determineCommPattern", &
                        "Deallocation error for nodeHalo1st, &
                        &cellHalo1st and cellHalo2nd")

       ! Deallocate the memory of periodicGlobal.

       deallocate(periodicGlobal, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("determineCommPattern", &
                        "Deallocation error for periodicGlobal")

       end subroutine determineCommPattern

       !-----------------------------------------------------------------

       subroutine deallocatePointersHaloList(entityHalo, nHalo)
!
!      ******************************************************************
!      *                                                                * 
!      * deallocatePointersHaloList deallocates the memory of the       *
!      * pointer variables of entityHalo.                               *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)               :: nHalo
       type(haloListType), dimension(*), intent(inout) :: entityHalo
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of halo's and deallocate the memory
       ! of the pointer variables, if allocated.

       do i=1,nHalo
         if( associated(entityHalo(i)%interp) ) then
           deallocate(entityHalo(i)%interp, stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("deallocatePointersHaloList", &
                            "Deallocation failure for interp")
         endif

         if( associated(entityHalo(i)%periodicSubfaces) ) then
           deallocate(entityHalo(i)%periodicSubfaces, stat=ierr)
           if(ierr /= 0)                                  &
             call terminate("deallocatePointersHaloList", &
                            "Deallocation failure for periodicSubfaces")
         endif
       enddo

       end subroutine deallocatePointersHaloList
