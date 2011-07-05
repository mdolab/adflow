!
!      ******************************************************************
!      *                                                                *
!      * File:          allocMemHaloList.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-12-2003                                      *
!      * Last modified: 11-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine allocMemHaloList(level)
!
!      ******************************************************************
!      *                                                                *
!      * allocMemHaloList allocates the memory for the variables        *
!      * needed to construct the communication lists and the periodic   *
!      * information. These variables are located in the module         *
!      * haloList. Only the 1st level halo variables are allocated here *
!      * to avoid unnecessary memory usage. The 2nd level cell halo     *
!      * are allocated later on in init2ndLevelCellHalos.               *
!      *                                                                *
!      ******************************************************************
!
       use block
       use haloList
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i
       integer(kind=intType) :: ie, je, ke, ib, jb, kb
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the 1st level cell and node halo lists.

       allocate(cellHalo1st(nCellHalo1st), &
                nodeHalo1st(nNodeHalo1st), stat=ierr)
       if(ierr /= 0)                                               &
         call terminate("allocMemHaloList",                        &
                        "Memory allocation failure for cellHalo1st &
                        &and nodeHalo1st")

       ! Initialize the level of indirectness to 0. This will only be
       ! overwritten for indirect boundary halo's. Also initialize the
       ! number of periodid subfaces to 0.

       do i=1,nCellHalo1st
         cellHalo1st(i)%levOfInd          = 0
         cellHalo1st(i)%nPeriodicSubfaces = 0
         nullify(cellHalo1st(i)%periodicSubfaces)
         nullify(cellHalo1st(i)%interp)
       enddo

       do i=1,nNodeHalo1st
         nodeHalo1st(i)%levOfInd          = 0
         nodeHalo1st(i)%nPeriodicSubfaces = 0
         nullify(nodeHalo1st(i)%periodicSubfaces)
         nullify(nodeHalo1st(i)%interp)
       enddo

       ! Allocate the memory to store the short hand of the transformation
       ! matrix for both cell and nodal halo's. In principle this matrix
       ! is only needed for the face (i.e. direct) halo's. However the
       ! difference between the number of 1st level halo's and the cell
       ! halo's is not so large and at this time of the program you should
       ! not worry too much about memory, because the metrics as well as
       ! the solution variables have not been allocated yet.

       allocate(transformCell(nCellHalo1st,3), &
                transformNode(nNodeHalo1st,3), stat=ierr)
       if(ierr /= 0)                                                 &
         call terminate("allocMemHaloList"  ,                        &
                        "Memory allocation failure for transformCell &
                        &and transformNode")

       ! Allocate the memory for nodeIndex and cellIndex, which will
       ! store the indices per block in the lists above.

       allocate(nodeIndex(nDom), cellIndex(nDom), stat=ierr)
       if(ierr /= 0)                                             &
         call terminate("allocMemHaloList",                      &
                        "Memory allocation failure for nodeIndex &
                        &and cellIndex")

       ! Loop over the number of blocks to allocate and initialize
       ! the elements of nodeIndex and cellIndex.

       do i=1,nDom

         ! Store the upper indices in the allocation a bit easier.

         ie = flowDoms(i,level,1)%ie
         je = flowDoms(i,level,1)%je
         ke = flowDoms(i,level,1)%ke

         ib = flowDoms(i,level,1)%ib
         jb = flowDoms(i,level,1)%jb
         kb = flowDoms(i,level,1)%kb

         ! Allocate the memory for entryList.

         allocate(nodeIndex(i)%entryList(0:ie,0:je,0:ke), &
                  cellIndex(i)%entryList(0:ib,0:jb,0:kb), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("allocMemHaloList", &
                          "Memory allocation failure for entryList")

         ! Initialize entryList to zero. This serves as a check later
         ! on. Cell halo's are uniquely defined via the 1 to 1 block
         ! connectivity, but for node halo's (on the boundary of a
         ! subface) several possibilities exist.

         nodeIndex(i)%entryList = 0
         cellIndex(i)%entryList = 0

       enddo

       end subroutine allocMemHaloList
