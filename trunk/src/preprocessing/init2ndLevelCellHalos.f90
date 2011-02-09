!
!      ******************************************************************
!      *                                                                *
!      * File:          init2ndLevelCellHalos.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-04-2003                                      *
!      * Last modified: 11-29-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine init2ndLevelCellHalos
!
!      ******************************************************************
!      *                                                                *
!      * init2ndLevelCellHalos initializes the 2nd level cell halo      *
!      * list. Basically the 1st level cell halo list is copied and the *
!      * counter iicell2nd is set to nCellHalo1st. This means that      *
!      * the 2nd level cell halo's are appended to the first level      *
!      * halo's. They are stored in a separate list, because the        *
!      * communication pattern of the 2nd level halo's is separate from *
!      * the 1st level halo's. Efficiency is the reason to do this; it  *
!      * is more efficient to send one big message than two smaller     *
!      * ones.                                                          *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, jj
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the 2nd level halo list.

       allocate(cellHalo2nd(nCellHalo2nd), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("init2ndLevelCellHalos", &
                        "Memory allocation failure for cellHalo2nd")

       ! Initialize iicell2nd to nCellHalo1st.

       iiCell2nd = nCellHalo1st

       ! Copy the information from the 1st level cell halo's into the
       ! second level cell halo list. Make sure to make a deep copy.

       do i=1,nCellHalo1st
         cellHalo2nd(i)%myBlock    = cellHalo1st(i)%myBlock
         cellHalo2nd(i)%myI        = cellHalo1st(i)%myI
         cellHalo2nd(i)%myJ        = cellHalo1st(i)%myJ
         cellHalo2nd(i)%myK        = cellHalo1st(i)%myK
         cellHalo2nd(i)%donorProc  = cellHalo1st(i)%donorProc
         cellHalo2nd(i)%donorBlock = cellHalo1st(i)%donorBlock
         cellHalo2nd(i)%dI         = cellHalo1st(i)%dI
         cellHalo2nd(i)%dJ         = cellHalo1st(i)%dJ
         cellHalo2nd(i)%dK         = cellHalo1st(i)%dK
         cellHalo2nd(i)%levOfInd   = cellHalo1st(i)%levOfInd

         nullify(cellHalo2nd(i)%interp)

         cellHalo2nd(i)%nPeriodicSubfaces = &
                               cellHalo1st(i)%nPeriodicSubfaces

         if(cellHalo2nd(i)%nPeriodicSubfaces > 0) then
           jj = cellHalo2nd(i)%nPeriodicSubfaces
           allocate(cellHalo2nd(i)%periodicSubfaces(jj), stat=ierr)
           if(ierr /= 0)                             &
             call terminate("init2ndLevelCellHalos", &
                            "Memory allocation failure for &
                            &periodicSubfaces")
           do j=1,jj
             cellHalo2nd(i)%periodicSubfaces(j) = &
                           cellHalo1st(i)%periodicSubfaces(j)
           enddo
         else
           nullify(cellHalo2nd(i)%periodicSubfaces)
         endif
       enddo

       ! Initialize the level of indirectness for the rest of the list
       ! to 0 and initialize the periodic data to 0 as well.

       do i=(nCellHalo1st+1), nCellHalo2nd
         cellHalo2nd(i)%levOfInd          = 0
         cellHalo2nd(i)%nPeriodicSubfaces = 0
         nullify(cellHalo2nd(i)%periodicSubfaces)
         nullify(cellHalo2nd(i)%interp)
       enddo

       end subroutine init2ndLevelCellHalos
