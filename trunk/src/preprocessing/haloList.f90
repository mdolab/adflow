!
!      ******************************************************************
!      *                                                                *
!      * File:          haloList.f90                                    *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 01-22-2003                                      *
!      * Last modified: 05-24-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module haloList
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains temporary variables to create the   *
!      * list of halo cells and nodes.                                  *
!      *                                                                *
!      ******************************************************************
!
       use precision
       implicit none
       save

       public
       private :: lessEqualHaloListType
       private :: lessHaloListType
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the variables for the 3 lists.               *
!      *                                                                *
!      ******************************************************************
!
       type haloListType

         ! myBlock:           local block ID of the halo.
         ! myI, myJ, myK:     i,j,k indices of the halo.
         ! donorProc :        processor where donor is stored. In case
         !                    the halo is a boundary halo, donorProc
         !                    is set to -1.
         ! donorBlock:        block ID of the donor. In case the halo
         !                    is a boundary halo donorBlock is set
         !                    to the corresponding boundary condition.
         ! dI, dJ, dK:        i,j,k indices of the donor.
         ! levOfInd:          level of indirectness.
         ! interp(..):        interpolants for the donor stencil; only
         !                    allocated for lists requiring this info.
         ! nPeriodicSubfaces: Number of periodic subfaces that are
         !                    crossed when going from the halo to the
         !                    donor. This is at most the level of
         !                    indirectness of the halo.
         ! periodicSubfaces:  The corresponding subfaces ID's according
         !                    to the sequence defined in periodicGlobal.

         integer(kind=intType) :: myBlock
         integer(kind=intType) :: myI, myJ, myK
         integer(kind=intType) :: donorProc, donorBlock
         integer(kind=intType) :: dI, dJ, dK
         integer(kind=intType) :: levOfInd

         real(kind=realType), dimension(:), pointer :: interp

         integer(kind=intType) :: nPeriodicSubfaces
         integer(kind=intType), dimension(:), pointer :: periodicSubfaces

       end type haloListType

       ! Interface for the extension of the operators <= and <.
       ! These are needed for the sorting of haloListType.

       interface operator(<=)
         module procedure lessEqualHaloListType
       end interface

       interface operator(<)
         module procedure lessHaloListType
       end interface

       ! nCellHalo1st: # of 1st level cell halo's
       ! nCellHalo2nd: # of 2nd level cell halo's
       ! nNodeHalo1st: # of 1st level node halo's

       integer(kind=intType) :: nCellHalo1st, nCellHalo2nd
       integer(kind=intType) :: nNodeHalo1st

       ! iiCell1st: Counter variable for the 1st level cell halo's
       ! iiCell2nd: Counter variable for the 2nd level cell halo's
       ! iiNode1st: Counter variable for the 1st level node halo's

       integer(kind=intType) :: iiCell1st, iiCell2nd, iiNode1st

       ! cellHalo1st(nCellHalo1st) :: List of halo info for 1st
       !                              level cell halo's.
       ! cellHalo2nd(nCellHalo2nd) :: Idem for 2nd level cell halo's.
       ! nodeHalo1st(nNodeHalo1st) :: Idem for 1st level node halo's.

       type(haloListType), dimension(:), allocatable :: cellHalo1st
       type(haloListType), dimension(:), allocatable :: cellHalo2nd
       type(haloListType), dimension(:), allocatable :: nodeHalo1st

       ! transformCell(nCellHalo1st,3) :: Short hand for the transformation
       !                                  matrix between the halo and
       !                                  the donor for cell based halo's.
       !                                  In principle the size equals the
       !                                  number of faces (i.e. direct)
       !                                  halo's, but the difference is
       !                                  not so large.
       ! transformNode(nNodeHalo1st,3) :: Idem for the nodes.

       integer(kind=intType), dimension(:,:), allocatable :: transformCell
       integer(kind=intType), dimension(:,:), allocatable :: transformNode
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the index variables, which store for each    *
!      * i,j,k in the block the index in the corresponding list.        *
!      * I know I'm wasting memory here (because only the halo's are    *
!      * relevant), but that's not too much of a problem. The reason is *
!      * that neither the metrics nor the variables have been allocated *
!      * yet. So later on, much more memory is needed than the single   *
!      * integer for each cell/node used here.                          *
!      *                                                                *
!      ******************************************************************
!
       type indexListType

         ! entryList(:,:,:)  :: Corresponding entry in the list.
         !                      Dimensions are either 0:ie,0:je,0:ke for
         !                      the node or 0:ib,0:jb,0:kb for the cell
         !                      based halo's. The latter is then suited
         !                      for the 2nd level halo's.

         integer(kind=intType), dimension(:,:,:), pointer :: entryList

       end type indexListType

       ! nodeIndex(nDom) :: The node indices for every block.
       ! cellIndex(nDom) :: Idem for the cells.

       type(indexListType), allocatable, dimension(:) :: nodeIndex
       type(indexListType), allocatable, dimension(:) :: cellIndex

       contains
!
!        ****************************************************************
!        *                                                              *
!        * Functions to simulate the operators <= and <.                *
!        *                                                              *
!        ****************************************************************
!
         logical function lessEqualHaloListType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * lessEqual returns .true. if g1 <= g2 and .false. otherwise.  *
!        * The comparison is firstly based on the processor ID of the   *
!        * donor. After that it depends whether the halo is a boundary  *
!        * halo or not. Note that boundary halo's have a donor processor*
!        * if of -1, such that they are always first in the list.       *
!        *                                                              *
!        ****************************************************************
!
         use BCTypes
         implicit none
!
!        Function arguments.
!
         type(haloListType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the donor processors first. If not equal,
         ! set lessEqual appropriately and return.

         if(g1%donorProc < g2%donorProc) then
           lessEqualHaloListType = .true.
           return
         else if(g1%donorProc > g2%donorProc) then
           lessEqualHaloListType = .false.
           return
         endif

         ! Donor processors are identical. Now it depends whether we are
         ! dealing with boundary halo's or not.

         boundary: if(g1%donorProc == -1) then ! And thus
                                               ! g2%donorProc == -1

           ! Both halo's are boundary halo's. Compare the block ID of
           ! the halo's.

           if(g1%myBlock < g2%myBlock) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myBlock > g2%myBlock) then
             lessEqualHaloListType = .false.
             return
           endif

           ! Compare the boundary conditions, which are stored in
           ! donorBlock. Note that the sequence in BCTypes is such that
           ! the most important BC has the highest number.

           if(g1%donorBlock < g2%donorBlock) then
             lessEqualHaloListType = .true.
             return
           else if(g1%donorBlock > g2%donorBlock) then
             lessEqualHaloListType = .false.
             return
           endif

           ! As it is possible that indirect halo's need donor info from
           ! direct halo's or even indirect halo's with a smaller level
           ! of indirectness, compare the level of indirectness.

           if(g1%levOfInd < g2%levOfInd) then
             lessEqualHaloListType = .true.
             return
           else if(g1%levOfInd > g2%levOfInd) then
             lessEqualHaloListType = .false.
             return
           endif

           ! Compare the indices of the halo. First k, then j and
           ! finally i.

           if(g1%myK < g2%myK) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myK > g2%myK) then
             lessEqualHaloListType = .false.
             return
           endif

           if(g1%myJ < g2%myJ) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myJ > g2%myJ) then
             lessEqualHaloListType = .false.
             return
           endif

           if(g1%myI < g2%myI) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myI > g2%myI) then
             lessEqualHaloListType = .false.
             return
           endif

           ! No need to compare anything else; g1 == g2.

         else boundary

           ! Both halo's are internal halo's, whose donor is stored on
           ! the same processor. Compare the donor blocks.

           if(g1%donorBlock < g2%donorBlock) then
             lessEqualHaloListType = .true.
             return
           else if(g1%donorBlock > g2%donorBlock) then
             lessEqualHaloListType = .false.
             return
           endif

           ! Also the blocks are identical. Compare the donor indices.
           ! First the k index.

           if(g1%dK < g2%dK) then
             lessEqualHaloListType = .true.
             return
           else if(g1%dK > g2%dK) then
             lessEqualHaloListType = .false.
             return
           endif

           ! The j index.

           if(g1%dJ < g2%dJ) then
             lessEqualHaloListType = .true.
             return
           else if(g1%dJ > g2%dJ) then
             lessEqualHaloListType = .false.
             return
           endif

           ! And the i index.

           if(g1%dI < g2%dI) then
             lessEqualHaloListType = .true.
             return
           else if(g1%dI > g2%dI) then
             lessEqualHaloListType = .false.
             return
           endif

           ! The donors are identical. Compare the halo's.
           ! First the block id.

           if(g1%myBlock < g2%myBlock) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myBlock > g2%myBlock) then
             lessEqualHaloListType = .false.
             return
           endif

           ! Halo blocks are also identical. Finally compare the
           ! halo indices. Start with k.

           if(g1%myK < g2%myK) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myK > g2%myK) then
             lessEqualHaloListType = .false.
             return
           endif

           ! The j index.

           if(g1%myJ < g2%myJ) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myJ > g2%myJ) then
             lessEqualHaloListType = .false.
             return
           endif

           ! The i index.

           if(g1%myI < g2%myI) then
             lessEqualHaloListType = .true.
             return
           else if(g1%myI > g2%myI) then
             lessEqualHaloListType = .false.
             return
           endif

         endif boundary

         ! Both entities are identical. So set lessEqual to .true.

         lessEqualHaloListType = .true.

         end function lessEqualHaloListType

!        ================================================================

         logical function lessHaloListType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 < g2 and .false.          *
!        * otherwise. It is basically the same as the lessEqual         *
!        * function, except that the equality is now considered as      *
!        * .false.                                                      *
!        *                                                              *
!        ****************************************************************
!
         use BCTypes
         implicit none
!
!        Function arguments.
!
         type(haloListType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Compare the donor processors first. If not equal,
         ! set the function appropriately and return.

         if(g1%donorProc < g2%donorProc) then
           lessHaloListType = .true.
           return
         else if(g1%donorProc > g2%donorProc) then
           lessHaloListType = .false.
           return
         endif

         ! Donor processors are identical. Now it depends whether we are
         ! dealing with boundary halo's or not.

         boundary: if(g1%donorProc == -1) then ! And thus
                                               ! g2%donorProc == -1

           ! Both halo's are boundary halo's. Compare the block ID of
           ! the halo's.

           if(g1%myBlock < g2%myBlock) then
             lessHaloListType = .true.
             return
           else if(g1%myBlock > g2%myBlock) then
             lessHaloListType = .false.
             return
           endif

           ! Compare the boundary conditions, which are stored in
           ! donorBlock. Note that the sequence in BCTypes is such that
           ! the most important bc has the highest number.

           if(g1%donorBlock < g2%donorBlock) then
             lessHaloListType = .true.
             return
           else if(g1%donorBlock > g2%donorBlock) then
             lessHaloListType = .false.
             return
           endif

           ! As it is possible that indirect halo's need donor info from
           ! direct halo's or even indirect halo's with a smaller level
           ! of indirectness, compare the level of indirectness.

           if(g1%levOfInd < g2%levOfInd) then
             lessHaloListType = .true.
             return
           else if(g1%levOfInd > g2%levOfInd) then
             lessHaloListType = .false.
             return
           endif

           ! Compare the indices of the halo. First k, then j and
           ! finally i.

           if(g1%myK < g2%myK) then
             lessHaloListType = .true.
             return
           else if(g1%myK > g2%myK) then
             lessHaloListType = .false.
             return
           endif

           if(g1%myJ < g2%myJ) then
             lessHaloListType = .true.
             return
           else if(g1%myJ > g2%myJ) then
             lessHaloListType = .false.
             return
           endif

           if(g1%myI < g2%myI) then
             lessHaloListType = .true.
             return
           else if(g1%myI > g2%myI) then
             lessHaloListType = .false.
             return
           endif

           ! No need to compare anything else. G1 == g2.

         else boundary

           ! Both halo's are internal halo's, whose donor is stored on
           ! the same processor. Compare the donor blocks.

           if(g1%donorBlock < g2%donorBlock) then
             lessHaloListType = .true.
             return
           else if(g1%donorBlock > g2%donorBlock) then
             lessHaloListType = .false.
             return
           endif

           ! Also the blocks are identical. Compare the donor indices.
           ! First the k index.

           if(g1%dK < g2%dK) then
             lessHaloListType = .true.
             return
           else if(g1%dK > g2%dK) then
             lessHaloListType = .false.
             return
           endif

           ! The j index.

           if(g1%dJ < g2%dJ) then
             lessHaloListType = .true.
             return
           else if(g1%dJ > g2%dJ) then
             lessHaloListType = .false.
             return
           endif

           ! And the i index.

           if(g1%dI < g2%dI) then
             lessHaloListType = .true.
             return
           else if(g1%dI > g2%dI) then
             lessHaloListType = .false.
             return
           endif

           ! The donors are identical. Compare the halo's.
           ! First the block id.

           if(g1%myBlock < g2%myBlock) then
             lessHaloListType = .true.
             return
           else if(g1%myBlock > g2%myBlock) then
             lessHaloListType = .false.
             return
           endif

           ! Halo blocks are also identical. Finally compare the
           ! halo indices. Start with k.

           if(g1%myK < g2%myK) then
             lessHaloListType = .true.
             return
           else if(g1%myK > g2%myK) then
             lessHaloListType = .false.
             return
           endif

           ! The j index.

           if(g1%myJ < g2%myJ) then
             lessHaloListType = .true.
             return
           else if(g1%myJ > g2%myJ) then
             lessHaloListType = .false.
             return
           endif

           ! The i index.

           if(g1%myI < g2%myI) then
             lessHaloListType = .true.
             return
           else if(g1%myI > g2%myI) then
             lessHaloListType = .false.
             return
           endif

         endif boundary

         ! Both entities are identical.
         ! So set lessHaloListType to .false.

         lessHaloListType = .false.

         end function lessHaloListType

       end module haloList
