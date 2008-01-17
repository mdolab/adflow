!
!      ******************************************************************
!      *                                                                *
!      * File:          partitionMod.f90                                *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 02-05-2003                                      *
!      * Last modified: 11-07-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       module partitionMod
!
!      ******************************************************************
!      *                                                                *
!      * This local module contains definitions of derived datatypes    *
!      * as well as variables used in the partitioning directory.       *
!      *                                                                *
!      ******************************************************************
!
       use constants
       implicit none
       save
!
!      ******************************************************************
!      *                                                                *
!      * The definition of the derived datatype distributionBlockType   *
!      *                                                                *
!      ******************************************************************
!
       type distributionBlockType
!
!        ****************************************************************
!        *                                                              *
!        * Block dimensions and local block ID.                         *
!        *                                                              *
!        ****************************************************************
!
         !  nx, ny, nz - block integer dimensions for no halo cell based
         !               quantities.
         !  il, jl, kl - block integer dimensions for no halo node based
         !               quantities.
         !  blockID    - local block ID on the processor this block is
         !               stored.

         integer(kind=intType) :: nx, ny, nz, &
                                  il, jl, kl
         integer(kind=intType) :: blockID
!
!        ****************************************************************
!        *                                                              *
!        * Total number cells and faces inside the block. In the number *
!        * faces the work for nonmatching block boundaries is included, *
!        * such that the load balance is still guaranteed.              *
!        *                                                              *
!        ****************************************************************
!
         ! Ncell     : total number of cells in this block.
         ! Nface     : total number of faces in this block.

         integer(kind=intType) :: ncell, nface
!
!        ****************************************************************
!        *                                                              *
!        * Block boundary conditions.                                   *
!        *                                                              *
!        ****************************************************************
!
         !  nSubface             - Number of subfaces on this block.
         !  n1to1                - Number of 1 to 1 block boundaries.
         !  nBocos               - Number of physical boundary subfaces.
         !  BCType(:)            - Boundary condition type for each
         !                         subface. See the module BCTypes for
         !                         the possibilities.
         !  BCFaceID(:)          - Block face location of each subface.
         !                         Possible values are: iMin, iMax, jMin,
         !                         jMax, kMin, kMax.
         !  cgnsSubface(:)       - The subface in the corresponding cgns
         !                         block. As cgns distinguishes between
         !                         boundary and internal boundaries, the
         !                         BCType of the subface is needed to
         !                         know which one to take. A zero
         !                         indicates that this face was obtained
         !                         by splitting a cgns block.
         !  inBeg(:), inEnd(:)   - Lower and upper limits for the nodes
         !  jnBeg(:), jnEnd(:)     in each of the index directions on a
         !  knBeg(:), knEnd(:)     given subface. Note that one of these
         !                         indices will not change since we will
         !                         be moving on a face.
         !  dinBeg(:), dinEnd(:) - Lower and upper limits for the nodes
         !  djnBeg(:), djnEnd(:)   in the each of the index directions
         !  dknBeg(:), dknEnd(:)   of the donor subface for this
         !                         particular subface. Note that one of
         !                         these indices will not change since we
         !                         will be moving on a face.
         !  neighBlock(:)        - Block number to which this subface
         !                         connects. This value is set to zero if
         !                         this subface is not connected to
         !                         another block.
         !  l1(:), l2(:),        - Short hand for the transformation
         !  l3(:)                  matrix between this subface and the
         !                         neighbor block. These value are set to
         !                         zero if this subface is not connected
         !                         to another block.
         !  groupNum(:)          - Group number to which this subface
         !                         belongs. If this subface does not
         !                         belong to any group, the corresponding
         !                         entry in this array is zeroed out.

         integer(kind=intType) :: nSubface, n1to1, nBocos

         integer(kind=intType), dimension(:), pointer :: BCType
         integer(kind=intType), dimension(:), pointer :: BCFaceID
         integer(kind=intType), dimension(:), pointer :: cgnsSubface

         integer(kind=intType), dimension(:), pointer :: inBeg, inEnd
         integer(kind=intType), dimension(:), pointer :: jnBeg, jnEnd
         integer(kind=intType), dimension(:), pointer :: knBeg, knEnd

         integer(kind=intType), dimension(:), pointer :: dinBeg, dinEnd
         integer(kind=intType), dimension(:), pointer :: djnBeg, djnEnd
         integer(kind=intType), dimension(:), pointer :: dknBeg, dknEnd

         integer(kind=intType), dimension(:), pointer :: neighBlock
         integer(kind=intType), dimension(:), pointer :: l1, l2, l3
         integer(kind=intType), dimension(:), pointer :: groupNum
!
!        ****************************************************************
!        *                                                              *
!        * Overset boundary (fringe) cells and blanked cells.           *
!        *                                                              *
!        * Note that the local and donor indices and the interpolants   *
!        * are not stored here to conserve memory since the blocks have *
!        * not been partitioned yet. The info here is used to copy the  *
!        * necessary indices and interpolants into the flowDoms.        *
!        *                                                              *
!        ****************************************************************
!
         !  nCellsOverset       - total number of overset cells.
         !  cgnsOver(:)         - index into connOver for cgnsDoms.
         !  ipntOver(:)         - point number which is the index into
         !                         all data for each cell.
         !  neighOver(:)        - block number to which each donor
         !                         cell belongs.
         !  overComm(:,:)       - amount of overset communication to
         !                         to this block. The 1st dimension is
         !                         the block number of the donor and the
         !                         2nd is the number of cells
         !                         interpolating from that donor.

         integer(kind=intType) :: nCellsOverset

         integer(kind=intType), dimension(:), pointer :: cgnsOver
         integer(kind=intType), dimension(:), pointer :: ipntOver
         integer(kind=intType), dimension(:), pointer :: neighOver

         integer(kind=intType), dimension(:,:), pointer :: overComm
!
!        ****************************************************************
!        *                                                              *
!        * Relation to the original cgns grid.                          *
!        *                                                              *
!        ****************************************************************
!
         ! cgnsBlockID    - block/zone number of the cgns grid to which
         !                  this block is related.
         ! iBegOr, iEndOr - range of points of this block in the
         ! jBegOr, jEndOr   corresponding cgns block, i.e. for this block
         ! kBegOr, kEndOr   iBegOr <= i <= iEndOr, jBegOr <= j <= jEndOr,
         !                  kBegOr <= k <= kEndOr.
         !                  It is of course possible that the entire
         !                  block is stored.

         integer(kind=intType) :: cgnsBlockID
         integer(kind=intType) :: iBegOr, iEndOr, jBegOr, jEndOr
         integer(kind=intType) :: kBegOr, kEndOr

       end type distributionBlockType

       ! nBlocks        : Number of blocks to be distributed.
       ! blocks(nBlocks): The array with the block info.

       integer(kind=intType) :: nBlocks
       type(distributionBlockType), dimension(:), allocatable :: blocks
!
!      ******************************************************************
!      *                                                                *
!      * Type definition to store the way the original cgns blocks are  *
!      * split for load balancing reasons.                              *
!      *                                                                *
!      ******************************************************************
!
       type splitCGNSType

         ! nSubBlocks:             # of subblocks into which the original
         !                         block is split.
         ! ranges(nsubblocks,3,2): The nodal range in the original
         !                         block of these subblocks.

         integer(kind=intType) :: nSubBlocks
         integer(kind=intType), dimension(:,:,:), pointer :: ranges

       end type splitCGNSType
!
!      ******************************************************************
!      *                                                                *
!      * Type definition needed to determine the processor ID's and     *
!      * nodal ranges of the subblocks for every CGNS block.            *
!      *                                                                *
!      ******************************************************************
!
       type subblocksOfCGNSType

         ! cgnsBlockID    - block/zone number of the cgns grid to which
         !                  this subblock is related.
         ! procID         - Processor ID on which this subblock is
         !                  stored.
         ! blockID        - local block ID on the processor this block is
         !                  stored.
         ! iBegOr, iEndOr - range of points of this block in the
         ! jBegOr, jEndOr   corresponding cgns block, i.e. for this block
         ! kBegOr, kEndOr   iBegOr <= i <= iEndOr, jBegOr <= j <= jEndOr,
         !                  kBegOr <= k <= kEndOr.
         !                  It is of course possible that the entire
         !                  block is stored.

         integer :: cgnsBlockID
         integer :: procID
         integer :: blockID
         integer :: iBegOr, iEndOr, jBegOr, jEndOr, kBegOr, kEndOr

       end type subblocksOfCGNSType

       ! Interface for the extension of the operators <= and <.
       ! These are needed for the sorting of subblocksOfCGNSType. Note
       ! that the = operator does not need to be defined, because
       ! subblocksOfCGNSType only contains primitive types.

       interface operator(<=)
         module procedure lessEqualSubblocksOfCGNSType
       end interface

       interface operator(<)
         module procedure lessSubblocksOfCGNSType
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Type definition needed to determine the number of distinct     *
!      * non-matching abutting subfaces in the CGNS file.               *
!      *                                                                *
!      ******************************************************************
!
       type subfaceNonMatchType

         ! iBeg, iEnd - Nodal subface range om the CGNS block, i-direction
         ! jBeg, jEnd - Idem in the j-direction
         ! kBeg, kEnd - Idem in the k-direction
         ! connID     - The cgns connectivity ID.

         integer :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd
         integer :: connID

       end type subfaceNonMatchType

       ! Interface for the extension of the operators <= and <.
       ! These are needed for the sorting of subfaceNonMatchType. Note
       ! that the = operator does not need to be defined, because
       ! subfaceNonMatchType only contains primitive types.

       interface operator(<=)
         module procedure lessEqualSubfaceNonMatchType
       end interface

       interface operator(<)
         module procedure lessSubfaceNonMatchType
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Variable to store the partition number (processor ID) of the   *
!      * computational blocks.                                          *
!      *                                                                *
!      ******************************************************************
!
       ! ubvec(2):      Tolerance for the constraints.
       ! part(nBlocks): The processor ID for each block, starting at 0.

       real, dimension(2) :: ubvec

       integer(kind=intType), dimension(:), allocatable :: part
!
!      ******************************************************************
!      *                                                                *
!      * Variables needed for the reading of the grid files.            *
!      *                                                                *
!      ******************************************************************
!
       ! nGridsRead:            Number of grids to read.
       ! fileIDs(nGridsRead):   The file ID's.
       ! gridFiles(nGridsRead): Names of the grid files to read.
       ! interpolSpectral:      Whether or not to interpolate the 
       !                        coordinates for the time spectral mode.

       integer(kind=intType) :: nGridsRead
       logical ::               interpolSpectral

       integer, dimension(:), allocatable :: fileIDs

       character(len=maxStringLen), dimension(:), allocatable :: gridFiles
!
!      ******************************************************************
!      *                                                                *
!      * Variables only needed for the reading of Plot3D files.         *
!      *                                                                *
!      ******************************************************************
!
       ! byteSwapGrids(nGridsRead):    Whether or not byte swapping must
       !                               be applied for the grids to be
       !                               read. Note that different 
       !                               endianness of the grid files is
       !                               allowed.
       ! blockFormatGrids(nGridsRead): What type of grid is stored, 
       !                               single block or multi-block.
       !                               Again every grid can be different.

       integer(kind=intType), dimension(:), allocatable :: &
                                                        blockFormatGrids
       logical, dimension(:), allocatable :: byteSwapGrids

       contains
!
!        ****************************************************************
!        *                                                              *
!        * Functions to simulate the operators <= and < for the derived *
!        * datatypes subblocksOfCGNSType and subfaceNonMatchType.       *
!        *                                                              *
!        ****************************************************************
!
         logical function lessEqualSubblocksOfCGNSType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 <= g2 and .false.         *
!        * otherwise. The comparison is firstly based on the CGNS block *
!        * ID, then the processor ID and finally the local block ID.    *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(subblocksOfCGNSType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Comparison of the CGNS block ID. If not equal, set
         ! lessEqualSubblocksOfCGNSType appropriately and return.

         if(g1%cgnsBlockID < g2%cgnsBlockID) then
           lessEqualSubblocksOfCGNSType = .true.
           return
         else if(g1%cgnsBlockID > g2%cgnsBlockID) then
           lessEqualSubblocksOfCGNSType = .false.
           return
         endif

         ! Compare the processor ID's.

         if(g1%procID < g2%procID) then
           lessEqualSubblocksOfCGNSType = .true.
           return
         else if(g1%procID > g2%procID) then
           lessEqualSubblocksOfCGNSType = .false.
           return
         endif

         ! Compare the local block ID's.

         if(g1%blockID < g2%blockID) then
           lessEqualSubblocksOfCGNSType = .true.
           return
         else if(g1%blockID > g2%blockID) then
           lessEqualSubblocksOfCGNSType = .false.
           return
         endif

         ! It does not make sense to compare the subranges, because
         ! these are identical if we came so far. Both entities are
         ! identical. So set lessEqualSubblocksOfCGNSType to .true.

         lessEqualSubblocksOfCGNSType = .true.

         end function lessEqualSubblocksOfCGNSType

         !===============================================================

         logical function lessSubblocksOfCGNSType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 < g2 and .false.          *
!        * otherwise. The comparison is firstly based on the CGNS block *
!        * ID, then the processor ID and finally the local blockID.     *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(subblocksOfCGNSType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Comparison of the CGNS block ID. If not equal, set
         ! lessSubblocksOfCGNSType appropriately and return.

         if(g1%cgnsBlockID < g2%cgnsBlockID) then
           lessSubblocksOfCGNSType = .true.
           return
         else if(g1%cgnsBlockID > g2%cgnsBlockID) then
           lessSubblocksOfCGNSType = .false.
           return
         endif

         ! Compare the processor ID's.

         if(g1%procID < g2%procID) then
           lessSubblocksOfCGNSType = .true.
           return
         else if(g1%procID > g2%procID) then
           lessSubblocksOfCGNSType = .false.
           return
         endif

         ! Compare the local block ID's.

         if(g1%blockID < g2%blockID) then
           lessSubblocksOfCGNSType = .true.
           return
         else if(g1%blockID > g2%blockID) then
           lessSubblocksOfCGNSType = .false.
           return
         endif

         ! It does not make sense to compare the subranges, because
         ! these are identical if we came so far. Both entities are
         ! identical. So set lessSubblocksOfCGNSType to .false.

         lessSubblocksOfCGNSType = .false.

         end function lessSubblocksOfCGNSType

         !===============================================================

         logical function lessEqualSubfaceNonMatchType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 <= g2 and .false.         *
!        * otherwise. The comparison is firstly based on the i-range,   *
!        * followed by the j-range and k-range. If these are all the    *
!        * same the connectivity ID is compared.                        *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(subfaceNonMatchType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Comparison of the iBeg value. If different set 
         ! lessEqualSubfaceNonMatchType appropriately and return.

         if(g1%iBeg < g2%iBeg) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%iBeg > g2%iBeg) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The iEnd value.

         if(g1%iEnd < g2%iEnd) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%iEnd > g2%iEnd) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The jBeg value.

         if(g1%jBeg < g2%jBeg) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%jBeg > g2%jBeg) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The jEnd value.

         if(g1%jEnd < g2%jEnd) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%jEnd > g2%jEnd) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The kBeg value.

         if(g1%kBeg < g2%kBeg) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%kBeg > g2%kBeg) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The kEnd value.

         if(g1%kEnd < g2%kEnd) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%kEnd > g2%kEnd) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! The subrange is identical. Compare the connectivity ID.

         if(g1%connID < g2%connID) then
           lessEqualSubfaceNonMatchType = .true.
           return
         else if(g1%connID > g2%connID) then
           lessEqualSubfaceNonMatchType= .false.
           return
         endif

         ! g1 and g2 are identical. Return .true.

         lessEqualSubfaceNonMatchType = .true.

         end function lessEqualSubfaceNonMatchType

         !===============================================================

         logical function lessSubfaceNonMatchType(g1, g2)
!
!        ****************************************************************
!        *                                                              *
!        * This function returns .true. if g1 < g2 and .false.          *
!        * otherwise. The comparison is firstly based on the i-range,   *
!        * followed by the j-range and k-range. If these are all the    *
!        * same the connectivity ID is compared.                        *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Function arguments.
!
         type(subfaceNonMatchType), intent(in) :: g1, g2
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
         ! Comparison of the iBeg value. If different set 
         ! lessEqualSubfaceNonMatchType appropriately and return.

         if(g1%iBeg < g2%iBeg) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%iBeg > g2%iBeg) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The iEnd value.

         if(g1%iEnd < g2%iEnd) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%iEnd > g2%iEnd) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The jBeg value.

         if(g1%jBeg < g2%jBeg) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%jBeg > g2%jBeg) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The jEnd value.

         if(g1%jEnd < g2%jEnd) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%jEnd > g2%jEnd) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The kBeg value.

         if(g1%kBeg < g2%kBeg) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%kBeg > g2%kBeg) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The kEnd value.

         if(g1%kEnd < g2%kEnd) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%kEnd > g2%kEnd) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! The subrange is identical. Compare the connectivity ID.

         if(g1%connID < g2%connID) then
           lessSubfaceNonMatchType = .true.
           return
         else if(g1%connID > g2%connID) then
           lessSubfaceNonMatchType= .false.
           return
         endif

         ! g1 and g2 are identical. Return .false.

         lessSubfaceNonMatchType = .false.

         end function lessSubfaceNonMatchType

       end module partitionMod
