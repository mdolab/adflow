module partitionMod
  !
  !       This local module contains definitions of derived datatypes    
  !       as well as variables used in the partitioning directory.       
  !
  use constants
  implicit none
  save
  !
  !       The definition of the derived datatype distributionBlockType   
  !
  type distributionBlockType
     !
     !         Block dimensions and local block ID.                         
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
     !         Total number cells and faces inside the block. In the number 
     !         faces the work for nonmatching block boundaries is included, 
     !         such that the load balance is still guaranteed.              
     !
     ! Ncell     : total number of cells in this block.
     ! Nface     : total number of faces in this block.

     integer(kind=intType) :: ncell, nface
     !
     !         Block boundary conditions.                                   
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
     !         Relation to the original cgns grid.                          
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
  !       Type definition to store the way the original cgns blocks are  
  !       split for load balancing reasons.                              
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
  !       Type definition needed to determine the processor ID's and     
  !       nodal ranges of the subblocks for every CGNS block.            
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
  end interface operator(<=)

  interface operator(<)
     module procedure lessSubblocksOfCGNSType
  end interface operator(<)
  !
  !       Type definition needed to determine the number of distinct     
  !       non-matching abutting subfaces in the CGNS file.               
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
  end interface operator(<=)

  interface operator(<)
     module procedure lessSubfaceNonMatchType
  end interface operator(<)

  type sortSubRangeType

     ! iMin: minimum i-index in the subrange.
     ! jMin: minimum j-index in the subrange.
     ! kMin: minimum k-index in the subrange.
     ! iMax: maximum i-index in the subrange.
     ! jMax: maximum j-index in the subrange.
     ! kMax: maximum k-index in the subrange.

     integer(kind=intType) :: iMin, jMin, kMin
     integer(kind=intType) :: iMax, jMax, kMax

  end type sortSubRangeType

  ! Interfaces for the definitions of the operators <=, < and /=.
  ! These are needed for the sorting of this derived data type.
  ! Note that the = operator does not need to be defined, because
  ! sortSubRangeType only contains primitive types.

  interface operator(<=)
     module procedure lessEqualSortSubRangeType
  end interface operator(<=)

  interface operator(<)
     module procedure lessSortSubRangeType
  end interface operator(<)

  ! Definition of the derived data type.

  type fourIntPlusRealType
     integer(kind=intType) :: n1, n2, n3, n4
#ifdef USE_COMPLEX
     complex(kind=realType)   :: dist
#else
     real(kind=realType)   :: dist
#endif
  end type fourIntPlusRealType

  ! Interfaces for the definitions of the operators <=, < and /=.
  ! These are needed for the sorting of this derived data type.
  ! Note that the = operator does not need to be defined, because
  ! fourIntPlusRealType only contains primitive types.

  interface operator(<=)
     module procedure lessEqualFourIntPlusRealType
  end interface operator(<=)

  interface operator(<)
     module procedure lessFourIntPlusRealType
  end interface operator(<)

  interface operator(/=)
     module procedure notEqualFourIntPlusRealType
  end interface operator(/=)

  ! ==========================================================================
  !
  !       Variable to store the partition number (processor ID) of the   
  !       computational blocks.                                          
  !
  ! ubvec(2):      Tolerance for the constraints.
  ! part(nBlocks): The processor ID for each block, starting at 0.

  real, dimension(2) :: ubvec

  integer(kind=intType), dimension(:), allocatable :: part
  !
  !       Variables needed for the reading of the grid files.            
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
  ! ==========================================================================


contains
  !
  !         Functions to simulate the operators <= and < for the derived 
  !         datatypes subblocksOfCGNSType and subfaceNonMatchType.       
  !
  logical function lessEqualSubblocksOfCGNSType(g1, g2)
    !
    !         This function returns .true. if g1 <= g2 and .false.         
    !         otherwise. The comparison is firstly based on the CGNS block 
    !         ID, then the processor ID and finally the local block ID.    
    !
    implicit none
    !
    !        Function arguments.
    !
    type(subblocksOfCGNSType), intent(in) :: g1, g2

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
    !         This function returns .true. if g1 < g2 and .false.          
    !         otherwise. The comparison is firstly based on the CGNS block 
    !         ID, then the processor ID and finally the local blockID.     
    !
    implicit none
    !
    !        Function arguments.
    !
    type(subblocksOfCGNSType), intent(in) :: g1, g2

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
    !         This function returns .true. if g1 <= g2 and .false.         
    !         otherwise. The comparison is firstly based on the i-range,   
    !         followed by the j-range and k-range. If these are all the    
    !         same the connectivity ID is compared.                        
    !
    implicit none
    !
    !        Function arguments.
    !
    type(subfaceNonMatchType), intent(in) :: g1, g2

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
    !         This function returns .true. if g1 < g2 and .false.          
    !         otherwise. The comparison is firstly based on the i-range,   
    !         followed by the j-range and k-range. If these are all the    
    !         same the connectivity ID is compared.                        
    !
    implicit none
    !
    !        Function arguments.
    !
    type(subfaceNonMatchType), intent(in) :: g1, g2

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

  subroutine qsortSubblocksOfCGNSType(arr, nn)
    !
    !       qsortSubblocksOfCGNSType sorts the array of the derived        
    !       datatype subblocksOfCGNSType in increasing order based on the  
    !       <= operator for this derived data type.                        
    !
    use constants
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(subblocksOfCGNSType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(subblocksOfCGNSType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortSubblocksOfCGNSType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubblocksOfCGNSType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubblocksOfCGNSType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubblocksOfCGNSType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubblocksOfCGNSType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortSubblocksOfCGNSType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                        &
               call terminate("qsortSubblocksOfCGNSType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortSubblocksOfCGNSType

  subroutine qsortSubfaceNonMatchType(arr, nn)
    !
    !       qsortSubfaceNonMatchType sorts the array of the derived        
    !       datatype subfaceNonMatchType in increasing order based on the  
    !       <= operator for this derived data type.                        
    !
    use constants
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(subfaceNonMatchType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(subfaceNonMatchType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortSubfaceNonMatchType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubfaceNonMatchType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubfaceNonMatchType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubfaceNonMatchType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortSubfaceNonMatchType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortSubfaceNonMatchType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                        &
               call terminate("qsortSubfaceNonMatchType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortSubfaceNonMatchType

  logical function lessEqualSortSubRangeType(g1, g2)
    !
    !         lessEqualSortSubRangeType defines the operator <= for the    
    !         derived datatype sortSubRangeType. The comparison is first   
    !         based on kMin, followed by jMin and finally iMin.            
    !         The comparison is therefore not based on the max values.     
    !
    implicit none
    !
    !        Function arguments.
    !
    type(sortSubRangeType), intent(in) :: g1, g2
    !
    !         Begin executation.                                           
    !
    ! Compare the kMin index and return .true. or .false. if they
    ! differ.

    if(g1%kMin < g2%kMin) then
       lessEqualSortSubRangeType = .true.
       return
    else if(g1%kMin > g2%kMin) then
       lessEqualSortSubRangeType = .false.
       return
    endif

    ! kMin indices are equal. Compare the jMin's.

    if(g1%jMin < g2%jMin) then
       lessEqualSortSubRangeType = .true.
       return
    else if(g1%jMin > g2%jMin) then
       lessEqualSortSubRangeType = .false.
       return
    endif

    ! Also the jMin's are equal. Compare iMin's.

    if(g1%iMin < g2%iMin) then
       lessEqualSortSubRangeType = .true.
       return
    else if(g1%iMin > g2%iMin) then
       lessEqualSortSubRangeType = .false.
       return
    endif

    ! g1 equals g2. Return .true.

    lessEqualSortSubRangeType = .true.

  end function lessEqualSortSubRangeType

  !===============================================================

  logical function lessSortSubRangeType(g1, g2)
    !
    !         lessSortSubRangeType defines the operator < for the derived  
    !         datatype sortSubRangeType. The comparison is first based on  
    !         kMin, followed by jMin and finally iMin.                     
    !         The comparison is therefore not based on the max values.     
    !
    implicit none
    !
    !        Function arguments.
    !
    type(sortSubRangeType), intent(in) :: g1, g2
    !
    !         Begin executation.                                           
    !
    ! Compare the kMin index and return .true. or .false. if they
    ! differ.

    if(g1%kMin < g2%kMin) then
       lessSortSubRangeType = .true.
       return
    else if(g1%kMin > g2%kMin) then
       lessSortSubRangeType = .false.
       return
    endif

    ! kMin indices are equal. Compare the jMin's.

    if(g1%jMin < g2%jMin) then
       lessSortSubRangeType = .true.
       return
    else if(g1%jMin > g2%jMin) then
       lessSortSubRangeType = .false.
       return
    endif

    ! Also the jMin's are equal. Compare iMin's.

    if(g1%iMin < g2%iMin) then
       lessSortSubRangeType = .true.
       return
    else if(g1%iMin > g2%iMin) then
       lessSortSubRangeType = .false.
       return
    endif

    ! g1 equals g2. Return .false.

    lessSortSubRangeType = .false.

  end function lessSortSubRangeType

  !      ==================================================================

  subroutine sortRangesSplitInfo(splitInfo)
    !
    !       sortRangesSplitInfo sort the ranges of the given subblocks in  
    !       increasing order such that a unique ordering is obtained,      
    !       independent of the history of the splitting.                   
    !
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    type(splitCGNSType), intent(inout) :: splitInfo
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, nSubBlocks

    type(sortSubRangeType), dimension(splitInfo%nSubBlocks) :: subRanges

    ! Copy the subface range from splitInfo into subRanges.

    nSubBlocks = splitInfo%nSubBlocks

    do i=1,nSubBlocks
       subRanges(i)%iMin = splitInfo%ranges(i,1,1)
       subRanges(i)%jMin = splitInfo%ranges(i,2,1)
       subRanges(i)%kMin = splitInfo%ranges(i,3,1)

       subRanges(i)%iMax = splitInfo%ranges(i,1,2)
       subRanges(i)%jMax = splitInfo%ranges(i,2,2)
       subRanges(i)%kMax = splitInfo%ranges(i,3,2)
    enddo

    ! Sort subRanges in increasing order.

    call qsortSortSubRangeType(subRanges, nSubBlocks)

    ! Copy the data back into splitInfo.

    do i=1,nSubBlocks
       splitInfo%ranges(i,1,1) = subRanges(i)%iMin
       splitInfo%ranges(i,2,1) = subRanges(i)%jMin
       splitInfo%ranges(i,3,1) = subRanges(i)%kMin

       splitInfo%ranges(i,1,2) = subRanges(i)%iMax
       splitInfo%ranges(i,2,2) = subRanges(i)%jMax
       splitInfo%ranges(i,3,2) = subRanges(i)%kMax
    enddo

  end subroutine sortRangesSplitInfo

  !      ==================================================================

  subroutine qsortSortSubRangeType(arr, nn)
    !
    !       qsortSortSubRangeType sorts the given number of halo's in      
    !       increasing order based on the <= operator for this derived     
    !       data type.                                                     
    !
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(sortSubRangeType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(sortSubRangeType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                             &
         call terminate("qsortSortSubRangeType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                             &
                  call terminate("qsortSortSubRangeType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                             &
                  call terminate("qsortSortSubRangeType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                             &
                  call terminate("qsortSortSubRangeType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                             &
                  call terminate("qsortSortSubRangeType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                             &
         call terminate("qsortSortSubRangeType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                     &
               call terminate("qsortSortSubRangeType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortSortSubRangeType

  !
  !         Functions to define the operators <, <= and /=.              
  !         Note that the comparison is only based on the integers.      
  !         The real contains additional info, the maximum deviation,    
  !         which is normally different even if the subfaces are         
  !         identical.                                                   
  !
  logical function lessEqualFourIntPlusRealType(g1, g2)
    implicit none
    type(fourIntPlusRealType), intent(in) :: g1, g2

    ! Compare the first element.

    if(g1%n1 < g2%n1) then
       lessEqualFourIntPlusRealType = .true.
       return
    else if(g1%n1 > g2%n1) then
       lessEqualFourIntPlusRealType = .false.
       return
    endif

    ! Compare the second element.

    if(g1%n2 < g2%n2) then
       lessEqualFourIntPlusRealType = .true.
       return
    else if(g1%n2 > g2%n2) then
       lessEqualFourIntPlusRealType = .false.
       return
    endif

    ! Compare the third element.

    if(g1%n3 < g2%n3) then
       lessEqualFourIntPlusRealType = .true.
       return
    else if(g1%n3 > g2%n3) then
       lessEqualFourIntPlusRealType = .false.
       return
    endif

    ! Compare the fourth element.

    if(g1%n4 < g2%n4) then
       lessEqualFourIntPlusRealType = .true.
       return
    else if(g1%n4 > g2%n4) then
       lessEqualFourIntPlusRealType = .false.
       return
    endif

    ! g1 equals g2. Return .true.

    lessEqualFourIntPlusRealType = .true.

  end function lessEqualFourIntPlusRealType

  !===============================================================

  logical function lessFourIntPlusRealType(g1, g2)
    implicit none
    type(fourIntPlusRealType), intent(in) :: g1, g2

    ! Compare the first element.

    if(g1%n1 < g2%n1) then
       lessFourIntPlusRealType = .true.
       return
    else if(g1%n1 > g2%n1) then
       lessFourIntPlusRealType = .false.
       return
    endif

    ! Compare the second element.

    if(g1%n2 < g2%n2) then
       lessFourIntPlusRealType = .true.
       return
    else if(g1%n2 > g2%n2) then
       lessFourIntPlusRealType = .false.
       return
    endif

    ! Compare the third element.

    if(g1%n3 < g2%n3) then
       lessFourIntPlusRealType = .true.
       return
    else if(g1%n3 > g2%n3) then
       lessFourIntPlusRealType = .false.
       return
    endif

    ! Compare the fourth element.

    if(g1%n4 < g2%n4) then
       lessFourIntPlusRealType = .true.
       return
    else if(g1%n4 > g2%n4) then
       lessFourIntPlusRealType = .false.
       return
    endif

    ! g1 equals g2. Return .false.

    lessFourIntPlusRealType = .false.

  end function lessFourIntPlusRealType

  !===============================================================

  logical function notEqualFourIntPlusRealType(g1, g2)
    implicit none
    type(fourIntPlusRealType), intent(in) :: g1, g2

    notEqualFourIntPlusRealType = .true.
    if(g1%n1 == g2%n1 .and. g1%n2 == g2%n2 .and. &
         g1%n3 == g2%n3 .and. g1%n4 == g2%n4) &
         notEqualFourIntPlusRealType = .false.

  end function notEqualFourIntPlusRealType


  !      ==================================================================

  subroutine sortBadEntities(nEntities, entities, dist, sortDist)
    !
    !       sortBadEntities sorts the given number of entities in          
    !       increasing order and gets rid of the multiple entries.         
    !
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(inout) :: nEntities
    integer(kind=intType), dimension(4,*), intent(inout) :: entities
#ifdef USE_COMPLEX
    complex(kind=realType),   dimension(*),   intent(inout) :: dist
#else
    real(kind=realType),   dimension(*),   intent(inout) :: dist
#endif

    logical, intent(in) :: sortDist
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm

    type(fourIntPlusRealType), dimension(nEntities) :: tmp

    ! Return immediately if there are no entities to be sorted.

    if(nEntities == 0) return

    ! Copy the info into tmp. If the distances must be sorted as
    ! well, copy the info. Otherwise simply put zero.

    do nn=1,nEntities
       tmp(nn)%n1   = entities(1,nn)
       tmp(nn)%n2   = entities(2,nn)
       tmp(nn)%n3   = entities(3,nn)
       tmp(nn)%n4   = entities(4,nn)
       if( sortDist ) then
          tmp(nn)%dist = dist(nn)
       else
          tmp(nn)%dist = zero
       endif
    enddo

    ! Sort tmp in increasing order.

    call qsortFourIntPlusRealType(tmp, nEntities)

    ! Get rid of the multiple entries. Note that the exceptional
    ! case of zero entities does not to be considered, because in
    ! that case this part of the subroutine is not executed.
    ! If multiple entries are present the distance is taken as
    ! the maximum of the two.

    mm = 1
    do nn=2,nEntities
       if(tmp(nn) /= tmp(mm)) then
          mm = mm + 1
          tmp(mm) = tmp(nn)
       else
          tmp(mm)%dist = max(tmp(mm)%dist, tmp(nn)%dist)
       endif
    enddo

    ! Copy the data back info entities and dist. The latter
    ! only if the distances should be sorted as well.

    nEntities = mm

    do nn=1,nEntities
       entities(1,nn) = tmp(nn)%n1
       entities(2,nn) = tmp(nn)%n2
       entities(3,nn) = tmp(nn)%n3
       entities(4,nn) = tmp(nn)%n4
       if( sortDist ) dist(nn) = tmp(nn)%dist
    enddo

  end subroutine sortBadEntities

  !      ==================================================================

  subroutine qsortFourIntPlusRealType(arr, nn)
    !
    !       qsortFourIntPlusRealType sorts the given number of halo's in   
    !       increasing order based on the <= operator for this derived     
    !       data type.                                                     
    !
    use constants
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nn

    type(fourIntPlusRealType), dimension(*), intent(inout) :: arr
    !
    !      Local variables.
    !
    integer(kind=intType), parameter :: m = 7

    integer(kind=intType) :: nStack
    integer(kind=intType) :: i, j, k, r, l, jStack, ii

    integer :: ierr

    type(fourIntPlusRealType) :: a, tmp

    integer(kind=intType), allocatable, dimension(:) :: stack
    integer(kind=intType), allocatable, dimension(:) :: tmpStack

    ! Allocate the memory for stack.

    nStack = 100
    allocate(stack(nStack), stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortFourIntPlusRealType", &
         "Memory allocation failure for stack")

    ! Initialize the variables that control the sorting.

    jStack = 0
    l      = 1
    r      = nn

    ! Start of the algorithm

    do

       ! Check for the size of the subarray.

       if((r-l) < m) then

          ! Perform insertion sort

          do j=l+1,r
             a = arr(j)
             do i=(j-1),l,-1
                if(arr(i) <= a) exit
                arr(i+1) = arr(i)
             enddo
             arr(i+1) = a
          enddo

          ! In case there are no more elements on the stack, exit from
          ! the outermost do-loop. Algorithm has finished.

          if(jStack == 0) exit

          ! Pop stack and begin a new round of partitioning.

          r = stack(jStack)
          l = stack(jStack-1)
          jStack = jStack - 2

       else

          ! Subarray is larger than the threshold for a linear sort.
          ! Choose median of left, center and right elements as
          ! partitioning element a.
          ! Also rearrange so that (l) <= (l+1) <= (r).

          k = (l+r)/2
          tmp      = arr(k)             ! Swap the elements
          arr(k)   = arr(l+1)           ! k and l+1.
          arr(l+1) = tmp

          if(arr(r) < arr(l)) then
             tmp    = arr(l)             ! Swap the elements
             arr(l) = arr(r)             ! r and l.
             arr(r) = tmp
          endif

          if(arr(r) < arr(l+1)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(r)           ! r and l+1.
             arr(r)   = tmp
          endif

          if(arr(l+1) < arr(l)) then
             tmp      = arr(l+1)         ! Swap the elements
             arr(l+1) = arr(l)           ! l and l+1.
             arr(l)   = tmp
          endif

          ! Initialize the pointers for partitioning.

          i = l+1
          j = r
          a = arr(l+1)

          ! The innermost loop

          do

             ! Scan up to find element >= a.
             do
                i = i+1
                if(a <= arr(i)) exit
             enddo

             ! Scan down to find element <= a.
             do
                j = j-1
                if(arr(j) <= a) exit
             enddo

             ! Exit the loop in case the pointers i and j crossed.

             if(j < i) exit

             ! Swap the element i and j.

             tmp    = arr(i)
             arr(i) = arr(j)
             arr(j) = tmp
          enddo

          ! Swap the entries j and l+1. Remember that a equals
          ! arr(l+1).

          arr(l+1) = arr(j)
          arr(j)   = a

          ! Push pointers to larger subarray on stack,
          ! process smaller subarray immediately.

          jStack = jStack + 2
          if(jStack > nStack) then

             ! Storage of the stack is too small. Reallocate.

             allocate(tmpStack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortFourIntPlusRealType", &
                  "Memory allocation error for tmpStack")
             tmpStack = stack

             ! Free the memory of stack, store the old value of nStack
             ! in tmp and increase nStack.

             deallocate(stack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortFourIntPlusRealType", &
                  "Deallocation error for stack")
             ii = nStack
             nStack = nStack + 100

             ! Allocate the memory for stack and copy the old values
             ! from tmpStack.

             allocate(stack(nStack), stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortFourIntPlusRealType", &
                  "Memory reallocation error for stack")
             stack(1:ii) = tmpStack(1:ii)

             ! And finally release the memory of tmpStack.

             deallocate(tmpStack, stat=ierr)
             if(ierr /= 0)                                &
                  call terminate("qsortFourIntPlusRealType", &
                  "Deallocation error for tmpStack")
          endif

          if((r-i+1) >= (j-l)) then
             stack(jStack)   = r
             r               = j-1
             stack(jStack-1) = j
          else
             stack(jStack)   = j-1
             stack(jStack-1) = l
             l               = j
          endif

       endif
    enddo

    ! Release the memory of stack.

    deallocate(stack, stat=ierr)
    if(ierr /= 0)                                &
         call terminate("qsortFourIntPlusRealType", &
         "Deallocation error for stack")

    ! Check in debug mode whether the array is really sorted.

    if( debug ) then
       do i=1,(nn-1)
          if(arr(i+1) < arr(i))                        &
               call terminate("qsortFourIntPlusRealType", &
               "Array is not sorted correctly")
       enddo
    endif

  end subroutine qsortFourIntPlusRealType

end module partitionMod
