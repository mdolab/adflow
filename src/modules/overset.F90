module overset

  use precision
  use adtData
  use block
  implicit none
  
  ! Helper dataType for communicated overset grid points. This data
  ! structure mirrros the blockType structure in block.F90, but only
  ! contains minimum amount of information required for computing
  ! overset connectivities. 

  ! Store the coordinates from a block that will need to be searched.

  ! A simple generic sparse matrix storage container for storing the
  ! (sparse) overlap structure of an overset mesh
  type CSRMatrix
     integer(kind=intType) :: nRow, nCol, nnz, nnzLocal
     integer(kind=intType), dimension(:), pointer :: colInd, rowPtr
     real(kind=realType), dimension(:), pointer :: data
     integer(Kind=intType), dimension(:), pointer :: assignedProc

     ! Flag if this matrix is allocated
     logical :: allocated=.False.

  end type CSRMatrix

  ! This derived type contains sufficient information to perfom ADT
  ! donor searches. The information is com
  type oversetBlock

     ! Sizes for the block
     integer(kind=intType) :: ib, jb, kb
     integer(kind=intType) :: ie, je, ke
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: nx, ny, nz
     
     ! This is the cell volume of the donor
     real(kind=realType), dimension(:, :), pointer :: qualDonor

     ! Connectivity for the ADT
     integer(kind=intType), dimension(:, :), pointer :: hexaConn

     ! Coordinates for the ADT
     real(kind=realType), dimension(:, :), pointer :: xADT

     ! Copy of global cell
     integer(kind=intType), dimension(:, :, :), pointer :: globalCell

     ! Whether or not a cell is "near" a wall
     integer(kind=intType), dimension(:, :, :), pointer :: nearWall
     
     ! Minimum volume for this block
     real(kind=realType) :: minVol

     ! The ADT for this block
     type(adtType) :: ADT

     ! The processor for this block
     integer(kind=intType) :: proc
     
     ! And the local block index
     integer(kind=intType) :: block

     ! Buffer space for sending/receiving the block
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! Flag if this block got allocated
     logical :: allocated=.False.

  end type oversetBlock

  type fringeListType
     type(fringeType), dimension(:), allocatable :: arr
  end type fringeListType

  ! These are the two main lists of derived types used for overset
  ! assembly
  type(oversetBlock), dimension(:), allocatable :: oBlocks
  type(fringeListType), dimension(:), allocatable :: fringeList

  ! This is the flattened list of the fringes next to the wall that we
  !  have actually found donors for.
  ! tmpFringePtr is only used if we need to realloc. 
  type(fringeType), dimension(:), pointer :: localWallFringes, wallFringes, tmpFringePtr
  integer(kind=intType) :: nLocalWallFringe, nWallFringe

  ! A receive buffer for fringes
  type(fringeType), dimension(:), allocatable :: tmpFringes

  ! These are the master overlap matrices
  type(CSRMatrix), dimension(:, :), allocatable, target :: overlapMatrix

  ! Some additional helper stuff
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType) :: nDomTotal

  ! The new fringe datatype
  integer :: oversetMPIFringe

end module overset
