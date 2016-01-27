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
     integer(kind=intType) :: il, jl, kl

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

     ! Whether or not a cell is not possible to be a donor. Ie. a forceRecv
     integer(kind=intType), dimension(:, :, :), pointer :: invalidDonor

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
     logical :: allocated = .False.

     ! Flag if the real/int Buffers are ready after receiving info
     logical :: realBufferReady = .False. 
     logical :: intBufferReady = .False. 

  end type oversetBlock

  type oversetFringe

     ! Sizes
     integer(kind=intType) :: il, jl ,kl

     ! Buffer space for sending/receiving the fringes
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! These are the coordinate of what we are searching
     real(kind=realType), dimension(:, :), allocatable :: x

     ! qualaity is the best quality that has been found from a
     ! DONOR cell. It is initialized to large. 
     real(kind=realType), dimension(:), allocatable :: quality

     ! This is the information regarding where the cell came from. 
     integer(kind=intType), dimension(:), allocatable :: myBlock
     integer(kind=intType), dimension(:), allocatable :: myIndex

     ! This is the information about the donor that was found. Note we
     ! use dI, dJ, dK, short for donorI, etc.
     integer(kind=intType), dimension(:), allocatable :: donorProc
     integer(kind=intType), dimension(:), allocatable :: donorBlock
     integer(kind=intType), dimension(:), allocatable :: dI, dJ, dK

     real(kind=realType), dimension(:, :), allocatable  :: donorFrac

     ! gInd are the global indices of the donor cells. We will need
     ! these for forming the PC for the Newton Krylov solver
     integer(kind=intType), dimension(:, :), allocatable :: gInd

     ! Flag specifying if this cell is next to a wall
     integer(kind=intType), dimension(:), allocatable  :: isWall

     ! Flag if this set of fringes got allocated
     logical :: allocated = .False.

     ! Flag if the real/int Buffers are ready after receiving info
     logical :: realBufferReady = .False. 
     logical :: intBufferReady = .False. 

     ! The number of actual fringes that need to communicated. 
     integer(kind=intType) :: fringeReturnSize = 0

  end type oversetFringe

  type oversetWall

     ! Sizes
     integer(kind=intType) :: il, jl, kl, nNodes, nCells

     ! Buffer space for sending/receiving the fringes
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! Surface nodes used to build the tree:
     real(kind=realType), dimension(:, :), pointer :: x

     ! Connectivity for the surface
     integer(kind=intType), dimension(:, :), pointer:: conn

     ! The ADT for this block's wall(s)
     type(adtType) :: ADT

     ! Flag if the real/int Buffers are ready after receiving info
     logical :: realBufferReady = .False. 
     logical :: intBufferReady = .False. 

     ! Flag if this wall got allocated
     logical :: allocated = .False.

  end type oversetWall


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

  ! Several different MPI data types depending on the data we need to send

  ! This contains only the data needed to search for coordinates
  integer :: oversetMPISearchCoord

  ! This contains all MPI fringe data
  integer :: oversetMPIFringe

end module overset
