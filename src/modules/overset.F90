module overset

  use precision
  use adtData
  use block
  use kdtree2_module
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
     integer(kind=intType) :: il, jl, kl
     integer(kind=intType) :: nNodes=0
     integer(kind=intType) :: nCells=0
     integer(kind=intType) :: maxCells=0
     integer(kind=intType) :: cluster=0

     ! Buffer space for sending/receiving the fringes
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! Surface nodes used to build the tree:
     real(kind=realType), dimension(:, :), pointer :: x

     ! Connectivity for the surface
     integer(kind=intType), dimension(:, :), pointer:: conn
     
     ! Blanking values for Nodes
     integer(kind=intType), dimension(:), allocatable :: iBlank
     integer(kind=intType), dimension(:), allocatable :: cellPtr

     ! Node to element array
     integer(kind=intType), dimension(:, :), allocatable :: nte

     ! The ADT for this block's wall(s)
     type(adtType) :: ADT

     ! This KDTree for this block's wall
     type(kdtree2), pointer :: tree

     ! Flag if the real/int Buffers are ready after receiving info
     logical :: realBufferReady = .False. 
     logical :: intBufferReady = .False. 

     ! Flag if this wall got allocated
     logical :: allocated = .False.

  end type oversetWall

  type oversetString
     ! This is a generic type defining a string list. It may be used
     ! as both a "parent" or a "child". 

     ! Sizes
     integer(kind=intType) :: nNodes, nElems

     ! My String's index
     integer(kind=intType) :: myID

     ! Whether or no this string is periodic
     logical :: isPeriodic=.False.

     ! The actual physical node locations. Size (3, nNodes)
     real(kind=realType), dimension(:, :), pointer :: x

     ! The (unit, averaged) surface normal for nodes. Same size as x
     real(kind=realType), dimension(:, :), pointer :: norm

     ! The orignal nodal index. Size nNodes.
     integer(kind=intType), dimension(:), pointer :: ind

     ! The connectivity of the nodes forming 1D bar elements. Size (2, nElems)
     integer(kind=intType), dimension(:, :), pointer :: conn

     ! The index of my node numbers on my parent
     integer(kind=intType), dimension(:), pointer :: pNodes

     ! The index of my elements numbers on my parent
     integer(kind=intType), dimension(:), pointer :: pElems

     ! The string ID and index of my nodes on a split substing
     integer(kind=intType), dimension(:, :), pointer :: cNodes

     ! The cloest string ID of each node *AND* the node index on the
     ! other string. Size (2, nNodes)
     integer(kind=intType), dimension(:, :), pointer :: otherID

     ! Node location on the other string. 
     real(kind=realType), dimension(:, :), pointer :: otherX

     ! The inverse of the connectivity node to elem array. Size (5,
     ! nNodes). First index is the number of elements, other 4 entries
     ! are the up to 4 possible element neighbours. 
     integer(kind=intType), dimension(:, :), pointer :: nte

     ! Two buffer used for storing element indices while creating
     ! chains. Size (2, nElem)
     integer(kind=intType), dimension(:, :), pointer :: subStr

     ! The sizes of the two substrings
     integer(kind=intType), dimension(2) :: NsubStr

     ! A array to keep track of the number of elements
     ! "consumed" during chain searches or during zipping.
     integer(kind=intType), dimension(:), pointer :: elemUsed

     ! The KD tree for this string for performing fast seaches. 
     !type(tree_master_record), pointer :: tree
     type(kdtree2), pointer :: tree

     ! Pointer to the parent string
     type(oversetString), pointer :: p

     ! Pointer for next string for a linked list
     type(oversetString), pointer :: next

     ! List of all all directed edges. 
     type(oversetEdge), pointer, dimension(:) :: edges

     ! List of all computed triangles
     integer(kind=intType), dimension(:, :), pointer :: tris

     ! Number of trianges
     integer(kind=intType) :: nTris
     
  end type oversetString

  type oversetEdge
     ! Simple data structure representing a directed edge from n1->n2
     integer(kind=intType) :: n1, n2
  end type oversetEdge

  ! This is the flattened list of the fringes next to the wall that we
  !  have actually found donors for.
  ! tmpFringePtr is only used if we need to realloc. 
  type(fringeType), dimension(:), pointer :: localWallFringes, wallFringes, tmpFringePtr
  integer(kind=intType) :: nLocalWallFringe, nWallFringe

  ! These are the master overlap matrices
  type(CSRMatrix), dimension(:, :), allocatable, target :: overlapMatrix

  ! Some additional helper stuff
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType) :: nDomTotal
  real(kind=realType), dimension(:), allocatable :: clusterAreas

  type XPlane
     real(kind=realType), dimension(:, :, :), pointer :: xx
     integer(kind=intType), dimension(:, :), pointer :: nearWall
  end type XPlane

end module overset
