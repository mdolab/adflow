module oversetData

  use constants
  use adtData, only : adtType
  use block, only : fringeType
  use kdtree2_module, only : kdtree2
#ifndef USE_TAPENADE
#include <petsc/finclude/petsc.h>
  use petsc
  implicit none

#endif
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
  ! donor searches.
  type oversetBlock

     ! Sizes for the block
     integer(kind=intType) :: il, jl, kl

     ! Cluster this block belongs to
     integer(kind=intType) :: cluster

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

     ! The processor where this set of fringes came from
     integer(kind=intType) :: proc

     ! The block number of these fringes on processor 'proc'
     integer(kind=intType) :: block

     ! Sizes
     integer(kind=intType) :: il, jl, kl, nx, ny, nz

     ! Cluster this set of fringes belongs to
     integer(kind=intType) :: cluster

     ! Buffer space for sending/receiving the fringes
     real(kind=realType), dimension(:), allocatable :: rBuffer
     integer(kind=intType), dimension(:), allocatable :: iBuffer

     ! These are the coordinates of the points we are searching for
     real(kind=realType), dimension(:, :), allocatable :: x

     ! These are the coordinate of its wall surface pt if applicable
     real(kind=realType), dimension(:, :), allocatable :: xSeed

     ! These are the indices of the wall surfaces
     integer(kind=intType), dimension(:), allocatable :: wallInd

     ! Flag specifying if this cell is next to a wall or not. 1 for
     ! next to wall, 0 otherwise.
     integer(kind=intType), dimension(:), allocatable  :: isWall

     ! This is where we will store all the potential donors that have
     ! been found for this set of fringes
     integer(kind=intType) :: nDonor=0
     integer(kind=intType), dimension(:,:), pointer :: fringeIntBuffer=>null()
     real(kind=realType), dimension(:,:), pointer :: fringeRealBuffer=>null()

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
     real(kind=realType), dimension(:, :), pointer :: x => null()

     ! Only primal mesh cell centers for dual mesh. Only allocated as
     ! needed.
     real(kind=realType), dimension(:, :),  pointer :: xPrimalCen => null()

     ! Surface nonal used for determining if point is "underneath" the
     ! surface.
     real(kind=realType), dimension(:, :), pointer :: norm => null()

     ! Local estimate of surface error
     real(kind=realType), dimension(:), pointer :: delta => null()

     ! Connectivity for the surface
     integer(kind=intType), dimension(:, :), pointer:: conn => null()

     ! ind: Global node index for nodes
     integer(kind=intType), dimension(:), pointer :: ind => null()

     ! indPrimal: Global node index. Only temporarly used to store
     ! index for strictly primal cells on dual mesh.
     integer(kind=intType), dimension(:), pointer :: indPrimal => null()

     ! indCell: Global cell index for wall cells
     integer(kind=intType), dimension(:), pointer :: indCell

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

     ! Whether or no this string is a pocket
     logical :: isPocket=.False.

     ! --------------------------------------------------------------------
     ! Node Data: The actual physical node locations, unit surface normal,
     ! perpNormal and mesh size. x is from index 1:3, normal from 4:6,
     ! perpNormal form 7:9 and h is index 10. This pointer gets allocated.
     real(kind=realType), dimension(:, :), pointer :: nodeData => null()

     ! Pointer for physical node location. Points to nodeData
     real(kind=realType), dimension(:, :), pointer :: x => null()

     ! Pointer for nodal unit normal. Points to nodeData
     real(kind=realType), dimension(:, :), pointer :: norm => null()

     ! Pointer for nodal unit perpendicual in-plane normal. Points to nodeData
     real(kind=realType), dimension(:, :), pointer :: perpNorm => null()

     ! Pointer for nodal element size. Points to nodeData
     real(kind=realType), dimension(:), pointer :: h => null()

     ! --------------------------------------------------------------------
     ! Integer Node Data: This stores the global node index (into
     ! xVec), the cluster ID of the node, and the family ID of the node
     integer(kind=intType), dimension(:, :), pointer :: intNodeData => null()

     ! The orignal nodal index. Size nNodes. Pointer into intNodeData
     integer(kind=intType), dimension(:), pointer :: ind => null()

     ! The cluster the node came from. Pointer into intNodeData
     integer(kind=intType), dimension(:), pointer :: cluster => null()

     ! The family the node came from. Pointer into intNodeData
     integer(kind=intType), dimension(:), pointer :: family => null()
     ! --------------------------------------------------------------------

     ! The connectivity of the nodes forming 1D bar elements. Size (2, nElems)
     integer(kind=intType), dimension(:, :), pointer :: conn => null()

     ! The index of my node numbers on my parent
     integer(kind=intType), dimension(:), pointer :: pNodes => null()

     ! The index of my elements numbers on my parent
     integer(kind=intType), dimension(:), pointer :: pElems => null()

     ! The string ID and index of my nodes on a split substing
     integer(kind=intType), dimension(:, :), pointer :: cNodes => null()

     ! The cloest string ID of each node *AND* the node index on the
     ! other string. Size (2, nNodes)
     integer(kind=intType), dimension(:, :), pointer :: otherID => null()

     ! The inverse of the connectivity node to elem array. Size (5,
     ! nNodes). First index is the number of elements, other 4 entries
     ! are the up to 4 possible element neighbours.
     integer(kind=intType), dimension(:, :), pointer :: nte => null()

     ! Two buffer used for storing element indices while creating
     ! chains. Size (2, nElem)
     integer(kind=intType), dimension(:, :), pointer :: subStr => null()

     ! The sizes of the two substrings
     integer(kind=intType), dimension(2) :: NsubStr

     ! A array to keep track of the number of elements
     ! "consumed" during chain searches or during zipping.
     integer(kind=intType), dimension(:), pointer :: elemUsed => null()

     ! Array to keep track of nodes used to contruct string pairs for
     ! crossZipping.
     integer(kind=intType), dimension(:), pointer :: XzipNodeUsed => null()

     ! The KD tree for this string for performing fast seaches.
     !type(tree_master_record), pointer :: tree
     type(kdtree2), pointer :: tree => null()

     ! Pointer to the parent string
     type(oversetString), pointer :: p => null()

     ! Pointer for next string for a linked list
     type(oversetString), pointer :: next => null()

     ! List of all all directed edges.
     type(oversetEdge), pointer, dimension(:) :: edges => null()

     ! nEdges: The number of new edges added due to triangles.
     integer(kind=intTYpe) :: nEdges

     ! List of all computed triangles
     integer(kind=intType), dimension(:, :), pointer :: tris => null()

     ! Number of trianges
     integer(kind=intType) :: nTris

     ! surfCellID(1:nTris)
     ! Global cellID of the primal cell containing the triangle centroid
     integer(kind=intType), dimension(:), pointer :: surfCellID


  end type oversetString

  type oversetEdge
     ! Simple data structure representing a directed edge from n1->n2
     integer(kind=intType) :: n1, n2
  end type oversetEdge

  interface operator(<=)
     module procedure lessEqualEdgeType
  end interface operator(<=)

  interface operator(<)
     module procedure lessEdgeType
  end interface operator(<)

  type pocketEdge
     ! Simple data structure representing a directed edge from n1->n2
     ! Similar to oversetEdge, but introduced to do different type
     ! of sort on edges
     integer(kind=intType) :: n1, n2
  end type pocketEdge

  interface operator(<=)
     module procedure lessEqualPocketEdgeN2
  end interface operator(<=)

  interface operator(<)
     module procedure lessPocketEdgeN2
  end interface operator(<)

  interface operator(==)
     module procedure EqualPocketEdgeN2
  end interface operator(==)

  type zipperMesh

     ! Data required for zipper mesh surface integer for a particular BCGroup
     integer(kind=intType), dimension(:, :), allocatable :: conn
     integer(kind=intType), dimension(:), allocatable :: fam, indices
     logical :: allocated=.False.
#ifndef USE_TAPENADE
     VecScatter :: scatter
     Vec :: localVal
#endif
  end type zipperMesh

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
  integer(kind=intType) :: nClusters
  integer(kind=intType), dimension(:), allocatable :: clusters
  real(kind=realType), dimension(:), allocatable :: clusterAreas
  real(kind=realType), dimension(:), allocatable :: clusterMarchDist
  type(oversetWall), dimension(:), allocatable, target :: clusterWalls

  ! Flag specifying if overset is present in mesh
  logical :: oversetPresent

  ! Zipper meshes
  type(zipperMesh), dimension(nFamExchange), target :: zipperMeshes

  ! Static arrays for doing timings
  real(kind=realType), dimension(iTotal) :: tStart
  real(kind=realType), dimension(iTotal) :: oversetTimes

  contains
  ! ==============================
  ! Operator overloading functions
  ! ==============================
  logical function lessEqualEdgeType(e1, e2)
    !
    !   lessEqualEdgeType returns .True. if e1<=e2 and .False. otherwise.
    !   Compared on the directed edge node indices n1 and n2.
    !   First compare wrt averaged node indices data. If same averaged
    !   node data, then compare wrt increasing or decreasing node indices.

    implicit none

    ! Input
    type(oversetEdge), intent(in) :: e1, e2

    ! Local variables
    integer(kind=intType) :: nsum1, nsum2, ndiff1, ndiff2

    ! Compare the averaged (or just sum of) node indices values.
    ! Positive sign for increasing order of node indices.
    nsum1 = e1%n1 + e1%n2
    nsum2 = e2%n1 + e2%n2
    ndiff1 = e1%n2 - e1%n1
    ndiff2 = e2%n2 - e2%n1

    ! Compare based on averaged node indices values
    if (abs(nsum1) < abs(nsum2)) then
       lessEqualEdgeType = .True.
       return
    else if (abs(nsum1) > abs(nsum2)) then
       lessEqualEdgeType = .False.
       return
    end if

    if (abs(nsum1) /= abs(nsum2)) &
       STOP ' *** Error in lessEqualEdgeType ***'

    ! Compare based on edge nodes difference
    if (abs(ndiff1) < abs(ndiff2)) then
       lessEqualEdgeType = .True.
       return
    else if (abs(ndiff1) > abs(ndiff2)) then
       lessEqualEdgeType = .False.
       return
    end if

    ! here abs(ndiff1) == abs(ndiff2) and
    ! abs(nsum1)== abs(nsum2), so same edge
    if (ndiff1 < ndiff2) then
       lessEqualEdgeType = .True.
       return
    else if (ndiff1 > ndiff2) then
       lessEqualEdgeType = .False.
       return
    end if

    ! here ndiff1 == ndiff2, hence .True.

    lessEqualEdgeType = .True.

  end function lessEqualEdgeType

  ! ---------------------------------
  logical function lessEdgeType(e1, e2)
    !
    !   lessEdgeType returns .True. if e1<e2 and .False. otherwise.
    !   Compared on the directed edge node indices n1 and n2.
    !   First compare wrt averaged node indices data. If same averaged
    !   node data, then compare wrt increasing or decreasing node indices.

    implicit none

    ! Input
    type(oversetEdge), intent(in) :: e1, e2

    ! Local variables
    integer(kind=intType) :: nsum1, nsum2, ndiff1, ndiff2

    ! Compare the averaged (or just sum of) node indices values.
    ! Positive sign for increasing order of node indices.
    nsum1 = e1%n1 + e1%n2
    nsum2 = e2%n1 + e2%n2
    ndiff1 = e1%n2 - e1%n1
    ndiff2 = e2%n2 - e2%n1

    ! Compare based on averaged node indices values
    if (abs(nsum1) < abs(nsum2)) then
       lessEdgeType = .True.
       return
    else if (abs(nsum1) > abs(nsum2)) then
       lessEdgeType = .False.
       return
    end if

    if (abs(nsum1) /= abs(nsum2)) &
       STOP ' *** Error in lessEdgeType ***'

    ! Compare based on edge nodes difference
    if (abs(ndiff1) < abs(ndiff2)) then
       lessEdgeType = .True.
       return
    else if (abs(ndiff1) > abs(ndiff2)) then
       lessEdgeType = .False.
       return
    end if

    ! here abs(ndiff1) == abs(ndiff2) and
    ! abs(nsum1)== abs(nsum2), so same edge
    if (ndiff1 < ndiff2) then
       lessEdgeType = .True.
       return
    else if (ndiff1 > ndiff2) then
       lessEdgeType = .False.
       return
    end if

    ! here ndiff1 == ndiff2, hence .False.

    lessEdgeType = .False.

  end function lessEdgeType


  logical function lessEqualPocketEdgeN2(e1, e2)
    !
    !   lessEqualPocketEdgeN2 returns .True. if e1%n2<=e2%n2 and .False.
    !   otherwise. Compared on the directed edge node indices n1 and n2.

    implicit none

    ! Input
    type(pocketEdge), intent(in) :: e1, e2

    if (e1%n2 < e2%n2) then
       lessEqualPocketEdgeN2 = .True.
       return
    else if (e1%n2 > e2%n2) then
       lessEqualPocketEdgen2 = .False.
       return
    end if

    ! Here e1%n2==e2%n2, so edges are equal, hence .True.
    lessEqualPocketEdgeN2 = .True.

  end function lessEqualPocketEdgeN2

  logical function lessPocketEdgeN2(e1, e2)
    !
    !   lessPocketEdgeN2 returns .True. if e1%N2<e2%N2 and .False.
    !   otherwise. Compared on the directed edge node indices n1 and n2.

    implicit none

    ! Input
    type(pocketEdge), intent(in) :: e1, e2

    if (e1%N2 < e2%N2) then
       lessPocketEdgeN2 = .True.
       return
    else if (e1%N2 > e2%N2) then
       lessPocketEdgeN2 = .False.
       return
    end if

    ! Here e1%N2==e2%n2, so edges are equal, hence .False.
    lessPocketEdgeN2 = .False.

  end function lessPocketEdgeN2

  logical function EqualPocketEdgeN2(e1, e2)
    !
    !   EqualPocketEdgeN2 returns .True. if e1%N2==e2%N2 and .False.
    !   otherwise. Compared on the directed edge node indices n1 and n2.

    implicit none

    ! Input
    type(pocketEdge), intent(in) :: e1, e2

    if (e1%N2 == e2%N2) then
       EqualPocketEdgeN2 = .True.
       return
    else
       EqualPocketEdgeN2 = .False.
       return
    end if

  end function EqualPocketEdgeN2
end module oversetData
