module liftDistributionData
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This module contains the definition of the derived data type   *
  !      * for slice which is the basic building block for lift           *
  !      * distribuions. It also contains the data required for           *
  !      * communicaiting force data. 
  !      *                                                                *
  !      ******************************************************************
  !

  use constants
  implicit none
  save

  type slice
     
     ! nNodes : Number of nodes for this slice
     ! ind(2, nNodes) : Indices of the nodes in the global node list on either 
     !                  side of node 
     ! w(2, nNodes) : Weights used to multiply the two global nodes defined in 
     !                ind to get compute nodal values (positions, forces etc)
     ! pL, vL, pD, vD : Pressure and viscous lift, pressure and viscous drag
     ! CLp, CLv, CDp, CDv : Coefficients of pressure and viscous lift and drag
     ! chord: chord of section

     character(len=maxStringLen) :: sliceName
     integer(kind=intType) :: nNodes
     integer(kind=intType), dimension(:,:), pointer :: ind
     real(kind=realType), dimension(:, :), pointer :: w
     real(kind=realType) :: pL, vL, pD, vD, CLp, CLv, CDp, CDv
     real(kind=realType) :: chord, twist, thickness
     real(kind=realType), dimension(3) :: pt, dir
  end type slice

  type liftDist
     ! nSegments: Number of nodes to use for distribution
     ! dir_ind: Index of direction..1 for x, 2 for y, 3 for z
     ! dir: Slice direction
     ! distName: Name of lift distribution
     ! slices: The list of slices this distribution will use
     ! delta: The current delta spacing for the distribution
     ! slicePoints: The list of points where the slices are taken
     character(len=maxStringLen) :: distName
     integer(kind=intType) :: nSegments, dir_ind
     real(kind=realType) :: dir(3)
     type(slice), dimension(:), allocatable :: slices
     real(kind=realType) :: delta
     real(kind=realType), dimension(:,:), allocatable :: slicePts
  end type liftDist

  ! Data used to define the global reduced FE surface mesh and data
  ! used to produce 

  ! nNodesLocal : Number of nodes on wall surfaces on this processor
  ! nCellsLocal : Number of cells (quads) on wall surface on this processor
  ! nNodesTotal : Total number of nodes on wall surfaces (sum of nNodesLocal across procs)
  ! nCellsTotal : Total number of cells on wall surfaces (sum of nCellsLocal across procs)
  ! nUnique : Total number of unique nodes on the global mesh

  ! nNodesProc(nProc) : Number of nodes on each proc, stored on each processor
  ! nCellsProc(nProc) : NUmber of cells on each proc, stroed on each processor
  ! cumNodesProc(0:nProc) : Cumulative number of nodes on each processor
  ! cumCellsProc(0:nProc) : Cululative number of cells on each procssor
  ! localCells(4, nCellsLocal) : Arrary used to store the connectivity of the local mesh
  ! allCells(4, nCellsTotal) : Array used to store connectivity of the global mesh. Only
  !                            allocated on root processor
  ! link(nNodesTotal) : The index array that for each node in allNodes, points to the 
  !                     cooresponding node in uniqueNodes
  ! localNodes(3, nNodesLocal) : Array for storing the nodes on the local wall surfaces
  ! allNodes(3, nNodesTotal)   : Array for storing the nodes on the global mesh. Only 
  !                              allocated on the root processor
  ! localForcesP(3, nNodesLocal) : Pressure forces on the local wall surfaces
  ! localForcesV(3, nNodesLocal) : Viscous forces on the local wall surfaces
  ! allForcesP(3, nNodesTotal) : Pressure forces on the global wall surface
  ! allForcesV(3, nNodesTotal) : Viscous forces on the global wall surface
  ! uniqueNodes(3, nNodesTotal): The list of unique nodes defining the global surface mesh
  !                              While this is allocated to size nNodesTotal, only the 
  !                              first nNodesUnique are meaningful. Only allocated on 
  !                              root proc
  ! uniqueTractionsP(3, nUnique) : Unique pressure tractions global mesh
  ! uniqueTractionsV(3, nUnique) : Unique viscous tractions on global mesh
  ! dualAreas(nUnqiue) : Dual areas (area surrounding nodes) for the global surface mesh
  ! fc(nUnique) : Array storing the current function value used for performing the 
  !               slicing operation
  ! liftDistInitialzed = Flag to specify if the lift distribution data has been allocated
  ! msCon1 : Lookup table for marching squares
  ! msCon2 :Get Nodes from Edges

  integer(kind=intType) :: nNodesLocal, nCellsLocal, nNodesTotal, nCellsTotal, nUnique
  integer(kind=intType), allocatable, dimension(:)    :: nNodesProc, nCellsProc
  integer(kind=intType), allocatable, dimension(:)    :: cumNodesProc, cumCellsProc
  integer(kind=intType), allocatable, dimension(:,:)  :: localCells, allCells
  integer(kind=intType), allocatable, dimension(:)    :: link
  real(kind=realType),   allocatable, dimension(:, :) :: localNodes, allNodes
  real(kind=realType),   allocatable, dimension(:, :) :: localForcesP, localForcesV
  real(kind=realType),   allocatable, dimension(:, :) :: allForcesP, allForcesV
  real(kind=realType),   allocatable, dimension(:, :) :: uniqueNodes
  real(kind=realType),   allocatable, dimension(:, :) :: uniqueTractionsP, uniqueTractionsV
  real(kind=realType),   allocatable, dimension(:)    :: dualAreas, fc
  logical :: liftDistInitialized = .False.
  integer(kind=intType) :: msCon1(16, 5), msCon2(4, 2)

  ! Data for the user supplied slices:
  integer(kind=intType), parameter :: nSliceMax=10000
  integer(kind=intType) :: nParaSlices, nAbsSlices
  type(slice), dimension(nSliceMax) :: paraSlices, absSlices

  ! Data for the user supplied lift distributions
  integer(kind=intType), parameter :: nLiftDistMax=100
  integer(kind=intType) :: nLiftDists
  type(liftDist), dimension(nLiftDistMax), target :: liftDists

  ! Tecplot Variable names of the data in the lift distribution data file:
  character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistName
  integer(kind=intType), parameter :: nLiftDistVar=18

end module liftDistributionData
