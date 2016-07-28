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
  use surfacefamilies, only : familyExchange
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
     integer(kind=intType) :: sps
     integer(kind=intType), dimension(:,:), allocatable :: ind, conn
     real(kind=realType), dimension(:, :), allocatable :: w, vars
     integer(kind=intType) :: nNodes
     real(kind=realType) :: pL, vL, pD, vD, CLp, CLv, CDp, CDv
     real(kind=realType) :: chord, twist, thickness
     real(kind=realType), dimension(3) :: le, te
     real(kind=realType), dimension(3) :: pt, dir
     integer(kind=intType), allocatable, dimension(:) :: famList
     type(familyExchange), pointer :: exch
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
     integer(kind=intType), dimension(:), allocatable :: famList
     real(kind=realType) :: dir(3)
     real(kind=realType) :: delta
     real(kind=realType), dimension(:,:), allocatable :: slicePts
     type(familyExchange), pointer :: exch
     
  end type liftDist

  logical :: liftDistInitialized = .False.
  integer(kind=intType) :: msCon1(16, 5), msCon2(4, 2)

 
  ! Data for the user supplied slices:
  integer(kind=intType), parameter :: nSliceMax=1000
  integer(kind=intType) :: nParaSlices=0
  integer(kind=intType) :: nAbsSlices=0
  type(slice), dimension(:, :), allocatable :: paraSlices, absSlices

  ! Data for the user supplied lift distributions
  integer(kind=intType), parameter :: nLiftDistMax=100
  integer(kind=intType) :: nLiftDists
  type(liftDist), dimension(nLiftDistMax), target :: liftDists

  ! Tecplot Variable names of the data in the lift distribution data file:
  character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistName
  integer(kind=intType), parameter :: nLiftDistVar=18

  
end module liftDistributionData

