module userSurfaceIntegrationData
  use constants

  type userSurfCommType
     ! Data required on each proc:

     ! nDonor: The number of donor points the proc will provide
     ! frac (3, nDonor) : The uvw coordinates of the interpolation point
     ! donorInfo(4, nDonor) : Donor information. 1 is the local block ID and 2-4 is the 
     !    starting i,j,k indices for the interpolation. 
     ! procSizes(0:nProc-1) : The number of donors on each proc
     ! procDisps(0:nProc) : Cumulative form of procSizes

     ! inv(nConn) : Array allocated only on root processor used to
     ! reorder the nodes or elements back to the original order. 

     integer(kind=intType) :: nDonor
     real(kind=realType), dimension(:,:), allocatable :: frac
     integer(kind=intType), dimension(:, :), allocatable :: donorInfo
     integer(kind=intTYpe), dimension(:), allocatable :: procSizes, procDisps
     integer(kind=intTYpe), dimension(:), allocatable :: inv
     logical, dimension(:), allocatable :: valid

  end type userSurfCommType

  type userIntSurf

     character(len=maxStringLen) :: famName
     integer(Kind=intType) :: famID
     real(kind=realType), dimension(:, :), allocatable :: pts
     integer(kind=intType), dimension(:, :), allocatable :: conn

     ! Two separate commes: One for the nodes (based on the primal
     ! mesh) and one for the variables (based on the dual mesh)
     type(userSurfCommType) :: nodeComm, faceComm

  end type userIntSurf

  integer(kind=intType), parameter :: nUserIntSurfsMax=25
  type(userIntSurf), dimension(nUserIntSurfsMax), target :: userIntSurfs
  integer(kind=intTYpe) :: nUserIntSurfs=0
end module userSurfaceIntegrationData
