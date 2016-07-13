subroutine initializeTractionScatter

  use blockPointers
  use ADjointPETSc
  use communication
  implicit none

  ! Working param
  integer(kind=intType) :: i,j, ierr, sps, nNodesLocal, nNodesTotal, nCellsLocal
  integer(kind=intType) :: nUnique
  real(kind=realType), dimension(:, :), allocatable :: localNodes, allNodes
  real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
  integer(kind=intType), dimension(:), allocatable :: link, localIndices
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  real(kind=realType) :: tol

  ! Set sps=1 since the communication will be same independent of the
  ! sps instance. 
  sps = 1

  ! Note that for overset, the reduction must be done by cluster. 
  
  ! Get the nodes for the surface
  call getForceSize(nNodesLocal, nCellsLocal)
  allocate(localNodes(3, nNodesLocal), nNodesProc(nProc), cumNodesProc(0:nProc))
  call getForcePoints(localNodes, nNodesLocal, sps)

  ! Determine the total number of nodes on each proc
  call mpi_allgather(nNodesLocal, 1, sumb_integer, nNodesProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Determine cumulative version
  cumNodesProc(0) = 0_intType
  nNodesTotal = 0
  do i=1, nProc
     nNodesTotal = nNodesTotal + nNodesProc(i)
     cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)
  end do

  ! Send all the nodes to everyone
  allocate(allNodes(3, nNodesTotal))
  call mpi_allgatherv(localNodes, nNodesLocal*3, sumb_real, allNodes, &
       nNodesProc*3, cumNodesProc*3, sumb_real, sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Local nodes is no longer necessary
  deallocate(localNodes)

  ! Now point reduce
  allocate(uniqueNodes(3, nNodestotal), link(nNodestotal))
  tol = 1e-12
 
  call pointReduce(allNodes, nNodesTotal, tol, uniqueNodes, link, nUnique)

  ! We can immediately discard everything but link since we are only
  ! doing logical operations here:
  deallocate(uniqueNodes, allNodes)

  ! Now back out the global indices for our local points
  allocate(localIndices(nNodesLocal))
  do i=1, nNodesLocal
     ! The -1 is to convert to 0-based ordering for petsc
     localIndices(i) = link(cumNodesProc(myid) + i)-1

  end do

  ! Create the basic (scalar) local vector
  call VecCreateMPI(SUMB_COMM_WORLD, nNodesLocal, PETSC_DETERMINE, &
       nodeValLocal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Create the basic global vector. For now it is evenly distributed,
  ! we could partition this more intelligently so that there is less
  ! communication. 
  call VecCreateMPI(SUMB_COMM_WORLD, PETSC_DETERMINE, nUnique, &
       nodeValGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(nodeValGlobal, sumGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now create the scatter that goes from the local vec to the global
  ! vec.

  ! Indices for the local vector is just a stride, starting at the
  ! offset

  call ISCreateStride(SUMB_COMM_WORLD, nNodesLocal, cumNodesProc(myid), &
       1, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Indices for the global vector are the "localIndices" we previously
  ! computed. 
  call ISCreateGeneral(sumb_comm_world, nNodesLocal, localIndices, &
       PETSC_COPY_VALUES, IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterCreate(nodeValLocal, IS1, nodeValGlobal, IS2, &
       tracScatter, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And dont' forget to destroy the index sets
  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

end subroutine initializeTractionScatter
