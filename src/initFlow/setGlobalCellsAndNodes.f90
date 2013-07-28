!
!     ******************************************************************
!     *                                                                *
!     * File:          setGlobalCellsAndNodes.f90                      *
!     * Author:        C.A.(Sandy) Mader, Gaetan Kenway                *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 12-28-2012                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine setGlobalCellsAndNodes(level)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Determine the global node numbering that is used to assemble   *
  !     * the adjoint system of equations. It take cares of all the halo *
  !     * nodes between the blocks.                                      *
  !     *                                                                *
  !     * The nodes are numbered according to the following sequence:    *
  !     *                                                                *
  !     * loop processor = 1, nProc                                      *
  !     *   loop domain = 1, nDom                                        *
  !     *     loop k = 2, kl                                             *
  !     *       loop j = 2, jl                                           *
  !     *         loop i = 2, il                                         *
  !     *                                                                *
  !     * Only the onwned nodes are numbered, meaning i/j/k span from 2  *
  !     * to il/jl/kl. The halo nodes receive the numbering from the     *
  !     * neighboring block that owns them.                              *
  !     *                                                                *
  !     * These variables are the same for all spectral modes, therefore *
  !     * only the 1st mode needs to be communicated.                    * 
  !     *                                                                *
  !     * This function will also set FMPointer which is only defined    *
  !     * on wall boundary conditions and points to the correct index    *
  !     * for the vectors that are of shape nsurface nodes               *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointVars 
  use BCTypes
  use blockpointers
  use communication
  use inputTimeSpectral

  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: level

  ! Local variables
  integer(kind=intType) :: nn, i, j, k, sps
  integer(kind=intType) :: ierr, istart
  logical :: commPressure, commViscous, commGamma
  integer(kind=intType), dimension(nProc) :: nNodes, nCells, nCellOffset, nNodeOffset
  integer(kind=intType), dimension(nDom) :: nCellBLockOffset,nNodeBLockOffset
  integer(kind=intType) :: npts, nts
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  integer(kind=intTYpe), dimension(:), allocatable :: nCellsProc, cumCellsProc
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii, jj,mm
  ! Determine the number of nodes and cells owned by each processor
  ! by looping over the local block domains.
  nCellsLocal(level) = 0
  nNodesLocal(level) = 0
  do nn=1,nDom
     ! Set to first spectral instance since we only need sizes
     call setPointers(nn, level, 1_intType)
     nCellsLocal(level) = nCellsLocal(level) + nx*ny*nz
     nNodesLocal(level) = nNodesLocal(level) + il*jl*kl
  enddo

  ! Reduce the number of cells in all processors: add up nCellsLocal
  ! into nCellsGlobal and sends the result to all processors.
  ! (use mpi sum operation)

  call mpi_allreduce(nCellsLocal(level), nCellsGlobal(level), 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)

  ! Gather the number of Cells per processor in the root processor.
  call mpi_gather(nCellsLocal(level), 1, sumb_integer, nCells, 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Repeat for the number of nodes.
  ! (use mpi sum operation)
  call mpi_allreduce(nNodesLocal(level), nNodesGlobal(level), 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)

  ! Gather the number of nodes per processor in the root processor.
  call mpi_gather(nNodesLocal(level), 1, sumb_integer, nNodes, 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Determine the global cell number offset for each processor.
  rootProc: if( myID==0) then
     nCellOffset(1) = 0
     nNodeOffset(1) = 0
     do nn=2,nProc
        nCellOffset(nn) = nCellOffset(nn-1) + nCells(nn-1)
        nNodeOffset(nn) = nNodeOffset(nn-1) + nNodes(nn-1)
     enddo
  endif rootProc

  ! Scatter the global cell number offset per processor.
  call mpi_scatter(nCellOffset, 1, sumb_integer, nCellOffsetLocal(level), 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Determine the global cell number offset for each local block.
  nCellBlockOffset(1) = nCellOffsetLocal(level)
  do nn=2,nDom
     call setPointers(nn-1, level, 1)
     nCellBlockOffset(nn) = nCellBlockOffset(nn-1)          &
          + nx*ny*nz
  enddo

  ! Repeat for nodes.
  call mpi_scatter(nNodeOffset, 1, sumb_integer, nNodeOffsetLocal(level), 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Determine the global node number offset for each local block.
  nNodeBlockOffset(1) = nNodeOffsetLocal(level)
  do nn=2,nDom
     call setPointers(nn-1, level, 1)
     nNodeBlockOffset(nn) = nNodeBLockOffset(nn-1) + il*jl*kl
  enddo

  ! Determine the global block row index for each (i,j,k) cell in
  ! each local block.

  do nn=1, nDom
     do sps=1, nTimeIntervalsSpectral
        call setPointers(nn, level, sps)
        do k=2, kl
           do j=2, jl
              do i=2, il
                 ! modified Timespectral indexing. Put all time
                 ! instances of a give block adjacent to each other in
                 ! the matrix
                 globalCell(i, j, k) = &
                      nCellBLockOffset(nn)*nTimeIntervalsSpectral+nx*ny*nz*(sps-1)+&
                      (i-2) +(j-2)*nx +(k-2)*nx*ny
              enddo
           enddo
        enddo
     enddo
  end do

  ! Determine the global block row index for each (i,j,k) node in
  ! each local block.
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=1, kl
           do j=1, jl
              do i=1, il
                 !modified Timespectral indexing. Put all time
                 !instances of a give block adjacent to each other in
                 !the matrix
                 globalNode(i, j, k) = &
                      nNodeBLockOffset(nn)*nTimeIntervalsSpectral + &
                      il*jl*kl*(sps-1) + (i-1)+(j-1)*il + (k-1)*il*jl
              end do
           end do
        end do
     end do
  end do

  ! The above procedure has uniquely numbered all cells and nodes
  ! owned on each processor. However we must also determine the
  ! indices of the halo cells/nodes from other processors. Do do this,
  ! we will cheat slightly. We will copy globalCell into the Pressure
  ! variable. Then do a double halo exchange, and copy them back out
  ! into globalCel. Then we will do the same thing for the globalNode
  ! values. Pressure values are saved in dw just in case.

  ! ------------ globalCell ----------------
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 dw(i ,j, k, 1) = P(i, j, k)
                 P(i,j,k) = real(globalCell(i,j,k))
              end do
           end do
        end do
     end do
  end do

  istart = 0
  iend   = -1
  commPressure = .True.
  commGamma    = .False.
  commViscous  = .False.
  call wHalo2(level, istart, iend, commPressure, commGamma, commViscous)

  ! Copy back out
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0, kb
           do j=0, jb
              do i=0, ib
                 globalCell(i, j, k) = int(P(i, j, k)) 
              end do
           end do
        end do
     end do
  end do

  ! ------------ globalNode ----------------
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0,ke
           do j=0,je
              do i=0,ie
                 dw(i, j, k, 2) = X(i, j, k, 1)
                 X(i, j, k, 1) = globalNode(i, j, k)
              end do
           end do
        end do
     end do
  end do

  call exchangeCoor(level)
  
  ! Copy back out
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0, ke
           do j=0, je
              do i=0, ie
                 globalNode(i, j, k) = X(i, j, k, 1)
              end do
           end do
        end do
     end do
  end do

  ! Reset Pressure 
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0,kb
           do j=0,jb
              do i=0,ib
                 P(i, j, k) = dw(i, j, k, 1)
                 dw(i, j, k, 1) = zero

              end do
           end do
        end do
     end do
  end do
  
  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=0,ke
           do j=0,je
              do i=0,ie
                 X(i, j, k, 1) = dw(i, j, k, 2)
                 dw(i, j, k, 2) = zero
              end do
           end do
        end do
     end do
  end do
  
  ! Now we will do indexing for FM. 
  
  ! First get the number of points on each proc and communicate
  call getForceSize(npts, ncells)

  allocate(nNodesProc(nProc), cumNodesProc(0:nProc))
  allocate(nCellsProc(nProc), cumCellsProc(0:nProc))
  nNodesProc(:) = 0_intType
  nCellsProc(:) = 0_intType

  call mpi_allgather(npts*nTS, 1, sumb_integer, nNodesProc, 1, sumb_integer, &
       sumb_comm_world, ierr)

  call mpi_allgather(ncells*nTS, 1, sumb_integer, nCellsProc, 1, sumb_integer, &
       sumb_comm_world, ierr)

  ! Sum and Allocate receive displ offsets
  cumNodesProc(0) = 0_intType
  cumCellsProc(0) = 0_intType
  do i=1, nProc
     cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)
     cumCellsProc(i) = cumCellsProc(i-1) + nCellsProc(i)
  end do
  
  ! Now we know the offset for the start of each processor. We can
  ! loop through in the desired order and just increment.
  ii = cumNodesProc(myid)
  jj = cumCellsProc(myid)

  do sps=1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn, 1_intType, sps)
        bocos: do mm=1,nBocos
           if(BCType(mm) == EulerWall.or.BCType(mm) == NSWallAdiabatic .or.&
                BCType(mm) == NSWallIsothermal) then
      
              jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    bcData(mm)%FMNodeIndex(i,j) = ii
                    ii = ii + 1
                 end do
              end do
              do j=jBeg+1, jEnd
                 do i=iBeg+1, iEnd
                    bcData(mm)%FMCellIndex(i,j) = jj
                    jj = jj + 1
                 end do
              end do
           end if
        end do bocos
     end do
  end do
  deallocate(nNodesProc, cumNodesProc)
  deallocate(nCellsProc, cumCellsProc)
end subroutine setGlobalCellsAndNodes
