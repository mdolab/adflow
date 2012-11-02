!
!     ******************************************************************
!     *                                                                *
!     * File:          setGlobalCellsAndNodes.f90                      *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
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
  !     ******************************************************************
  !
  use ADjointVars ! nNodesGlobal, nNodesLocal, nOffsetLocal
  use blockpointers
  use communication
  use inputTimeSpectral !nTimeIntervalsSpectral

  implicit none
  !
  !     Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !     Local variables.
  !
  integer(kind=intType) :: nn, i, j, k, sps,modFamID!, il, jl, kl, ie, je, ke
  integer :: ierr, nHalo

  integer(kind=intType), dimension(nProc) :: nNodes, nCells, nCellOffset, nNodeOffset

  integer(kind=intType), dimension(nDom) :: nCellBLockOffset,nNodeBLockOffset


  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution                                                *
  !     *                                                                *
  !     ******************************************************************
  !
  !set spectral level to 1

  sps=1

  ! Determine the number of nodes and cells owned by each processor
  ! by looping over the local block domains.

  nCellsLocal = 0
  nNodesLocal = 0

  do nn=1,nDom
     call setPointers(nn,level,sps)

     nCellsLocal = nCellsLocal + nx*ny*nz
     nNodesLocal = nNodesLocal + il  *jl * kl
  enddo

  ! Reduce the number of cells in all processors: add up nCellsLocal
  ! into nCellsGlobal and sends the result to all processors.
  ! (use mpi sum operation)

  call mpi_allreduce(nCellsLocal, nCellsGlobal, 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)

  ! Gather the number of Cells per processor in the root processor.
  call mpi_gather(nCellsLocal, 1, sumb_integer, nCells, 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Repeat for the number of nodes.
  ! (use mpi sum operation)
  call mpi_allreduce(nNodesLocal, nNodesGlobal, 1, sumb_integer, &
       mpi_sum, SUmb_comm_world, ierr)

  ! Gather the number of nodes per processor in the root processor.
  call mpi_gather(nNodesLocal, 1, sumb_integer, nNodes, 1, &
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
  call mpi_scatter(nCellOffset, 1, sumb_integer, nCellOffsetLocal, 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)


  ! Determine the global cell number offset for each local block.
  nCellBlockOffset(1) = nCellOffsetLocal
  do nn=2,nDom
     call setPointers(nn-1,level,sps)
     nCellBlockOffset(nn) = nCellBlockOffset(nn-1)          &
          + nx*ny*nz
  enddo


  ! Repeat for nodes.
  call mpi_scatter(nNodeOffset, 1, sumb_integer, nNodeOffsetLocal, 1, &
       sumb_integer, 0, SUmb_comm_world, ierr)

  ! Determine the global node number offset for each local block.

  nNodeBlockOffset(1) = nNodeOffsetLocal
  do nn=2,nDom
     call setPointers(nn-1,level,sps)
     nNodeBlockOffset(nn) = nNodeBLockOffset(nn-1)          &
          + il *jl * kl
  enddo

  ! Determine the global block row index for each (i,j,k) cell in
  ! each local block.

  do nn=1,nDom
     do sps = 1,nTimeIntervalsSpectral
        call setPointers(nn,level,sps)
        do k=2,kl
           do j=2,jl
              do i=2,il
                 !original steady indexing
                 !              flowDoms(nn,level,sps)%globalCell(i,j,k) &
                 !                = nCellBLockOffset(nn) +(i-2) +(j-2)*il +(k-2)*il*jl

                 ! modified Timespectral indexing. Put all time instances of a give block adjacent to each other in the matrix
                 flowDoms(nn,level,sps)%globalCell(i,j,k)=&
                      nCellBLockOffset(nn)*nTimeIntervalsSpectral+nx*ny*nz*(sps-1)+&
                      (i-2) +(j-2)*nx +(k-2)*nx*ny
              enddo
           enddo
        enddo
     enddo
  end do

  ! Determine the global block row index for each (i,j,k) node in
  ! each local block.

  do sps = 1,nTimeIntervalsSpectral
     do nn=1,nDom
        call setPointers(nn,level,sps)
        do k=1,kl
           do j=1,jl
              do i=1,il
                 !modified Timespectral indexing. Put all time instances of a give block adjacent to each other in the matrix
                 flowDoms(nn,level,sps)%globalNode(i,j,k) &
                      =nNodeBLockOffset(nn)*nTimeIntervalsSpectral + il *jl * kl*(sps-1)+&
                      (i-1) +(j-1)*il +(k-1)*il*jl
              enddo
           enddo
        enddo
     enddo
  enddo


  nHalo=2

  ! Create a new copy of the communication routines to update the 
  ! indices in the halo cells and nodes.

  call haloIndexCommunication(level, nHalo)

end subroutine setGlobalCellsAndNodes
