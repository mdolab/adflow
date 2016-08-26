!
!      File:          setGlobalCellsAndNodes.f90                      
!      Author:        C.A.(Sandy) Mader, Gaetan Kenway                
!      Starting date: 01-14-2008                                      
!      Last modified: 12-28-2012                                      
!
subroutine setGlobalCellsAndNodes(level)
  !
  !      Determine the global node numbering that is used to assemble   
  !      the adjoint system of equations. It take cares of all the halo 
  !      nodes between the blocks.                                      
  !      The nodes are numbered according to the following sequence:    
  !      loop processor = 1, nProc                                      
  !        loop domain = 1, nDom                                        
  !          loop k = 2, kl                                             
  !            loop j = 2, jl                                           
  !              loop i = 2, il                                         
  !      Only the onwned nodes are numbered, meaning i/j/k span from 2  
  !      to il/jl/kl. The halo nodes receive the numbering from the     
  !      neighboring block that owns them.                              
  !      These variables are the same for all spectral modes, therefore 
  !      only the 1st mode needs to be communicated.                    
  !      This function will also set FMPointer which is only defined    
  !      on wall boundary conditions and points to the correct index    
  !      for the vectors that are of shape nsurface nodes               
  !
  use ADjointVars 
  use blockpointers
  use communication
  use inputTimeSpectral
  use utils, only: setPointers, terminate 
  use haloExchange, only : whalo1to1intgeneric
  implicit none

  ! Input variables
  integer(kind=intType), intent(in) :: level

  ! Local variables
  integer(kind=intType) :: nn, i, j, k, sps, iDim
  integer(kind=intType) :: ierr, istart
  logical :: commPressure, commLamVis, commEddyVis, commGamma
  integer(kind=intType), dimension(nProc) :: nNodes, nCells, nCellOffset, nNodeOffset
  integer(kind=intType), dimension(nDom) :: nCellBLockOffset,nNodeBLockOffset
  integer(kind=intType) :: npts, nCell, nNode
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  integer(kind=intTYpe), dimension(:), allocatable :: nCellsProc, cumCellsProc
  integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii, jj,mm

  do sps=1, nTimeIntervalsSpectral
     do nn=1, nDom
        call setPointers(nn, level, sps)
        ! Allocate memory for the cell and node indexing...only on sps=1
        allocate(flowDoms(nn,level,sps)%globalCell(0:ib,0:jb,0:kb), &
             flowDoms(nn,level,sps)%globalNode(0:ie,0:je,0:ke), stat=ierr)
        if (ierr /=0) then
           call terminate("setGlobalCellsAndNodes", "Allocation failure for globalCell/Node")
        end if
        ! Assign a 'magic number' of -5 to globalCell and global Node:
        flowDoms(nn,level,sps)%globalCell = -5
        flowDoms(nn,level,sps)%globalNode = -5
     end do
  end do

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
  ! indices of the halo cells/nodes from other processors. To do this
  ! we just run the specific halo exchanges for the cells and one for
  ! the nodes

  spectralModes: do sps=1,nTimeIntervalsSpectral
     domainLoop:do nn=1, nDom
        flowDoms(nn, level, sps)%intCommVars(1)%var => &
             flowDoms(nn, level, sps)%globalNode(:, :, :)
     end do domainLoop

     ! Run the generic integer exchange
     call wHalo1to1IntGeneric(1, level, sps, commPatternNode_1st, internalNode_1st)
  end do spectralModes

  spectralModes2: do sps=1,nTimeIntervalsSpectral
     domainLoop2:do nn=1, nDom
        flowDoms(nn, level, sps)%intCommVars(1)%var => &
             flowDoms(nn, level, sps)%globalCell(:, :, :)
     end do domainLoop2

     ! Run the generic integer exchange
     call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)
  end do spectralModes2

end subroutine setGlobalCellsAndNodes
