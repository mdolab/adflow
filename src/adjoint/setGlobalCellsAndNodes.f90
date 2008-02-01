!
!     ******************************************************************
!     *                                                                *
!     * File:          setGlobalCellsAndNodes.f90                              *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine setGlobalCellsAndNodes(level,sps)
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
      use block
      use communication
      implicit none
!
!     Subroutine arguments.
!
      integer(kind=intType), intent(in) :: level, sps
!
!     Local variables.
!
      integer(kind=intType) :: nn, i, j, k, il, jl, kl, ie, je, ke
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
      ! Determine the number of nodes and cells owned by each processor
      ! by looping over the local block domains.

      nCellsLocal = 0
      nNodesLocal = 0
      do nn=1,nDom
        nCellsLocal = nCellsLocal + flowDoms(nn,level,sps)%il &
                                  * flowDoms(nn,level,sps)%jl &
                                  * flowDoms(nn,level,sps)%kl
        nNodesLocal = nNodesLocal + flowDoms(nn,level,sps)%ie &
                                  * flowDoms(nn,level,sps)%je &
                                  * flowDoms(nn,level,sps)%ke
      enddo
      print *,'local cells counted'

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

      print *,' local nodes gathered'

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
      
      print *,'scatter'

      ! Determine the global cell number offset for each local block.

      nCellBlockOffset(1) = nCellOffsetLocal
      do nn=2,nDom
        nCellBlockOffset(nn) = nCellBlockOffset(nn-1)          &
                         + flowDoms(nn-1,level,sps)%il &
                         * flowDoms(nn-1,level,sps)%jl &
                         * flowDoms(nn-1,level,sps)%kl
      enddo
      print *,'global cell offsets determined'

      ! Repeat for nodes.

      call mpi_scatter(nNodeOffset, 1, sumb_integer, nNodeOffsetLocal, 1, &
                      sumb_integer, 0, SUmb_comm_world, ierr)

      ! Determine the global cell number offset for each local block.

      nNodeBlockOffset(1) = nNodeOffsetLocal
      do nn=2,nDom
        nNodeBlockOffset(nn) = nNodeBLockOffset(nn-1)          &
                         + flowDoms(nn-1,level,sps)%il &
                         * flowDoms(nn-1,level,sps)%jl &
                         * flowDoms(nn-1,level,sps)%kl
      enddo
      print *,'global node offsets determined'

     
      ! Determine the global block row index for each (i,j,k) cell in
      ! each local block.

      do nn=1,nDom
         call setPointers(nn,level,sps)
!        il = flowDoms(nn,level,sps)%il
!        jl = flowDoms(nn,level,sps)%jl
!        kl = flowDoms(nn,level,sps)%kl
        do k=2,kl
          do j=2,jl
            do i=2,il
              flowDoms(nn,level,sps)%globalCell(i,j,k) &
                = nCellBLockOffset(nn) +(i-2) +(j-2)*il +(k-2)*il*jl
            enddo
          enddo
        enddo
      enddo
      print *,'global cell bock row determined'
     

      ! Determine the global block row index for each (i,j,k) node in
      ! each local block.

      do nn=1,nDom
         call setPointers(nn,level,sps)
!        ie = flowDoms(nn,level,sps)%ie
!        je = flowDoms(nn,level,sps)%je
!        ke = flowDoms(nn,level,sps)%ke
        do k=1,ke
          do j=1,je
            do i=1,ie
              flowDoms(nn,level,sps)%globalNode(i,j,k) &
                = nNodeBLockOffset(nn) +(i-1) +(j-1)*il +(k-1)*il*jl
            enddo
          enddo
        enddo
      enddo
      print *,'end'

      ! Synchronize the processors, just to be sure.

      call mpi_barrier(SUmb_comm_world, ierr)

      ! Determine the number of halo's to be exchanged
      ! and communicate the global node indices for the halos.
      
      nHalo=2

      ! Create a new copy of the communication routines to update the 
      ! indices in the halo cells and nodes.

      call haloIndexCommunication(level, nHalo)

     

    end subroutine setGlobalCellsAndNodes
