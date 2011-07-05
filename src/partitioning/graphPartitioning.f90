!
!      ******************************************************************
!      *                                                                *
!      * File:          graphPartitioning.f90                           *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-22-2004                                      *
!      * Last modified: 11-22-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine graphPartitioning(emptyPartitions, commNeglected)
!
!      ******************************************************************
!      *                                                                *
!      * graphPartitioning partitions the corresponding graph of the    *
!      * computational blocks such that both the number of cells and    *
!      * number of faces is about equal on all processors.              *
!      *                                                                *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use partitionMod
       use inputParallel
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(out) :: emptyPartitions, commNeglected
!
!      Variables to store the graph in metis format.
!
       ! nVertex        : Number of vertices in the graph, equals nBlocks.
       ! nCon           : Number of contraints, 2.
       ! xadj(0:nVertex): Number of edges per vertex, cumulative storage
       !                  format.
       ! adjncy(:)      : End vertex of the edge; the size of adjncy is
       !                  xadj(nVertex).
       ! vwgt(:,:)      : Vertex weights, size equals nCon,nVertex. The
       !                  vertex weights are stored contiguously.
       ! adjwgt(:)      : Edge weights, size equals xadj(nVertex). Note
       !                  that the edge weights of edge i-j can be
       !                  different from the weight of edge j-i.
       ! wgtflag        : Whether or not to use weights on edges. Here
       !                  wgtflag should always be 1 to indicate that
       !                  edge weights are used.
       ! numflag        : Flag to indicate the numbering convention,
       !                  starting from 0 or 1. Here we start from 0.
       ! nParts         : Number of parts to split the graph. This is
       !                  nProc.
       ! ubvec(2)       : Tolerance for the constraints. Stored in the
       !                  module dpartitionMod.
       ! options(5)     : Option array; normally the default is used
       !                  indicated by options(1) = 0.
       ! edgecut        : On return it contains the edge cut of the
       !                  distributed graph.
       ! part(nVertex)  : On return the processor ID for each block.
       !                  It will be returned in fortran numbering,
       !                  i.e. starting at 1.  Stored in the module
       !                  distributionMod.

       integer :: nVertex, nCon, wgtflag, numflag, nParts, edgecut
       integer, dimension(5) :: options

       integer(kind=intType), dimension(:),   allocatable :: xadj, adjncy
       integer(kind=intType), dimension(:),   allocatable :: adjwgt
       integer(kind=intType), dimension(:,:), allocatable :: vwgt
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j
       integer(kind=intType) :: nEdges, nEdgesMax, ii, jj, kk

       integer(kind=intType), dimension(0:nProc-1) :: nBlockPerProc

       integer(kind=intType), dimension(:), allocatable :: tmp

       integer(kind=8) :: nCellsTotal    ! 8 byte integers to avoid
       integer(kind=8) :: nFacesTotal    ! overflow.
!
!      Function definition
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check whether part is allocated from a previous call. If so,
       ! release the memory.

       if( allocated(part) ) then
         deallocate(part, stat=ierr)
         if(ierr /= 0)                       &
         call terminate("graphPartitioning", &
                        "Deallocation failure for part")
       endif

       ! Determine the number of edges in the graph and the maximum
       ! number for a vertex in the graph.

       nEdges = 0
       nEdgesMax = 0
       do i=1,nBlocks
         ii        = blocks(i)%n1to1 + ubound(blocks(i)%overComm,1)
         nEdges    = nedges + ii
         nEdgesMax = max(nedgesMax, ii)
       enddo

       ! Initialize some values for the graph.

       nVertex  = nBlocks
       nCon     = 2
       wgtflag  = 1
       numflag  = 0
       nParts   = nProc
       ubvec(1) = one + loadImbalance
       ubvec(2) = one + loadImbalance
       options  = 0

       !  Allocate the memory to store/build the graph.

       allocate(xadj(0:nVertex), vwgt(nCon,nVertex), adjncy(nEdges), &
                adjwgt(nEdges), part(nVertex), tmp(nEdgesMax), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("graphPartitioning", &
                        "Memory allocation failure for graph variables")

       ! Initialize xadj(0) to 0.
       ! Furthermore initialize adjwgt to 0, as these values are
       ! accumulated due to multiple subfaces between blocks.

       xadj(0) = 0
       adjwgt  = 0

       ! Loop over the number of blocks to build the graph.

       graphVertex: do i=1,nBlocks

         ! Store the both vertex weights.

         vwgt(1,i) = blocks(i)%nCell
         vwgt(2,i) = blocks(i)%nFace

         ! Sort the neighbors in increasing order and neglect the
         ! communication to myself, i.e. do not allow an edge to myself.
         ! The sorting is necessary, because a block might have several
         ! subfaces with another block

         nEdges = 0
         do j=1,blocks(i)%n1to1
           ii = blocks(i)%nBocos + j
           if(blocks(i)%neighBlock(ii) /= i) then
             nEdges = nEdges +1
             tmp(nEdges) = blocks(i)%neighBlock(ii)
           endif
         enddo

         do j = 1,ubound(blocks(i)%overComm,1)
           if(blocks(i)%overComm(j,1) /= i) then
             nedges = nedges +1
             tmp(nedges) = blocks(i)%overComm(j,1)
           endif
         enddo

         ! Sort tmp in increasing order and get rid of the possible
         ! multiple entries.

         call qsortIntegers(tmp, nEdges)

         ii = min(nEdges,1_intType)  ! Be aware of nEdges == 0
         do j=2,nEdges
           if(tmp(j) /= tmp(ii)) then
             ii = ii + 1
             tmp(ii) = tmp(j)
           endif
         enddo

         ! Set nEdges to ii and update xadj(i).

         nEdges  = ii
         xadj(i) = xadj(i-1) + nEdges

         ! Repeat the loop over the subfaces, but now store
         ! the edge info.

         Edges1to1: do j=1,blocks(i)%n1to1

           ii = blocks(i)%nBocos + j
           if(blocks(i)%neighBlock(ii) /= i) then

             ! Search for the block ID and add the offset of xadj(i-1)
             ! to obtain the correct index to store the edge info.
             ! The -1 to the adjncy is present because C-numbering
             ! is used when calling Metis.

             jj = xadj(i-1) &
                + bsearchIntegers(blocks(i)%neighBlock(ii), tmp, nEdges)
             adjncy(jj) = blocks(i)%neighBlock(ii) - 1

             ! The weight equals the number of 1st and 2nd level halo
             ! cells to be communicated between the blocks. The weights
             ! are accumulated, as multiple subfaces between blocks are
             ! possible.

             kk = 1
             adjwgt(jj) = adjwgt(jj) + 2 * &
               ( max(abs(blocks(i)%inEnd(ii) - blocks(i)%inBeg(ii)), kk) &
               * max(abs(blocks(i)%jnEnd(ii) - blocks(i)%jnBeg(ii)), kk) &
               * max(abs(blocks(i)%knEnd(ii) - blocks(i)%knBeg(ii)), kk) )
           endif

         enddo Edges1to1

         ! Repeat the loop over the overset edges.

         EdgesOverset: do j = 1, ubound(blocks(i)%overComm,1)

           if(blocks(i)%overComm(j,1) /= i) then

             ! Search for the block ID and add the offset of xadj(i-1)
             ! to obtain the correct index to store the edge info.
             ! The -1 to the adjncy is present because C-numbering
             ! is used when calling Metis.

             jj = xadj(i-1) &
                + bsearchIntegers(blocks(i)%overComm(j,1), tmp, nedges)
             adjncy(jj) = blocks(i)%overComm(j,1) - 1

             ! The weight equals the number overset cells being
             ! communicated to this block. The weights are
             ! accumulated in case of repeats.

             adjwgt(jj) = adjwgt(jj) + blocks(i)%overComm(j,2)
           endif

         enddo EdgesOverset

       enddo graphVertex

       ! Metis has problems when the total number of cells or faces
       ! used in the weights exceeds 2Gb. Therefore the sum of these
       ! values is determined and an appropriate weight factor is 
       ! determined. Note that the type of nCellsTotal and nFacesTotal
       ! is integer*8.

       nCellsTotal = 0
       nFacesTotal = 0
       do i=1,nBlocks
         nCellsTotal = nCellsTotal + vwgt(1,i)
         nFacesTotal = nFacesTotal + vwgt(2,i)
       enddo

       if(nCellsTotal > 2147483647 .or. nFacesTotal > 2147483647) then
         nCellsTotal = nCellsTotal/2147483647 + 1
         nFacesTotal = nFacesTotal/2147483647 + 1

         do i=1,nBlocks
           vwgt(1,i) = vwgt(1,i)/nCellsTotal
           vwgt(2,i) = vwgt(2,i)/nFacesTotal
         enddo
       endif

       ! Loop over the number of attempts to partition the graph.
       ! In the first attempt the communication is taken into account.
       ! If not successful, i.e. empty partitions present, the metis
       ! routine is called once more, but now with zero adjwgt. This
       ! means that the communication cost is neglected and metis
       ! normally gives a valid partitioning.
       ! Initialize commNeglected to .false. This will change if in
       ! the loop below the first call to metis is not successful.

       commNeglected = .false.
       attemptLoop: do ii=1,2

         ! Call the graph partitioner.

         call metisInterface(nVertex, nCon, xadj, adjncy, vwgt, &
                             adjwgt, wgtflag, numflag, nParts,  &
                             ubvec, options, edgecut, part)

         ! Determine the number of blocks per processor.

         nBlockPerProc = 0
         do i=1,nBlocks
           nBlockPerProc(part(i)) = nBlockPerProc(part(i)) + 1
         enddo

         ! Check for empty partitions.

         emptyPartitions = .false.
         do i=0,nProc-1
           if(nBlockPerProc(i) == 0) emptyPartitions = .true.
         enddo

         ! Exit the loop if no empty partitions are present or if
         ! this is the second time this loop is executed.

         if(ii == 2 .or. (.not. emptyPartitions)) exit attemptLoop

         ! The first call to metis resulted in empty partitions.
         ! Ignore the communication, i.e. set the number of
         ! neighbors to 0, and try again.

         commNeglected = .true.
         xadj          = 0

       enddo attemptLoop

       ! Deallocate the memory for the graph except part.

       deallocate(xadj, vwgt, adjncy, adjwgt, tmp, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("graphPartitioning", &
                        "Deallocation failure for graph variables")

       end subroutine graphPartitioning
