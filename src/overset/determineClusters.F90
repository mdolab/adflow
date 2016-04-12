!
!      ******************************************************************
!      *                                                                *
!      * determineClusters determines which blocks are connected with   *
!      * 1to1 cgns connections. Essentially what we are doing is        *
!      * identifying the consitutive multiblock meshes that make up     *
!      * an overset mesh. There should be precisely 1 cluster for a     *
!      * face mached mesh and 2 or more for a overset mesh              *
!      *                                                                *
!      ******************************************************************

subroutine determineClusters(clusters, N, cumDomsProc, clusterID)

  use constants
  use cgnsGrid
  use blockPointers
  use communication
  implicit none

  ! Input/output variables
  integer(kind=intType), intent(in) :: N
  integer(kind=intType), dimension(N), intent(out) :: clusters
  integer(kind=intType), dimension(0:nProc), intent(in) :: cumDomsProc
  integer(kind=intType), intent(out) :: clusterID

  ! Working variables
  integer(kind=intType) :: numBlocks, blockID, cgnsBlk, ierr
  integer(kind=intType) :: i, nn
  integer(kind=intType), dimension(N) :: clustersLocal
  logical :: blocksAvailable

  ! Initialize the cluster of each of the CGNSDoms to 0
  do i=1, cgnsNDom
     cgnsDoms(i)%cluster = 0
  end do

  ! Initialize cluster counter
  clusterID = 0

  ! Initialize counter of classified blocks
  blockID = 0

  ! Initialize variable to state that we have unclassified blocks
  blocksAvailable = .True. 

  ! Loop until all blocks are checked
  do while (blocksAvailable)

     ! Find position of the available block
     blocksAvailable = .false.

     do while ((.not. blocksAvailable) .and. (blockID .lt. cgnsnDom))
        blockID = blockID + 1 ! Increment counter
        if (cgnsDoms(blockID)%cluster .eq. 0) then
           blocksAvailable = .true.
        end if
     end do

     ! If we have blocks available, we start the search
     if (blocksAvailable) then
        clusterID = clusterID + 1 ! Increment the running cluser counter
        cgnsDoms(blockID)%cluster = clusterID
        call clusterSearch(blockID)
     end if

  end do

  ! Set the clusters to 0 so we can just all reduce
  clustersLocal = 0

  ! Set the cluster ID for all my blocks:
  do nn=1,nDom
     cgnsBlk = flowDoms(nn, 1, 1)%cgnsBlockID
     clustersLocal(cumDomsProc(myid) + nn) = cgnsDoms(cgnsBlk)%cluster
  end do
  call MPI_Allreduce(clustersLocal, clusters, N, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)

contains

  recursive subroutine clusterSearch(blockID)

    ! This is the recursive part of cluster search
    implicit none

    ! Subroutine inputs
    integer(kind=intType), intent(in) :: blockID

    ! Working variables
    integer(kind=intTYpe) :: clusterID, connID, connBlock

    ! Get the cluster ID from the reference block
    clusterID = cgnsDoms(blockID)%cluster

    ! Loop over all connections of this block
    do connID = 1, cgnsDoms(blockID)%n1to1

       connBlock = cgnsDoms(blockID)%conn1to1(connID)%donorBlock

       ! Check if connected block is already classified
       if (cgnsDoms(connBlock)%cluster == 0) then
          cgnsDoms(connBlock)%cluster = clusterID ! Assign block to the same cluster
          call clusterSearch(connBlock) ! Start search on the new block

       else if (cgnsDoms(connBlock)%cluster .ne. clusterID) then ! Check symmetry
          print *,'Non-symmetric connection between CGNS blocks:', blockID, ' and', connBlock
          stop
       end if
    end do

  end subroutine clusterSearch
end subroutine determineClusters
