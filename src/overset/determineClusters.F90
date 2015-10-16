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

subroutine determineClusters

  use constants
  use cgnsGrid

  implicit none

  ! Working variables
  integer(kind=intType) :: numNodes, blockID, clusterID, i
  logical :: nodesAvailable

  ! Initialize the cluster of each of the CGNSDoms to 0
  do i=1, cgnsNDom
     cgnsDoms(i)%cluster = 0
  end do

  ! Initialize cluster counter
  clusterID = 0

  ! Initialize counter of classified nodes
  blockID = 0

  ! Initialize variable to state that we have unclassified nodes
  nodesAvailable = .True. 

  ! Loop until all nodes are checked
  do while (nodesAvailable)

     ! Find position of the available node
     nodesAvailable = .false.

     do while ((.not. nodesAvailable) .and. (blockID .lt. cgnsnDom))
        blockID = blockID + 1 ! Increment counter
        if (cgnsDoms(blockID)%cluster .eq. 0) then
           nodesAvailable = .true.
        end if
     end do

     ! If we have nodes available, we start the search
     if (nodesAvailable) then
        clusterID = clusterID + 1 ! Increment the running cluser counter
        cgnsDoms(blockID)%cluster = clusterID
        call clusterSearch(blockID)
     end if

  end do

contains

  recursive subroutine clusterSearch(blockID)

    ! This is the recursive part of cluster search
    implicit none

    ! Subroutine inputs
    integer(kind=intType), intent(in) :: blockID

    ! Working variables
    integer(kind=intTYpe) :: clusterID, connID, connBlock

    ! Get the cluster ID from the reference node
    clusterID = cgnsDoms(blockID)%cluster

    ! Loop over all connections of this node
    do connID = 1, cgnsDoms(blockID)%n1to1

       connBlock = cgnsDoms(blockID)%conn1to1(connID)%donorBlock

       ! Check if connected node is already classified
       if (cgnsDoms(connBlock)%cluster == 0) then
          cgnsDoms(connBlock)%cluster = clusterID ! Assign node to the same cluster
          call clusterSearch(connBlock) ! Start search on the new node

       else if (cgnsDoms(connBlock)%cluster .ne. clusterID) then ! Check symmetry
          print *,'Non-symmetric connection between CGNS blocks:', blockID, ' and', connBlock
          stop
       end if
    end do

  end subroutine clusterSearch
end subroutine determineClusters
