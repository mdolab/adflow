module gapBoundaries
contains
  subroutine makeGapBoundaryStrings(zipperFamList, level, sps, master)

    use constants
    use adtBuild, only : buildSerialQuad
    use blockPointers, only : x, globalCell, globalNOde, BCData, nBocos, &
         il, jl, kl, nDom, rightHanded, BCFaceID, BCType
    use communication, only : adflow_comm_world, myid, nProc
    use oversetData, only : oversetString, oversetWall, nClusters, clusters, cumDomProc
    use stringOps
    use utils, only : setPointers, EChk, myNorm2, cross_prod
    use oversetPackingRoutines, only : getWallSize
    use sorting, only : famInList
    implicit none

    ! Input Params
    integer(kind=intType), intent(in), dimension(:) :: zipperFamLIst
    integer(kind=intType), intent(in) :: level, sps

    ! Working
    integer(kind=intType) :: i, j, k, nn, mm, ii, jj, kk, c, e,  idx
    integer(kind=intType) :: i1, i2, j1, j2, iBeg, iEnd, jBeg, jEnd
    integer(kind=intType) :: i3, i4, j3, j4
    integer(kind=intType) :: iStart, iSize, ierr, iProc, firstElem, curElem
    integer(kind=intType) :: below, above, left, right, nNodes, nElems
    integer(kind=intType) :: patchNodeCounter, nZipped, gc
    integer(kind=intType), dimension(:), allocatable :: nElemsProc, nNodesProc
    integer(kind=intType), dimension(:, :), pointer :: gcp
    real(kind=realType), dimension(:, :, :), pointer :: xx
    real(kind=realType), dimension(3) :: s1, s2, s3, s4, v1, v2, v3, v4, x0
    real(kind=realType) ::  fact, timeA, minNorm

    real(kind=realType), dimension(:, :, :), allocatable :: patchNormals
    real(kind=realType), dimension(:, :), allocatable :: patchH
    integer(kind=intType), dimension(:), allocatable :: epc, surfaceSeeds, inverse
    logical,  dimension(:), allocatable :: badString
    type(oversetString), dimension(:), allocatable, target :: localStrings
    type(oversetString), dimension(:), allocatable, target :: globalStrings
    type(oversetString) :: master, pocketMaster
    type(oversetString), pointer :: stringsLL, str
    type(oversetString), dimension(:), allocatable, target :: strings
    integer(kind=intType) :: nFullStrings, nUnique, famID
    logical :: regularOrdering
    integer mpiStatus(MPI_STATUS_SIZE)

    ! Wall search related
    integer(kind=intType) :: ncells
    type(oversetWall), dimension(:), allocatable, target :: walls
    type(oversetWall),  target :: fullWall
    character(80) :: fileName

    ! Set small number to avoid division by zero when computing normal vectors
    minNorm = 1.0e-14

    ! Loop over the wall faces counting up the edges that stradle a
    ! compute cell and a blanked (or interpolated) cell.

    allocate(epc(nClusters)) ! epc = elementsPerCluster
    epc = 0

    ! Get a (very) large overestimate of the total number of edges in a
    ! cluster: as twice the  number of nodes.
    domainLoop: do nn=1, nDom
       call setPointers(nn, level, sps)
       call getWallSize(zipperFamList, nNodes, nElems, .False.)
       c = clusters(cumDomProc(myid) + nn)
       epc(c) = epc(c) + 2*nNodes
    end do domainLoop


    ! Allocate the space we need in the local strings.
    allocate(localStrings(nClusters))
    do c=1, nClusters
       call nullifyString(localStrings(c))
    end do

    do c=1, nClusters
       allocate(&
            localStrings(c)%conn(2, epc(c)), localStrings(c)%nodeData(10, 2*epc(c)), &
            localStrings(c)%intNodeData(3, 2*epc(c)))
       localStrings(c)%nodeData = zero
       localStrings(c)%nNodes = 0
       localStrings(c)%nElems = 0

       ! Assign string pointers immediately after allocation
       call setStringPointers(localStrings(c))

    end do
    deallocate(epc)
    ! And now loop back through the walls and add in the
    ! elems/nodes/normals/indices for each edge.

    ! Reset the elems per cluster. We will count up the actual number
    ! now.

    patchNodeCounter = 0
    domainLoop2: do nn=1, nDom
       call setPointers(nn, level, sps)
       ! The current cluster is 'c'
       c = clusters(cumDomProc(myid) + nn)

       bocoLoop: do mm=1, nBocos
          famID = BCData(mm)%famID
          if (famInList(famID, zipperFamList)) then

             select case (BCFaceID(mm))
             case (iMin)
                xx => x(1, :, :, :)
                gcp => globalCell(2, :, :)
                fact = one
                regularOrdering = .True.
             case (iMax)
                xx => x(il, :, :, :)
                gcp => globalCell(il, :, :)
                fact = -one
                regularOrdering = .False.
             case (jMin)
                xx => x(:, 1, :, :)
                gcp => globalCell(:, 2, :)
                fact = -one
                regularOrdering = .False.
             case (jMax)
                xx => x(:, jl, :, :)
                gcp => globalCell(:, jl, :)
                fact = one
                regularOrdering = .True.
             case (kMin)
                xx => x(:, :, 1, :)
                gcp => globalCell(:, :, 2)
                fact = one
                regularOrdering = .True.
             case (kMax)
                xx => x(:, :, kl, :)
                gcp => globalCell(:, :, kl)
                fact = -one
                regularOrdering = .False.
             end select

             ! Need to reverse once more for a left-handed block
             if (.not. rightHanded) then
                fact = -fact
                regularOrdering = .not. (regularOrdering)
             end if

             ! Before we go through and find the actual elems,
             ! precompute the patch numbering and node-based averaged
             ! unit normals
             jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd

             allocate(patchNormals(3, iBeg:iEnd, jBeg:jEnd), &
                  patchH(iBeg:iEnd, jBeg:jEnd))

             do j=jBeg, jEnd
                do i=iBeg, iEnd
                   patchNodeCounter = patchNodeCounter + 1
                   x0 = xx(i+1, j+1, :)

                   ! Normalized normals for each surrounding face.
                   v1 = xx(i+2, j+1, :) - x0
                   v2 = xx(i+1, j+2, :) - x0
                   v3 = xx(i  , j+1, :) - x0
                   v4 = xx(i+1, j  , :) - x0

                   call cross_prod(v1, v2, s1)
                   call cross_prod(v2, v3, s2)
                   call cross_prod(v3, v4, s3)
                   call cross_prod(v4, v1, s4)

                   ! When we have an 0-grid node, two of the v vectors will be the same.
                   ! Therefore, one of the s vectors will be zero. So we define minNorm
                   ! to avoid a division by zero. This will not affect the averaged normal
                   ! vector, since we will normalize it anyway in the end.

                   s1 = s1/max(minNorm,mynorm2(s1))
                   s2 = s2/max(minNorm,mynorm2(s2))
                   s3 = s3/max(minNorm,mynorm2(s3))
                   s4 = s4/max(minNorm,mynorm2(s4))

                   ! Average and do final normalization including
                   ! correcting for inward normals.
                   s1 = fourth*(s1 + s2 + s3 + s4)
                   patchNormals(:, i, j) = s1/mynorm2(s1)*fact

                   ! Get the maximum edge length for this node. Use the
                   ! 4 diagonal nodes:
                   v1 = xx(i+2, j+2, :) - x0
                   v2 = xx(i  , j+2, :) - x0
                   v3 = xx(i  , j  , :) - x0
                   v4 = xx(i+2, j  , :) - x0

                   patchH(i, j) = max(mynorm2(v1), mynorm2(v2), mynorm2(v3), mynorm2(v4))

                end do
             end do

             ! ------------------
             ! Check the i-edges
             ! ------------------
             do j=jBeg, jEnd       ! <------- Node loop
                do i=iBeg+1, iEnd  ! <------- Face Loop
                   if (gcp(i+1, j+1) >= 0 .and. gcp(i+1, j+2) >= 0) then
                      below = max(BCData(mm)%iBlank(i, j), 0)
                      above = max(BCData(mm)%iBlank(i, j+1), 0)

                      if ((below == 0 .and. above == 1) .or. (below == 1 .and. above == 0)) then
                         localStrings(c)%nNodes = localStrings(c)%nNodes + 2
                         localStrings(c)%nElems = localStrings(c)%nElems + 1
                         e = localStrings(c)%nElems

                         ! Make sure the real cell is on the LEFT
                         ! of the directed edge

                         if (below == 0) then
                            i1 = i-1; j1 = j
                            i2 = i  ; j2 = j

                            i3 = i1; j3 = j + 1
                            i4 = i2; j4 = j + 1
                         else
                            i1 = i  ; j1 = j
                            i2 = i-1; j2 = j

                            i3 = i1; j3 = j - 1
                            i4 = i2; j4 = j - 1
                         end if

                         ! Don't forget pointer offset for xx
                         if (regularOrdering) then
                            localStrings(c)%nodeData(1:3, 2*e-1) = xx(i1+1, j1+1, :)
                            localStrings(c)%nodeData(1:3, 2*e  ) = xx(i2+1, j2+1, :)


                            ! Global index of node on reduced global surface.
                            localStrings(c)%intNodeData(1, 2*e-1) = BCData(mm)%surfIndex(i1, j1)
                            localStrings(c)%intNodeData(1, 2*e  ) = BCData(mm)%surfIndex(i2, j2)
                         else
                            localStrings(c)%nodeData(1:3, 2*e  ) = xx(i1+1, j1+1, :)
                            localStrings(c)%nodeData(1:3, 2*e-1) = xx(i2+1, j2+1, :)

                            localStrings(c)%intNodeData(1, 2*e  )= BCData(mm)%surfIndex(i1, j1)
                            localStrings(c)%intNodeData(1, 2*e-1) = BCData(mm)%surfIndex(i2, j2)

                         end if
                         v1 = xx(i1+1, j1+1, :) - xx(i3+1, j3+1, :)
                         v1 = v1 / mynorm2(v1)

                         v2 = xx(i2+1, j2+1, :) - xx(i4+1, j4+1, :)
                         v2 = v2 / mynorm2(v2)

                         ! Perpendicular vector
                         localStrings(c)%nodeData(7:9, 2*e-1) = v1
                         localStrings(c)%nodeData(7:9, 2*e  ) = v2

                         ! Averaged node normal
                         localStrings(c)%nodeData(4:6, 2*e-1) = patchNormals(:, i1, j1)
                         localStrings(c)%nodeData(4:6, 2*e  ) = patchNormals(:, i2, j2)

                         ! Surface deviation estimation
                         localStrings(c)%nodeData(10, 2*e-1) = patchH(i1, j1)
                         localStrings(c)%nodeData(10, 2*e  ) = patchH(i2, j2)


                         ! Cluster of the node
                         localStrings(c)%intNodeData(2, 2*e-1) = c
                         localStrings(c)%intNodeData(2, 2*e  ) = c

                         ! Family ID of node. Not implemented yet.
                         localStrings(c)%intNodeData(3, 2*e-1) = famID
                         localStrings(c)%intNodeData(3, 2*e  ) = famID

                         ! Connectivity
                         localStrings(c)%conn(:, e) = (/2*e-1, 2*e/)
                      end if
                   end if
                end do
             end do

             ! -----------------
             ! Check the j-edges
             ! -----------------
             do j=jBeg+1, jEnd   ! <------- Face loop
                do i=iBeg, iEnd  ! <------- Node Loop
                   if (gcp(i+1, j+1) >= 0 .and. gcp(i+2, j+1) >= 0)then
                      left = max(BCData(mm)%iBlank(i, j), 0)
                      right = max(BCData(mm)%iBlank(i+1,  j), 0)

                      if ((left == 0 .and. right == 1) .or. (left == 1 .and. right == 0)) then
                         localStrings(c)%nNodes = localStrings(c)%nNodes + 2
                         localStrings(c)%nElems = localStrings(c)%nElems + 1

                         e = localStrings(c)%nElems

                         ! Again, make sure the real cell is on the LEFT
                         ! of the directed edge
                         if (left == 0) then
                            i1 = i  ; j1 = j
                            i2 = i  ; j2 = j-1

                            i3 = i1+1; j3 = j1
                            i4 = i2+1; j4 = j2

                         else
                            i1 = i  ; j1 = j-1
                            i2 = i  ; j2 = j

                            i3 = i1-1; j3 = j1
                            i4 = i2-1; j4 = j2
                         end if

                         ! Don't forget pointer offset xx
                         if (regularOrdering)  then
                            localStrings(c)%nodeData(1:3, 2*e-1) = xx(i1+1, j1+1, :)
                            localStrings(c)%nodeData(1:3, 2*e  ) = xx(i2+1, j2+1, :)

                            ! Index of global node
                            localStrings(c)%intNodeData(1, 2*e-1) = BCData(mm)%surfIndex(i1, j1)
                            localStrings(c)%intNodeData(1, 2*e  ) = BCData(mm)%surfIndex(i2, j2)

                         else
                            localStrings(c)%nodeData(1:3, 2*e  ) = xx(i1+1, j1+1, :)
                            localStrings(c)%nodeData(1:3, 2*e-1) = xx(i2+1, j2+1, :)

                            ! Index of global node
                            localStrings(c)%intNodeData(1, 2*e  ) = BCData(mm)%surfIndex(i1, j1)
                            localStrings(c)%intNodeData(1, 2*e-1) = BCData(mm)%surfIndex(i2, j2)

                         end if

                         v1 = xx(i1+1, j1+1, :) - xx(i3+1, j3+1, :)
                         v1 = v1 / mynorm2(v1)

                         v2 = xx(i2+1, j2+1, :) - xx(i4+1, j4+1, :)
                         v2 = v2 / mynorm2(v2)

                         ! Perpendicular vector
                         localStrings(c)%nodeData(7:9, 2*e-1) = v1
                         localStrings(c)%nodeData(7:9, 2*e  ) = v2

                         ! Averaged node normal
                         localStrings(c)%nodeData(4:6, 2*e-1) = patchNormals(:, i1, j1)
                         localStrings(c)%nodeData(4:6, 2*e  ) = patchNormals(:, i2, j2)

                         ! Surface deviation estimation
                         localStrings(c)%nodeData(10, 2*e-1) = patchH(i1, j1)
                         localStrings(c)%nodeData(10, 2*e  ) = patchH(i2, j2)



                         ! Cluster of the node
                         localStrings(c)%intNodeData(2, 2*e-1) = c
                         localStrings(c)%intNodeData(2, 2*e  ) = c

                         ! Family ID of node. Not implemented yet.
                         localStrings(c)%intNodeData(3, 2*e-1) = famID
                         localStrings(c)%intNodeData(3, 2*e  ) = famID

                         ! Connectivity
                         localStrings(c)%conn(:, e) = (/2*e-1, 2*e/)
                      end if
                   end if
                end do
             end do
             deallocate(patchNormals, patchH)
          end if
       end do bocoLoop
    end do domainLoop2

    ! Before we send the gap strings to the root proc, reduce them so
    ! the root proc has a little less work to do.
    do c=1, nClusters
       call reduceGapString(localStrings(c))
    end do


    ! Allocate the global list of strings on the root proc
    if (myid == 0) then
       allocate(globalStrings(nClusters))
       do c=1, nClusters
          call nullifyString(globalStrings(c))
       end do
    end if

    ! Next for each each cluster, gather to the root the gap boundary strings

    allocate(nNodesProc(0:nProc), nElemsProc(0:nProc))

    do c=1, nClusters
       ! Now let the root processor know how many nodes/elements my
       ! processor will be sending:
       nElemsProc(0) = 0
       nNodesProc(0) = 0

       call MPI_Gather(localStrings(c)%nElems, 1, adflow_integer, nElemsProc(1:nProc), 1, adflow_integer, 0, &
            adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       call MPI_Gather(localStrings(c)%nNodes, 1, adflow_integer, nNodesProc(1:nProc), 1, adflow_integer, 0, &
            adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       if (myid == 0) then

          ! Before we can receive stuff, we need to determine the node
          ! off-sets such that the conn from the strings on each processor
          ! don't overlap.

          do i=2, nProc
             ! The 0 and 1st entry of the nEdgeProc and nNodeProc arrays are already correct:
             nNodesProc(i) = nNodesProc(i) + nNodesProc(i-1)
             nElemsProc(i) = nElemsProc(i) + nElemsProc(i-1)
          end do

          allocate(globalStrings(c)%nodeData(10, nNodesProc(nProc)), &
               globalStrings(c)%intNodeData(3, nNodesProc(nProc)), &
               globalStrings(c)%conn(2, nElemsProc(nProc)))

          ! Always set the pointers immediately after allocation
          call setStringPointers(globalStrings(c))

          ! Put proc 0's own nodes/normals/indices in the global list if we have any
          do i=1, localStrings(c)%nNodes
             globalStrings(c)%nodeData(:, i) = localStrings(c)%nodeData(:, i)
             globalStrings(c)%intNodeData(:, i) = localStrings(c)%intNodeData(:, i)
          end do

          ! Put proc 0's own elements in the global list if we have any
          do i=1, localStrings(c)%nElems
             globalStrings(c)%conn(:, i) = localStrings(c)%conn(:, i)
          end do

          ! Set my total sizes
          globalStrings(c)%nNodes = nNodesProc(nProc)
          globalStrings(c)%nElems = nElemsProc(nProc)

          ! Now receive from each of the other procs.
          do iProc=1, nProc-1
             ! Check if this proc actually has anything to send:
             if ((nElemsProc(iProc+1) - nElemsProc(iProc)) > 0) then
                iStart = nNodesProc(iProc) + 1
                iEnd =   nNodesProc(iProc+1)
                iSize = iEnd - iStart + 1

                ! ----------- Node sized arrays -------------
                call MPI_Recv(globalStrings(c)%nodeData(:, iStart:iEnd), iSize*10, adflow_real, iProc, iProc, &
                     adflow_comm_world, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                call MPI_Recv(globalStrings(c)%intNodeData(:, iStart:iEnd), iSize*3, adflow_integer, iProc, iProc, &
                     adflow_comm_world, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                ! ----------- Element sized arrays -------------
                iStart = nElemsProc(iProc) + 1
                iEnd =   nElemsProc(iProc+1)
                iSize = iEnd - iStart + 1
                call MPI_Recv(globalStrings(c)%conn(:, iStart:iEnd), iSize*2, adflow_integer, iProc, iProc, &
                     adflow_comm_world, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                ! Increment the conn we just received by the node offset:
                do i=iStart, iEnd
                   globalStrings(c)%conn(:, i) = globalStrings(c)%conn(:, i) + nNodesProc(iProc)
                end do
             end if
          end do
       else
          ! Not root proc so send my stuff if we have anything:
          if (localStrings(c)%nElems > 0) then

             ! ----------- Node sized arrays -------------
             call MPI_Send(localStrings(c)%nodeData, 10*localStrings(c)%nNodes, adflow_real, 0, myid, &
                  adflow_comm_world, ierr)
             call ECHK(ierr, __FILE__, __LINE__)

             call MPI_Send(localStrings(c)%intNodeData, 3*localStrings(c)%nNodes, adflow_integer, 0, myid, &
                  adflow_comm_world, ierr)
             call ECHK(ierr, __FILE__, __LINE__)

             ! ----------- Element sized arrays -------------
             call MPI_Send(localStrings(c)%conn, 2*localStrings(c)%nElems, adflow_integer, 0, myid, &
                  adflow_comm_world, ierr)
             call ECHK(ierr, __FILE__, __LINE__)

          end if
       end if
    end do

    ! Everyone is now done with the local strings
    do c=1, nClusters
       call deallocateString(localStrings(c))
    end do
    deallocate(localStrings)

    ! Before we perform serial code implementations, surface wall info
    ! need to be communicated to root proc (zero) too. This will be used
    ! later to identify zipper triangle containment search for identifying
    ! quad surface cell info for force integration. Search will be
    ! performed on dual surface cells.

    ! ---------- Begin wall data accumulation -------------
    allocate(walls(nClusters))
    ! Build primal quad walls ADT
    call buildClusterWalls(level, sps, .False., walls, zipperFamLIst, size(zipperFamList))


    ! Finally build up a "full wall" that is made up of all the cluster
    ! walls.

    nNodes = 0
    nCells = 0
    do i=1, nClusters
       nNodes = nNodes+ walls(i)%nNodes
       nCells = nCells + walls(i)%nCells
    end do

    allocate(fullWall%x(3, nNodes))
    allocate(fullWall%conn(4, nCells))
    allocate(fullWall%ind(nNodes))
    allocate(fullWall%indCell(nCells))

    nNodes = 0
    nCells = 0
    ii = 0
    do i=1, nClusters

       ! Add in the nodes/elements from this cluster

       do j=1, walls(i)%nNodes
          nNodes = nNodes + 1
          fullWall%x(:, nNodes) = walls(i)%x(:, j)
          fullWall%ind(nNodes) = walls(i)%ind(j)
       end do

       do j=1, walls(i)%nCells
          nCells = nCells + 1
          fullWall%conn(:, nCells) = walls(i)%conn(:, j) + ii
          fullWall%indCell(nCells) = walls(i)%indCell(j)
       end do

       ! Increment the node offset
       ii = ii + walls(i)%nNodes
    end do

    ! Finish the setup of the full wall.
    fullWall%nCells = nCells
    fullWall%nNodes = nNodes
    call buildSerialQuad(nCells, nNodes, fullWall%x, fullWall%conn, fullWall%ADT)

    ! Now all the procs have fullWall info.  Note: this is overkill,
    ! since only proc 0 needs them for zipper triangle containment
    ! search.

    ! =================================================================
    !                   Serial code from here on out
    ! =================================================================


    if (myid == 0) then
       timea = mpi_wtime()

       ! First thing we do is reduce each of the global cluster gap
       ! strings
       do c=1, nClusters
          call reduceGapString(globalStrings(c))
       end do

       ! Combine all global strings together into a masterString. First
       ! count up the sizes
       nElems = 0
       nNodes = 0
       do c=1, nClusters
          nElems = nElems + globalStrings(c)%nElems
          nNodes = nNodes + globalStrings(C)%nNodes
       end do

       call nullifyString(master)
       master%nNodes = nNodes
       master%nElems = nElems
       allocate(master%nodeData(10, nNodes), master%conn(2, nElems), &
            master%intNodeData(3, nNodes))

       ! Set the string pointers to the individual arrays
       call setStringPointers(master)

       nNodes = 0 ! This is our running counter for offseting nodes
       ii = 0
       jj = 0

       do c=1, nClusters
          do i=1, globalStrings(c)%nNodes
             ii = ii + 1
             master%nodeData(:, ii) = globalStrings(c)%nodeData(:, i)
             master%intNodeData(:, ii) = globalStrings(c)%intNodeData(:, i)
          end do

          do i=1, globalStrings(c)%nElems
             jj = jj + 1
             master%conn(:, jj) = globalStrings(c)%conn(:, i) + nNodes
          end do
          nNodes =ii
       end do

       ! Now the root is done with the global strings so deallocate that
       ! too.
       do c=1, nClusters
          call deallocateString(globalStrings(c))
       end do
       deallocate(globalStrings)
    end if

  end subroutine makeGapBoundaryStrings
end module gapBoundaries
