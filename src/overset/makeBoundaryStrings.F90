subroutine makeGapBoundaryStrings(level, sps)

  use adtAPI
  use blockPointers
  use bctypes
  use communication
  use overset
  use stringOps
  use kdtree2_module
  use adjointvars
  use wallDistanceData, only : xVolumeVec, IS1
  implicit none

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps

  ! Working
  integer(kind=intType) :: i, j, k, nn, mm, ii, jj, kk, c, e,  idx
  integer(kind=intType) :: i1, i2, j1, j2, iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: i3, i4, j3, j4
  integer(kind=intType) :: iStart, iSize, ierr, iProc, firstElem, curElem
  integer(kind=intType) :: below, above, left, right, nNodes, nElems
  integer(kind=intType) :: patchNodeCounter, nZipped, gc
  integer(kind=intType), dimension(:), allocatable :: nElemsProc, nNodesProc
  integer(kind=intType), dimension(:, :), pointer :: gcp, gnp
  real(kind=realType), dimension(:, :, :), pointer :: xx
  real(kind=realType), dimension(3) :: s1, s2, s3, s4, v1, v2, v3, v4, x0
  real(kind=realType) ::  fact, dStar, curDist, minDist, edgeLength
  logical :: isWallType

  real(kind=realType), dimension(:, :, :), allocatable :: patchNormals
  real(kind=realType), dimension(:, :), allocatable :: patchH
  integer(kind=intType), dimension(:), allocatable :: epc, surfaceSeeds, inverse
  logical,  dimension(:), allocatable :: badString
  type(oversetString), dimension(:), allocatable :: localStrings
  type(oversetString), dimension(:), allocatable :: globalStrings
  type(oversetString) :: master, pocketMaster
  type(oversetString), pointer :: stringsLL, str
  type(oversetString), dimension(:), allocatable, target :: strings

  integer(kind=intType) :: nFullStrings, nALloc, nUnique, nSearch
  type(kdtree2_result), allocatable, dimension(:) :: results
  logical :: checkLeft, checkRight, concave, nodeInFrontOfEdges
  logical :: checkLeft2, checkRight2, concave2
  logical :: leftOK, rightOK, overlappedEdges, overlappedEdges2
  real(kind=realType) :: timeA,  pt(3),   v(3), cosTheta,  cutOff, dist, maxH, ratio
  real(kind=realType), dimension(3) :: xj, xjp1, xjm1, normj
  real(kind=realType), dimension(3) :: xk, xkp1, xkm1, normk
  real(kind=realType), dimension(3) :: myPt, otherPt, eNorm
  integer(kind=intTYpe) :: otherID, otherIndex, closestOtherIndex, closestOtherString
  integer(kind=intType) :: id, index
  integer(kind=intType), dimension(:), allocatable :: nodeIndices, cellIndices
  integer status(MPI_STATUS_SIZE) 
  integer :: nSelfZipTris

  ! Wall search related
  integer(kind=intType) :: ncells
  type(oversetWall), dimension(:), allocatable, target :: walls
  type(oversetWall),  target :: fullWall
  character(80) :: fileName

  ! Loop over the wall faces counting up the edges that stradle a
  ! compute cell and a blanked (or interpolated) cell. 

  allocate(epc(nClusters)) ! epc = elementsPerCluster
  epc = 0

  ! Get a (very) large overestimate of the total number of edges in a
  ! cluster: as twice the  number of nodes. 
  domainLoop: do nn=1, nDom
     call setPointers(nn, level, sps)
     call getWallSize(nNodes, nElems, .False.)
     c = clusters(cumDomProc(myid) + nn)
     epc(c) = epc(c) + 2*nNodes
  end do domainLoop


  ! Allocate the spae we need in the local strings. 
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
        if (isWallType(BCType(mm))) then 

           select case (BCFaceID(mm))
           case (iMin)
              xx => x(1, :, :, :)
              gnp => globalNode(1, :, :)
              gcp => globalCell(2, :, :)
              fact = one
           case (iMax)
              xx => x(il, :, :, :)
              gnp => globalNode(il, :, :)
              gcp => globalCell(il, :, :)
              fact = -one
           case (jMin)
              xx => x(:, 1, :, :)
              gnp => globalNode(:, 1, :)
              gcp => globalCell(:, 2, :)
              fact = -one
           case (jMax)
              xx => x(:, jl, :, :)
              gnp => globalNode(:, jl, :)
              gcp => globalCell(:, jl, :)
              fact = one
           case (kMin)
              xx => x(:, :, 1, :)          
              gnp => globalNode(:, :, 1)
              gcp => globalCell(:, :, 2)
              fact = one
           case (kMax)
              xx => x(:, :, kl, :)
              gnp => globalNode(:, :, kl)
              gcp => globalCell(:, :, kl)
              fact = -one
           end select
           
           ! Need to reverse once more for a left-handed block
           if (.not. rightHanded) then 
              fact = -fact
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

                 s1 = s1/norm2(s1)
                 s2 = s2/norm2(s2)
                 s3 = s3/norm2(s3)
                 s4 = s4/norm2(s4)

                 ! Average and do final normalization including
                 ! correcting for inward normals. 
                 s1 = fourth*(s1 + s2 + s3 + s4)
                 patchNormals(:, i, j) = s1/norm2(s1)*fact
                 
                 ! Get the maximum edge length for this node. Use the
                 ! 4 diagonal nodes:
                 v1 = xx(i+2, j+2, :) - x0
                 v2 = xx(i  , j+2, :) - x0
                 v3 = xx(i  , j  , :) - x0
                 v4 = xx(i+2, j  , :) - x0

                 patchH(i, j) = max(norm2(v1), norm2(v2), norm2(v3), norm2(v4))

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
                       localStrings(c)%nodeData(1:3, 2*e-1) = xx(i1+1, j1+1, :)
                       localStrings(c)%nodeData(1:3, 2*e  ) = xx(i2+1, j2+1, :)

                       v1 = xx(i1+1, j1+1, :) - xx(i3+1, j3+1, :)
                       v1 = v1 / norm2(v1)

                       v2 = xx(i2+1, j2+1, :) - xx(i4+1, j4+1, :)
                       v2 = v2 / norm2(v2)

                       ! Perpendicular vector
                       localStrings(c)%nodeData(7:9, 2*e-1) = v1
                       localStrings(c)%nodeData(7:9, 2*e  ) = v2

                       ! Averaged node normal
                       localStrings(c)%nodeData(4:6, 2*e-1) = patchNormals(:, i1, j1)
                       localStrings(c)%nodeData(4:6, 2*e  ) = patchNormals(:, i2, j2)

                       ! Surface deviation estimation
                       localStrings(c)%nodeData(10, 2*e-1) = patchH(i1, j1)
                       localStrings(c)%nodeData(10, 2*e  ) = patchH(i2, j2)

                       ! Index of global node
                       localStrings(c)%intNodeData(1, 2*e-1) = gnp(i1+1, j1+1)
                       localStrings(c)%intNodeData(1, 2*e  ) = gnp(i2+1, j2+1)

                       ! Cluster of the node
                       localStrings(c)%intNodeData(2, 2*e-1) = c
                       localStrings(c)%intNodeData(2, 2*e  ) = c

                       ! Family ID of node. Not implemented yet. 
                       localStrings(c)%intNodeData(3, 2*e-1) = -1
                       localStrings(c)%intNodeData(3, 2*e  ) = -1

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
                       localStrings(c)%nodeData(1:3, 2*e-1) = xx(i1+1, j1+1, :)
                       localStrings(c)%nodeData(1:3, 2*e  ) = xx(i2+1, j2+1, :)

                       v1 = xx(i1+1, j1+1, :) - xx(i3+1, j3+1, :)
                       v1 = v1 / norm2(v1)

                       v2 = xx(i2+1, j2+1, :) - xx(i4+1, j4+1, :)
                       v2 = v2 / norm2(v2)

                       ! Perpendicular vector
                       localStrings(c)%nodeData(7:9, 2*e-1) = v1
                       localStrings(c)%nodeData(7:9, 2*e  ) = v2

                       ! Averaged node normal
                       localStrings(c)%nodeData(4:6, 2*e-1) = patchNormals(:, i1, j1)
                       localStrings(c)%nodeData(4:6, 2*e  ) = patchNormals(:, i2, j2)

                       ! Surface deviation estimation
                       localStrings(c)%nodeData(10, 2*e-1) = patchH(i1, j1)
                       localStrings(c)%nodeData(10, 2*e  ) = patchH(i2, j2)

                       ! Index of global node
                       localStrings(c)%intNodeData(1, 2*e-1) = gnp(i1+1, j1+1)
                       localStrings(c)%intNodeData(1, 2*e  ) = gnp(i2+1, j2+1)

                       ! Cluster of the node
                       localStrings(c)%intNodeData(2, 2*e-1) = c
                       localStrings(c)%intNodeData(2, 2*e  ) = c

                       ! Family ID of node. Not implemented yet.
                       localStrings(c)%intNodeData(3, 2*e-1) = -1
                       localStrings(c)%intNodeData(3, 2*e  ) = -1

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
  
  ! Allocate the gloabl list of strings on the root proc
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

     call MPI_Gather(localStrings(c)%nElems, 1, sumb_integer, nElemsProc(1:nProc), 1, sumb_integer, 0, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     call MPI_Gather(localStrings(c)%nNodes, 1, sumb_integer, nNodesProc(1:nProc), 1, sumb_integer, 0, &
          sumb_comm_world, ierr)
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
              call MPI_Recv(globalStrings(c)%nodeData(:, iStart:iEnd), iSize*10, sumb_real, iProc, iProc, &
                   sumb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              call MPI_Recv(globalStrings(c)%intNodeData(:, iStart:iEnd), iSize*3, sumb_integer, iProc, iProc, &
                   sumb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              ! ----------- Element sized arrays -------------
              iStart = nElemsProc(iProc) + 1
              iEnd =   nElemsProc(iProc+1)
              iSize = iEnd - iStart + 1
              call MPI_Recv(globalStrings(c)%conn(:, iStart:iEnd), iSize*2, sumb_integer, iProc, iProc, &
                   sumb_comm_world, status, ierr)
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
           call MPI_Send(localStrings(c)%nodeData, 10*localStrings(c)%nNodes, sumb_real, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
           
           call MPI_Send(localStrings(c)%intNodeData, 3*localStrings(c)%nNodes, sumb_integer, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           ! ----------- Element sized arrays -------------
           call MPI_Send(localStrings(c)%conn, 2*localStrings(c)%nElems, sumb_integer, 0, myid, &
                sumb_comm_world, ierr)
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
  call buildClusterWalls(level, sps, .False., walls)

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

  ! Now all the procs have fullWall info. 
  ! Note: this is overkill, since only proc 0 needs them for zipper 
  !       triangle containment search.

  !! Debug fullWalls
  !! ----------------------
  !write (fileName,"(a,I2.2,a)") "fullwall_", myid, ".dat"
  !open(unit=101,file=trim(fileName),form='formatted')
  !write(101,*) 'TITLE = "mywalls"'
  !write(101,*) 'Variables = "X", "Y", "Z"'
  !write(101,*) "Zone T=fullwall"
  !write (101,*) "Nodes = ", fullWall%nNodes, " Elements= ", fullWall%nCells, " ZONETYPE=FEQUADRILATERAL"
  !write(101, *) "DATAPACKING=POINT"
  !do i=1, fullWall%nNodes
  !   write(101, '(3(E20.12,x))')fullWall%x(1:3,i)
  !end do
  !do i=1, fullWall%nCells
  !   write(101, '(4(I5,x))') fullWall%conn(1:4, i)
  !end do
  !close(101)
  !! ----------------------
  

  ! ---------- End wall data accumulation in root proc -------------
  
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

     ! We are now left with just the single "master" string. Create
     ! the node to element data structure for master.

     call  createNodeToElem(master)
     
     ! The next step is to create ordered strings based on the
     ! connectivity. This is a purely logical operation. We don't know
     ! how many actual strings we will need so we will use a linked
     ! list as we go. 
     
     ! Allocate some additional arrays we need for doing the chain
     ! searches. 
     nElems = master%nElems
     nNodes = master%nNodes
     allocate(master%elemUsed(nElems), master%subStr(2, nElems), &
          master%cNodes(2, nNodes))

      master%cNodes = 0
      master%elemUsed = 0
      curElem = 1
      nFullStrings = 0
      do while (curElem < master%nElems)

        ! Arbitrarily get the first node for my element:
        iStart = master%conn(1, curElem)
        nElems = master%nte(1, iStart)

        ! ----------------------
        ! First side of chain:
        ! ----------------------
        firstElem = master%nte(2, iStart)
        master%subStr(1, 1) = firstElem
        call doChain(master, iStart, 1)

        ! ----------------------
        ! Second side of chain:
        ! ----------------------
        if (nElems > 1) then 
           firstElem = master%nte(3, iStart)

           ! Make sure the second one wasn't wrapped around on a
           ! periodic chain
           if (master%elemUsed(firstElem) == 0) then 

              master%subStr(2, 1) = firstElem
              call doChain(master, iStart, 2)
              call combineChainBuffers(master)
           end if
        end if

        ! We now have a boundary string stored in master%subString(1,
        ! :nSubStr(1)). These are actually the element numbers of the
        ! master that form a continuous chain.
        
        ! Create or add a new string to our linked list
        ! "stringsLL".
        if (nFullStrings == 0) then 
           allocate(stringsLL)
           nFullStrings = 1
           stringsLL%next => stringsLL
           str => stringsLL
        else
           allocate(str%next)
           str%next%next => stringsLL
           str => str%next
           nFullStrings = nFullStrings + 1 
        end if

        ! Create a substring from master based on the elements we
        ! have in the buffer
        call createSubStringFromElems(master, str, nFullStrings)

        ! Scan through until we find the next unused element:
        do while(master%elemUsed(curElem) == 1 .and. curElem < master%nElems) 
           curElem = curElem + 1
        end do
     end do

     ! Put the strings into an regular array which will be easier to
     ! manipulate.
     allocate(strings(nFullStrings))
     str => stringsLL
     i = 0 
     do while (i < nFullStrings)
        i = i + 1
        strings(i) = str ! This is derived type assigment. 
        call nullifyString(str)
        str => str%next
     end do

     ! =================== Initial Self Zipping -=================
     master%myID = 99

     ! Allocate space for the maximum number of directed edges. This
     ! is equal to the initial number of edges (nElems) plus 3 times
     ! the number of triangles we will add, which is also nElems. Now,
     ! we will probably not actualy have that many since we save a
     ! triangle and 3 edges for every self zip that is
     ! applied. Therefore we know this will always be enough 
     allocate(master%edges(4*master%nElems))

     master%nEdges = 0

     do i=1, nFullStrings
        str => strings(i)
        do j=1, str%nElems
           master%nEdges = master%nEdges + 1
           master%edges(master%nEdges)%n1 = str%p%conn(1, str%pElems(j)) !<-- first node
           master%edges(master%nEdges)%n2 = str%p%conn(2, str%pElems(j)) !<-- second node
        end do
     end do 

    ! ! Reset master%elemused so that it can be used while zipping.
    !  master%elemUsed = 0

     ! Allocate space for the triangles. Again, this can be at most,
     ! nElems, but can be smaller due to self zipping. 
     allocate(master%tris(3, master%nElems))
     master%nTris = 0

     ! Build the master tree
     master%tree => kdtree2_create(master%x, sort=.True.)

     ! Loop over the full strings and try to self zip. 
     do i=1, nFullStrings
        zipperLoop: do j=1, 5
           if (j== 1) then
              cutOff = 120_realType
           else
              cutOff = 90_realType
           end if
           call selfZip(strings(i), cutOff, nZipped)
           if (nZipped == 0) then 
              exit zipperLoop
           end if
        end do zipperLoop
     end do
     nSelfZipTris = master%nTris
     print*,' nSelfZipTris ', nSelfZipTris

     ! Now that we have self zipped the edges, allocated space for
     ! elemUsed. The remaining edges on each string must be either
     ! cross zipped or pocket zipped
     do i=1, nFullStrings
        allocate(strings(i)%elemUsed(strings(i)%nElems))
     end do
     
     ! =============================================================

     ! Now make we determine the nearest point on another substring
     ! for each point. 
     nAlloc = 50
     allocate(results(nAlloc))

     ! Loop over the fullStrings
     do i=1, nFullStrings
        str => strings(i) ! Easier readability

        ! Allocate space for otherID as it is not done yet
        allocate(str%otherID(2, str%nNodes))
        str%otherID = -1

        ! Loop over my nodes and search for it in master tree
        nodeLoop:do j=1, str%nNodes

           ! Reinitialize initial maximum number of neighbours
           nSearch = 50
           
           ! We have to be careful since single-sided chains have only
           ! 1 neighbour at each end. 

           call getNodeInfo(str, j, checkLeft, checkRight, concave, &
                xj, xjm1, xjp1,  normj)
     
           outerLoop: do
              minDist = large
              closestOtherIndex = -1
              call kdtree2_n_nearest(master%tree, xj, nSearch, results)

              ! Only check edges connected to nodes within the
              ! distance the maximum element size of my self or the
              ! closest node. We put in a fudge factor of 1.5. 

              innerLoop: do k=1, nSearch

                 ! Since we know the results are sorted, if the
                 ! distance(k) > than our current minDist, we can stop
                 ! since there is no possible way that any of the
                 ! remaining points can be closer given that the modified
                 ! D* is always larger than the original D
                 
                 ! Extract current information to make things a little
                 ! easier to read
                 curDist = sqrt(results(k)%dis)
                 idx = results(k)%idx
                 pt = master%x(:, idx)

                 ! --------------------------------------------- 
                 ! Exit Condition: We can stop the loop if the current
                 ! uncorrected distance is larger than our current
                 ! minimum. This guarantees the minimum corrected
                 ! distance is found.
                 ! ---------------------------------------------

                 if (curDist > minDist) then 
                    exit outerLoop
                 end if

                 ! ---------------------------------------------
                 ! Check 1: If the node we found isn't on our
                 ! substring. we don't need to do anything
                 ! ---------------------------------------------

                 if (master%cNodes(1, idx) == str%myID) then 
                    cycle innerLoop 
                 end if
                 
                 ! ---------------------------------------------

                 ! Check 1b: If the node we found has been removed due
                 ! to self zipping, we can just keep going
                 ! --------------------------------------------
                 if (master%cNodes(2, idx) == 0) then 
                    cycle innerLoop 
                 end if

                 ! The first time we make it here, idx will be the
                 ! index of the closest node on another string that
                 ! isn't me. 
                 if (closestOtherIndex == -1) then 
                    closestOtherString = master%cNodes(1, idx)
                    closestOtherIndex = master%cNodes(2, idx)
                 end if

                 ! --------------------------------------------- 
                 ! Check 2: Check if the node we found violates the
                 ! the "in front" test. For a concave corner TWO
                 ! triangle areas formed by the point and the two
                 ! edges must be positive. For a convex corner only
                 ! one of the triangle areas needs to be positive.
                 ! ---------------------------------------------
                 if (.not. nodeInFrontOfEdges(pt, concave, checkLeft, checkRight, &
                      xj, xjm1, xjp1, normj)) then 
                    cycle innerLoop
                 end if

                 ! --------------------------------------------- 
                 ! Check 3: This is the *reverse* of check 2: Is the
                 ! node we're searching for visible from the potential
                 ! closest other node. 
                 ! ---------------------------------------------
                 otherID = master%cNodes(1, idx)
                 otherIndex = master%cNodes(2, idx)

                 call getNodeInfo(strings(otherID), otherIndex, checkLeft2, &
                      checkRight2, concave2, xk, xkm1, xkp1,  normk)
     
                 if (.not. nodeInFrontOfEdges(xj, concave2, checkLeft2, &
                      checkRight2, xk, xkm1, xkp1, normk)) then 
                    cycle innerLoop
                 end if
                 
                 ! --------------------------------------------- 
                 ! Check 4a: Check if the potential node intersects
                 ! itself.
                 ! ---------------------------------------------
                 if (overlappedEdges(str, j, pt)) then 
                    cycle
                 end if

                 ! --------------------------------------------- 
                 ! Check 4b: OR if the other node would have to
                 ! intersect *ITSELF* to get back to me. This is used 
                 ! to catch closest points crossing over thin strips. 
                 ! ---------------------------------------------
                 
                 if (overlappedEdges(strings(otherID), otherIndex, xj)) then 
                    cycle
                 end if

                 ! --------------------------------------------- 
                 ! Check 4c: Make sure it doesn't inersect the closest
                 ! string if that happens to be different from the
                 ! cloest one.  string. This should only check very
                 ! rare cases the other checks miss.
                 ! ---------------------------------------------
                 
                 if (otherID /= closestOtherString) then 
                    if (overlappedEdges2(&
                         strings(closestOtherString), xj, normj, pt)) then 
                       cycle
                    end if
                 end if
                    
                 ! --------------------------------------------- 
                 ! Check 5: Now that the point has passed the previous
                 ! checks, we can compute the agumented distance
                 ! function and see if it better than the exisitng min
                 ! distance.
                 ! ---------------------------------------------

                 ! Now calculate our new distance
                 v =  pt - xj
                 v = v/norm2(v)

                 ! Recompute the distance function
                 cosTheta = abs(dot_product(normj, v))
                 
                 ! Update distFunction 
                 dStar = curDist / (max(1-cosTheta, 1e-6))
                 
                 if (dStar < minDist) then 
                    ! Save the string ID and the index.
                    minDist = dStar
                    str%otherID(:, j) = master%cNodes(:, idx)
                 end if
              end do innerLoop

              ! If we have already searched the max, we have to quit the loop
              if (nSearch == master%Nnodes) then 
                 exit outerLoop
              end if

              ! We are not 100% sure that we found the minium
              ! yet. Make nAlloc twice as big and start over. 
              nSearch = nSearch * 2
              nSearch = min(nSearch, master%nNodes)
              if (nSearch > nAlloc) then 
                 deallocate(results)
                 nAlloc = nAlloc*2
                 allocate(results(nAlloc))
              end if
           end do outerLoop
        end do nodeLoop
     end do

     print *,'search time:', mpi_wtime()-timea

     ! ---------------------------------------------------------------
     ! Xzip 2: Call the actual Xzip by providing all gap strings data.
     ! ---------------------------------------------------------------
     call makeCrossZip(master, strings, nFullStrings)
  
     ! =============== DEBUGGING =================
     call writeOversetTriangles(master, "fullTriangulation.dat")

     print *, 'nFullStrings:', nFullStrings
     open(unit=101, file="fullGapStrings.dat", form='formatted')
     write(101,*) 'TITLE = "Gap Strings Data" '
     write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
          "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
     do i=1, nFullStrings  
      call writeOversetString(strings(i), strings, nFullStrings, 101)
     end do
     close(101)
     ! ===========================================

     ! ===============================================================
     ! Do pocket zipping
     ! ---------------------------------------------------------------
     ! Sort through zipped triangle edges and the edges which have not
     ! been used twice (orphan edges) will be ultimately gathered to 
     ! form polygon pockets to be zipped.
     call makePocketZip(master, strings, nFullStrings, pocketMaster)


     ! -------------------------------------------------------------
     ! Perform comm data preparation for force integration on zipper
     ! triangles. Do containment search for zipper triangle
     ! containement in primal wall quad cells. Here 'surfCellID'
     ! store the indices of the global cells containing zipper
     ! triangle's cell center.
     ! (Look at ../wallDistance/determineWallAssociation.F90).

     call determineZipperWallAssociation(master, pocketMaster, fullWall)
     ! -------------------------------------------------------------

     ! (3 nodes per triangle and 3 DOF per ndoe)
     allocate(nodeIndices(3*3*(master%ntris + pocketMaster%ntris)))
     print*,'Triangles: (master,pocket,total) ',master%nTris, pocketMaster%nTris, master%nTris+pocketMaster%nTris
     
     ii = 0
     do i=1, master%nTris
        do j=1, 3 ! Triangle node loop
           jj = master%tris(j, i)
           do k=1, 3 ! DOF Loop 
              ! Zero-based ordering for petsc
              nodeIndices(9*ii+3*(j-1)+k) = 3*master%ind(jj) + k-1 
           end do
        end do
        ii = ii + 1
     end do

     do i=1, pocketMaster%nTris
        do j=1, 3 ! Triangle node loop
           jj = pocketMaster%tris(j, i)
           do k=1, 3 ! DOF Loop 
              ! Zero-based ordering for petsc
              nodeIndices(9*ii+3*(j-1)+k) = 3*pocketMaster%ind(jj) + k-1 
           end do
        end do
        ii = ii + 1
     end do

     ! Save indices of the primal quad cells containing cell center
     ! of each zipper triangle. Comes from containment search 
     ! routine. surfCellID are global cell IDs of real primal cells.

     ! Now the indices of the "global traction vector"
     allocate(cellIndices(3*(master%nTris + pocketMaster%nTris)))

     ii = 0
     do i=1, master%nTris
        do k=1, 3 ! DOF Loop 
           ! Zero-based ordering for petsc
           cellIndices(3*ii + k) = 3*master%surfCellID(i) + k-1 
        end do
        ii = ii + 1
     end do

     do i=1, pocketMaster%nTris
        do k=1, 3 ! DOF Loop 
           ! Zero-based ordering for petsc
           cellIndices(3*ii + k) = 3*pocketMaster%surfCellID(i) + k-1 
        end do
        ii = ii + 1
     end do
  else

     ! Other procs don't get any triangles :-(
     allocate(nodeIndices(0))
     allocate(cellIndices(0))
     
  end if
  
  ! Do not need walls and fullWall, deallocate them.
  ! ------------------------------------------------
  do i=1, nClusters
     deallocate(walls(i)%x, walls(i)%conn, walls(i)%ind)
     call destroySerialQuad(walls(i)%ADT)
  end do
  deallocate(walls)

  deallocate(fullWall%x, fullWall%conn, fullWall%ind)
  call destroySerialQuad(fullWall%ADT)
  ! ------------------------------------------------

  ! This is the vector we will scatter the nodes into. 
  call VecCreateMPI(SUMB_COMM_WORLD, size(nodeIndices), PETSC_DETERMINE, &
       zipperNodes, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPI(SUMB_COMM_WORLD, 3*nNodesLocal(1), PETSC_DETERMINE, &
       globalNodes, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISCreateGeneral(sumb_comm_world, size(nodeIndices), &
       nodeIndices, PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterCreate(globalNodes, IS1, zipperNodes, PETSC_NULL_OBJECT, &
       nodeZipperScatter, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! For the tractions we need to create the global vector as well as
  ! the local one

  call VecCreateMPI(SUMB_COMM_WORLD, 3*nCellsLocal(1), PETSC_DETERMINE, &
       globalPressureTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecSetBlockSize(globalPressureTractions, 3, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecDuplicate(globalPressureTractions, globalViscousTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPI(SUMB_COMM_WORLD, size(cellIndices), PETSC_DETERMINE, &
       zipperPressureTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecSetBlockSize(zipperPressureTractions, 3, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecDuplicate(zipperPressureTractions, zipperViscousTractions, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISCreateGeneral(sumb_comm_world, size(cellIndices), &
       cellIndices, PETSC_COPY_VALUES, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterCreate(globalPressureTractions, IS1, zipperPressureTractions,& 
       PETSC_NULL_OBJECT, tractionZipperScatter, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  !! Test output
  !! --------------
  !call PetscViewerASCIIOpen(sumb_comm_world,'indices.dat',viewer,ierr)
  !call ISView(IS1, viewer, ierr)
  !call PetscViewerDestroy(viewer,ierr)
  !! --------------

  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! Free the remaining memory
  deallocate(nElemsProc, nNodesProc, nodeIndices, cellIndices)
end subroutine makeGapBoundaryStrings


subroutine pointInTriangle(x1, x2, x3, pt, inTri) 

  use constants
  implicit none
  real(kind=realType), dimension(3), intent(in) :: x1, x2, x3, pt
  logical, intent(out) :: inTri
  
  if (sameSide(pt,x1, x2,x3) .and. sameSide(pt,x2, x1,x3) .and. sameSide(pt,x3, x1,x2)) then 
     inTri = .True. 
  else
     inTri = .false.
  end if

contains 
  function sameSide(p1, p2, a, b)
    
    implicit none
    logical :: sameSide
    real(kind=realType), dimension(3) ::p1, p2, a, b, cp1, cp2

    sameSide = .False.
    call cross_prod(b-a, p1-a, cp1)
    call cross_prod(b-a, p2-a, cp2)
    if (dot_product(cp1, cp2) >= zero) then 
       sameSide = .true. 
    end if
  end function SameSide
end subroutine pointInTriangle

function positiveTriArea(p1, p2, p3, norm)

  use constants
  implicit none
  real(kind=realType), intent(in), dimension(3) :: p1, p2, p3, norm
  real(kind=realType), dimension(3) :: n
  logical :: positiveTriArea
  
  call cross_prod(p2-p1, p3-p1, n)
  if (dot_product(n, norm) > zero) then 
     positiveTriArea = .True. 
  else
     positiveTriArea = .False.
  end if
end function positiveTriArea

subroutine getNodeInfo(str, j, checkLeft, checkRight, concave, xj, xjm1, xjp1, normj)

  use overset
  implicit none

  type(oversetString) :: str
  integer(kind=intType) :: j
  logical ::checkLeft, checkRight, concave 
  real(kind=realType), dimension(3) :: xj, xjm1, xjp1, normj
  real(kind=realType), dimension(3) :: v
  checkLeft = .True. 
  checkRight = .True. 
  xj = str%x(:, j)
  normj = str%norm(:, j)
  concave = .False.
  if (str%isPeriodic) then 
     if (j > 1 .and. j < str%nNodes) then 
        xjm1 = str%x(:, j-1)
        xjp1 = str%x(:, j+1)
     else if (j == 1) then 
        xjm1 = str%x(:, str%nNodes)
        xjp1 = str%x(:, j+1)
     else if (j == str%nNodes) then 
        xjm1 = str%x(:, j-1)
        xjp1 = str%x(:, 1)
     end if
  else
     ! Not periodic. Assume the ends are concave. This will
     ! forces checking if both the left and right are ok,
     ! which since the leftOK and rightOK's default to
     ! .True., it just checks the one triangle which is what
     ! we want.
     if (j == 1) then 
        checkLeft = .False.
        concave = .True. 
     end if
     
     if (j == str%nNodes) then 
        checkRight = .False.
        concave = .True. 
     end if

     if (checkLeft)  &
          xjm1 = str%x(:, j-1)
     if (checkRight) & 
          xjp1 = str%x(:, j+1)
  end if

  if (checkLeft .and. checkRight) then 
              
     ! Determine if the point is convex or concave provided
     !  we have both neighbours.
     call cross_prod(xjm1 - xj, xjp1 - xj, v)

     if (dot_product(v, normj) > zero) then 
        concave = .True. 
     end if
  end if
  
end subroutine getNodeInfo

function nodeInFrontOfEdges(pt, concave, checkLeft, checkRight, xj, xjm1, xjp1, normj)

  use constants
  implicit none

  real(kind=realType), dimension(3), intent(in) :: pt, xj, xjm1, xjp1, normj
  logical, intent(in) :: concave, checkLeft, checkRight
  logical :: positiveTriArea, nodeInFrontOfEdges
  logical :: leftOK, rightoK
  
  nodeInFrontOfEdges = .True.
  leftOK = .True. 
  rightOK = .True. 
  if (checkLeft .and. .not. positiveTriArea(xj, xjm1, pt, normj)) then
     leftOK = .False.
  end if

  if (checkRight .and. .not. positiveTriArea(xjp1, xj, pt, normj)) then 
     rightOK = .False.
  end if
  
  if (concave) then 
     if (.not. (leftOK .and. rightOK)) then 
        nodeInFrontofEdges = .False. 
     end if
  else
     if (.not. (leftOK .or. rightOK)) then 
        nodeInFrontOfEdges = .False.
     end if
  end if
end function nodeInFrontOfEdges


function overlappedEdges(str, j, pt)

  use overset
  implicit none

  ! Input/output
  real(kind=realType), dimension(3), intent(in) :: pt
  type(oversetString) , intent(in) :: str
  integer(kind=intType), intent(in) :: j
  logical :: overlappedEdges

  ! Working
  integer(kind=intType) :: i
  real(kind=realType), dimension(3) :: v, p1, p2, u, normA,  normB, x0, norm
  real(kind=realType) :: uNrm, x1, x2, x3, x4, y1, y2, y3, y4, idet, Px, Py
  real(kind=realType) :: u1, u2, v1, v2, w1, w2
  real(kind=realType) :: s1, s2, tmp, line(2), vec(2), tol
  overlappedEdges = .False. 
  tol = 1e-6
  ! We will conver this completely into a 2D problem by projecting
  ! everything onto the plane defined by norm. x0 is at the origin of
  ! the 2D system and the xaxis point from x0 to pt
 
  x0 = str%x(:, j)
  norm = str%norm(:, j)

  u = pt - x0
  uNrm = norm2(u)
  u = u/uNrm

  call cross_prod(norm, u, v)
  v = v /norm2(v)
  
  ! Now u,v,norm is an orthogonal coordinate system
  x1 = zero
  y1 = zero
  x2 = uNrm
  y2 = zero
  overLappedEdges = .False.

  ! Loop over the number of edges on my string
  elemLoop: do i=1, str%nElems
     ! Don't check the ones right next to me, since they will
     ! "overlap" exactly at x0

     if (str%conn(1, i) == j .or. str%conn(2, i) == j) then 
        cycle
     end if

     ! Project the two points into the plane
     p1 = str%x(:, str%conn(1, i))
     p2 = str%x(:, str%conn(2, i))

     normA = str%norm(:, str%conn(1, i))
     normB = str%norm(:, str%conn(2, i))

     ! Make sure the edges are on the same plane, otherwise this is
     ! meaningless
     if (dot_product(normA, norm) < half .or. dot_product(normb, norm) < half) then 
        cycle
     end if
     ! Project the two points onto the plane
     p1 = p1 - norm*dot_product(p1 - x0, norm)
     p2 = p2 - norm*dot_product(p2 - x0, norm)

     ! Now get the 2D coordinates
     x3 = dot_product(p1-x0, u)
     y3 = dot_product(p1-x0, v)
     x4 = dot_product(p2-x0, u)
     y4 = dot_product(p2-x0, v)

     u1 = x2 - x1
     y2 = y2 - y1

     v1 = x4 - x3
     v2 = y4 - y3
     
     w1 = x1- x3
     w2 = y1- y3
     
     s1 = (v2*w1 - v1*w2)/(v1*u2 - v2*u1)
     s2 = (u1*w2 - u2*w1)/(u1*v2 - u2*v1)

     if (s1 > tol .and. s1 < one - tol .and. s2 > tol .and. s2 < one - tol) then 
        overlappedEdges = .True. 
        exit elemLoop 
     end if
  end do elemLoop

end function overlappedEdges

function overlappedEdges2(str, pt1, norm, pt2)

  use overset
  implicit none

  ! Input/output
  real(kind=realType), dimension(3), intent(in) :: pt1, pt2, norm
  type(oversetString) , intent(in) :: str
  logical :: overlappedEdges2

  ! Working
  integer(kind=intType) :: i
  real(kind=realType), dimension(3) :: v, p1, p2, u, normA,  normB, x0
  real(kind=realType) :: uNrm, x1, x2, x3, x4, y1, y2, y3, y4, idet, Px, Py
  real(kind=realType) :: u1, u2, v1, v2, w1, w2
  real(kind=realType) :: s1, s2, tmp, line(2), vec(2), tol

  tol = 1e-6
  ! We will conver this completely into a 2D problem by projecting
  ! everything onto the plane defined by norm. x0 is at the origin of
  ! the 2D system and the xaxis point from x0 to pt
 
  x0 = pt1

  u = pt2 - x0
  uNrm = norm2(u)
  u = u/uNrm

  call cross_prod(norm, u, v)
  v = v /norm2(v)
  
  ! Now u,v,norm is an orthogonal coordinate system
  x1 = zero
  y1 = zero
  x2 = uNrm
  y2 = zero
  overLappedEdges2 = .False.

  ! Loop over the number of edges on my string
  elemLoop: do i=1, str%nElems
  
     ! Project the two points into the plane
     p1 = str%x(:, str%conn(1, i))
     p2 = str%x(:, str%conn(2, i))

     normA = str%norm(:, str%conn(1, i))
     normB = str%norm(:, str%conn(2, i))

     ! Make sure the edges are on the same plane, otherwise this is
     ! meaningless
     if (dot_product(normA, norm) < half .or. dot_product(normb, norm) < half) then 
        cycle
     end if
     ! Project the two points onto the plane
     p1 = p1 - norm*dot_product(p1 - x0, norm)
     p2 = p2 - norm*dot_product(p2 - x0, norm)

     ! Now get the 2D coordinates
     x3 = dot_product(p1-x0, u)
     y3 = dot_product(p1-x0, v)
     x4 = dot_product(p2-x0, u)
     y4 = dot_product(p2-x0, v)

     u1 = x2 - x1
     y2 = y2 - y1

     v1 = x4 - x3
     v2 = y4 - y3
     
     w1 = x1- x3
     w2 = y1- y3
     
     s1 = (v2*w1 - v1*w2)/(v1*u2 - v2*u1)
     s2 = (u1*w2 - u2*w1)/(u1*v2 - u2*v1)

     if (s1 > tol .and. s1 < one - tol .and. s2 > tol .and. s2 < one - tol) then 
        overlappedEdges2 = .True. 
        exit elemLoop 
     end if
  end do elemLoop

end function overlappedEdges2
