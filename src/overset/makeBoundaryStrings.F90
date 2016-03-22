subroutine makeGapBoundaryStrings(level, sps, clusters)

  use blockPointers
  use bctypes
  use communication
  use overset
  use stringOps
  use kdtree2_module
  implicit none

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps
  integer(kind=intType), intent(in) :: clusters(nDomTotal)

  ! Working
  integer(kind=intType) :: i, j, k, nn, mm, ii, jj, kk, c, e,  idx, nClusters
  integer(kind=intType) :: i1, i2, j1, j2, iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: iStart, iSize, ierr, iProc, firstElem, curElem
  integer(kind=intType) :: below, above, left, right, nNodes, nElems
  integer(kind=intType) :: patchNodeCounter, eBlank(4), nZipped
  integer(kind=intType), dimension(:), allocatable :: nElemsProc, nNodesProc
  integer(kind=intType), dimension(:, :), pointer :: gcp
  real(kind=realType), dimension(:, :, :), pointer :: xx
  real(kind=realType), dimension(3) :: s1, s2, s3, s4, v1, v2, v3, v4, x0
  real(kind=realType) ::  fact, dStar, curDist, minDist
  logical :: isWallType

  integer(kind=intType), dimension(:, :), allocatable :: patchNums
  real(kind=realType), dimension(:, :, :), allocatable :: patchNormals
  integer(kind=intType), dimension(:), allocatable :: epc
  integer(kind=intType), dimension(:), allocatable :: elemUsed
  type(oversetString), dimension(:), allocatable :: localStrings
  type(oversetString), dimension(:), allocatable :: globalStrings
  type(oversetString) :: master
  type(oversetString), pointer :: fullStrings, strings, str
  integer(kind=intType) :: nFullStrings, nALloc, nUnique
  type(kdtree2_result), allocatable, dimension(:) :: results
  integer(kind=intType), allocatable, dimension(:) :: nearEdgeList
  logical :: checkLeft, checkRight, concave, positiveTriArea
  logical :: leftOK, rightOK, overlappedEdges
  real(kind=realType) :: timeA,  pt(3),   v(3), cosTheta,  cutOff
  real(kind=realType), dimension(3) :: xj, xjp1, xjm1, normj
  integer status(MPI_STATUS_SIZE) 

  ! Loop over the wall faces counting up the edges that stradle a
  ! compute cell and a blanked (or interpolated) cell. 
  nClusters = maxval(clusters)
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

  ! We also need to gather all the node counts such that we can
  ! produce offsets and label every node on each patch. We will have
  ! to carry this index information forward as it will determine which
  ! of node tractions are needed for the triangle integration. 

  allocate(nNodesProc(0:nProc), nElemsProc(0:nProc))

  nNodesProc(0) = 0
  nElemsProc(0) = 0

  ! The number of nodes my proc has is sum(epc)/2 (all clusters and get rid of the *2 above)
  call MPI_allGather(sum(epc)/2, 1, sumb_integer, nNodesProc(1:nProc), 1, sumb_integer, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Finish the cumulatie form
  do i=2, nProc
     ! The 0 and 1st entry of the nEdgeProc and nNodeProc arrays are already correct:
     nNodesProc(i) = nNodesProc(i) + nNodesProc(i-1)
  end do

  ! Allocate the spae we need in the local strings. 
  allocate(localStrings(nClusters))
  do c=1, nClusters
     call nullifyString(localStrings(c))
  end do

  do c=1, nClusters
     allocate(&
          localStrings(c)%conn(2, epc(c)), localStrings(c)%x(3, 2*epc(c)), &
          localStrings(c)%norm(3, 2*epc(c)), localStrings(c)%ind(2*epc(c)))
     localStrings(c)%x = zero
     localStrings(c)%nNodes = 0
     localStrings(c)%nElems = 0
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
              gcp => globalCell(2, :, :)
              fact = one
           case (iMax)
              xx => x(il, :, :, :)
              gcp => globalCell(il, :, :)
              fact = -one
           case (jMin)
              xx => x(:, 1, :, :)
              gcp => globalCell(:, 2, :)
              fact = -one
           case (jMax)
              xx => x(:, jl, :, :)
              gcp => globalCell(:, jl, :)
              fact = one
           case (kMin)
              xx => x(:, :, 1, :)
              gcp => globalCell(:, :, 2)
              fact = one
           case (kMax)
              xx => x(:, :, kl, :)
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
           
           allocate(patchNums(iBeg:iEnd, jBeg:jEnd), &
                patchNormals(3, iBeg:iEnd, jBeg:jEnd))

           do j=jBeg, jEnd
              do i=iBeg, iEnd
                 patchNodeCounter = patchNodeCounter + 1
                 patchNums(i, j) = patchNodeCounter + nNodesProc(myid)
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
                 
              end do
           end do

           ! ------------------
           ! Check the i-edges
           ! ------------------
           do j=jBeg, jEnd       ! <------- Node loop
              do i=iBeg+1, iEnd  ! <------- Face Loop
                 if (gcp(i+1, j+1) > 0 .and. gcp(i+1, j+2) > 0) then 
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
                       else
                          i1 = i  ; j1 = j  
                          i2 = i-1; j2 = j  
                       end if

                       ! Don't forget pointer offset for xx
                       localStrings(c)%x(:, 2*e-1) = xx(i1+1, j1+1, :)
                       localStrings(c)%x(:, 2*e  ) = xx(i2+1, j2+1, :)

                       localStrings(c)%norm(:, 2*e-1) = patchNormals(:, i1, j1)
                       localStrings(c)%norm(:, 2*e  ) = patchNormals(:, i2, j2)

                       localStrings(c)%ind(2*e-1) = patchNums(i1, j1)
                       localStrings(c)%ind(2*e  ) = patchNums(i2, j2)

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
                 if (gcp(i+1, j+1) > 0 .and. gcp(i+2, j+1)> 0)then 
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
                       else
                          i1 = i  ; j1 = j-1
                          i2 = i  ; j2 = j  
                       end if

                       ! Don't forget pointer offset xx
                       localStrings(c)%x(:, 2*e-1) = xx(i1+1, j1+1, :)
                       localStrings(c)%x(:, 2*e  ) = xx(i2+1, j2+1, :)

                       localStrings(c)%norm(:, 2*e-1) = patchNormals(:, i1, j1)
                       localStrings(c)%norm(:, 2*e  ) = patchNormals(:, i2, j2)

                       localStrings(c)%ind(2*e-1) = patchNums(i1, j1)
                       localStrings(c)%ind(2*e  ) = patchNums(i2, j2)

                       localStrings(c)%conn(:, e) = (/2*e-1, 2*e/)

                    end if
                 end if
              end do
           end do
           deallocate(patchNums, patchNormals)
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

        allocate(globalStrings(c)%x(3, nNodesProc(nProc)), &
             globalStrings(c)%norm(3, nNodesProc(nProc)), &
             globalStrings(c)%ind(nNodesProc(nProc)), &
             globalStrings(c)%conn(2, nElemsProc(nProc)))

        ! Put proc 0's own nodes/normals/indices in the global list if we have any
        do i=1, localStrings(c)%nNodes
           globalStrings(c)%x(:, i) = localStrings(c)%x(:, i)
           globalStrings(c)%norm(:, i) = localStrings(c)%norm(:, i)
           globalStrings(c)%ind(i) = localStrings(c)%ind(i)
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
              call MPI_Recv(globalStrings(c)%x(:, iStart:iEnd), iSize*3, sumb_real, iProc, iProc, &
                   sumb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              call MPI_Recv(globalStrings(c)%norm(:, iStart:iEnd), iSize*3, sumb_real, iProc, iProc, &
                   sumb_comm_world, status, ierr)
              call ECHK(ierr, __FILE__, __LINE__)

              call MPI_Recv(globalStrings(c)%ind(iStart:iEnd), iSize, sumb_integer, iProc, iProc, &
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
                 globalStrings(c)%conn(:, i) = globalStrings(c)%conn(:, i) + nNodesProc(iProc-1)
              end do
           end if
        end do
     else
        ! Not root proc so send my stuff if we have anything:
        if (localStrings(c)%nElems > 0) then 
 
          ! ----------- Node sized arrays -------------
           call MPI_Send(localStrings(c)%x, 3*localStrings(c)%nNodes, sumb_real, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           call MPI_Send(localStrings(c)%norm, 3*localStrings(c)%nNodes, sumb_real, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
           
           call MPI_Send(localStrings(c)%ind, localStrings(c)%nNodes, sumb_integer, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           ! ----------- Element sized arrays -------------
           call MPI_Send(localStrings(c)%conn, 2*localStrings(c)%nElems, sumb_integer, 0, myid, &
                sumb_comm_world, ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end if
  end do

  ! Everyone is now down with the local strings
  do c=1, nClusters
     call deallocateString(localStrings(c))
  end do
  deallocate(localStrings)
  
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
     allocate(master%x(3, nNodes), master%conn(2, nElems), master%norm(3, nNodes), &
          master%ind(nNodes))

     nNodes = 0 ! This is our running counter for offseting nodes
     ii = 0
     jj = 0
     do c=1, nClusters
        do i=1, globalStrings(c)%nNodes
           ii = ii + 1
           master%x(:, ii) = globalStrings(c)%x(:, i)
           master%norm(:, ii) = globalStrings(c)%norm(:, i)
           master%ind(ii) = globalStrings(c)%ind(i)
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
        ! "fullStrings".
        if (nFullStrings == 0) then 
           allocate(fullStrings)
           nFullStrings = 1
           fullStrings%next => fullStrings
           str => fullStrings
        else
           allocate(str%next)
           str%next%next => fullStrings
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

     ! =================== Initial Self Zipping -=================
     master%myID = 99

     ! Allocate space for the maximum number of directed edges. This
     ! is equal to the initial number of edges (nElems) plus 3 times
     ! the number of triangles we will add, which is also nElems. Now,
     ! we will probably not actualy have that many since we save a
     ! triangle and 3 edges for every self zip that is
     ! applied. Therefore we know this will always be enough 
     allocate(master%edges(4*master%nElems))

     ! Allocate space for the triangles. Again, this can be at most,
     ! nElems, but can be smaller due to self zipping. 
     allocate(master%tris(3, master%nElems))
     master%nTris = 0

     ! Build the master tree
     master%tree => kdtree2_create(master%x, sort=.True.)

     ! Loop over the full strings and try to self zip. 
     str => fullStrings
     i = 0 
     do while (i < nFullStrings)
        i = i + 1
        zipperLoop: do j=1, 5
           if (j== 1) then
              cutOff = 120_realType
           else
              cutOff = 90_realType
           end if
           call selfZip(str, cutOff, nZipped)
           if (nZipped == 0) then 
              exit zipperLoop
           end if
        end do zipperLoop
        str => str%next
     end do

     ! =============================================================

     ! Now make we determine the nearest point on another substring
     ! for each point. 
     nAlloc = 50
     allocate(results(nAlloc))
     allocate(nearEdgeList(100))
     ! Loop over the fullStrings
     str => fullStrings
     i = 0 
     do while (i < nFullStrings)
        i = i + 1
       
        ! Allocate space for otherID as it is not done yet
        allocate(str%otherID(2, str%nNodes))
        allocate(str%otherX(3, str%nNodes))

        ! Loop over my nodes and search for it in master tree
        nodeLoop:do j=1, str%nNodes

           ! Reinitialize initial maximum number of neighbours
           nAlloc = 50

           ! Set the elements to not check, the two right next to me
           eBlank(1:2) = 0
           do jj=1, master%nte(1, str%pNodes(j))
              eBlank(jj) = master%nte(1+jj, str%pNodes(j))
           end do

           ! We have to be careful since single-sided chains have only
           ! 1 neighbour at each end. 
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

           outerLoop: do
              minDist = large
              kk = 0
              call kdtree2_n_nearest(master%tree, xj, nAlloc, results)
              
              ! Add the edges connected to the closest 50 nodes
              do k=1, 50
                 idx = results(k)%idx
                 do jj=1, master%nte(1, idx)
                    kk = kk + 1
                    nearEdgeList(kk) = master%nte(1+jj, idx)
                 end do
              end do
             
              innerLoop: do k=1, nAlloc

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
                 ! Check 2: Check if the node we found violates the
                 ! the "in front" test. For a concave corner TWO
                 ! triangle areas formed by the point and the two
                 ! edges must be positive. For a convex corner only
                 ! one of the triangle areas needs to be positive.
                 ! ---------------------------------------------
                 
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
                       cycle innerLoop
                    end if
                 else
                    if (.not. (leftOK .or. rightOK)) then 
                       cycle innerLoop
                    end if
                 end if
                 
                 ! --------------------------------------------- 
                 ! Check 3: Check if the vector to our point cross
                 ! over any edges in the near edge list. That will
                 ! also invalidate the canidate point. But only do it
                 ! if the edges are pointing in the same direction
                 ! ---------------------------------------------

                 ! Keep track of the 1 or 2 elements connected to the
                 ! current search node so we know not to test them
                 eBlank(3:4) = 0
                 do jj=1, master%nte(1, idx)
                    eBlank(jj+2) = master%nte(1+jj, idx)
                 end do

                 if (dot_product(normj, master%norm(:, idx)) > 0.5) then 
                    if (overlappedEdges(pt, xj, normj, master, eBlank,nearEdgeList, kk)) then 
                       cycle innerLoop
                    end if
                 end if
                 ! --------------------------------------------- 
                 ! Check 4: Now that the point has passed the previous
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
                    str%otherX(:, j) = pt
                 end if
              end do innerLoop

              ! We are not 100% sure that we found the minium
              ! yet. Make nAlloc twice as big and start over. 
              nAlloc = nAlloc * 2
              if (nAlloc > size(results)) then 
                 deallocate(results)
                 allocate(results(nAlloc))
              end if

           end do outerLoop
        end do nodeLoop

        ! Set pointer to next substring
        str => str%next
     end do

     print *,'search time:', mpi_wtime()-timea

     ! =============== DEBUGGING =================

     call writeOversetTriangles(master, "fullTriangulation.dat")

     print *, 'nFullStrings:', nFullStrings
     open(unit=101, file="fullGapStrings.dat", form='formatted')
     write(101,*) 'TITLE = "Gap Strings Data" '
     write(101,*) 'Variables = "X", "Y", "Z", "Nx", "Ny", "Nz", "ind" "gapID" "gapIndex" "otherID" "otherIndex"'
     i = 0 
     str => fullStrings
     do while (i < nFullStrings)
        i = i + 1
        call writeOversetString(str, 101)
        str => str%next
     end do
     close(101)
     ! ===========================================
  end if

  ! Free the remaining memory
  deallocate(nElemsProc, nNodesProc)
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

function overlappedEdges(pt, x0, norm, s, eBlank, nearEdgeList, n)

  use overset
  implicit none

  ! Input/output
  real(kind=realType), dimension(3), intent(in) :: pt, x0, norm
  type(oversetString) , intent(in) :: s
  integer(kind=intType), intent(in) :: n, eBlank(4)
  integer(kind=intType), dimension(n), intent(in) :: nearEdgeList
  logical :: overlappedEdges

  ! Working
  integer(kind=intType) :: i
  real(kind=realType), dimension(3) :: v, p1, p2, u, normA,  normB
  real(kind=realType) :: uNrm, x1, x2, x3, x4, y1, y2, y3, y4, idet, Px, Py
  real(kind=realType) :: u1, u2, v1, v2, w1, w2
  real(kind=realType) :: s1, s2, tmp, line(2), vec(2), tol
  overlappedEdges = .False. 
  tol = 1e-6
  ! We will conver this completely into a 2D problem by projecting
  ! everything onto the plane defined by norm. x0 is at the origin of
  ! the 2D system and the xaxis point from x0 to pt
 
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

  ! Loop over the number of edges in the list 
  elemLoop: do i=1, n
     if (eBlank(1) == nearEdgeList(i) .or. &
          eBlank(2) == nearEdgeList(i) .or. &
          eBlank(3) == nearEdgeList(i) .or. &
          eBlank(4) == nearEdgeList(i)) then 
        cycle
     end if

     ! Project the two points into the plane
     p1 = s%x(:, s%conn(1, nearEdgeList(i)))
     p2 = s%x(:, s%conn(2, nearEdgeList(i)))

     normA = s%norm(:, s%conn(1, nearEdgeList(i)))
     normB = s%norm(:, s%conn(2, nearEdgeList(i)))

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
