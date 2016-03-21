subroutine makeGapBoundaryStrings(level, sps, clusters)

  use blockPointers
  use bctypes
  use communication
  use overset
  use kdtree2_module
  implicit none

  ! Input Params
  integer(kind=intType), intent(in) :: level, sps
  integer(kind=intType), intent(in) :: clusters(nDomTotal)

  ! Working
  integer(kind=intType) :: i, j, k, nn, mm, ii, jj, c, e,  idx, nClusters
  integer(kind=intType) :: i1, i2, j1, j2, iBeg, iEnd, jBeg, jEnd
  integer(kind=intType) :: iStart, iSize, ierr, iProc, firstElem, curElem
  integer(kind=intType) :: below, above, left, right, nNodes, nElems
  integer(kind=intType) :: patchNodeCounter
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
  integer(kind=intType) :: nFullStrings, nALloc
  type(kdtree2_result), allocatable, dimension(:) :: results
  logical :: checkLeft, checkRight, isConvex
  real(kind=realType) :: timeA,  pt(3),   v(3), nDist, cosTheta
  real(kind=realType), dimension(3) :: xj, xjp1, xjm1
  integer status(MPI_STATUS_SIZE) 
  interface
     subroutine createSubString(p, s, id)
       use overset
       implicit none
       type(oversetString), target, intent(in) :: p
       type(oversetString), intent(out) :: s
       integer(kind=intType), intent(in) :: id
     end subroutine createSubString
  end interface

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

  allocate(nNodesProc(0:nProc))
  nNodesProc(0) = 0

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
  allocate(nElemsProc(0:nProc))

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
     ! too.  with global strings too
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
     
     ! ! ================= DEBUG ===============
     ! call writeOversetString(master, 99)
     ! ! =======================================

     ! Allocate some additional arrays we need for doing the chain
     ! searches. 
     nElems = master%nElems
     nNodes = master%nNodes
     allocate(master%elemUsed(nElems), master%subStr(2, nElems), &
          master%cNodes(2, nNodes), master%cElems(2, nElems))
      master%cElems = 0
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
        ! have in buffer
        call createSubString(master, str, nFullStrings)

        ! Build the tree for this string. 
         str%tree => kdtree2_create(str%x, sort=.True., rearrange=.True.)

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
        call selfZip(str)
        str => str%next
     end do

     ! =============================================================


     ! Now make we determine the nearest point on another substring
     ! for each point. 
     nAlloc = 50
     allocate(results(nAlloc))
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

           ! Reinitialize the max number of neighbours
           nAlloc = 50

           ! We have to be careful since single-sided chains have only
           ! 1 neighbour at each end. 
           checkLeft = .True. 
           checkRight = .True. 
           xj = str%x(:, j)

           ! if (str%isPeriodic) then 
           !    if (j > 1 .and. j < str%nNodes) then 
           !       xjm1 = str%x(:, j-1)
           !       xjp1 = str%x(:, j+1)
           !    else if (j == 1) then 
           !       xjm1 = str%x(:, str%nNodes)
           !       xjp1 = str%x(:, j+1)
           !    else if (j == str%nNodes) then 
           !       xjm1 = str%x(:, j-1)
           !       xjp1 = str%x(:, 1)
           !    end if
           ! else
           !    ! Not periodic
           !    if (j == 1) &
           !         checkLeft = .False.
           !    if (j == str%nNodes) &
           !         checkRight = .False.

           !    if (checkLeft) &
           !         xjm1 = str(:, j-1)
           !    if (checkRight) & 
           !         xjp1 = str(:, j+1)
           ! else if 

              



           outerLoop: do
              minDist = large

              call kdtree2_n_nearest(master%tree, xj, nAlloc, results)

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

                 ! Check if the node we found isn't on our substring
                 if (master%cNodes(1, idx) == str%myID) then 
                    cycle innerLoop
                 end if

                 ! Check if the node we found violates the the "in
                 ! front" test. For a concave corner TWO triangle
                 ! areas formed by the point and the two edges must be
                 ! positive. For a convex corner only one of the
                 ! triangle areas needs to be positive. 




                 if (curDist > minDist) then 
                    ! We can stop the loop if the current uncorrected
                    ! distance is larger than our current
                    ! minimum. This guarantees the minimum corrected
                    ! distance is found. 
                    exit outerLoop
                 end if

                 ! Now calculate our new distance
                 v =  master%x(:, idx) - xj
                 nDist = dot_product(v, str%norm(:, j))
                 v = v/norm2(v)
                 pt = master%x(:, idx) - str%norm(:, j)*nDist

                 ! Recompute the distance function
                 cosTheta = abs(dot_product(str%norm(:, j), v))
                 
                 ! Update distFunction 
                 dStar = curDist / (max(1-cosTheta, 1e-6))
                 
                 if (dStar < minDist) then 
                    ! Save the string ID and the index.
                    minDist = dStar
                    str%otherID(:, j) = master%cNodes(:, idx)
                    str%otherX(:, j) = master%x(:, idx)
                 end if
              end do innerLoop

              ! We are not 100% sure that we found the minium yet. 
              nAlloc = nAlloc *2
              print *, 'reallocing:', i, j, nAlloc
              if (nAlloc > size(results)) then 
                 deallocate(results)
                 allocate(results(nAlloc))
              end if
           end do outerLoop

        end do nodeLoop


        ! Set pointer to next substring
        str => str%next
     end do



     ! Now make we determine the nearest point on another substring
     ! for each point. 
     
     ! ! Loop over the fullStrings
     ! str => fullStrings
     ! i = 0 
     ! do while (i < nFullStrings)
     !    i = i + 1
        
     !    ! Current min distance for each node. 
     !    allocate(minDist(str%nNodes))
     !    minDist = large

     !    ! Allocate space for otherID as it is not done yet
     !    allocate(str%otherID(2, str%nNodes))
     !    allocate(str%otherX(3, str%nNodes))
     !    ! Loop over all *OTHER* substring
     !    str2 => fullStrings
     !    j = 0
     !    do while(j < nFullStrings)
     !       j = j + 1

     !       ! Obviously don't check my own string. 
     !       if (j /= i ) then 
              
     !          ! Loop over my nodes and search in the other tree.
     !          do k=1, str%nNodes
     !             !call n_nearest_to_brute_force(str2%tree, str%x(:, k), 1, indexes, distances)
     !             call kdtree2_n_nearest(str2%tree, str%x(:, k), 1, results)

     !             ! We can't use the *actual* distance. Insteady we
     !             ! use the following: D* = dist_to_projected_pt +
     !             ! (K*normal_dist)**2. K can usually be quite large. 
             
     !             v =  str2%x(:, results(1)%idx) - str%x(:, k)
     !             nDist = dot_product(v, str%norm(:, k))
     !             v = v/norm2(v)
     !             pt =str2%x(:, results(1)%idx) - str%norm(:, k)*nDist

     !             ! Recompute the distance function
     !             cosTheta = abs(dot_product(str%norm(:, k), v))

     !             ! Update distFunction 
     !             results(1)%dis = results(1)%dis / (max(1-cosTheta, 1e-6))
                 
     !             if (sqrt(results(1)%dis) < minDist(k)) then 
     !                ! Save the string ID and the index.
     !                minDist(k) = sqrt(results(1)%dis)
     !                str%otherID(:, k) = (/j, results(1)%idx/)
     !                str%otherX(:, k) = str2%x(:, results(1)%idx)
     !             end if

     !          end do
     !       end if
     !       str2 => str2%next
     !    end do

     !    deallocate(minDist)
     !    str => str%next
     ! end do

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

subroutine reduceGapString(string)

  use overset

  ! Generic routine for removing duplicate nodes on the given
  ! string. The string is returned with the nodes and connectivities
  ! adjusted accordingly.
  implicit none

  ! Input/Ouput
  type(oversetString), intent(inout) :: string

  ! Working:
  real(kind=realType) :: minEdge
  integer(kind=intType) :: nUnqiue, i, n1, n2, nUnique
  integer(kind=intType), dimension(:), allocatable :: link
  real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
  real(kind=realType), dimension(:, :), pointer :: xptr, normPtr
  integer(kind=intType) , dimension(:), pointer :: indPtr
  ! interface
  !    subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
  !      use precision
  !      use kdtree2_module
  !      implicit none
  !      integer(kind=intType), intent(in) :: N
  !      real(kind=realType), intent(in), dimension(3, N) :: pts
  !      real(kind=realType), intent(in) :: tol
  !      real(kind=realType), intent(out), dimension(3, N) :: uniquePts
  !      integer(kind=intType), intent(out), dimension(N) :: link
  !      integer(kind=intType), intent(out) :: nUnique
  !    end subroutine pointReduce
  ! end interface

  ! We will do a sort of adaptive tolernace here: Get the minium edge
  ! length and base the tolerance on that:

  minEdge = huge(one)

  do i=1, string%nElems
     n1 = string%conn(1, i)
     n2 = string%conn(2, i)
     minEdge = min(minEdge, norm2(string%x(:, n1) - string%x(:, n2)))
  end do

  allocate(link(string%nNodes), uniqueNodes(3, string%nNodes))

  call pointReduce(string%x, string%nNodes, minEdge/1000.0, uniqueNodes, link, nUnique)
  
  ! Update the connectivity to use the new set of nodes
  do i=1, string%nElems
     string%conn(1, i) = link(string%conn(1, i))
     string%conn(2, i) = link(string%conn(2, i))
  end do

  ! Reallocate the node based data to the correct size. Set pointers
  ! to original data first. 
  xPtr => string%x
  normPtr => string%norm
  indPtr => string%ind
  allocate(string%x(3, nUnique), string%norm(3, nUnique), string%ind(nUnique))

  do i=1, string%nNodes
     string%x(:, link(i)) = xPtr(:, i)
     string%norm(:, link(i)) = normPtr(:, i)
     string%ind(link(i)) = indPtr(i)
  end do
  string%nNodes = nUnique

  ! deallocate the pointer data which is the original data
  deallocate(xPtr, normPtr, indPtr, link, uniqueNodes)

end subroutine reduceGapString

recursive subroutine createNodeToElem(string)

  ! Next thing we have to do is produce the inverse of the
  ! connectivity...the nodeToElem array. Normally, each node should
  ! point to 1 element (at a boundary) or two elements for a normal
  ! part of a chain. There are however, some pathalogical cases
  ! that will result in more: 
  
  ! 1. The same edge may be added twice if it is on a boundary
  ! between two block and each add the same edge. 
  
  ! 2. Chains that are "just" touching such as middle node of the
  ! following: (double lines are chains)
  !
  !     +----------++------------+
  !     +          ||            |
  !     +  ib=1    /\   ib=0     |
  !     +          ||            |
  !     +===>>=====++=====<<=====+
  !     +          ||            |
  !     +  ib=0    \/   ib=1     |
  !     +          ||            |
  !     +----------++------------+
  ! 
  ! The center node has 4 edges eminating from it. It is part of
  ! the chain that enters from the left and cointinues to the top
  ! as well as the chain that enters from the right and continues
  ! to the bottom. 
  
  ! This algorithm should take linear time, but has lots of
  ! indirect addressing (welcome to unstructured grids). It looks
  ! like it might be quadratic, but there should never be more than
  ! 4 edges attached to each node. 

  use overset
  implicit none

  ! Input/Output
  type(oversetString) :: string

  ! Working
  integer(kind=intType) :: i, j, ii, jj, n(2), m(2), curElem, nDup
  integer(kind=intType), dimension(string%nElems) :: duplicated
  integer(kind=intType), dimension(:, :), pointer :: tmpConn
  logical :: duplicateElement

  allocate(string%nte(3, string%nNodes))
  string%nte = 0 
  duplicated = 0

  do i=1, string%nElems
     ! Node numbers we're working with:
     n = string%conn(:, i)

     ! For each node check which elements (if any) are already
     ! connected. We need to check them again the node numbers n1 and n2

     duplicateElement = .False.
     do jj=1, 2

        do j=1, string%nte(1, n(jj)) ! Loop over the element numbers already here:
           curElem = string%nte(j+1, n(jj))

           m = string%conn(:, curElem)

           if (m(1) == n(1) .and. m(2) == n(2)) then 
              duplicateElement = .True.

           else if(m(1) == n(2) .and. m(2) == n(1)) then 
              ! Element exists, but it is the wrong order...don't
              ! know what to do with this, probably an error or
              ! maybe a corner case I haven't thought of.
              call terminate("makeBoundaryString", "Inconsistnet duplicate edge.")
           end if
        end do
     end do

     if (.not. duplicateElement) then 
        do jj=1, 2
           string%nte(1, n(jj)) = string%nte(1, n(jj)) + 1
           ii = string%nte(1, n(jj))
           string%nte(ii+1, n(jj)) = i
        end do
     else
        ! Well, we've figured out that this element is actually a
        ! duplicate so we'll make a note of that
        duplicated(i) = 1
     end if
  end do

  ! If we have duplicated elements, modify the conn to adjust for this. 
  nDup = sum(duplicated) 
  if (nDup > 0) then 
     tmpConn => string%conn
     allocate(string%conn(2, string%nElems - nDup))

     j = 0
     do i=1, string%nElems
        if (duplicated(i) == 0) then 
           j = j + 1
           string%conn(:, j) = tmpConn(:, i)
        end if
     end do
     
     ! Set the new number of elements
     string%nElems = string%nElems - nDup

     ! Don't forget to deallocate the tmpConn pointer which is
     ! actually the original conn data.
     deallocate(tmpConn)

     ! Destrory nte and call myself again to get the final correct nte
     ! without the duplicates.
     deallocate(string%nte)
     call createNodeToElem(string)
  end if
end subroutine createNodeToElem

subroutine nullifyString(string)

  use overset
  implicit none
  type(oversetString) :: string

  nullify(string%x, string%norm, string%ind, string%conn, &
       string%otherID, string%nte, string%subStr, string%elemUsed)
end subroutine nullifyString

subroutine deallocateString(string)

  use overset
  implicit none
  type(oversetString) :: string
  integer(kind=intType) :: i

  if (associated(string%x)) & 
       deallocate(string%x)

 if (associated(string%norm)) & 
       deallocate(string%norm)

 if (associated(string%ind)) & 
       deallocate(string%ind)

 if (associated(string%conn)) & 
       deallocate(string%conn)

 if (associated(string%otherID)) & 
       deallocate(string%otherID)

 if (associated(string%nte)) & 
       deallocate(string%nte)

 if (associated(string%subStr)) & 
       deallocate(string%subStr)

 if (associated(string%elemUsed)) &
       deallocate(string%elemUsed)

end subroutine deallocateString

subroutine doChain(master, iStart, iSub)

  use overset
  implicit none
  ! Input/OUtput
  type(oversetString) :: master
  integer(kind=intType), intent(in) :: iStart, iSub

  ! Working
  integer(Kind=intType) :: i, j, jj, c, n1, n2, curNode, nextNode
  integer(Kind=intType) :: jElem, elem1, elem2, curElem, nextElem
  logical :: finished
  real(kind=realType), dimension(3) :: v1, v2, v3
  real(kind=realType), dimension(4) :: testDP
  integer(kind=intType) :: N

  ! The number of elements in this substring
  N = 1

  curNode = iStart

  chainLoop: do

     ! Get the currnet element
     curElem = master%subStr(iSub, N)

     ! Flag the element as used:
     master%elemUsed(curElem) = 1

     ! Get the two nodes for the current element:
     n1 = master%conn(1, curElem)
     n2 = master%conn(2, curElem)

     if (n1 == curNode) then 
        nextNode = n2
     else
        nextNode = n1
     end if

     ! Exit condition 1: Next node was our starting node:
     if (nextNode == iStart) then 
        print *,'Peroidic exit'
        exit chainLoop
     end if

     ! Exit condition 2: The next node has only 1 element, (the one
     ! we're currently on) so that means the the chain is finished
     c = master%nte(1, nextNode)
     
     if (c == 1) then 
        print *,'Single chain exit'
        exit chainLoop

     else if (c == 2) then 
        ! With c=2 this easy, just extract the two elements
        elem1 = master%nte(2, nextNode)
        elem2 = master%nte(3, nextNode)
        
        if (elem1 == curElem) then 
           nextElem = elem2
        else
           nextElem = elem1
        end if
     ! else
     !    ! Uggh, we have a bow-tie topology. These are pretty
     !    ! annonying. There are 4 elements attached, 1 of which is our
     !    ! current element. We need to find the out of the three to
     !    ! take. This is done using the directed edge and the stored
     !    ! normals. 

     !    ! v1 is the vector along the curent element
     !    v1 = master%x(:, n2) - master%x(:, n1)

     !    ! Test each of the 4 elements
     !    testDP = zero
     !    do jj=1, 4
     !       jElem = master%nte(1+jj, nextNode)
     !       if (jElem /= curElem) then 
     !          ! Not our current element. 
     !          n1 = master%conn(1, jElem)
     !          n2 = master%conn(2, jElem)
     !          v2 = master%x(:, n2) - master%x(:, n1)

     !          call cross_prod(v1, v2, v3)
     !          v3 = v3 / norm2(v3)

     !          ! Dot product with "nextNode"'s normal
     !          testDP(jj) = dot_product(master%norm(:, nextNode), v3)
     !       end if
     !    end do

     !    ! The next element is the one with the largest DP, it should
     !    ! be positive.

     !    j = maxloc(testDP, dim=1)
     !    nextElem = master%nte(1+j, nextNode)
     end if

     ! Now add the "nextElem" to our chain:
     N = N + 1
     master%subStr(iSub, N) = nextElem
     
     ! Flag this elemet as being used
     master%elemUsed(nextElem) = 1 
     
     ! Finally set the nextNode back to the current node for the next
     ! iteration
     curNode = nextNode
  end do chainLoop
  master%nSubStr(iSub) = N
end subroutine doChain

subroutine createSubString(p, s, id)

  use overset
  implicit none

  ! Input/output
  type(oversetString), target, intent(in) :: p
  type(oversetString), intent(out) :: s
  integer(kind=intType), intent(in) :: id

  ! Working 
  integer(kind=intType) :: i, j, n1, n2, k
  integer(kind=intType), dimension(:), allocatable :: nodeUsed

  ! Firstly we can set the number of elements, since we know precisely
  ! what this is:
  s%nElems = p%nSubStr(1)
  s%myID = id

  ! Next determine the number of nodes. This is done by flagging the
  ! nodes in the parent that are used by 's'
  allocate(nodeUsed(p%nNodes))
  nodeUsed = 0
  k = 0
  do i=1, s%nElems
     n1 = p%conn(1, p%subStr(1, i))
     n2 = p%conn(2, p%subStr(1, i))
     if (nodeUsed(n1) == 0) then 
        k = k + 1
        nodeUsed(n1) = k
     end if

     if (nodeUsed(n2) == 0) then 
        k = k + 1
        nodeUsed(n2) = k
     end if
     
     ! Set the child relationship information for the elements
     p%cElems(:, p%subStr(1, i)) = (/s%myID, i/)

  end do

  ! We can now set the number of nodes the substring has
  s%nNodes = k
  
  ! The number of nodes will equal the number of elements iff the
  ! string is period. Otherwise we will have 1 more node than element.
  
  n1 = p%subStr(1, 1)
  n2 = p%subStr(1, s%nElems)
  if (n1 == n2 .and. s%nNodes ==  s%nElems) then 
     s%isPeriodic = .True. 
  end if

  ! Allocate and set the node and element parent information
  allocate(s%pElems(s%nElems), s%pNodes(s%nNodes))

  do i=1, s%nElems
     s%pElems(i) = p%subStr(1, i)
  end do

  ! Now create the pNodes ("link") array such that pNodes(i) points to
  ! the node index in the parent
  j = 0
  do i=1, p%nNodes
     if (nodeUsed(i) /= 0) then 
       s%pNodes(nodeUsed(i)) = i
     end if
  end do

  ! Set the parent's cNode to point to my nodes
  do i=1, s%nNodes
     p%cNodes(:, s%pNodes(i)) = (/s%myID, i/)
  end do

  ! Now that we know the mapping between by local nodes-based
  ! quantities and the parent, we can allocate and set all the
  ! node-based quantities.
  
  allocate(s%x(3, s%nNodes), s%norm(3, s%nNodes), s%ind(s%nNodes))
  do i=1, s%nNodes
     s%x(:, i) = p%x(:, s%pNodes(i))
     s%norm(:, i) = p%norm(:, s%pNodes(i))
     s%ind(i) = p%ind(s%pNodes(i))
  end do

  ! We can now create the local conn too, *USING THE LOCAL NODE NUMBERS*
  allocate(s%conn(2, s%nElems))

  do i=1, s%nElems
     s%conn(1, i) = nodeUsed(p%conn(1, s%pElems(i)))
     s%conn(2, i) = nodeUsed(p%conn(2, s%pElems(i)))
  end do

  ! Set the pointer to my parent. 
  s%p => p

  deallocate(nodeUsed)

  ! Last thing we can do is create the nodeToElem for the substring. 
  call  createNodeToElem(s)

end subroutine createSubString

subroutine combineChainBuffers(s)

  use overset
  implicit none
  type(oversetString), intent(inout) :: s
  integer(kind=intType) :: N1, N2

  N1 = s%nSubStr(1)
  N2 = s%nSubStr(2)

  ! First reverse the direction of string 2 of the nodes we found
  s%subStr(2, 1:N2) = s%subStr(2, N2:1:-1)
              
  ! Now String 1 can be tacked on the end of string2
  s%subStr(2, N2+1:N2+N1) = s%subStr(1, 1:N1)
              
  ! And finally copied back to string1
  s%subStr(1, 1:N1+N2) = s%subStr(2, 1:N1+N2)
  s%nSubStr(1) = N1 + n2
              
end subroutine combineChainBuffers

         

subroutine selfZip(s)

  use overset
  use kdtree2_module
  implicit none
  type(oversetString), intent(inout) :: s
  
  ! Working
  integer(kind=intType) :: i, k,  N, ii, im1, ip1, nalloc, idx, nFound
  logical :: lastNodeZipper, inTri, overlapFound
  real(kind=realType), dimension(3) :: v1, v2, norm, c
  real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm
  type(kdtree2_result), dimension(:), allocatable  :: results

  ! Perform self zipping on the supplied string. The string at this
  ! point should be either peroidic or since sinded --- no multiple
  ! loops should be left. Therefore, we can count on the nodes being
  ! in order.

  cosCutoff = cos(120*pi/180)

  if (s%isPeriodic) then 
     im1 = s%nNodes
     ii = 1
     ip1 = 2
     N = s%nElems
  else
     im1 = 1
     ii = 2
     ip1 = 3
     N = s%nNodes - 2
  end if

  nAlloc = 25
  allocate(results(nAlloc))

  lastNodeZipper = .False. 
  do while (ii < N)
     
     ! Determine the anlge between the vectors
     v1 = s%x(:, ip1) - s%x(:, ii)
     v2 = s%x(:, im1) - s%x(:, ii)
     v1nrm = norm2(v1)
     v2nrm = norm2(v2)
     call cross_prod(v2, v1, norm)
     norm = norm / norm2(norm)
     
     if (dot_product(norm, s%norm(:, ii)) > zero) then 

        costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)
     
        if (costheta > cosCutoff) then 

           ! We may have a valid triangle. We need to make sure we
           ! don't overlap anyone else. 
           !
           !     +
           !     | \
           !     |   \
           !     |     c
           !     |       \
           !     |         \
           !     +----------+
           !
           ! We do a ball search based at 'c' which is just the
           ! (average of xip1 and xim1) using a radius defined as the
           ! maximum of (the distance between 'c' and 'xi', half
           ! length of xip1 to xim1)
           ! 
           c = half*(s%x(:, ip1) + s%x(:, im1))
           r2 = (c(1) - s%x(1, ii))**2 +  (c(2) - s%x(2, ii))**2 +  (c(3) - s%x(3, ii))**2

           r2 = max(r2, (s%x(1, ip1) - s%x(1, im1))**2 + (s%x(2, ip1) - s%x(2, im1))**2 + &
                (s%x(3, ip1) - s%x(3, im1))**2)

           nFound = 0
           outerLoop: do 

              call kdtree2_r_nearest(s%p%tree, c, r2, nfound, nalloc, results) 
              if (nFound < nAlloc) then 
                 exit outerLoop
              end if

              ! Allocate more space and keep going
              deallocate(results)
              nAlloc = nAlloc * 2
              allocate(results(nAlloc))
           end do outerLoop
           
           ! We can now be sure that we have all the points inside our
           ! ball. Next we proceed to systematically check them. 
           overlapFound = .False.
           nodeFoundLoop: do k=1, nFound
              ! Note that we do check nodes from our own string,
              ! except for the the three nodes we're dealing
              ! with. Remember that we are working in our parent's
              ! ording here.
              idx = results(k)%idx 

              notPartofTriangle: if (idx /= s%pNodes(im1) .and. &
                   idx /= s%pNodes(ii) .and. idx /= s%pNodes(ip1)) then 
                 
                 ! Only check if the node normal of the point we're
                 ! checking is in the same direction as the triangle. 
                 if (dot_product(s%norm(:, ii), s%p%norm(:, idx)) > zero) then 

                
                    ! Finally do the actual trianlge test
                    call pointInTriangle(s%x(:, ip1), s%x(:, ii), s%x(:, im1), &
                         s%p%x(:, idx), inTri)
                    if (inTri) then 
                       ! As soon as 1 is in the triangle, we know the
                       ! gap string is no good. 
                       overlapFound = .True. 
                       exit nodeFoundLoop
                    end if
                 end if
              end if notPartofTriangle
           end do nodeFoundLoop

           if (.not. overlapFound) then 
              ! This trialge is good!
              print *,'adding triangle:', s%myID, (/s%pNodes(im1), s%pNodes(ii),s%pNodes(ip1)/)
              s%p%nTris = s%p%nTris+ 1
              s%p%tris(:, s%p%nTris) = (/s%pNodes(im1), s%pNodes(ii),s%pNodes(ip1)/)
           end if
        end if
     end if

     ! Just shuffle along
     ii = ii + 1
     ip1 = ii + 1
     im1 = ii -1
     
  end do


end subroutine selfZip

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
