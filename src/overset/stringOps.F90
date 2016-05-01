module stringOps

  use overset

contains
  subroutine nullifyString(string)

    use overset
    implicit none
    type(oversetString) :: string

    nullify(string%nodeData, &
         string%x, &
         string%norm, &
         string%perpNorm, &
         string%h, &
         string%intNodeData, &
         string%ind, &
         string%cluster, &
         string%family, & 
         string%conn, &
         string%pNodes, &
         string%pElems, &
         string%cNodes, &
         string%otherID, &
         string%nte, &
         string%subStr, &
         string%elemUsed, &
         string%XzipNodeUsed, &
         string%tris)

  end subroutine nullifyString

  subroutine deallocateString(string)

    use overset
    implicit none
    type(oversetString) :: string
    integer(kind=intType) :: i

    if (associated(string%nodeData)) & 
         deallocate(string%nodeData)

    if (associated(string%intNodeData)) & 
         deallocate(string%intNodeData)

    if (associated(string%conn)) & 
         deallocate(string%conn)

    if (associated(string%pNodes)) &
         deallocate(string%pNodes)

    if (associated(string%pElems)) &
         deallocate(string%pElems)

    if (associated(string%cNodes)) &
         deallocate(string%cNodes)

    if (associated(string%otherID)) & 
         deallocate(string%otherID)

    if (associated(string%nte)) & 
         deallocate(string%nte)

    if (associated(string%subStr)) & 
         deallocate(string%subStr)

    if (associated(string%elemUsed)) &
         deallocate(string%elemUsed)

    if (associated(string%xZipNodeUsed)) &
         deallocate(string%xZipNodeUsed)

    if (associated(string%tris)) &
         deallocate(string%tris)

    call nullifyString(string)

  end subroutine deallocateString

  subroutine setStringPointers(string)

    use overset
    implicit none
    type(oversetString) :: string
    string%x => string%nodeData(1:3, :)
    string%norm => string%nodeData(4:6, :)
    string%perpNorm => string%nodeData(7:9, :)
    string%h => string%nodeData(10, :)

    string%ind => string%intNodeData(1, :)
    string%cluster => string%intNodeData(2, :)
    string%family => string%intNodeData(3, :)

  end subroutine setStringPointers

  subroutine reduceGapString(string)

    ! Generic routine for removing duplicate nodes on the given
    ! string. The string is returned with the nodes and connectivities
    ! adjusted accordingly.

    use overset
    implicit none

    ! Input/Ouput
    type(oversetString), intent(inout) :: string

    ! Working:
    real(kind=realType) :: minEdge
    integer(kind=intType) :: nUnqiue, i, n1, n2, nUnique
    integer(kind=intType), dimension(:), allocatable :: link
    real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
    real(kind=realType), dimension(:, :), pointer :: nodeDataPtr
    integer(kind=intType) , dimension(:, :), pointer :: intNodeDataPtr

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
    nodeDataPtr => string%nodeData
    intNodeDataPtr => string%intNodeData
    allocate(string%nodeData(10, nUnique), string%intNodeData(3, nUnique))

    ! Reset the pointers 
    call setStringPointers(string)

    do i=1, string%nNodes
       string%nodeData(:, link(i)) = nodeDataPtr(:, i)
       string%intNodeData(:, link(i)) = intNodeDataPtr(:, i)
    end do
    string%nNodes = nUnique

    ! deallocate the pointer data which is actually the original data
    deallocate(nodeDataPtr, intNodeDataPtr, link, uniqueNodes)

  end subroutine reduceGapString

  recursive subroutine createNodeToElem(string)

    ! Produce the inverse of the connectivity...the nodeToElem
    ! array. Each node should point to 1 element (at a
    ! boundary) or two elements for a normal part of a chain.

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

       ! Destroy nte and call myself again to get the final correct nte
       ! without the duplicates.
       deallocate(string%nte)
       call createNodeToElem(string)
    end if
  end subroutine createNodeToElem

  subroutine doChain(master, iStart, iSub)

    use overset
    implicit none
    ! Input/OUtput
    type(oversetString) :: master
    integer(kind=intType), intent(in) :: iStart, iSub

    ! Working
    integer(Kind=intType) :: i, j, jj, c, n1, n2, curNode, nextNode
    integer(Kind=intType) ::  elem1, elem2, curElem, nextElem
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
          exit chainLoop
       end if

       ! Exit condition 2: The next node has only 1 element, (the one
       ! we're currently on) so that means the the chain is finished
       c = master%nte(1, nextNode)

       if (c == 1) then 
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

  subroutine createSubStringFromElems(p, s, id)

    use overset
    implicit none

    ! Input/output
    type(oversetString), target, intent(in) :: p
    type(oversetString), intent(out) :: s
    integer(kind=intType), intent(in) :: id

    ! Working 
    integer(kind=intType) :: i, j, n1, n2, k
    integer(kind=intType), dimension(:), allocatable :: nodeUsed

    ! First thing we always have to do with a new string is to nullify
    ! all the poitners
    call nullifyString(s)

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
    end do

    ! We can now set the number of nodes the substring has
    s%nNodes = k

    ! The number of nodes will equal the number of elements iff the
    ! string is period. Otherwise we will have 1 more node than element.

    if (s%nNodes ==  s%nElems) then 
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

    allocate(s%nodeData(10, s%nNodes), s%intNodeData(3, s%nNodes))

    ! Set the string pointers
    call setStringPointers(s)

    do i=1, s%nNodes
       s%nodeData(:, i) = p%nodeData(:, s%pNodes(i))
       s%intNodeData(:, i) = p%intNodeData(:, s%pNodes(i))
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
  end subroutine createSubStringFromElems

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

  subroutine selfZip(s, cutOff, nZipped)

    use overset
    use kdtree2_module
    implicit none
    type(oversetString), intent(inout), target :: s
    integer(Kind=intType), intent(out) :: nZipped
    real(kind=realType), intent(in) :: cutOff
    ! Working
    integer(kind=intType) :: i, j, k,  N, ii, im1, ip1, nalloc, idx, nFound
    integer(kind=intType) :: nNodes, nElems, elem1, elem2
    logical :: lastNodeZipper, inTri, overlapFound
    real(kind=realType), dimension(3) :: v1, v2, norm, c
    real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm
    integer(Kind=intType), dimension(:), allocatable :: nodeMap, elemMap
    type(kdtree2_result), dimension(:), allocatable  :: results
    real(kind=realType), dimension(:, :), pointer :: nodeDataTmp
    integer(kind=intType), dimension(:, :), pointer :: connTmp, intNodeDataTmp
    integer(kind=intType), dimension(:), pointer :: pNodesTmp
    ! Perform self zipping on the supplied string. The string at this
    ! point should be either peroidic or since sinded --- no multiple
    ! loops should be left. Therefore, we can count on the nodes being
    ! in order.

    cosCutoff = cos(cutOff*pi/180)
    nzipped = 0

    ! Peroidic string starts at node 1, and uses node 'N' as the previous
    ! node. Single chains start at node 2 and only go to the N-1 node. 
    if (s%isPeriodic) then 
       im1 = s%nNodes
       ii = 1
       ip1 = 2
       N = s%nNodes
    else
       im1 = 1
       ii = 2
       ip1 = 3
       N = s%nNodes - 1
    end if

    nAlloc = 25
    allocate(results(nAlloc))
    allocate(nodeMap(s%nNodes), elemMap(s%nElems))
    nodeMap = 1
    elemMap = 1

    do while (ii <= N)

       ! Peroidic string at end...loop around
       if (s%isPeriodic .and. ii == N) then 
          ip1 = 1
       end if

       lastNodeZipper = .False. 

       ! Determine the anlge between the vectors
       v1 = s%x(:, ip1) - s%x(:, ii)
       v2 = s%x(:, im1) - s%x(:, ii)
       v1nrm = norm2(v1)
       v2nrm = norm2(v2)
       call cross_prod(v2, v1, norm)
       norm = norm / norm2(norm)

       if (dot_product(norm, s%norm(:, ii)) > zero) then 

          ! the dot product of the im1 and ip1 nodes have to be close
          if (dot_product(s%norm(:, ip1), s%norm(:, im1)) > 0.80) then 

             costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)

             if (costheta > cosCutoff) then 

                ! We may have a valid triangle. We need to make sure we
                ! don't overlap anyone else. 
                !
                ! xim1 +
                !      | \
                !      |   \
                !      |     c
                !      |       \
                !      |         \
                !      +----------+
                !      xi         xip1
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
                   ! This triangle is good!
                   s%p%nTris = s%p%nTris+ 1
                   s%p%tris(:, s%p%nTris) = (/s%pNodes(ip1), s%pNodes(ii),s%pNodes(im1)/)
                   lastNodeZipper = .True.
                   nZipped = nZipped + 1

                   ! Flag this node as gone
                   nodeMap(ii) = 0

                   ! Flag the two edges on either side of this node as
                   ! also being gone
                   elemMap(s%nte(2, ii)) = 0
                   elemMap(s%nte(3, ii)) = 0

                   ! Flag the two edges that got removed as being used
                   ! in the parent.
                   elem1 = s%p%nte(2, s%pNodes(ii))
                   elem2 = s%p%nte(3, s%pNodes(ii))
                   s%p%elemUsed(elem1) = 1
                   s%p%elemUsed(elem2) = 1

                   ! Add these three edges to oversetEdge's 
                   ! Edge 1:
                   s%p%nEdges = s%p%nEdges + 1
                   s%p%edges(s%p%nEdges)%n1 = s%pNodes(ip1)
                   s%p%edges(s%p%nEdges)%n2 = s%pNodes(ii)
                   ! Edge 2:
                   s%p%nEdges = s%p%nEdges + 1
                   s%p%edges(s%p%nEdges)%n1 = s%pNodes(ii)
                   s%p%edges(s%p%nEdges)%n2 = s%pNodes(im1)
                   ! Edge 3:
                   s%p%nEdges = s%p%nEdges + 1
                   s%p%edges(s%p%nEdges)%n1 = s%pNodes(im1)
                   s%p%edges(s%p%nEdges)%n2 = s%pNodes(ip1)
                end if
             end if
          end if
       end if
       if (lastNodeZipper) then 
          ! Skip the next node...we'll get it on the next pass
          ii = ii + 2
          im1 = ii -1
          ip1 = ii + 1
       else
          ! Just shuffle along
          ii = ii + 1
          ip1 = ii + 1
          im1 = ii -1
       end if
    end do

    ! Now we will modify our string to remove the elements and nodes
    ! that got knocked off due to self zipping. This way the calling
    ! process still sees the same string, it just gets a little
    ! shorter. 

    ! Save pointers to existing data
    nNodes = s%nNodes
    nElems = s%nElems
    nodeDataTmp => s%nodeData
    intNodeDataTmp => s%intNodeData
    connTmp => s%conn
    pNodesTmp => s%pNodes

    ! Convert the nodeMap which currently contains a one if the node
    ! still exists and 0 if it doesn't. This will convert it to the new
    ! node numbers. Ie nodeMap(i) gives the new node index of the
    ! shorted chain. If nodeMap(i) = 0, it is no longer part of the
    ! chain.
    j = 0
    do i=1, s%nNodes
       if (nodeMap(i) == 1) then 
          j = j + 1
          nodeMap(i) = j
       end if
    end do

    !  Update the cNodes in the parent so they point to the updated node
    ! numbers. Note that the nodes that have been eliminated, have cNode
    ! = 0, which will identify that it no longer has a child node. 
    do i=1, s%nNodes
       s%p%cNodes(:, s%pNodes(i)) = (/s%myID, nodeMap(i)/)
    end do

    ! Update the number of nodes/elems in our shorted chain. Every
    ! zipper reduces the number of nodes and number of elems by 1
    s%nNodes = s%nNodes - nZipped
    s%nElems = s%nElems - nZipped

    allocate(s%nodeData(10, s%nNodes), s%intNodeData(3, s%nNodes), &
         s%pNodes(s%nNodes), s%conn(2, s%nElems))

    ! Set the pointers for the new string
    call setStringPointers(s)

    do i=1, nNodes
       if (nodeMap(i) /= 0) then 
          s%nodeData(:, nodeMap(i)) = nodeDataTmp(:, i)
          s%intNodeData(:, nodeMap(i)) = intNodeDataTmp(:, i)
          s%pNodes(nodeMap(i)) = pNodesTmp(i)
       end if
    end do

    ! Since we know the string was in order, we can simply redo the connectivity
    do i=1, s%nElems
       s%conn(:, i) = (/i, i+1/)
    end do
    if (s%isPeriodic) then 
       s%conn(2, s%nElems) = 1
    end if

    ! Dellocate the existing memory
    deallocate(nodeDataTmp, intNodeDataTmp, connTmp, pNodesTmp)
    deallocate(nodeMap, elemMap)

    ! Recrate the node to elem
    call createNodeToElem(s)

  end subroutine selfZip

  subroutine crossZip(str1, N1, N2, str2, N3, N4)

    implicit none
    type(oversetString), intent(inout) :: str1, str2
    integer(kind=intType) :: N1, N2, N3, N4

    ! Working
    integer(kind=intType) :: stepsA, stepsB, nStepsA, nStepsB
    integer(kind=intType) :: nTriToAdd, ii, i, j, k, A, B, Ap, Bp
    real(kind=realType), dimension(3) :: ptA, ptB, ptAp, ptBp
    real(kind=realType), dimension(3) :: Aoff, Boff, ApOff, BpOff
    real(kind=realType), dimension(3) :: normA, normB, normAp, normBp
    real(kind=realType), dimension(3) :: perpA, perpB, perpAp, perpBp
    real(kind=realType), dimension(3) :: triNorm1, quadNorm1
    real(kind=realType), dimension(3) :: triNorm2, quadNorm2
    logical :: aValid, bValid, advanceA, aPreferred, area1, area2
    logical :: advanceB
    logical :: positiveTriArea, changeA, changeB
    real(kind=realType) ::  sum1, sum2, h, dpa, dpb
    !am real(kind=realType), parameter :: cutOff = 0.95*3
    real(kind=realType), parameter :: cutOff = 0.85*3
    ! First determine the the total number of triangles we will add
    ! total. It is equal to the total number of triangles on each
    ! string. This will form the index on the do loop.

    ! Str1 goes forward
    if (N2 > N1) then 
       nStepsA = N2 - N1
    else if (N2 < N1) then 
       nStepsA = N2 + str1%nNodes - N1
    else ! N1 == N2
       nStepsA = str1%nElems
    end if

    ! Str2 goes backwards
    if (N3 < N4) then 
       nStepsB = N3 + str2%nNodes - N4
    else if (N3 > N4) then 
       nStepsB = N3 - N4
    else ! N3 == N4
       nStepsB =  str2%nElems
    end if

    ! Initialize the front: 
    A = N1
    B = N3
    ptA = str1%x(:, A)
    ptB = str2%x(:, B)

    normA = str1%norm(:, A)
    normB = str2%norm(:, B)

    perpA = str1%perpNorm(:, A)
    perpB = str2%perpNorm(:, B)

    Ap = nextNode(str1, A, .True.)
    Bp = nextNode(str2, B, .False.)
    ptAp = str1%x(:, Ap)
    ptBp = str2%x(:, Bp)
    normAp = str1%norm(:, Ap)
    normBp = str2%norm(:, Bp)
    perpAp = str1%perpNorm(:, Ap)
    perpBp = str2%perpNorm(:, Bp)

    ! The number of steps we've performed in each edge
    stepsA = 0
    stepsB = 0

    ! Cross zip nodes N1 to N2 on str1 to nodes N3 to N4 on str2
    ii = 0
    do while (ii < nStepsA + nStepsB)

       aValid = .True. 
       bValid = .True. 
       ! ---------------------------------------------------------------
       ! Check 1: Point-in-Triangle test: This test considers the
       ! triangle ABA+ and determines if any of the neighbouring points
       ! on either of the two strings is contained inside the
       ! triangle. If the test is positive, A+ must be rejected. The
       ! same test is repeated for B+. 
       ! ---------------------------------------------------------------

       if (triOverlap(ptA, ptB, ptAp, str1, A, Ap) .or. &
            triOverlap(ptA, ptB, ptAp, str2, B, B)) then
          aValid = .False.
       end if

       if (triOverlap(ptA, ptB, ptBp, str1, A, A) .or. & 
            triOverlap(ptA, ptB, ptBp, str2, B, Bp)) then
          bValid = .False.
       end if

       ! ---------------------------------------------------------------
       ! Check 2: Convex quadrilaterl test: This test considers the
       ! quadrilateral ABB+A+ and determines if it is convex. For
       ! connection to point A+ to be valid, the vector areas of
       ! triangles ABA+ and BB+A+ should have the same size. For
       ! connection to B+ to be valid, the vector areas of trianges ABB+
       ! and AB+A+ should be the same sign.  NOTE THAT THIS TEST DOES NOT
       ! ACTUALLY WORK. IT IS 100% INCORRECT!!! THERE ARE CASES WHERE
       ! THE SIGN OF BOTH AREAS ARE OPPOSITE! IT CANNOT BE SAFELY USED. 
       ! ---------------------------------------------------------------

       ! area1 = positiveTriArea(ptA, ptB, ptAp, normB)
       ! area2 = positiveTriArea(ptB, ptBp, ptAp, normB)

       ! if (area1 .neqv. area2) then 
       !    aValid = .False. 
       ! end if

       ! area1 = positiveTriArea(ptA, ptB, ptBp, normA)
       ! area2 = positiveTriArea(ptAp, ptBp, ptA, normA)

       ! if (area1 .neqv. area2) then 
       !    bValid = .False. 
       ! end if

       ! Instead, check if the triangle we're going to add has a
       ! positive or negative vector area

       area1 = positiveTriArea(ptA, ptB, ptAp, normA)
       if (area1 .eqv. .False.) then 
          aValid = .False. 
       end if

       area2 = positiveTriArea(ptA, ptB, ptBp, normB)
       if (area2 .eqv. .False.) then 
          bValid = .False. 
       end if

       ! ---------------------------------------------------------------
       ! Check 3: Prism volume test: Using the surface normals,
       ! "extrude" a prisim in the direction of each surface normal and
       ! find it's volume. It is is not positive, reject the
       ! triangle. Since we don't have the node off wall, we will have
       ! to make do with the normal vectors and average cell size. We
       ! average the cell size and divide by 1000 to give an approximate
       ! offwall distance. Then we use the norm veectors to offset in
       ! that distance to produce the "off" points. 
       ! ---------------------------------------------------------------

       ! h = quarter*(str1%h(A) + str1%h(Ap) + str2%h(B) + str2%h(Bp)) / 1000
       ! AOff = ptA + normA * h
       ! BOff = ptB + Bnorm * h
       ! ApOff = ptAp + normAp * h
       ! BpOff = ptBp + normBp * h

       ! if (prismVol(A, B, Ap, Aoff, Boff, ApOff) < zero) then 
       !    aValid = .False. 
       ! end if

       ! if (prismVol(A, B, Bp, Aoff, Boff, BpOff) < zero) then 
       !    bValid = .False. 
       ! end if

       ! ---------------------------------------------------------------
       ! Check 4: Interpolation stencil test: This one isn't implemented
       ! ---------------------------------------------------------------

       ! ---------------------------------------------------------------
       ! Check 5: Surface normal compatibility test. The surface normal
       ! from the triangle should be pointing (mostly) in the same
       ! direction as the normal of the quad that this triangle shares
       ! and edge with. THIS ALSO DOES NOT WORK! What we have to do
       ! instead, is check the normal tri normal against the node
       ! normals it would be using. This is simplier and is vastly
       ! superior. 
       ! ---------------------------------------------------------------

       call cross_prod(ptB-ptA, ptAp-ptA, triNorm1)
       triNorm1 = triNorm1 / norm2(triNorm1)

       ! Compute the sum of the dot product of the nodal norms with the triNorm
       sum1 = dot_product(triNorm1, normA) + dot_product(triNorm1, normB) + &
            dot_product(triNorm1, normAp)

       call cross_prod(ptB-ptA, ptBp-ptA, triNorm2)
       triNorm2 = triNorm2 / norm2(triNorm2)

       sum2 = dot_product(triNorm2, normA) + dot_product(triNorm2, normB) + &
            dot_product(triNorm2, normBp)

       if (aValid .and. bValid) then 
          ! Only use this to help pick one if both are still valid:

          if (sum1 < cutoff .and. sum2 > cutoff) then 
             aValid = .False.

          else if(sum2 < cutoff .and. sum1 > cutoff) then 
             bValid = .False.

          else if (sum1 < cutoff .and. sum2 < cutoff) then 
             ! Both bad. Take the least bad one
             if (sum1 > sum2) then 
                bValid = .False. 
             else
                aValid = .False. 
             end if
          end if
       end if

       ! ---------------------------------------------------------------
       ! Check 6: Front angle test: Try to keep the front as close as
       ! possible to the gap edges. 
       ! ---------------------------------------------------------------

       ! Triangle ABA+. Original implemnetation
       !am sum1 = vecAngle(ptA-ptAp, ptB-ptAp) + vecAngle(ptBp-ptB, ptAp-ptB) 
       !am sum2 = vecAngle(ptB-ptBp, ptA-ptBp) + vecAngle(ptAp-ptB, ptBp-ptB)
       sum1 = abs(vecAngle(ptA-ptAp, ptB-ptAp)) + abs(vecAngle(ptBp-ptB, ptAp-ptB)) 
       sum2 = abs(vecAngle(ptA-ptBp, ptB-ptBp)) + abs(vecAngle(ptBp-ptA, ptAp-ptA))

       if (sum1 > sum2) then 
          aPreferred = .True. 
       else
          aPreferred = .False. 
       end if

       ! ---------------------------------------------------------------
       ! Check 7: End of string test
       ! ---------------------------------------------------------------

       if (A == Ap) then 
          aValid = .False. 
          bValid = .True. 
       end if

       if (B == Bp) then 
          bValid = .False. 
          aValid = .True. 
       end if

       ! ---------------------------------------------------------------
       ! Decide on the triangle we want to take. 
       ! ---------------------------------------------------------------

       if (aValid .and. .not. bValid) then 

          ! We have no choice but to take A+

          call addTri(A, str1, B, str2, Ap, str1)
          advanceA = .True. 
          advanceB = .False. 

       else if (bValid .and. .not. aValid) then 

          ! We have no choice but to take B+

          call addTri(A, str1, B, str2, Bp, str2)
          advanceA = .False. 
          advanceB = .True. 

       else if (aValid .and. bValid) then 

          ! We could take either. Use the preferred triangle. 
          if (aPreferred) then 

             call addTri(A, str1, B, str2, Ap, str1)
             advanceA = .True.
             advanceB = .False.

          else

             call addTri(A, str1, B, str2, Bp, str2)
             advanceA = .False. 
             advanceB = .True. 

          end if

       else 

          ! Ewww. neither triangle is valid. Do not add any triangle,
          ! leave it for pocket zipping. Just move forward to Ap and Bp.
          print *,' ****** eww, skipping both A and B:', A, B, aValid, bValid, aPreferred, ptA, ptB

          advanceA = .False.
          advanceB = .False.

       end if

       ! Now we have to shuffle along the string. 
       if (advanceA .and. .not.advanceB) then 

          stepsA = stepsA + 1

          ! Copy the Ap to A
          A = Ap
          ptA = ptAp
          normA = normAp
          perpA = perpAp

          ! And get the new data for Ap
          Ap = nextNode(str1, A, .True.)
          ptAp = str1%x(:, Ap)
          normAp = str1%norm(:, Ap)
          perpAp = str1%perpNorm(:, Ap)

       else if (advanceB .and. .not.advanceA) then

          stepsB = stepsB + 1

          ! Copy the Bp to B
          B = Bp
          ptB = ptBp
          normB = normBp
          perpB = perpBp

          ! And get the new data for Bp
          Bp = nextNode(str2, B, .False.)
          ptBp = str2%x(:, Bp)
          normBp = str2%norm(:, Bp)
          perpBp = str2%perpNorm(:, Bp)

       else if (.not.advanceA .and. .not.advanceB) then

          ! Move A
          ! -------------------
          stepsA = stepsA + 1

          ! Copy the Ap to A
          A = Ap
          ptA = ptAp
          normA = normAp
          perpA = perpAp

          ! And get the new data for Ap
          Ap = nextNode(str1, A, .True.)
          ptAp = str1%x(:, Ap)
          normAp = str1%norm(:, Ap)
          perpAp = str1%perpNorm(:, Ap)
          ! -------------------

          ! Move B
          ! -------------------
          stepsB = stepsB + 1

          ! Copy the Bp to B
          B = Bp
          ptB = ptBp
          normB = normBp
          perpB = perpBp

          ! And get the new data for Bp
          Bp = nextNode(str2, B, .False.)
          ptBp = str2%x(:, Bp)
          normBp = str2%norm(:, Bp)
          perpBp = str2%perpNorm(:, Bp)
          ! -------------------
       else
          Stop ' *** Error: Can not advance both A and B ***'
       end if

       ! Finally increment the number of triangles we've used so far. 
       ii = ii + 1

       ! Account for two skipped triangles (i.e. one extra count)
       ! if both A and B are skipped and advanced to Ap and Bp.
       if (.not.advanceA .and. .not.advanceB) ii = ii + 1
    end do

  contains

    function nextNode(str, i, pos)

      implicit none
      type(oversetString), intent(iN) :: str
      integer(kind=intType), intent(in) :: i
      logical, intent(in) :: pos
      integer(kind=intType) :: nextNode

      if (pos) then 
         if (stepsA == nStepsA) then 
            nextNode = i
         else
            nextNode = i + 1
            if (nextNode > str%nNodes) then 
               if (str%isPeriodic) then 
                  ! Loop back around
                  nextNode = 1
               else
                  ! Leave it at the same node
                  nextNode = i
               end if
            end if
         end if
      else
         if (stepsB == nStepsB) then 
            nextNode = i
         else
            nextNode = i - 1
            if (nextNode < 1) then 
               if (str%isPeriodic) then 
                  ! Loop back around
                  nextNode = str%nNodes
               else
                  ! Leave it at the same node
                  nextNode = i
               end if
            end if
         end if
      end if
    end function nextNode

    function vecAngle(vec1, vec2)

      implicit none

      ! Input/Output
      real(kind=realType), dimension(3), intent(in) :: vec1, vec2
      real(kind=realType) :: vecAngle

      ! Working
      real(kind=realType), dimension(3) :: vecA, vecB

      vecA = vec1 / norm2(vec1)
      vecB = vec2 / norm2(vec2)

      vecAngle = acos(dot_product(vecA, vecB))

    end function vecAngle

    function elemBetweenNodes(str, a, b)
      implicit none

      ! Input/Output
      type(oversetString), intent(in) :: str
      integer(kind=intType), intent(in) :: a,b
      integer(kind=intType) :: elemBetweenNodes

      ! Working
      integer(kind=intType) :: e1, e2, e3, e4

      if (str%nte(1, a) == 1) then 
         e1 = str%nte(2, a)
         e2 = e1
      else
         e1 = str%nte(2, a)
         e2 = str%nte(3, a)
      end if

      if (str%nte(1, b) == 1) then 
         e3 = str%nte(2, b)
         e4 = e3
      else
         e3 = str%nte(2, b)
         e4 = str%nte(3, b)
      end if

      ! Two of the edges are the same. And this is the one that must
      ! be between the two nodes. 

      if (e1 == e3 .or. e1 == e4) then 
         elemBetweenNodes = e1
      else
         elemBetweenNodes = e2
      end if

    end function elemBetweenNodes

    function triArea(pt1, pt2, pt3)

      implicit none

      ! Input/Output
      real(kind=realType), intent(in), dimension(3) :: pt1, pt2, pt3
      real(kind=realType) :: triArea

      ! Working
      real(kind=realType), dimension(3) :: norm

      call cross_prod(pt2-pt1, pt3-pt1, norm)
      triArea = half * norm2(norm)

    end function triArea

    subroutine addTri(A, sA, B, sB, C, sC)

      ! Form a triangle from index 'A' on string 'sA' , index 'B' on
      ! string 'sB' and index 'C' on string 'sC'

      implicit none

      ! Input/Output
      integer(kind=intType), intent(in) :: A, B, C
      type(oversetString), intent(in) :: sA, sB, sC

      ! Working 
      type(oversetString), pointer :: p
      integer(kind=intType) :: mn1, mn2, mn3
      p => sA%p

      p%nTris = p%nTris+ 1

      ! mn = master node
      mn1 = sA%pNodes(A)
      mn2 = sB%pNodes(B)
      mn3 = sC%pNodes(C)

      p%tris(:, p%nTris) = (/mn1, mn2, mn3/)

      ! Add these three edges to master list of edges

      ! Edge 1:
      p%nEdges = p%nEdges + 1
      p%edges(p%nEdges)%n1 = mn1
      p%edges(p%nEdges)%n2 = mn2

      ! Edge 2:
      p%nEdges = p%nEdges + 1
      p%edges(p%nEdges)%n1 = mn2
      p%edges(p%nEdges)%n2 = mn3

      ! Edge 3:
      p%nEdges = p%nEdges + 1
      p%edges(p%nEdges)%n1 = mn3
      p%edges(p%nEdges)%n2 = mn1

      ! Flag the edge that got used on master
      if (sA%myID == sC%myID) then 
         p%elemUsed(sA%pElems(elemBetweenNodes(sA, A, C))) = 1
      else if (sB%myID == sC%myID) then 
         p%elemUsed(sB%pElems(elemBetweenNodes(sB, B, C))) = 1
      end if

    end subroutine addTri

  end subroutine crossZip

  subroutine makeCrossZip(p, strings, nStrings)
    use overset
    implicit none

    ! Input/output
    integer(kind=intType), intent(in) :: nStrings
    type(oversetString), intent(inout) :: p, strings(nStrings)

    ! Local variables
    integer(kind=intType) :: i, j, inode, jnode, iDir, inzip, jnzip, jp1, jm1
    integer(kind=intType) :: inodeS, inodeE, jnodeS, jnodeE, jsym, isym, jbeg
    integer(kind=intType) :: oID, oIDx, nsplits, oidJ, oidxJ, iSpl, iSubstr
    integer(kind=intType) :: inodetmp, inodep, inodem, inodeStmp, inodeEtmp
    integer(kind=intType) :: jnodep
    integer(kind=intType), allocatable, dimension(:) :: Ins, Ine, Jns, Jne
    logical :: jEndfound(2), checkNodeUsed
  !   ! ------------------------------------------------

  !   ! Should the node used be checked to avoid string pairings?
  !   checkNodeUsed = .True.

  !   ! Allocate arrays to keep track of nodes to avoid same substring pairs for
  !   ! crossZipping
  !   do i=1, nstrings
  !      allocate(strings(i)%XzipNodeUsed(1:strings(i)%nNodes))
  !      strings(i)%XzipNodeUsed = 0 ! Not used, 1 = used
  !   end do

  !   ! For debugging sub-strings pairings
  !   ! ----------------------------------
  !   open(unit=101, file="subGapStrings.dat", form='formatted')
  !   write(101,*) 'TITLE = "SubGap Strings Data" '
  !   write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
  !        "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
  !   ! ----------------------------------

  !   iSubstr = 0
  !   do i=1, nStrings

  !      ! Find total splits
  !      nSplits = 0
  !      inode = 1
  !      do while (inode <= strings(i)%nNodes-1) 
  !         if ( (strings(i)%otherID(1, inode) /= &
  !              strings(i)%otherID(1, inode+1)) .and. &
  !              inode > 1 ) then
  !            nsplits = nsplits + 1
  !         end if
  !         inode = inode + 1
  !      end do

  !      ifNsplits: if (nsplits == 0) then
  !         iSpl = 0
  !         ! No splits found, cross zip with the full other gap string.
  !         inodeS = 1
  !         inodeE = strings(i)%nNodes

  !         J = strings(i)%otherID(1, inodeS)

  !         jnodeS = 1
  !         jnodeE = strings(j)%nNodes

  !         ! If I is periodic J should be periodic too. Fix end nodes of J
  !         ! strings wrt otherID of inodeS/inodeE
  !         if (strings(i)%isPeriodic) then

  !            if (.not.strings(j)%isPeriodic) stop ' This is wrong '

  !            ! Find the first symmetric pair with J strings
  !            jnode = 1
  !            loopPerJsym: do while (jnode <= strings(j)%nNodes) 

  !               oidJ = strings(j)%otherID(1, jnode) ! string I
  !               oidxJ = strings(j)%otherID(2, jnode) ! index on string I

  !               if (oidJ /= strings(i)%myID) stop ' Should not be here, string I'

  !               if (strings(i)%otherID(1, oidxJ) /= strings(j)%myID) &
  !                    stop ' Something is wrong '

  !               if (strings(i)%otherID(2, oidxJ) == jnode) then
  !                  ! Found symmetric point
  !                  inodeE = oidxJ
  !                  jnodeS = jnode
  !                  exit loopPerJsym
  !               end if
  !               jnodep = jnode + 1
  !               if (jnode == strings(j)%nNodes) then
  !                  exit loopPerJsym
  !               end if

  !               jnode = jnodep
  !            end do loopPerJsym
  !            if (jnode == strings(j)%nNodes) stop ' One cycle over '

  !            inodeS = inodeE + 1
  !            if (inodeE == strings(i)%nNodes) inodeS = 1
  !            jnodeE = jnodeS - 1 
  !            if (jnodeS == 1) jnodeE = strings(j)%nNodes
  !         end if

  !         ! Check if the nodes have been used already
  !         if (checkNodeUsed) then
  !            inodeStmp = inodeS
  !            inodeEtmp = inodeE
  !            if (inodeStmp <= inodeEtmp) then

  !               inodeS = inodeStmp
  !               inodeE = inodeStmp
  !               chkinodeE10:do inode=inodeStmp, inodeEtmp
  !                  if (strings(i)%XzipNodeUsed(inode) == 0) then
  !                     inodeE = inode
  !                  else
  !                     exit chkinodeE10
  !                  end if
  !               end do chkinodeE10
  !               if (inodeE == inodeS) exit ifNsplits
  !            else ! inodeStmp > inodeEtmp
  !               inodeS = inodeStmp
  !               inodeE = inodeStmp
  !               inode = inodeStmp
  !               chkinodeE20:do while (inode >= inodeStmp .or. &
  !                    inode <= inodeEtmp) 
  !                  if (strings(i)%XzipNodeUsed(inode) == 0) then
  !                     inodeE = inode

  !                     inodep = inode + 1
  !                     if (inode == strings(i)%nNodes) inodep = 1

  !                     if (inodep == inodeStmp) then
  !                        ! Come full cycle around the string
  !                        exit chkinodeE20
  !                     end if

  !                     inode = inodep
  !                  else
  !                     exit chkinodeE20
  !                  end if
  !               end do chkinodeE20
  !               if (inodeE == inodeS) exit ifNsplits
  !            end if
  !         end if !checkNodeUsed

  !         if (jnodeS == jnodeE) exit ifNsplits

  !         print '(A,6(I4,x),A,2(I4,x))',&
  !              '      I, J, inodeS, inodeE, jnodeS, jnodeE ', &
  !              I, J, inodeS, inodeE, jnodeS, jnodeE, &
  !              'SubStr: ',iSubstr+1, iSubStr+2

  !         ! ! Debug sub-strings pairings
  !         ! iSubStr = iSubStr + 1
  !         ! call writeOversetSubString(strings(i), inodeS, inodeE, iSubStr, &
  !         !      strings, nStrings, 101)
  !         ! iSubStr = iSubStr + 1
  !         ! call writeOversetSubString(strings(j), jnodeS, jnodeE, iSubStr, &
  !         !      strings, nStrings, 101)

  !         ! Remember I string nodes used
  !         if (inodeS <= inodeE) then
  !            strings(i)%XzipNodeUsed(inodeS:inodeE) = 1
  !         else ! inodeS > inodeE
  !            strings(i)%XzipNodeUsed(inodeS:strings(i)%nNodes) = 1
  !            strings(i)%XzipNodeUsed(1:inodeE) = 1
  !         end if

  !         ! Remember J string nodes used
  !         if (jnodeS <= jnodeE) then
  !            strings(j)%XzipNodeUsed(jnodeS:jnodeE) = 1
  !         else ! jnodeS > jnodeE
  !            strings(j)%XzipNodeUsed(jnodeS:strings(j)%nNodes) = 1
  !            strings(j)%XzipNodeUsed(1:jnodeE) = 1
  !         end if

  !         if (strings(i)%isPeriodic .and. strings(j)%isPeriodic) then
  !            ! For both periodic I and J strings, jnodeS <--> inodeE 
  !            ! are symmetric pairs.
  !            call crossZip(strings(i), inodeE, inodeE, strings(j), &
  !                 jnodeS, jnodeS)
  !         else 
  !            ! single sided strings, end points on symmetric planes
  !            call crossZip(strings(i), inodeS, inodeE, strings(j), &
  !                 jnodeE, jnodeS)
  !         end if

  !      else ifNsplits

  !         ! Define the node ranges for 
  !         inodeS = 1
  !         inode = 1
  !         iSpl = 0

  !         loopInode: do while (inode <=strings(i)%nNodes-1)

  !            ! Cycle if node has already been used in previous strings pairings
  !            if (checkNodeUsed .and. strings(i)%XzipNodeUsed(inode)==1) then
  !               inode = inode + 1

  !               ! Update inodeS 
  !               inodeS = inode

  !               cycle loopInode 
  !            end if

  !            splitIf: if ( (strings(i)%otherID(1, inode) /= &
  !                 strings(i)%otherID(1, inode+1)) .and. &
  !                 inode > inodeS ) then

  !               inodeE = inode
  !               ispl = ispl + 1

  !               ! Potential other string
  !               J = strings(i)%otherID(1, inodeE)

  !               ! If I strings is periodic, change inodeS and inodeE
  !               ! --------------------------------------------------
  !               if (strings(i)%isPeriodic) then

  !                  ! Change inodeS:
  !                  inodetmp = inodeS
  !                  perInodeS:do 
  !                     inodem = inodetmp - 1
  !                     if (inodetmp == 1) inodem = strings(i)%nNodes
  !                     if (strings(i)%otherID(1, inodetmp) == J .and. &
  !                          strings(i)%otherID(1, inodem) /= J) then

  !                        inodeS = inodetmp
  !                        exit perInodeS
  !                     end if
  !                     inodetmp = inodem
  !                     if (inodetmp == inodeS) then
  !                        ! Retain inodeS if came full circle
  !                        print*, ' Come full circle inodeS'
  !                        exit perInodeS
  !                     end if
  !                  end do perInodeS

  !                  ! Change inodeE:
  !                  inodetmp = inodeE
  !                  perInodeE:do 
  !                     inodep = inodetmp + 1
  !                     if (inodetmp == strings(i)%nNodes) inodep = 1
  !                     if (strings(i)%otherID(1, inodetmp) == J .and. &
  !                          strings(i)%otherID(1, inodep) /= J) then

  !                        inodeE = inodetmp
  !                        exit perInodeE
  !                     end if
  !                     inodetmp = inodep
  !                     if (inodetmp == inodeE) then
  !                        ! Retain inodeE if came full circle
  !                        print*, ' Come full circle inodeE '
  !                        exit perInodeE
  !                     end if
  !                  end do perInodeE

  !               end if
  !               ! --- end periodic strings inodeS/inodeE -----------

  !               ! Find symmetric point
  !               jnode = 1
  !               loopjsym: do while (jnode <= strings(j)%nNodes)

  !                  oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !                  oidxJ = strings(j)%otherID(2, jnode) ! hopeful strings(i) Idx

  !                  !am if (oidJ == strings(i)%myID .and.  &
  !                  !am     oidxJ >= inodeS .and. oidxJ <= inodeE) then
  !                  if (oidJ == strings(i)%myID) then
  !                     if ( (inodeS < inodeE .and. oidxJ >= inodeS .and. &
  !                          oidxJ <= inodeE) .or. &
  !                          (inodeS > inodeE .and. (oidxJ >= inodeS .or. &
  !                          oidxJ <= inodeE)) ) then

  !                        if (strings(i)%otherID(1, oidxJ) == strings(j)%myID .and. &
  !                             strings(i)%otherID(2, oidxJ) == jnode) then
  !                           ! Found symmetry point
  !                           jsym = jnode
  !                           isym = oidxJ
  !                           exit loopjsym
  !                        end if
  !                     end if
  !                  end if
  !                  jnode = jnode + 1 
  !               end do loopjsym

  !               ! Traverse in +ve jDir to find jnodeE
  !               jnode = jsym 
  !               loopJnodeE: do 

  !                  ! Already at the end then exit
  !                  if ( jnode == strings(j)%nNodes .and. &
  !                       .not.strings(j)%isPeriodic) then
  !                     exit loopJnodeE
  !                  end if

  !                  jp1 = jnode + 1
  !                  ! Treat different for periodic strings 
  !                  if (jnode == strings(j)%nNodes .and. strings(j)%isPeriodic) then
  !                     jp1 = 1
  !                  end if

  !                  if ( strings(j)%otherID(1, jp1) /= strings(i)%myID) then
  !                     oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !                     oidxJ = strings(j)%otherID(2, jnode) ! strings(i) index

  !                     !am if(strings(i)%isPeriodic) inodeS = oidxJ
  !                     jnodeE = jnode
  !                     exit loopJnodeE

  !                  else if (strings(j)%otherID(1, jp1) == strings(i)%myid .and.  &
  !                       jp1 == strings(j)%nNodes .and. &
  !                       .not.strings(j)%isPeriodic) then
  !                     jnodeE = jnode+1
  !                     exit loopJnodeE
  !                  end if

  !                  jnode = jp1
  !                  ! Exit if counter has cycled through all nodes in
  !                  ! periodic string J
  !                  if (jnode == jsym) exit loopJnodeE 
  !               end do loopJnodeE

  !               ! Traverse in -ve jDir to find jnodeE
  !               jnode = jsym 
  !               loopJnodeS: do 

  !                  ! Already at the end then exit
  !                  if ( jnode == 1 .and. &
  !                       .not.strings(j)%isPeriodic) then
  !                     jnodeS = jnode
  !                     exit loopJnodeS
  !                  end if

  !                  jm1 = jnode - 1
  !                  ! Treat different for periodic strings 
  !                  if (jnode == 1 .and. strings(j)%isPeriodic) then
  !                     jm1 = strings(j)%nNodes
  !                  end if

  !                  if ( strings(j)%otherID(1, jm1) /= strings(i)%myID) then
  !                     oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !                     oidxJ = strings(j)%otherID(2, jnode) ! strings(i) index

  !                     !am if(strings(i)%isPeriodic) inodeE = oidxJ
  !                     jnodeS = jnode
  !                     exit loopJnodeS

  !                  else if (strings(j)%otherID(1, jm1) == strings(i)%myid .and.  &
  !                       jm1 == 1 .and. &
  !                       .not.strings(j)%isPeriodic) then
  !                     jnodeS = jm1
  !                     exit loopJnodeS
  !                  end if

  !                  jnode = jm1
  !                  ! Exit if counter has cycled through all nodes in
  !                  ! periodic string J
  !                  if (jnode == jsym) exit loopJnodeS 
  !               end do loopJnodeS

  !               if (jnodeS == jnodeE) exit ifNsplits

  !               print '(A,6(I4,x),A,2(I4,x))',&
  !                    '      I, J, inodeS, inodeE, jnodeS, jnodeE ', &
  !                    I, J, inodeS, inodeE, jnodeS, jnodeE, &
  !                    'SubStr: ',iSubstr+1, iSubStr+2

  !               ! ! Debug sub-strings pairings
  !               ! iSubStr = iSubStr + 1
  !               ! call writeOversetSubString(strings(i), inodeS, inodeE, iSubstr,&
  !               !      strings, nStrings, 101)
  !               ! iSubStr = iSubStr + 1
  !               ! call writeOversetSubString(strings(j), jnodeS, jnodeE, iSubstr,&
  !               !      strings, nStrings, 101)

  !               ! Remember I string nodes used
  !               if (inodeS <= inodeE) then
  !                  strings(i)%XzipNodeUsed(inodeS:inodeE) = 1
  !               else ! inodeS > inodeE
  !                  strings(i)%XzipNodeUsed(inodeS:strings(i)%nNodes) = 1
  !                  strings(i)%XzipNodeUsed(1:inodeE) = 1
  !               end if

  !               ! Remember J string nodes used
  !               if (jnodeS <= jnodeE) then
  !                  strings(j)%XzipNodeUsed(jnodeS:jnodeE) = 1
  !               else ! jnodeS > jnodeE
  !                  strings(j)%XzipNodeUsed(jnodeS:strings(j)%nNodes) = 1
  !                  strings(j)%XzipNodeUsed(1:jnodeE) = 1
  !               end if

  !               ! Do crossZip
  !               call crossZip(strings(i), inodeS, inodeE, strings(j), &
  !                    jnodeE, jnodeS)

  !               inodeS = inodeE + 1
  !            end if splitIf

  !            inode = inode + 1
  !         end do loopInode


  !         iSpl = iSpl + 1

  !         ! ---------------------------------------------------------
  !         ! Now do the remaining nodes of string I
  !         ! ---------------------------------------------------------
  !         inodeE = strings(i)%nNodes

  !         ! Potential other string
  !         J = strings(i)%otherID(1, inodeE)

  !         ! If I strings is periodic, change inodeS and inodeE
  !         ! --------------------------------------------------
  !         if (strings(i)%isPeriodic) then

  !            ! Change inodeS:
  !            inodetmp = inodeS
  !            perInodeS1:do 
  !               inodem = inodetmp - 1
  !               if (inodetmp == 1) inodem = strings(i)%nNodes
  !               if (strings(i)%otherID(1, inodetmp) == J .and. &
  !                    strings(i)%otherID(1, inodem) /= J) then

  !                  inodeS = inodetmp
  !                  exit perInodeS1
  !               end if
  !               inodetmp = inodem
  !               if (inodetmp == inodeS) then
  !                  ! Retain inodeS if came full circle
  !                  print*, ' Come full circle inodeS'
  !                  exit perInodeS1
  !               end if
  !            end do perInodeS1

  !            ! Change inodeE:
  !            inodetmp = inodeE
  !            perInodeE1:do 
  !               inodep = inodetmp + 1
  !               if (inodetmp == strings(i)%nNodes) inodep = 1
  !               if (strings(i)%otherID(1, inodetmp) == J .and. &
  !                    strings(i)%otherID(1, inodep) /= J) then

  !                  inodeE = inodetmp
  !                  exit perInodeE1
  !               end if
  !               inodetmp = inodep
  !               if (inodetmp == inodeE) then
  !                  ! Retain inodeE if came full circle
  !                  print*, ' Come full circle inodeE '
  !                  exit perInodeE1
  !               end if
  !            end do perInodeE1

  !         end if
  !         ! --- end periodic strings inodeS/inodeE -----------


  !         ! Find symmetric point
  !         jnode = 1
  !         loopjsym1: do while (jnode <= strings(j)%nNodes)

  !            oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !            oidxJ = strings(j)%otherID(2, jnode) ! hopeful strings(i) Idx

  !            if (oidJ == strings(i)%myID) then
  !               if ( (inodeS < inodeE .and. oidxJ >= inodeS .and. &
  !                    oidxJ <= inodeE) .or. &
  !                    (inodeS > inodeE .and. (oidxJ >= inodeS .or. &
  !                    oidxJ <= inodeE)) ) then

  !                  if (strings(i)%otherID(1, oidxJ) == strings(j)%myID .and. &
  !                       strings(i)%otherID(2, oidxJ) == jnode) then
  !                     ! Found symmetry point
  !                     jsym = jnode
  !                     isym = oidxJ
  !                     exit loopjsym1
  !                  end if
  !               end if
  !            end if

  !            jnode = jnode + 1 
  !         end do loopjsym1

  !         ! Traverse in +ve jDir to find jnodeE
  !         jnode = jsym 
  !         loopJnodeE1: do 

  !            ! Already at the end then exit
  !            if ( jnode == strings(j)%nNodes .and. &
  !                 .not.strings(j)%isPeriodic) then
  !               jnodeE = jnode
  !               exit loopJnodeE1
  !            end if

  !            jp1 = jnode + 1
  !            ! Treat different for periodic strings 
  !            if (jnode == strings(j)%nNodes .and. strings(j)%isPeriodic) then
  !               jp1 = 1
  !            end if

  !            if ( strings(j)%otherID(1, jp1) /= strings(i)%myID) then
  !               oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !               oidxJ = strings(j)%otherID(2, jnode) ! strings(i) index

  !               !am if(strings(i)%isPeriodic) inodeS = oidxJ
  !               jnodeE = jnode
  !               exit loopJnodeE1

  !            else if (strings(j)%otherID(1, jp1) == strings(i)%myid .and.  &
  !                 jp1 == strings(j)%nNodes .and. &
  !                 .not.strings(j)%isPeriodic) then
  !               jnodeE = jnode+1
  !               exit loopJnodeE1
  !            end if

  !            jnode = jp1
  !            ! Exit if counter has cycled through all nodes in
  !            ! periodic string J
  !            if (jnode == jsym) exit loopJnodeE1
  !         end do loopJnodeE1

  !         ! Traverse in -ve jDir to find jnodeE
  !         jnode = jsym 
  !         loopJnodeS1: do 

  !            ! Already at the end then exit
  !            if ( jnode == 1 .and. &
  !                 .not.strings(j)%isPeriodic) then
  !               jnodeS = jnode
  !               exit loopJnodeS1
  !            end if

  !            jm1 = jnode - 1
  !            ! Treat different for periodic strings 
  !            if (jnode == 1 .and. strings(j)%isPeriodic) then
  !               jm1 = strings(j)%nNodes
  !            end if

  !            if ( strings(j)%otherID(1, jm1) /= strings(i)%myID) then
  !               oidJ = strings(j)%otherID(1, jnode) ! hopeful strings(i) ID
  !               oidxJ = strings(j)%otherID(2, jnode) ! strings(i) index

  !               !am if(strings(i)%isPeriodic) inodeE = oidxJ
  !               jnodeS = jnode
  !               exit loopJnodeS1

  !            else if (strings(j)%otherID(1, jm1) == strings(i)%myid .and.  &
  !                 jm1 == 1 .and. &
  !                 .not.strings(j)%isPeriodic) then
  !               jnodeS = jm1
  !               exit loopJnodeS1
  !            end if

  !            jnode = jm1
  !            ! Exit if counter has cycled through all nodes in
  !            ! periodic string J
  !            if (jnode == jsym) exit loopJnodeS1
  !         end do loopJnodeS1



  !         ! Check it has nodes that have not been used already
  !         if (checkNodeUsed) then
  !            inodeStmp = inodeS
  !            inodeEtmp = inodeE
  !            if (inodeStmp <= inodeEtmp) then

  !               inodeS = inodeStmp
  !               inodeE = inodeStmp
  !               chkinodeE1:do inode=inodeStmp, inodeEtmp
  !                  if (strings(i)%XzipNodeUsed(inode) == 0) then
  !                     inodeE = inode
  !                  else
  !                     exit chkinodeE1
  !                  end if
  !               end do chkinodeE1
  !               if (inodeE == inodeS) exit ifNsplits
  !            else ! inodeStmp > inodeEtmp
  !               inodeS = inodeStmp
  !               inodeE = inodeStmp
  !               inode = inodeStmp
  !               chkinodeE2:do while (inode >= inodeStmp .or. &
  !                    inode <= inodeEtmp) 
  !                  if (strings(i)%XzipNodeUsed(inode) == 0) then
  !                     inodeE = inode

  !                     inodep = inode + 1
  !                     if (inode == strings(i)%nNodes) inodep = 1

  !                     inode = inodep
  !                  else
  !                     exit chkinodeE2
  !                  end if
  !               end do chkinodeE2
  !               if (inodeE == inodeS) exit ifNsplits
  !            end if
  !         end if !checkNodeUsed

  !         if (jnodeS == jnodeE) exit ifNsplits

  !         print '(A,6(I4,x),A,2(I4,x))',&
  !              'last: I, J, inodeS, inodeE, jnodeS, jnodeE ', &
  !              I, J, inodeS, inodeE, jnodeS, jnodeE, &
  !              'SubStr: ',iSubstr+1, iSubStr+2

  !         ! ! Debug sub-strings pairings
  !         ! iSubStr = iSubStr + 1
  !         ! call writeOversetSubString(strings(i), inodeS, inodeE, iSubstr,&
  !         !      strings, nStrings, 101)
  !         ! iSubStr = iSubStr + 1
  !         ! call writeOversetSubString(strings(j), jnodeS, jnodeE, iSubstr,&
  !         !      strings, nStrings, 101)

  !         ! Remember I string nodes used
  !         if (inodeS <= inodeE) then
  !            strings(i)%XzipNodeUsed(inodeS:inodeE) = 1
  !         else ! inodeS > inodeE
  !            strings(i)%XzipNodeUsed(inodeS:strings(i)%nNodes) = 1
  !            strings(i)%XzipNodeUsed(1:inodeE) = 1
  !         end if

  !         ! Remember J string nodes used
  !         if (jnodeS <= jnodeE) then
  !            strings(j)%XzipNodeUsed(jnodeS:jnodeE) = 1
  !         else ! jnodeS > jnodeE
  !            strings(j)%XzipNodeUsed(jnodeS:strings(j)%nNodes) = 1
  !            strings(j)%XzipNodeUsed(1:jnodeE) = 1
  !         end if

  !         ! Do crossZip
  !         call crossZip(strings(i), inodeS, inodeE, strings(j), &
  !              jnodeE, jnodeS)

  !         ! ---------------------------------------------------------
  !         ! End remaining nodes treatment
  !         ! ---------------------------------------------------------

  !      end if ifNsplits

  !   end do ! nStrings

  !   ! Debug sub-strings pairings
  !   close(101)


  !   do i=1, nStrings
  !      deallocate(strings(i)%XzipNodeUsed)
  !   end do

end subroutine makeCrossZip

  subroutine makePocketZip(p, strings, nStrings, pocketMaster)
    use overset
    implicit none

    ! Input/output
    integer(kind=intType), intent(in) :: nStrings
    type(oversetString), intent(in) :: p, strings(nStrings)
    type(oversetString) :: pocketMaster
    ! Local variables
    integer(kind=intType) :: i, nsum1, nsum2, ndiff1, ndiff2, ipedge, icur
    integer(kind=intType) :: n1, n2, npolyEdges, npolyEdgestmp, nEdgeUsed
    integer(kind=intType) :: nNodes1, nNodes2, cn1, cn2, str1, str2, nends
    type(oversetEdge), pointer, dimension(:) :: polyEdges, polyEdgestmp
    !am type(pocketEdge), pointer, dimension(:) :: pocketEdges
    integer(kind=intType), allocatable, dimension(:) :: edgeMap, edgeMaptmp
    integer(kind=intType), allocatable, dimension(:) :: nodeList, nodeMap
    logical :: isEndEdge
    logical, allocatable, dimension(:) :: PocketEdgeUsed

    type(oversetString), pointer  :: pocketStrings, strPkt
    type(oversetString), allocatable, dimension(:), target :: tmpStrings
    type(oversetString), allocatable, dimension(:), target :: pocketStringsArr
    integer(kind=intType) :: npocketEdges, npocketStrings, nPktEdges, nUnique
    integer(kind=intType) :: ip, curElem, nElems, iStart, firstElem
    real(kind=realType) :: triArea
    ! ---------------------------------------------------------------


    ! ---------------------------------------------------------------
    ! PocketZip 1: 
    ! First sort the edges.
    ! ---------------------------------------------------------------
    call qsortEdgeType(p%Edges, p%nEdges)

    ! ---------------------------------------------------------------------
    ! PocketZip 2:
    ! Now cancel out edges counted twice in master%Edges. The edges are
    ! sorted such that the edges with same nodes are in consecutive order.
    ! ---------------------------------------------------------------------

    ! Over estimate of remaining pocket edges to zip
    allocate(polyEdges(p%nEdges)) 
    allocate(polyEdgestmp(p%nEdges)) 
    allocate(edgeMap(p%nEdges))
    allocate(edgeMaptmp(p%nEdges))
    edgeMap = -1
    edgeMaptmp = -1

    ! Initialize
    npolyEdgestmp = p%nEdges
    edgeMaptmp = -1
    do i=1, p%nEdges
       polyEdgesTmp(i)%n1 = p%Edges(i)%n1
       polyEdgesTmp(i)%n2 = p%Edges(i)%n2
       edgeMaptmp(i) = i
    end do

    ! Eliminate the edges going through the ordered edges. 
    ! The sorted opposite edges are canceled in pairs.
    npolyEdges = 0

    nends = 0

    loop_outer: do 

       npolyEdges = 0
       edgeMap = -1
       i = 1
       nends = 0
       loopEdge: do while (i < npolyEdgestmp) 

          nsum1 = (polyEdgesTmp(i)%n1   + polyEdgesTmp(i)%n2)
          nsum2 = (polyEdgesTmp(i+1)%n1 + polyEdgesTmp(i+1)%n2)

          ndiff1 = (polyEdgesTmp(i)%n2   - polyEdgesTmp(i)%n1)
          ndiff2 = (polyEdgesTmp(i+1)%n2 - polyEdgesTmp(i+1)%n1)

          ! If the edge joins the end nodes of single sided 
          ! fullStrings pair, eliminate it.
          ! -----------------------------------------------------------
          isEndEdge = .False.

          ! Parent nodes of this edge
          n1 = polyEdgesTmp(i)%n1
          n2 = polyEdgesTmp(i)%n2

          str1    = p%cNodes(1, n1) ! node1's child fullStrings ID
          cn1     = p%cNodes(2, n1) ! node1's child fullStrings node index
          !am nNodes1 = p%cNodes(3, n1) ! node1's child fullStrings nNodes size
          nNodes1 = strings(str1)%nNodes ! node1's child fullStrings nNodes size

          str2    = p%cNodes(1, n2) ! node2's child fullStrings ID
          cn2     = p%cNodes(2, n2) ! node2's child fullStrings node index
          !am nNodes2 = p%cNodes(3, n2) ! node2's child fullStrings nNodes size
          nNodes2 = strings(str2)%nNodes ! node2's child fullStrings nNodes size

          if (str1 /= str2 ) then
             if (.not.strings(str1)%isperiodic .and. &
                  .not.strings(str2)%isperiodic .and. &
                  (cn1==1 .or. cn1==nNodes1) .and. (cn2==1 .or. cn2==nNodes2)) then
                ! This is an end edge, eliminate this one too.
                nends = nends + 1

                i = i + 1
                cycle loopEdge
             end if
          end if
          ! --- End end nodes check -----------------------------------


          if (nsum1 == nsum2 .and. ndiff1 + ndiff2 == 0) then
             ! Found ordered edges pair. Eliminate these two edges.

             ! Jump to 2nd next edge
             i = i + 2
             cycle loopEdge
          else
             ! Add this edge to new edge list
             npolyEdges = npolyEdges + 1
             polyEdges(npolyEdges)%n1 = polyEdgesTmp(i)%n1
             polyEdges(npolyEdges)%n2 = polyEdgesTmp(i)%n2

             edgeMap(npolyEdges) = edgeMaptmp(i)

             if (i+1 == npolyEdgestmp) then
                ! Add this last edge to new edge list
                npolyEdges = npolyEdges + 1
                polyEdges(npolyEdges)%n1 = polyEdgesTmp(i+1)%n1
                polyEdges(npolyEdges)%n2 = polyEdgesTmp(i+1)%n2

                edgeMap(npolyEdges) = edgeMaptmp(i+1)

                exit loopEdge
             end if

             ! Loop to next edge
             i = i + 1
          end if
       end do loopEdge
       print*, '===========================> nends ',nends

       if (npolyEdgestmp - npolyEdges == 0) then
          print*, ' No more edges to cancel ', npolyEdges
          exit loop_outer
       end if

       ! Update polyedgesTmp  and edgeMaptmp before cycling loop
       ! -------------------------------------------------------
       ! First zero out
       do i=1, npolyEdgesTmp
          polyEdgesTmp(i)%n1 = -1
          polyEdgesTmp(i)%n2 = -1
       end do
       edgeMaptmp = -1

       ! New tmp sizes
       npolyEdgestmp = npolyEdges
       do i=1, npolyEdges
          polyEdgesTmp(i)%n1 = polyEdges(i)%n1
          polyEdgesTmp(i)%n2 = polyEdges(i)%n2
          edgeMaptmp(i) = edgeMap(i)
       end do
       ! -------------------------------------------------------

    end do loop_outer
    ! ---------------------------------------------------------------------
    ! End canceling out edges.
    ! ---------------------------------------------------------------------
    print*, 'npolyEdges ', npolyEdges
    do i=1, npolyEdges
       write(4001,*)i, edgeMap(i), polyEdges(i)%n1, polyEdges(i)%n2, &
            (polyedges(i)%n1+polyedges(i)%n2)*1.e5 + (polyedges(i)%n2-polyedges(i)%n1)
    end do

    ! Debug polygonEdges
    ! -----------------------------------------------------------
    open(unit=101, file="polygonEdges.dat", form='formatted')
    write(101,*) 'TITLE = "PolygonEdges Data" '
    write(101,*) 'Variables = "X", "Y", "Z"'
    write(101,*) "Zone T=Pockets"
    write (101,*) "Nodes = ", npolyEdges*2, " Elements= ", npolyEdges, " ZONETYPE=FELINESEG"
    write (101,*) "DATAPACKING=POINT"

    ! node data
    do i=1, npolyEdges
       n1 = polyEdges(i)%n1
       n2 = polyEdges(i)%n2
       ! node 1
       write(101,'(3(E20.12,x))')p%x(1, n1), p%x(2, n1), p%x(3, n1)
       ! node 2
       write(101,'(3(E20.12,x))')p%x(1, n2), p%x(2, n2), p%x(3, n2)
    end do

    ! Edge data
    do i=1, npolyEdges
       write(101,'(3(I5,x))')2*i-1, 2*i
    end do
    close(101)
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    ! PocketZip 3:
    ! Accumulate full pocket string edges.
    ! Do similar to how master edges were accumulated.
    ! 3.1: create pocketMaster elems and nodes
    ! 3.2: perform doChain on pocketMaster elems and create pocket strings.
    ! 3.3: selfZip pocketStrings
    !-------------------------------------------------------------

    ! ----------------------------------------
    ! 3.1: create pocketMaster elems and nodes
    ! ----------------------------------------

    ! First create unique nodeList
    ! Number of nodes is twice the size of edges
    allocate(nodeList(2*nPolyEdges), nodeMap(2*nPolyEdges))
    do i=1, nPolyEdges
       nodeList(2*i-1) = PolyEdges(i)%n1
       nodeList(2*i)   = PolyEdges(i)%n2
    end do

    call unique(nodeList, 2*nPolyEdges, nUnique, nodeMap)

    ! Define pocketMaster string
    call nullifyString(pocketMaster)

    pocketMaster%nNodes = nUnique
    pocketMaster%nElems = nPolyEdges

    allocate(pocketMaster%nodeData(10, pocketMaster%nNodes), &
         pocketMaster%conn(2, pocketMaster%nElems), &
         pocketMaster%intNodeData(3, pocketMaster%nNodes))
    allocate(pocketMaster%pNodes(pocketMaster%nNodes))

    allocate(pocketMaster%otherID(2, pocketMaster%nNodes))
    pocketMaster%otherID = -1

    ! Set the string pointers to the individual arrays
    call setStringPointers(pocketMaster)

    ! Node data
    do i=1, nUnique
       ip = nodeList(i)

       !Copy x, norm, perpNorm, and h data from global master data
       pocketMaster%nodeData(:, i) = p%nodeData(:, ip)  
       pocketMaster%intNodeData(:, i) = p%intNodeData(:, ip)

       ! Save the original master node index
       pocketMaster%pNodes(i) = ip
    end do

    ! Element data
    do i=1, nPolyEdges
       pocketMaster%conn(1, i) = nodeMap(2*i-1) !<-- map to the unique node index
       pocketMaster%conn(2, i) = nodeMap(2*i) !<-- map to the unique node index
    end do

    ! Debug pocketMaster
    ! --------------------
    pocketMaster%myID = 88
    open(unit=101, file="pocketMaster.dat", form='formatted')
    write(101,*) 'TITLE = "PocketMaster Data" '
    write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
         "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
    allocate(tmpStrings(1))
    tmpStrings(1) = pocketMaster ! Derived type assignment
    call writeOversetString(tmpStrings(1), tmpStrings, 1, 101)
    close(101)
    deallocate(tmpStrings)
    ! --------------------

    ! Create nte info
    call createNodeToElem(pocketMaster)

    ! End 3.1 create pocketMaster elems and nodes

    ! ---------------------------------------------------------------------
    ! 3.2: perform doChain on pocketMaster elems and create pocket strings.
    ! ---------------------------------------------------------------------

    ! Create ordered pocketStrings based on connectivity using
    ! linked list.

    ! Some additional arrays.
    allocate(pocketMaster%elemUsed(pocketMaster%nElems), &
         pocketMaster%subStr(2, pocketMaster%nElems), &
         pocketMaster%cNodes(2, pocketMaster%nNodes)) ! Third index saves the size of the substrings
    !am pocketMaster%cNodes(3, pocketMaster%nNodes)) ! Third index saves the size of the substrings

    ! Initialize
    pocketMaster%cNodes = 0
    pocketMaster%elemUsed = 0
    curElem = 1
    nPocketStrings = 0

    print*,' Create pocketStrings from pocketMaster >>>>>>>>> '
    ! Do similar to creation to fullStrings
    do while (curElem < pocketMaster%nElems)

       ! First node of the first unused element
       iStart = pocketMaster%conn(1, curElem)
       nElems = pocketMaster%nte(1, iStart)

       ! --------------------
       ! First side of chain:
       ! --------------------
       firstElem = pocketMaster%nte(2, iStart)
       pocketMaster%subStr(1, 1) = firstElem
       call doChain(pocketMaster, iStart, 1)

       ! ---------------------
       ! Second side of chain: 
       ! ---------------------
       ! Ideally, pocketStrings should be periodic. Should not require
       ! second side of chain. Do a sanity check anyway.

       if (nElems > 1) then
          firstElem = pocketMaster%nte(3, iStart)

          ! Check the second one is already end of the periodic chain
          ! done above.
          if (pocketMaster%elemUsed(firstElem) == 0) then
             print*, ' First side did not create a periodic chain'
             stop ' Error'
          end if
       end if

       ! Extract pocketStrings linked list from elements present in
       ! pocketMaster%subStr buffer

       ! Create or add a new string to 'pocketStrings' linked list
       if (nPocketStrings == 0) then
          allocate(pocketStrings) ! Create first linked list node
          nPocketStrings = 1
          pocketStrings%next => pocketStrings
          strPkt => pocketStrings ! Work with the first termporary linked list node
       else
          allocate(strPkt%next) ! Create new linked list node and link to it
          strPkt%next%next => pocketStrings ! point the end back to original pointer.
          strPkt => strPkt%next ! Work with the new temporary linked list node

          nPocketStrings = nPocketStrings + 1
       end if

       ! Create a pocketStrings substring from pocketMaster buffer
       call createSubStringFromElems(pocketMaster, strPkt, nPocketStrings)

       do while (pocketMaster%elemUsed(curElem) == 1 .and. &
            curElem < pocketMaster%nElems)
          curElem = curElem + 1
       end do
    end do ! main while loop

    ! Debug pocketStrings
    ! ---------------------------------------------
    print *, 'nPocketStrings:', nPocketStrings

    ! Temporary strings array for plotting and pocketZipping
    allocate(pocketStringsArr(nPocketStrings))
    strPkt => pocketStrings
    i = 0
    do while(i < nPocketStrings)
       i = i + 1
       pocketStringsArr(i) = strPkt ! Derived type assignment
       call nullifyString(strPkt)
       strPkt => strPkt%next
    end do

    open(unit=101, file="pocketStrings.dat", form='formatted')
    write(101,*) 'TITLE = "PocketStrings Data" '

    write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
         "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
    do i=1, nPocketStrings
       ! Temporarily allocate otherID
       allocate(pocketStringsArr(i)%otherID(2, pocketStringsArr(i)%nNodes))
       pocketStringsArr(i)%otherID = -1

       call writeOversetString(pocketStringsArr(i), pocketStringsArr, &
            nPocketStrings, 101)
    end do
    close(101)

    ! ---------------------------------------------
    ! ------------------------------------------------------------------------
    ! End 3.2, perform doChain on pocketMaster elems and create pocket strings.
    ! ------------------------------------------------------------------------
    print*, ' Perform pocketZipping ======>'

    ! --------------------------
    ! 3.3: selfZip pocketStrings
    ! --------------------------

    ! ======= Do actual pocketZip ================

    ! This is again going to be similar to selfZip of fullStrings

    pocketMaster%myID = 88

    ! Allocate space for pocket triangles.
    ! (n-sided polygon -> n-2 triangles)
    allocate(pocketMaster%tris(3, pocketMaster%nElems))
    pocketMaster%nTris = 0

    ! Build the pocketMaster tree
    pocketMaster%tree => kdtree2_create(pocketMaster%x, sort=.True.)

    ! Loop over pocketStrings and begin pocketZip starting 
    ! from smallest convex ear.
    do i=1, nPocketStrings
       pocketZiploop: do while (pocketStringsArr(i)%nNodes > 2) 
          ! Each pass zips one triangle. Keep zipping
          ! until last triangle is zipped in the pocket polygon.
          call pocketZip(pocketStringsArr(i))
       end do pocketZiploop
    end do

    call writeOversetTriangles(pocketMaster, "pocketTriangulation.dat")
    ! ------------------------------
    ! End 3.3, selfZip pocketStrings
    ! ------------------------------
    !am call computeTriSurfArea(pocketMaster, triArea)
    !am print*,'triArea ',triArea

  end subroutine makePocketZip

  !
  ! ======================================
  subroutine pocketZip(s)

    use overset
    use kdtree2_module
    implicit none

    ! Input parameters
    type(oversetString), intent(inout), target :: s

    ! Local variables
    integer(kind=intType) :: i, j, k, ii, im1, ip1, nalloc, idx, nFound, N
    integer(kind=intType) :: nNodes, nElems, elem1, elem2, imin
    logical :: lastNodeZipper, inTri, overlapFound
    real(kind=realType), dimension(3) :: v1, v2, norm, c
    real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm, costhetaMax

    integer(Kind=intType), dimension(:), allocatable :: nodeMap
    integer(Kind=intType), dimension(:), allocatable :: badNode
    type(kdtree2_result), dimension(:), allocatable  :: results
    !am real(kind=realType), dimension(:, :), pointer :: xTmp, normTmp
    real(kind=realType), dimension(:, :), pointer :: nodeDataTmp
    integer(kind=intType), dimension(:, :), pointer :: connTmp, intNodeDataTmp
    integer(kind=intType), dimension(:), pointer :: pNodesTmp
    logical :: foundZipNode

    ! ----------------------------------
    ! pocketStrings 's' should be all peridic
    if (.not.s%isPeriodic) then
       print*,' Non-periodic pocketStrings ID: ',s%myID
       stop
    end if

    N = s%nNodes

    nAlloc = 25
    allocate(results(nAlloc))
    allocate(nodeMap(s%nNodes))
    allocate(badNode(s%nNodes))
    nodeMap = 1
    badNode = 0 ! 1 if bad

    foundZipNode = .False.

    !=====================================
    outerZiploop: do while (.not.foundZipNode)
       !=====================================

       ! Find min angled ear
       costhetaMax = -Large 
       nodeloop: do ii=1, N

          if (badNode(ii) == 1) cycle nodeloop

          im1 = ii - 1
          ip1 = ii + 1

          ! Peroidic string at end...loop around
          if (ii == N) then 
             ip1 = 1
          else if (ii == 1) then
             im1 = s%nNodes
          end if

          ! Determine the angle between the vectors
          v1 = s%x(:, ip1) - s%x(:, ii)
          v2 = s%x(:, im1) - s%x(:, ii)
          v1nrm = norm2(v1)
          v2nrm = norm2(v2)
          call cross_prod(v2, v1, norm)
          norm = norm / norm2(norm)

          ! Interior node norm and potential node norm should point
          ! in the same direction.
          if (dot_product(norm, s%norm(:, ii)) > zero) then

             ! Dot product of im1 and ip1 nodes should be close
             if (dot_product(s%norm(:, ip1), s%norm(:, im1)) > 0.80) then

                costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)

                ! cos(theta) is the largest for smallest angle
                if (costhetaMax <= costheta) then
                   costhetaMax = costheta
                   imin = ii
                end if
             end if
          end if

       end do nodeloop !ii 

       ! Zip about imin node
       ii = imin

       im1 = ii - 1
       ip1 = ii + 1

       ! Peroidic string at end...loop around
       if (ii == N) then 
          ip1 = 1
       else if (ii == 1) then
          im1 = s%nNodes
       end if

       ! Check that this triangle does not contain any other 
       ! pocketStrings nodes

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
          ! This triangle is good!

          s%p%nTris = s%p%nTris+ 1
          s%p%tris(:, s%p%nTris) = (/s%pNodes(ip1), s%pNodes(ii),s%pNodes(im1)/)
          nodeMap(ii) = 0

          ! The two shorted string edges have been used for selfZip
          ! 
          elem1 = s%p%nte(2, s%pNodes(ii))
          elem2 = s%p%nte(3, s%pNodes(ii))
          s%p%elemUsed(elem1) = 1
          s%p%elemUsed(elem2) = 1

          foundZipNode = .True.
       else 
          ! Bad node. Need to cycle through rest of pocket nodes.
          ! Remember this bad node in next cycle.
          badNode(ii) = 1 
          cycle outerZiploop
       end if

       ! Modify the pocketStrings to remove the two elements and the node
       ! that got eliminated due to pocketZipping. 

       ! Save pointers to existing data
       nNodes = s%nNodes
       nElems = s%nElems
       nodeDataTmp => s%nodeData
       intNodeDataTmp => s%intNodeData
       connTmp => s%conn
       pNodesTmp => s%pNodes

       ! nodeMap(imin) = 0, and 1 for the rest nodes. Create new nodeMap
       ! by taking off node 'imin'.
       j = 0
       do i=1, s%nNodes
          if (nodeMap(i) == 1) then 
             j = j + 1
             nodeMap(i) = j
          end if
       end do

       ! Update the number of nodes/elems in our shorted chain. Every
       ! zipper reduces the number of nodes and number of elems by 1
       s%nNodes = s%nNodes - 1
       s%nElems = s%nElems - 1

       allocate(s%nodeData(10, s%nNodes), s%intNodeData(3, s%nNodes), &
            s%pNodes(s%nNodes), s%conn(2, s%nElems))

       ! Set the pointers for the new string
       call setStringPointers(s)

       do i=1, nNodes
          if (nodeMap(i) /= 0) then 
             s%nodeData(:, nodeMap(i)) = nodeDataTmp(:, i)
             s%intNodeData(:, nodeMap(i)) = intNodeDataTmp(:, i)
             s%pNodes(nodeMap(i)) = pNodesTmp(i)

             ! Update string's parent's child node data
             !am s%p%cNodes(:, s%pNodes(nodeMap(i))) = (/s%myID, nodeMap(i), s%nNodes/)
             s%p%cNodes(:, s%pNodes(nodeMap(i))) = (/s%myID, nodeMap(i) /)
          end if
       end do

       ! Since we know the string was in order, we can simply redo the connectivity
       do i=1, s%nElems
          s%conn(:, i) = (/i, i+1/)
       end do
       if (s%isPeriodic) then 
          s%conn(2, s%nElems) = 1
       end if

       ! Dellocate the existing memory
       deallocate(nodeDataTmp, intNodeDataTmp, connTmp, pNodesTmp)

    end do outerZiploop

    deallocate(nodeMap, badNode)

  end subroutine pocketZip

  subroutine computeTriSurfArea(master, area)

    ! Computes area sum of all triangles belonging to object master
    use overset
    implicit none

    ! Input parameters
    type(oversetString), intent(in) :: master
    real(kind=realType), intent(out) :: area

    ! Local variables
    integer(kind=intType) :: i, n1, n2, n3
    real(kind=realType), dimension(3) :: v1, v2, norm

    area = 0.0
    do i=1, master%nTris
       n1 = master%tris(1, i)
       n2 = master%tris(2, i)
       n3 = master%tris(3, i)

       v1 = master%x(:, n2) - master%x(:, n1)
       v2 = master%x(:, n3) - master%x(:, n1)
       call cross_prod(v1, v2, norm)
       area = area + half*norm2(norm)
    end do

  end subroutine computeTriSurfArea

  function triOverlap(pt1, pt2, pt3, str, i1, i2)

    implicit none
    ! Input/Output
    real(kind=realType), dimension(3), intent(in) :: pt1, pt2, pt3
    integer(kind=intType), intent(in) :: i1, i2
    type(oversetString), intent(in) :: str

    ! Working
    logical :: triOverlap, inTri
    integer(kind=intType) :: i
    real(kind=realType) :: triNorm(3)

    ! Note: This is a dumb loop. We need to do a spatial serch here to
    ! only check the nodes around the current point. 

    call cross_prod(pt2-pt1, pt3-pt1, triNorm)
    triNorm = triNorm / norm2(triNorm)

    triOverlap = .False. 
    do i=1, str%nNodes
       if (i /= i1 .and. i/= i2) then
          if (dot_product(str%norm(:, i), triNorm) > 0.8) then 
             call pointInTriangle(pt1, pt2, pt3, str%x(:, i), inTri)
             if (inTri) then 
                triOverlap = .true. 
                exit
             end if
          end if
       end if
    end do
  end function triOverlap


end module stringOps
