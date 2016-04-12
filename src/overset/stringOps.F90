module stringOps

  use overset

contains
  subroutine nullifyString(string)

    use overset
    implicit none
    type(oversetString) :: string

    nullify(string%x, string%norm, string%h, string%ind, string%conn, &
         string%gc, string%otherID, string%nte, string%subStr, string%elemUsed)
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

    if (associated(string%h)) & 
         deallocate(string%h)

    if (associated(string%ind)) & 
         deallocate(string%ind)

    if (associated(string%conn)) & 
         deallocate(string%conn)

    if (associated(string%gc)) & 
         deallocate(string%gc)

    if (associated(string%otherID)) & 
         deallocate(string%otherID)

    if (associated(string%nte)) & 
         deallocate(string%nte)

    if (associated(string%subStr)) & 
         deallocate(string%subStr)

    if (associated(string%elemUsed)) &
         deallocate(string%elemUsed)

    call nullifyString(string)

  end subroutine deallocateString


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
    real(kind=realType), dimension(:, :), pointer :: xptr, normPtr
    real(kind=realType), dimension(:), pointer :: hPtr
    integer(kind=intType) , dimension(:), pointer :: indPtr

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
    hPtr => string%h
    indPtr => string%ind
    allocate(string%x(3, nUnique), string%norm(3, nUnique), string%h(nUnique), string%ind(nUnique))

    do i=1, string%nNodes
       string%x(:, link(i)) = xPtr(:, i)
       string%norm(:, link(i)) = normPtr(:, i)
       string%h(link(i)) = hPtr(i)
       string%ind(link(i)) = indPtr(i)
    end do
    string%nNodes = nUnique

    ! deallocate the pointer data which is actually the original data
    deallocate(xPtr, normPtr, indPtr, link, uniqueNodes)

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
    integer(kind=intType), dimension(:), pointer :: tmpGC
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
       tmpGC => string%gc
       allocate(string%conn(2, string%nElems - nDup))
       allocate(string%gc(string%nElems - nDup))
       j = 0
       do i=1, string%nElems
          if (duplicated(i) == 0) then 
             j = j + 1
             string%conn(:, j) = tmpConn(:, i)
             string%gc(i) = tmpGC(i)
          end if
       end do

       ! Set the new number of elements
       string%nElems = string%nElems - nDup

       ! Don't forget to deallocate the tmpConn pointer which is
       ! actually the original conn data.
       deallocate(tmpConn, tmpGC)

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

    allocate(s%x(3, s%nNodes), s%norm(3, s%nNodes), s%h(s%nNodes), s%ind(s%nNodes))
    do i=1, s%nNodes
       s%x(:, i) = p%x(:, s%pNodes(i))
       s%norm(:, i) = p%norm(:, s%pNodes(i))
       s%h(i) = p%h(s%pNodes(i))
       s%ind(i) = p%ind(s%pNodes(i))
    end do

    ! We can now create the local conn too, *USING THE LOCAL NODE NUMBERS*
    allocate(s%conn(2, s%nElems), s%gc(s%nElems))

    do i=1, s%nElems
       s%conn(1, i) = nodeUsed(p%conn(1, s%pElems(i)))
       s%conn(2, i) = nodeUsed(p%conn(2, s%pElems(i)))
       s%gc(i) = p%gc(s%pElems(i))
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
  integer(kind=intType) :: nNodes, nElems
  logical :: lastNodeZipper, inTri, overlapFound
  real(kind=realType), dimension(3) :: v1, v2, norm, c
  real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm
  integer(Kind=intType), dimension(:), allocatable :: nodeMap
  type(kdtree2_result), dimension(:), allocatable  :: results
  real(kind=realType), dimension(:, :), pointer :: xTmp, normTmp
  real(kind=realType), dimension(:), pointer :: hTmp
  integer(kind=intType), dimension(:, :), pointer :: connTmp
  integer(kind=intType), dimension(:), pointer :: indTmp, pNodesTmp
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
  allocate(nodeMap(s%nNodes))
  nodeMap = 1

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
                 s%p%tris(:, s%p%nTris) = (/s%pNodes(im1), s%pNodes(ii),s%pNodes(ip1)/)
                 lastNodeZipper = .True.
                 nZipped = nZipped + 1
                 nodeMap(ii) = 0
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
  xTmp => s%x
  normTmp => s%norm
  hTmp => s%h
  indTmp => s%ind
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

  ! Update the number of nodes/elems in our shorted chain. Every
  ! zipper reduces the number of nodes and number of elems by 1
  s%nNodes = s%nNodes - nZipped
  s%nElems = s%nElems - nZipped

  allocate(s%x(3, s%nNodes), s%norm(3, s%nNodes), s%h(s%nNodes), &
       s%ind(s%nNodes), s%pNodes(s%nNodes), s%conn(2, s%nElems))

  do i=1, nNodes
     if (nodeMap(i) /= 0) then 
        s%x(:, nodeMap(i)) = xTmp(:, i)
        s%norm(:, nodeMap(i)) = normTmp(:, i)
        s%h(nodeMap(i)) = hTmp(i)
        s%ind(nodeMap(i)) = indTmp(i)
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
  deallocate(xTmp, normTmp, indTmp, connTmp, pNodesTmp)
end subroutine selfZip

end module stringOps
