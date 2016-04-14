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
             string%gc(j) = tmpGC(i)
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
  integer(kind=intType) :: nNodes, nElems, elem1, elem2
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
                 s%p%tris(:, s%p%nTris) = (/s%pNodes(ip1), s%pNodes(ii),s%pNodes(im1)/)
                 lastNodeZipper = .True.
                 nZipped = nZipped + 1
                 nodeMap(ii) = 0


                 ! The two shorted string edges have been used for selfZip
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

  ! subroutine makeCrossZip(p, strings, nStrings)
  !   use overset
  !   implicit none

  !   ! Input/output
  !   integer(kind=intType), intent(in) :: nStrings
  !   type(oversetString), intent(inout) :: p, strings(nStrings)

  !   ! Local variables
  !   integer(kind=intType) :: i, j, inode, jnode, iDir, inzip, jnzip
  !   integer(kind=intType) :: inodeS, inodeE, jnodeS, jnodeE, jsym, isym
  !   integer(kind=intType), allocatable, dimension(:) :: nsplits
  !   integer(kind=intType) :: oID, oIDx
  !   ! ------------------------------------------------

  !   allocate(nsplits(nStrings))
    
  !   nSplits = 0
  !   ! do i=1, nStrings
  !   do i=nStrings,1, -1
  !      ! Decide the ranges of this i'th string and its neighbouring 
  !      ! strings

  !      inodeS = 1 !<-- start index of this split string
  !      inodeE = 1 !<-- end index of this split string
  !      jnodeS = 1 !<-- start index of the other zip pairing string
  !      jnodeE = 1 !<-- end index of the other zip pairing string

  !      ! print*,''
  !      ! print*,'string, myid: ',i,strings(i)%myID

  !      inode = 1
  !      loopinode: do while (inode <= strings(i)%nNodes-1)

  !         ! Change in nearest neighbor fullstringID, i.e. otherID(1,:)
  !         ! Make sure the split string has at least two nodes.
  !         splitifI: if ( (strings(i)%otherID(1, inode) /= &
  !                         strings(i)%otherID(1, inode+1)) .and. &
  !                         (inode - inodeS + 1) >= 2 .and. &
  !                         inode <= strings(i)%nNodes-2 ) then

  !            ! End index of this gap string
  !            inodeE = inode
             
  !            ! Potential pairing strings id
  !            j = strings(i)%otherID(1, inode)

  !            ! Loop over other gap string's nodes to set the 
  !            ! starting symmetric otherID point.
  !            isym = -1
  !            jsym = -1
  !            loopjnode: do jnode=1, strings(j)%nNodes
  !               oID  = strings(j)%otherID(1, jnode)  ! hopefully myid of string I 
  !               oIDx = strings(j)%otherID(2, jnode)  ! index in string I
  !               if (oIDx >=inodeS .and. oIDx <=inodeE) then
  !                  if (oID == strings(i)%myID .and. &
  !                      strings(i)%otherID(2, oIDx) == jnode) then
 
  !                     ! Found symmetric otherID 
  !                     isym = oIDx  ! index of strings(i)
  !                     jsym = jnode ! index of strings(j)

  !                     exit loopjnode
  !                  end if 
  !               end if
  !            end do loopjnode

  !            ! Couldn't find any symmetric nearest neighbor. Should 
  !            ! move to the next inode.
  !            if (isym == -1 .and. jsym == -1) cycle loopinode

  !            nsplits(i) = nsplits(i) + 1
   
  !            print *,' nsplits, I, J, inodeS, inodeE, isym, jsym ', nsplits(i), I, J, inodeS, inodeE, isym, jsym
  !            ! Advancing front starting from (isym, jsym)

  !            ! Loop over split sections of the 'i' string 
  !            ! direction 1: isym--> inodeE (increasing)
  !            ! direction 2: isym--> inodeS (decreasing)

  !            ! direction-1 (+ve):
  !            ! ------------------------------------------------------
  !            ! Advancing front A--B is inzip--jnzip
  !            inzip = isym ! A
  !            jnzip = jsym ! B
  !            iDir = 1     ! increasing direction of 'i' strings

  !            !print*,''
  !            !print*,'I, J, begin iDir ',I,J, iDir
  !            !print*,'-----------------'
  !            call stepZip(p, strings(i), strings(j), Inzip, Jnzip, iDir, inodeS, inodeE)
  !            !print*,'-----------------'
  !            !print*,'I, J,   end iDir ',I,J, iDir
  !            !print*,''

  !            ! direction-2 (-ve):
  !            ! ------------------------------------------------------
  !            ! Advancing front A--B is inzip--jnzip
  !            inzip = isym ! A
  !            jnzip = jsym ! B
  !            iDir = -1    ! decreasing direction of 'i' strings

  !            !print*,''
  !            !print*,'I, J, begin iDir ',I,J, iDir
  !            !print*,'-----------------'
  !            call stepZip(p, strings(i), strings(j), Inzip, Jnzip, iDir, inodeS, inodeE)
  !            !print*,'-----------------'
  !            !print*,'I, J,   end iDir ',I,J, iDir
  !            !print*,''

  !            ! Next split string starts after current inodeE
  !            inodeS = inodeE + 1
  !         end if splitifI

  !         inode = inode + 1
  !      end do loopinode

  !      ! Zip the remaining nodes of strings(i) with its neighbor string.

  !      ! End index of this gap string
  !      inodeE = strings(i)%nNodes
       
  !      ! Potential pairing strings id (inode = strings(i)%nNodes-1)
  !      j = strings(i)%otherID(1, inodeS)

  !      ! Loop over other gap string's nodes to set the 
  !      ! starting symmetric otherID point.
  !      isym = -1
  !      jsym = -1
  !      loopjnode1: do jnode=1, strings(j)%nNodes
  !         oID  = strings(j)%otherID(1, jnode)  ! hopefully myid of string I 
  !         oIDx = strings(j)%otherID(2, jnode)  ! index in string I
  !         if (oIDx >=inodeS .and. oIDx <=inodeE) then
  !            if (oID == strings(i)%myID .and. &
  !                strings(i)%otherID(1, oIDx) == j .and. &
  !                strings(i)%otherID(2, oIDx) == jnode) then
 
  !               ! Found symmetric otherID 
  !               isym = oIDx  ! index of strings(i)
  !               jsym = jnode ! index of strings(j)

  !               exit loopjnode1
  !            end if 
  !         end if
  !      end do loopjnode1

  !      print *,' last:    I, J, inodeS, inodeE, isym, jsym             ', I,J, inodeS, inodeE, isym, jsym
  !      print *,''

  !      ! Couldn't find any symmetric nearest neighbor. *** PROBLEM ***
  !      if (isym == -1 .and. jsym == -1) STOP ' Problem finding symmetric nodes '

  !      ! Advancing front starting from (isym, jsym)

  !      ! Loop over split sections of 'i' string
  !      ! direction 1: isym--> inodeE (increasing)
  !      ! direction 2: isym--> inodeS (decreasing)

  !      ! direction-1 (+ve):
  !      ! ------------------------------------------------------
  !      ! Advancing front A--B is inzip--jnzip
  !      inzip = isym ! A
  !      jnzip = jsym ! B
  !      iDir = 1     ! increasing direction of 'i' strings

  !      !print*,''
  !      !print*,'last: I, J, begin iDir ',I,J, iDir
  !      !print*,'-----------------'
  !      call stepZip(p, strings(i), strings(j), Inzip, Jnzip, iDir, inodeS, inodeE)
  !      !print*,'-----------------'
  !      !print*,'last: I, J,   end iDir ',I,J, iDir
  !      !print*,''

  !      ! direction-2 (-ve):
  !      ! ------------------------------------------------------
  !      ! Advancing front A--B is inzip--jnzip
  !      inzip = isym ! A
  !      jnzip = jsym ! B
  !      iDir = -1    ! decreasing direction of 'i' strings

  !      !print*,''
  !      !print*,'last: I, J, begin iDir ',I,J, iDir
  !      !print*,'-----------------'
  !      call stepZip(p, strings(i), strings(j), Inzip, Jnzip, iDir, inodeS, inodeE)
  !      !print*,'-----------------'
  !      !print*,'last: I, J,   end iDir ',I,J, iDir
  !      !print*,''

  !   end do ! nStrings

  !   deallocate(nsplits)
  ! end subroutine makeCrossZip

  ! recursive subroutine stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
  !   !
  !   ! ---------------------------------------------------------------------
  !   ! Choose between B or B+ based on some tests.
  !   !
  !   !               A-B-B+     .or.        A-B-A+
  !   !           --------------         -------------
  !   !  String1:     A-----A+               A---A+
  !   !              / \                    /  . 
  !   !             /   \        .or.      / .  
  !   !  String2:  B-----B+               B-----B+
  !   !
  !   ! 
  !   ! Several tests need to be made to ensure the triangle formed:
  !   !   - Doesn't contain any other points of either of the gap strings pair
  !   !   - Doesn't create a concave triangle. 
  !   !   
  !   !     In other words, the above two can be ensured if none of the newly 
  !   !     created sides of the triangle intersects with any gap strings.
  !   ! ---------------------------------------------------------------------
  !   !
  !   use overset
  !   implicit none

  !   ! Input/output
  !   integer(kind=intType), intent(in)    :: inodeS, inodeE
  !   integer(kind=intType), intent(inout) :: Inzip, Jnzip, iDir
  !   type(oversetString),   intent(in)    :: stringI, stringJ
  !   type(oversetString),   intent(inout) :: p

  !   ! Local variables
  !   logical :: ipend, jpend, foundselfZipEdge
  !   integer(kind=intType) :: n1, n2, n3, InzipP, JnzipP, nnte, inte, elem, curElem
  !   integer(kind=intType) :: oID_i, oIDx_i, oID_j, oIDx_j, elem1, elem2, elemp1, elemp2
  !   real(kind=realType) :: triNorm1(3), triNorm2(3), Coor(1:3, 1:4), r1(3), r2(3)
  !   real(kind=realType) :: angsum1, angsum2, angtmp
  !   ! -----------------------------------------------------

  !   ! if A+ or B+ are not present
  !   ipend = .False.
  !   jpend = .False.

  !   ! Initialize 
  !   ! Coordinates of A  = coor(:, 1) 
  !   ! Coordinates of A+ = coor(:, 2) 
  !   ! Coordinates of B  = coor(:, 3) 
  !   ! Coordinates of B+ = coor(:, 4) 
  !   coor(1:3, 1:4) = -Large

  !   InzipP = Inzip + iDir
  !   JnzipP = Jnzip - iDir ! opposite of stringI node traverse direction
  !   if ( InzipP <inodeS .or. InzipP > inodeE) ipend = .True.
  !   if ( JnzipP <1 .or. JnzipP > stringJ%nNodes) jpend = .True.

  !   ! If neither A+ or B+ exists, exit routine.
  !   if (ipend .and. jpend) then
  !      ! print*,'return here : iDir, Inzip, Jnzip ', iDir, inzip, jnzip
  !      return
  !   end if
    
 
  !   ! Point A:
  !   coor(1:3, 1) = stringI%x(1:3, Inzip)
  !   ! Point A+:
  !   !am? if (.not.ipend) coor(1:3, 2) = stringI%x(1:3, InzipP)

  !   if (.not.ipend) then

  !      oID_i = stringI%otherID(1, InzipP)  ! stringJ
  !      oIDx_i = stringI%otherID(2, InzipP) ! index of stringJ

  !      !am? if (stringJ%otherID(1, oIDx_i) == stringI%myID) then
  !      if ( (abs(Jnzip-1) < abs(Jnzip-stringJ%nNodes) .and. &
  !            abs(oIDx_i-1) > abs(oIDx_i-stringJ%nNodes)) .or. &
  !           (abs(Jnzip-1) > abs(Jnzip-stringJ%nNodes) .and. &
  !            abs(oIDx_i-1) < abs(oIDx_i-stringJ%nNodes)) ) then
  !         ! If oIDx_i is closer to the other end of periodic stringJ
  !         ! then A+ is not valid.
  !         ipend = .True.
  !         ! print*,'======> not valid A+, I, InzipP ', stringI%myid, InzipP
  !      else
  !         coor(1:3, 2) = stringI%x(1:3, InzipP)
  !      end if
  !   end if
        
  !   ! Point B:
  !   coor(1:3, 3) = stringJ%x(1:3, Jnzip)
  !   ! Point B+:
  !   if (.not.jpend) then
  !      oID_j = stringJ%otherID(1, JnzipP)  ! stringI
  !      oIDx_j = stringJ%otherID(2, JnzipP) ! index of stringI

  !      !am? if (stringJ%otherID(1, JnzipP) == stringI%myID) then
  !      if (oID_j == stringI%myID .and. oIDx_j >= inodeS &
  !          .and. oIDx_j <= inodeE) then
  !         ! B+ is valid if its nearest neighbor is still stringI
  !         ! and within the current stringI index bounds.

  !         ! Special consideration for periodic stringI
  !         if ( (abs(Inzip-1) < abs(Inzip-stringI%nNodes) .and. &
  !               abs(oIDx_j-1) > abs(oIDx_j-stringI%nNodes)) .or. &
  !              (abs(Inzip-1) > abs(Inzip-stringI%nNodes) .and. &
  !               abs(oIDx_j-1) < abs(oIDx_j-stringI%nNodes)) ) then
  !            ! If oIDx_j is closer to the other end of this periodic string
  !            ! then B+ is not valid.
  !            jpend = .True.
  !            ! print*,'======> not valid B+, J, JnzipP', stringJ%myid, JnzipP
  !         else
  !            coor(1:3, 4) = stringJ%x(1:3, JnzipP)
  !         end if
  !      else
  !         jpend = .True.
  !      end if
  !   end if

  !   ! If neither A+ or B+ is valid node, exit routine.
  !   if (ipend .and. jpend) then
  !      return
  !   end if
 
  !   ! Add third point if it passes the following tests

  !   ! "Front angle test"
    
  !   ! Front angle test will be valid only when both A+ and B+ are present.
  !   angsum1 = 0.0
  !   angsum2 = 0.0

  !   ifall4nodes: if (.not.ipend .and. .not.jpend) then
  !      ! choose A+
  !      ! ------------------------------------------------
  !      angsum1 = 0.0
  !      ! A-A+-B 
  !      r1(1:3) = coor(1:3, 1) - coor(1:3, 2)
  !      r2(1:3) = coor(1:3, 3) - coor(1:3, 2)

  !      ! cos(theta) = (vec(r1) .dot vec(r2))/(abs(vec(r1))*abs(vec(r2)))
  !      angtmp = 0.0
  !      angtmp = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  !      angtmp = angtmp/(sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)*&
  !                         sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2))
  !      angsum1 = angsum1 + abs(acos(angtmp))

  !      ! A+-B-B+ 
  !      r1(1:3) = coor(1:3, 2) - coor(1:3, 3)
  !      r2(1:3) = coor(1:3, 4) - coor(1:3, 3)

  !      ! cos(theta) = (vec(r1) .dot vec(r2))/(abs(vec(r1))*abs(vec(r2)))
  !      angtmp = 0.0
  !      angtmp = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  !      angtmp = angtmp/(sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)*&
  !                         sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2))
  !      angsum1 = angsum1 + abs(acos(angtmp))
  !      ! --- end choose A+ ------------------------------

  !      ! choose B+
  !      ! ------------------------------------------------
  !      angsum2 = 0.0

  !      ! A+-A-B+
  !      r1(1:3) = coor(1:3, 2) - coor(1:3, 1)
  !      r2(1:3) = coor(1:3, 4) - coor(1:3, 1)

  !      ! cos(theta) = (vec(r1) .dot vec(r2))/(abs(vec(r1))*abs(vec(r2)))
  !      angtmp = 0.0
  !      angtmp = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  !      angtmp = angtmp/(sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)*&
  !                         sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2))
  !      angsum2 = angsum2 + abs(acos(angtmp))

  !      ! A-B+-B 
  !      r1(1:3) = coor(1:3, 1) - coor(1:3, 4)
  !      r2(1:3) = coor(1:3, 3) - coor(1:3, 4)

  !      ! cos(theta) = (vec(r1) .dot vec(r2))/(abs(vec(r1))*abs(vec(r2)))
  !      angtmp = 0.0
  !      angtmp = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  !      angtmp = angtmp/(sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)*&
  !                         sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2))
  !      angsum2 = angsum2 + abs(acos(angtmp))
  !      ! --- end choose B+ ------------------------------

  !   end if ifall4nodes

  !   ifangsum: if ( (.not.ipend .and. .not.jpend .and. angsum1 > angsum2) &
  !                 .or. jpend) then

  !      ! 1. Triangle 1 (A-A+-B) (1-2-3)
  !      !    Pick A+ if it passes the angle sum test or if B+ doesn't exist at all


  !      !! Make sure the triangle A-A+-B does not contain any gap nodes
  !      !! (To be done later).

  !      ! Identify A-A+ edge
  !      ! -------------------------------------------
  !      nnte = stringI%p%nte(1, stringI%pNodes(Inzip))

  !      curElem = -1

  !      foundselfZipEdge = .False.
  !      inteI: do inte=1, nnte
  !         elem = stringI%p%nte(1+inte, stringI%pNodes(Inzip))
  !         n1 = stringI%p%conn(1, elem)
  !         n2 = stringI%p%conn(2, elem)
  !         if ( (n1 == stringI%pNodes(Inzip) .and. &
  !               n2 == stringI%pNodes(InzipP)) .or. &
  !              (n2 == stringI%pNodes(Inzip) .and. &
  !               n1 == stringI%pNodes(InzipP)) ) then
  !            curElem = elem ! A-A+ edge
  !            exit inteI
  !         end if
  !      end do inteI
       
  !      if (curElem == -1) then

  !         ! Could not find any full string edge to either A or A+.
  !         ! Possibly A-A+ is a gap edge of a selfZip triangle
  !         elem1 = stringI%p%nte(2, stringI%pNodes(Inzip))
  !         elem2 = stringI%p%nte(3, stringI%pNodes(Inzip))
  !         elemp1 = stringI%p%nte(2, stringI%pNodes(InzipP))
  !         elemp2 = stringI%p%nte(3, stringI%pNodes(InzipP))
  !         if ( (stringI%p%elemUsed(elem1)==1 .or. &
  !               stringI%p%elemUsed(elem2)==1) .and. &
  !               (stringI%p%elemUsed(elemp1)==1 .or. &
  !               stringI%p%elemUsed(elemp2)==1)) then
  !            foundselfZipEdge = .True.
             
  !         end if
  !      end if

  !      if (curElem == -1 .and. .not.foundselfZipEdge) STOP ' Problem identifying curElem A-A+ in stepZip'
  !      ! -------------------------------------------
       
  !      ! Zip the triangle if the edge A-A+ has not been used
  !      if ( curElem /= -1 .or. foundselfZipEdge) then

  !         ! Skip A+ test 1: 
  !         ! Skip A-A+ if edge A-A+ is already used for zipping triangle.
  !         if (curElem /= -1) then 
  !            if (stringI%p%elemUsed(curElem) == 1) then

  !               ! Do not add A-A+ to any triangle. Move front to A+-B
  !               ! ---------------------------------------------------
  !               Inzip = InzipP

  !               ! Jnzip remains the same. 
  !               call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
  !               return 
  !            end if
  !         end if
 
  !         ! Skip A+ test 2:
  !         ! Skip A-A+ if A-A+ is a selfZip edge and these nodes have
  !         ! already been used for crossZip.
  !         if (foundselfZipEdge .and. &
  !             stringI%p%selfXzipNodeUsed(stringI%pNodes(Inzip)) == 1 .and. &
  !             stringI%p%selfXzipNodeUsed(stringI%pNodes(InzipP)) == 1) then

  !            ! Do not add A-A+ to any triangle. Move front to A+-B
  !            ! ---------------------------------------------------
  !            Inzip = InzipP

  !            ! Jnzip remains the same. 
  !            call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
  !            return 
  !         end if

  !         ! ------------------------------------------------------------------------
  !         ! Skip A+ test 3: 
  !         ! Skip A-A+ any master full string node is inside potential triangle A-A+-B
  !         ! TBD later.
  !         ! ------------------------------------------------------------------------

  !         stringI%p%nTris = stringI%p%nTris + 1
  
  !         ! Decide node order based on iDir. Triangle node direction should
  !         ! be opposite of fullStrings node order.
  !         if (Inzip < InzipP) then
  !            ! Inzip < InzipP : InzipP-->Inzip-->Jnzip (A+-->A-->B)
  !            stringI%p%tris(:, stringI%p%nTris) = &
  !               (/stringI%pNodes(InzipP), stringI%pNodes(Inzip), &
  !                 stringJ%pNodes(Jnzip) /)

  !            if (foundSelfZipEdge) &
  !            print*,'stringI, n1, n2 ',stringI%myid, stringI%pNodes(InzipP), stringI%pNodes(Inzip)

  !            ! Add master oversetEdge's
  !            ! Edge 1: InzipP-->Inzip: A+-->A 
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringI%pNodes(InzipP)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            ! Edge 2: Inzip-->Jnzip: A-->B
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            ! Edge 3: Jnzip-->InzipP: B-->A+
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringI%pNodes(InzipP)

  !            !am? ! Add master oversetEdge's
  !            !am? ! Edge 1: InzipP-->Inzip: A+-->A 
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(InzipP)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            !am? ! Edge 2: Inzip-->Jnzip: A-->B
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            !am? ! Edge 3: Jnzip-->InzipP: B-->A+
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(InzipP)

  !         else 
  !            ! Inzip > InzipP : Inzip-->InzipP-->Jnzip (A-->A+-->B)
  !            stringI%p%tris(:, stringI%p%nTris) = &
  !               (/stringI%pNodes(Inzip), stringI%pNodes(InzipP), &
  !                 stringJ%pNodes(Jnzip) /)

  !            if (foundSelfZipEdge) &
  !            print*,'stringI, n1, n2 ',stringI%myid, stringI%pNodes(Inzip), stringI%pNodes(InzipP)

  !            ! Add master oversetEdge's
  !            ! Edge 1: Inzip-->InzipP: A-->A+
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringI%pNodes(InzipP)
  !            ! Edge 2: InzipP-->Jnzip: A+-->B
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringI%pNodes(InzipP)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            ! Edge 3: Jnzip-->Inzip: B-->A
  !            stringI%p%nEdges = stringI%p%nEdges + 1
  !            stringI%p%edges(stringI%p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            stringI%p%edges(stringI%p%nEdges)%n2 = stringI%pNodes(Inzip)

  !            !am? ! Add master oversetEdge's
  !            !am? ! Edge 1: Inzip-->InzipP: A-->A+
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(InzipP)
  !            !am? ! Edge 2: InzipP-->Jnzip: A+-->B
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(InzipP)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            !am? ! Edge 3: Jnzip-->Inzip: B-->A
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(Inzip)
  !         end if

  !         ! If A-A+ is a master edge, remember it has already been used
  !         if (curElem /= -1 .and. .not.foundSelfZipEdge) then
  !            stringI%p%elemUsed(curElem) = 1
  !            stringI%p%elemUsedDebug(curElem) = stringI%p%elemUsedDebug(curElem) + 1
  !         end if

  !        ! If A-A+ is a selfZip edge, remember it was used
  !        if (foundSelfZipEdge) then
  !           stringI%p%selfXzipNodeUsed(stringI%pNodes(Inzip)) = 1
  !           stringI%p%selfXzipNodeUsed(stringI%pNodes(InzipP)) = 1
  !        end if
          
  !      end if
  !      ! print*,'I, J, Inzip, Jnzip ',stringI%myid, stringJ%myid, Inzip, Jnzip

  !      ! Move front 
  !      ! -------------------------
  !      ! The front is A+-B
  !      Inzip = InzipP

  !      ! Jnzip remains the same. 
  !      call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
       

  !   else if ((.not.ipend .and. .not.jpend .and. angsum1 < angsum2) &
  !            .or. ipend) then ifangsum

  !      ! 2. Triangle 2 (A-B-B+) (1-3-4)
  !      !    Pick B+ if it passes the angle sum test or if A+ doesn't exist at all

  !      !! Make sure the triangle A-B-B+ does not contain any gap nodes
  !      !! (To be done later).

  !      ! Identify B-B+ edge
  !      ! -------------------------------------------
  !      nnte = stringJ%p%nte(1, stringJ%pNodes(Jnzip))

  !      curElem = -1

  !      foundselfZipEdge = .False.
  !      inteJ: do inte=1, nnte
  !         elem = stringJ%p%nte(1+inte, stringJ%pNodes(Jnzip))
  !         n1 = stringJ%p%conn(1, elem)
  !         n2 = stringJ%p%conn(2, elem)

  !         if ( (n1 == stringJ%pNodes(Jnzip) .and. &
  !               n2 == stringJ%pNodes(JnzipP)) .or. &
  !              (n2 == stringJ%pNodes(Jnzip) .and. &
  !               n1 == stringJ%pNodes(JnzipP)) ) then
  !            curElem = elem ! B-B+ edge
  !            exit  inteJ
  !         end if

  !      end do inteJ

  !      if (curElem == -1) then

  !         ! Could not find any full string edge to either A or A+.
  !         ! Possibly A-A+ is a gap edge of a selfZip triangle
  !         elem1 = stringJ%p%nte(2, stringJ%pNodes(Jnzip))
  !         elem2 = stringJ%p%nte(3, stringJ%pNodes(Jnzip))
  !         elemp1 = stringJ%p%nte(2, stringJ%pNodes(JnzipP))
  !         elemp2 = stringJ%p%nte(3, stringJ%pNodes(JnzipP))
  !         if ( (stringJ%p%elemUsed(elem1)==1 .or. &
  !               stringJ%p%elemUsed(elem2)==1) .and. &
  !              (stringJ%p%elemUsed(elemp1)==1 .or. &
  !               stringJ%p%elemUsed(elemp2)==1)) then
  !            foundselfZipEdge = .True.
  !         end if
  !      end if

  !      if (curElem == -1 .and. .not.foundselfZipEdge) STOP ' Problem identifying curElem B-B+ in stepZip'
  !      ! -------------------------------------------

  !      ! Zip the triangle if the edge B-B+ has not been used
  !      if (curElem /= -1 .or. foundselfZipEdge) then
       
  !         ! Skip B+ test 1: 
  !         ! Skip B-B+ if edge B-B+ is already used for zipping triangle.
  !         if (curElem /= -1) then
  !            if (stringJ%p%elemUsed(curElem) == 1) then
                
  !               ! Do not use B-B+ for any triangle. Move front to A-B+
  !               ! ----------------------------------------------------
  !               Jnzip = JnzipP

  !               ! Inzip, iDir remain the same. 
  !               call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
  !               return 
  !            end if
  !         end if

  !         ! Skip B+ test 2:
  !         ! Skip B-B+ if B-B+ is a selfZip edge and these nodes have
  !         ! already been used for crossZip.
  !         if (foundselfZipEdge .and. &
  !             stringJ%p%selfXzipNodeUsed(stringJ%pNodes(Jnzip)) == 1 .and. &
  !             stringJ%p%selfXzipNodeUsed(stringJ%pNodes(JnzipP)) == 1) then

  !            ! Do not add B-B+ to any triangle. Move front to A-B+
  !            ! ---------------------------------------------------
  !            Jnzip = JnzipP

  !            ! Inzip remains the same. 
  !            call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)
  !            return 
  !         end if

  !         ! ------------------------------------------------------------------------
  !         ! Skip B+ test 3: 
  !         ! Skip B-B+ any master full string node is inside potential triangle A-B-B+
  !         ! TBD later.
  !         ! ------------------------------------------------------------------------

  !         stringJ%p%nTris = stringJ%p%nTris + 1

  !         ! Decide node order based on jDir. Triangle node direction should
  !         ! be opposite of fullStrings node order.
  !         if (Jnzip < JnzipP) then
  !            ! Jnzip < JnzipP : JnzipP-->Jnzip-->Inzip (B+-->B-->A)
  !            stringJ%p%tris(:, stringJ%p%nTris) = &
  !               (/stringJ%pNodes(JnzipP), stringJ%pNodes(Jnzip), &
  !                 stringI%pNodes(Inzip) /)

  !            if (foundSelfZipEdge) &
  !            print*,'stringJ, n1, n2 ',stringJ%myid, stringJ%pNodes(JnzipP), stringJ%pNodes(Jnzip)

  !            ! Add master oversetEdge's
  !            ! Edge 1: JnzipP-->Jnzip: B+-->B 
  !            stringJ%p%nEdges = stringJ%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringJ%pNodes(JnzipP)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            ! Edge 2: Jnzip-->Inzip: B-->A
  !            stringJ%p%nEdges = stringJ%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            ! Edge 3: Inzip-->JnzipP: A-->B+
  !            stringJ%p%nEdges = stringI%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringJ%pNodes(JnzipP)

  !            !am? ! Add master oversetEdge's
  !            !am? ! Edge 1: JnzipP-->Jnzip: B+-->B 
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(JnzipP)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(Jnzip)
  !            !am? ! Edge 2: Jnzip-->Inzip: B-->A
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            !am? ! Edge 3: Inzip-->JnzipP: A-->B+
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(JnzipP)

  !         else
  !            ! Jnzip > JnzipP : Jnzip-->JnzipP-->Inzip (B-->B+-->A)
  !            stringJ%p%tris(:, stringJ%p%nTris) = &
  !               (/stringJ%pNodes(Jnzip), stringJ%pNodes(JnzipP), &
  !                 stringI%pNodes(Inzip) /)

  !            if (foundSelfZipEdge) &
  !            print*,'stringJ, n1, n2 ',stringJ%myid, stringJ%pNodes(Jnzip), stringJ%pNodes(JnzipP)

  !            ! Add master oversetEdge's
  !            ! Edge 1: Jnzip-->JnzipP: B-->B+ 
  !            stringJ%p%nEdges = stringJ%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringJ%pNodes(JnzipP)
  !            ! Edge 2: JnzipP-->Inzip: B+-->A
  !            stringJ%p%nEdges = stringJ%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringJ%pNodes(JnzipP)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            ! Edge 3: Inzip-->Jnzip: A-->B
  !            stringJ%p%nEdges = stringI%p%nEdges + 1
  !            stringJ%p%edges(stringJ%p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            stringJ%p%edges(stringJ%p%nEdges)%n2 = stringJ%pNodes(Jnzip)

  !            !am? ! Add master oversetEdge's
  !            !am? ! Edge 1: Jnzip-->JnzipP: B-->B+ 
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(Jnzip)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(JnzipP)
  !            !am? ! Edge 2: JnzipP-->Inzip: B+-->A
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringJ%pNodes(JnzipP)
  !            !am? p%edges(p%nEdges)%n2 = stringI%pNodes(Inzip)
  !            !am? ! Edge 3: Inzip-->Jnzip: A-->B
  !            !am? p%nEdges = p%nEdges + 1
  !            !am? p%edges(p%nEdges)%n1 = stringI%pNodes(Inzip)
  !            !am? p%edges(p%nEdges)%n2 = stringJ%pNodes(Jnzip)

  !         end if

  !         if (curElem /= -1 .and. .not.foundSelfZipEdge) then
  !            stringJ%p%elemUsed(curElem) = 1
  !            stringJ%p%elemUsedDebug(curElem) = stringJ%p%elemUsedDebug(curElem) + 1
  !         end if

  !        ! If B-B+ is a selfZip edge, remember it was used
  !        if (foundSelfZipEdge) then
  !           stringJ%p%selfXzipNodeUsed(stringJ%pNodes(Jnzip)) = 1
  !           stringJ%p%selfXzipNodeUsed(stringJ%pNodes(JnzipP)) = 1
  !        end if
  !      end if
  !      ! print*,'I, J, Inzip, Jnzip ',stringI%myid, stringJ%myid, Inzip, Jnzip

  !      ! Move front 
  !      ! -------------------------
  !      ! The front is A-B+
  !      Jnzip = JnzipP

  !      ! Inzip, iDir remain the same. 
  !      call stepZip(p, stringI, stringJ, Inzip, Jnzip, iDir, inodeS, inodeE)

  !   else
  !      print*,' Problem picking up either A+ or B+ '
  !      Stop 

  !   end if ifangsum

  ! end subroutine stepZip

  ! subroutine makePocketZip(p, strings, nStrings)
  !   use overset
  !   implicit none

  !   ! Input/output
  !   integer(kind=intType), intent(in) :: nStrings
  !   type(oversetString), intent(in) :: p, strings(nStrings)

  !   ! Local variables
  !   integer(kind=intType) :: i, nsum1, nsum2, ndiff1, ndiff2, ipedge, icur
  !   integer(kind=intType) :: n1, n2, npolyEdges, npolyEdgestmp, nEdgeUsed
  !   integer(kind=intType) :: nNodes1, nNodes2, cn1, cn2, str1, str2, nends
  !   type(oversetEdge), pointer, dimension(:) :: polyEdges, polyEdgestmp
  !   type(pocketEdge), pointer, dimension(:) :: pocketEdges
  !   integer(kind=intType), allocatable, dimension(:) :: edgeMap, edgeMaptmp
  !   integer(kind=intType), allocatable, dimension(:) :: nodeList, nodeMap
  !   logical :: isEndEdge
  !   logical, allocatable, dimension(:) :: PocketEdgeUsed

  !   type(oversetString) :: pocketMaster
  !   type(oversetString), pointer  :: pocketStrings, strPkt
  !   integer(kind=intType) :: npocketEdges, npocketStrings, nPktEdges, nUnique
  !   integer(kind=intType) :: ip, curElem, nElems, iStart, firstElem
  !   real(kind=realType) :: triArea
  !   ! ---------------------------------------------------------------


  !   ! ---------------------------------------------------------------
  !   ! PocketZip 1: 
  !   ! First sort the edges.
  !   ! ---------------------------------------------------------------
  !   call qsortEdgeType(p%Edges, p%nEdges)

  !   ! Check order of edges
  !   ! do i=1, p%nEdges
  !   !    write(2001,*)i, p%Edges(i)%n1, p%Edges(i)%n2, &
  !   !                (p%edges(i)%n1+p%edges(i)%n2)*1.e5 + abs(p%edges(i)%n2-p%edges(i)%n1)
  !   ! end do

  !   ! ---------------------------------------------------------------------
  !   ! PocketZip 2:
  !   ! Now cancel out edges counted twice in master%Edges. The edges are
  !   ! sorted such that the edges with same nodes are in consequetive order.
  !   ! ---------------------------------------------------------------------

  !   ! Over estimate of remaining pocket edges to zip
  !   allocate(polyEdges(p%nEdges)) 
  !   allocate(polyEdgestmp(p%nEdges)) 
  !   allocate(edgeMap(p%nEdges))
  !   allocate(edgeMaptmp(p%nEdges))
  !   edgeMap = -1
  !   edgeMaptmp = -1

  !   ! Initialize
  !   npolyEdgestmp = p%nEdges
  !   edgeMaptmp = -1
  !   do i=1, p%nEdges
  !      polyEdgesTmp(i)%n1 = p%Edges(i)%n1
  !      polyEdgesTmp(i)%n2 = p%Edges(i)%n2
  !      edgeMaptmp(i) = i
  !   end do

  !   ! Eliminate the edges going through the ordered edges. 
  !   ! The sorted opposite edges are canceled in pairs.
  !   npolyEdges = 0

  !   nends = 0

  !   loop_outer: do 
       
  !      npolyEdges = 0
  !      edgeMap = -1
  !      i = 1
  !      nends = 0
  !      loopEdge: do while (i < npolyEdgestmp) 

  !         nsum1 = (polyEdgesTmp(i)%n1   + polyEdgesTmp(i)%n2)
  !         nsum2 = (polyEdgesTmp(i+1)%n1 + polyEdgesTmp(i+1)%n2)

  !         ndiff1 = (polyEdgesTmp(i)%n2   - polyEdgesTmp(i)%n1)
  !         ndiff2 = (polyEdgesTmp(i+1)%n2 - polyEdgesTmp(i+1)%n1)
         
  !         ! If the edge joins the end nodes of single sided 
  !         ! fullStrings pair, eliminate it.
  !         ! -----------------------------------------------------------
  !         isEndEdge = .False.
          
  !         ! Parent nodes of this edge
  !         n1 = polyEdgesTmp(i)%n1
  !         n2 = polyEdgesTmp(i)%n2

  !         str1    = p%cNodes(1, n1) ! node1's child fullStrings ID
  !         cn1     = p%cNodes(2, n1) ! node1's child fullStrings node index
  !         nNodes1 = p%cNodes(3, n1) ! node1's child fullStrings nNodes size

  !         str2    = p%cNodes(1, n2) ! node2's child fullStrings ID
  !         cn2     = p%cNodes(2, n2) ! node2's child fullStrings node index
  !         nNodes2 = p%cNodes(3, n2) ! node2's child fullStrings nNodes size

  !         if (str1 /= str2 ) then
  !            if (.not.strings(str1)%isperiodic .and. &
  !                .not.strings(str2)%isperiodic .and. &
  !                (cn1==1 .or. cn1==nNodes1) .and. (cn2==1 .or. cn2==nNodes2)) then
  !               ! This is an end edge, eliminate this one too.
  !               nends = nends + 1

  !               i = i + 1
  !               cycle loopEdge
  !            end if
  !         end if
  !         ! --- End end nodes check -----------------------------------
          

  !         if (nsum1 == nsum2 .and. ndiff1 + ndiff2 == 0) then
  !            ! Found ordered edges pair. Eliminate these two edges.

  !            ! Jump to 2nd next edge
  !            i = i + 2
  !            cycle loopEdge
  !         else
  !            ! Add this edge to new edge list
  !            npolyEdges = npolyEdges + 1
  !            polyEdges(npolyEdges)%n1 = polyEdgesTmp(i)%n1
  !            polyEdges(npolyEdges)%n2 = polyEdgesTmp(i)%n2
             
  !            edgeMap(npolyEdges) = edgeMaptmp(i)
           
  !            if (i+1 == npolyEdgestmp) then
  !               ! Add this last edge to new edge list
  !               npolyEdges = npolyEdges + 1
  !               polyEdges(npolyEdges)%n1 = polyEdgesTmp(i+1)%n1
  !               polyEdges(npolyEdges)%n2 = polyEdgesTmp(i+1)%n2
                
  !               edgeMap(npolyEdges) = edgeMaptmp(i+1)
 
  !               exit loopEdge
  !            end if   

  !            ! Loop to next edge
  !            i = i + 1
  !         end if
  !      end do loopEdge
  !      print*, '===========================> nends ',nends

  !      if (npolyEdgestmp - npolyEdges == 0) then
  !         print*, ' No more edges to cancel ', npolyEdges
  !         exit loop_outer
  !      end if

  !      ! Update polyedgesTmp  and edgeMaptmp before cycling loop
  !      ! -------------------------------------------------------
  !      ! First zero out
  !      do i=1, npolyEdgesTmp
  !         polyEdgesTmp(i)%n1 = -1
  !         polyEdgesTmp(i)%n2 = -1
  !      end do
  !      edgeMaptmp = -1

  !      ! New tmp sizes
  !      npolyEdgestmp = npolyEdges
  !      do i=1, npolyEdges
  !         polyEdgesTmp(i)%n1 = polyEdges(i)%n1
  !         polyEdgesTmp(i)%n2 = polyEdges(i)%n2
  !         edgeMaptmp(i) = edgeMap(i)
  !      end do
  !      ! -------------------------------------------------------

  !   end do loop_outer
  !   ! ---------------------------------------------------------------------
  !   ! End canceling out edges.
  !   ! ---------------------------------------------------------------------
  !   print*, 'npolyEdges ', npolyEdges
  !   do i=1, npolyEdges
  !      write(4001,*)i, edgeMap(i), polyEdges(i)%n1, polyEdges(i)%n2, &
  !                  (polyedges(i)%n1+polyedges(i)%n2)*1.e5 + (polyedges(i)%n2-polyedges(i)%n1)
  !   end do

  !   ! Debug polygonEdges
  !   ! -----------------------------------------------------------
  !   open(unit=101, file="polygonEdges.dat", form='formatted')
  !   write(101,*) 'TITLE = "PolygonEdges Data" '
  !   write(101,*) 'Variables = "X", "Y", "Z"'
  !   write(101,*) "Zone T=Pockets"
  !   write (101,*) "Nodes = ", npolyEdges*2, " Elements= ", npolyEdges, " ZONETYPE=FELINESEG"
  !   write (101,*) "DATAPACKING=POINT"

  !   ! node data
  !   do i=1, npolyEdges
  !      n1 = polyEdges(i)%n1
  !      n2 = polyEdges(i)%n2
  !      ! node 1
  !      write(101,'(3(E20.12,x))')p%x(1, n1), p%x(2, n1), p%x(3, n1)
  !      ! node 2
  !      write(101,'(3(E20.12,x))')p%x(1, n2), p%x(2, n2), p%x(3, n2)
  !   end do

  !   ! Edge data
  !   do i=1, npolyEdges
  !      write(101,'(3(I5,x))')2*i-1, 2*i
  !   end do
  !   close(101)
  !   !-------------------------------------------------------------

  !   !-------------------------------------------------------------
  !   ! PocketZip 3:
  !   ! Accumulate full pocket string edges.
  !   ! Do similar to how master edges were accumulated.
  !   ! 3.1: create pocketMaster elems and nodes
  !   ! 3.2: perform doChain on pocketMaster elems and create pocket strings.
  !   ! 3.3: selfZip pocketStrings
  !   !-------------------------------------------------------------

  !   ! ----------------------------------------
  !   ! 3.1: create pocketMaster elems and nodes
  !   ! ----------------------------------------
    
  !   ! First create unique nodeList
  !   ! Number of nodes is twice the size of edges
  !   allocate(nodeList(2*nPolyEdges), nodeMap(2*nPolyEdges))
  !   do i=1, nPolyEdges
  !      nodeList(2*i-1) = PolyEdges(i)%n1
  !      nodeList(2*i)   = PolyEdges(i)%n2
  !   end do

  !   call unique(nodeList, 2*nPolyEdges, nUnique, nodeMap)
  !   do i=1, nPolyEdges
  !      write(6001,*)2*i-1,nodeList(2*i-1), nodeMap(2*i-1), polyEdges(i)%n1
  !      write(6001,*)2*i,  nodeList(2*i),   nodeMap(2*i),   polyEdges(i)%n2
  !   end do

  !   ! Define pocketMaster string
  !   call nullifyString(pocketMaster)

  !   pocketMaster%nNodes = nUnique
  !   pocketMaster%nElems = nPolyEdges

  !   allocate(pocketMaster%x(3, pocketMaster%nNodes), &
  !            pocketMaster%norm(3, pocketMaster%nNodes), &
  !            pocketMaster%ind(pocketMaster%nNodes), &
  !            pocketMaster%pNodes(pocketMaster%nNodes), &
  !            pocketMaster%conn(3, pocketMaster%nElems))
             
  !   ! Node data
  !   do i=1, nUnique
  !      ip = nodeList(i)
  !      pocketMaster%x(:, i) = p%x(:, ip)
  !      pocketMaster%norm(:, i) = p%norm(:, ip)
  !      pocketMaster%ind(i) = p%ind(ip)

  !      ! Save the original master node index
  !      pocketMaster%pNodes(i) = ip
  !   end do

  !   ! Element data
  !   do i=1, nPolyEdges
  !      pocketMaster%conn(1, i) = nodeMap(2*i-1) !<-- map to the unique node index
  !      pocketMaster%conn(2, i) = nodeMap(2*i) !<-- map to the unique node index
  !     print*,'conn(:, i) ',nodeMap(2*i-1), nodeMap(2*i)
  !   end do

  !   ! Debug pocketMaster
  !   ! --------------------
  !   pocketMaster%myID = 88
  !   open(unit=101, file="pocketMaster.dat", form='formatted')
  !   write(101,*) 'TITLE = "PocketMaster Data" '
  !   write(101,*) 'Variables = "X", "Y", "Z", "Nx", "Ny", "Nz", "ind" "gapID" "gapIndex" "otherID" "otherIndex"'
  !   call writeOversetString(pocketMaster, 101)
  !   close(101)
  !   ! --------------------

  !   ! Create nte info
  !   call createNodeToElem(pocketMaster)

  !   ! End 3.1 create pocketMaster elems and nodes

  !   ! ---------------------------------------------------------------------
  !   ! 3.2: perform doChain on pocketMaster elems and create pocket strings.
  !   ! ---------------------------------------------------------------------

  !   ! Create ordered pocketStrings based on connectivity using
  !   ! linked list.

  !   ! Some additional arrays.
  !   allocate(pocketMaster%elemUsed(pocketMaster%nElems), &
  !            pocketMaster%subStr(2, pocketMaster%nElems), &
  !            pocketMaster%cNodes(3, pocketMaster%nNodes)) ! Third index saves the size of the substrings
   
  !   ! Initialize
  !   pocketMaster%cNodes = 0
  !   pocketMaster%elemUsed = 0
  !   curElem = 1
  !   nPocketStrings = 0

  !   ! Do similar to creation to fullStrings
  !   do while (curElem < pocketMaster%nElems)

  !      ! First node of the first unused element
  !      iStart = pocketMaster%conn(1, curElem)
  !      nElems = pocketMaster%nte(1, iStart)

  !      ! --------------------
  !      ! First side of chain:
  !      ! --------------------
  !      firstElem = pocketMaster%nte(2, iStart)
  !      pocketMaster%subStr(1, 1) = firstElem
  !      call doChain(pocketMaster, iStart, 1)

  !      ! ---------------------
  !      ! Second side of chain: 
  !      ! ---------------------
  !      ! Ideally, pocketStrings should be periodic. Should not require
  !      ! second side of chain. Do a sanity check anyway.

  !      if (nElems > 1) then
  !         firstElem = pocketMaster%nte(3, iStart)

  !         ! Check the second one is already end of the periodic chain
  !         ! done above.
  !         if (pocketMaster%elemUsed(firstElem) == 0) then
  !            print*, ' First side did not create a periodic chain'
  !            stop ' Error'
  !         end if
  !      end if

  !      ! Extract pocketStrings linked list from elements present in
  !      ! pocketMaster%subStr buffer

  !      ! Create or add a new string to 'pocketStrings' linked list
  !      if (nPocketStrings == 0) then
  !         allocate(pocketStrings) ! Create first linked list node
  !         nPocketStrings = 1
  !         pocketStrings%next => pocketStrings
  !         strPkt => pocketStrings ! Work with the first termporary linked list node
  !      else
  !         allocate(strPkt%next) ! Create new linked list node and link to it
  !         strPkt%next%next => pocketStrings ! point the end back to original pointer.
  !         strPkt => strPkt%next ! Work with the new temporary linked list node

  !         nPocketStrings = nPocketStrings + 1
  !      end if

  !      ! Create a pocketStrings substring from pocketMaster buffer
  !      call createSubStringFromElems(pocketMaster, strPkt, nPocketStrings)

  !      do while (pocketMaster%elemUsed(curElem) == 1 .and. &
  !                curElem < pocketMaster%nElems)
  !         curElem = curElem + 1
  !      end do
  !   end do ! main while loop

  !   ! Debug pocketStrings
  !   ! ---------------------------------------------
  !   print *, 'nPocketStrings:', nPocketStrings
  !   open(unit=101, file="pocketStrings.dat", form='formatted')
  !   write(101,*) 'TITLE = "PocketStrings Data" '
  !   write(101,*) 'Variables = "X", "Y", "Z", "Nx", "Ny", "Nz", "ind" "gapID" "gapIndex" "otherID" "otherIndex"'
  !   i = 0 
  !   strPkt => pocketStrings
  !   do while (i < nPocketStrings)
  !      i = i + 1
  !      call writeOversetString(strPkt, 101)
  !      strPkt => strPkt%next
  !   end do
  !   close(101)
  !   ! ---------------------------------------------
  !   ! ------------------------------------------------------------------------
  !   ! End 3.2, perform doChain on pocketMaster elems and create pocket strings.
  !   ! ------------------------------------------------------------------------

  !   ! --------------------------
  !   ! 3.3: selfZip pocketStrings
  !   ! --------------------------

  !   ! ======= Do actual pocketZip ================

  !   ! This is again going to be similar to selfZip of fullStrings

  !   pocketMaster%myID = 88

  !   ! Allocate space for pocket triangles.
  !   ! (n-sided polygon -> n-2 triangles)
  !   allocate(pocketMaster%tris(3, pocketMaster%nElems))
  !   pocketMaster%nTris = 0
   
  !   ! Build the pocketMaster tree
  !   pocketMaster%tree => kdtree2_create(pocketMaster%x, sort=.True.)

  !   ! Loop over pocketStrings and begin pocketZip starting 
  !   ! from smallest convex ear.
  !   strPkt => pocketStrings
  !   i = 0
  !   do while (i < nPocketStrings)
  !      i = i + 1
       
  !      pocketZiploop: do while (strPkt%nNodes > 2) 
  !         ! Each pass zips one triangle. Keep zipping
  !         ! until last triangle is zipped in the pocket polygon.
  !         call pocketZip(strPkt)
  !      end do pocketZiploop
       
  !      strPkt => strPkt%next 
  !   end do ! while loop

  !   call writeOversetTriangles(pocketMaster, "pocketTriangulation.dat")
  !   ! ------------------------------
  !   ! End 3.3, selfZip pocketStrings
  !   ! ------------------------------
  !   call computeTriSurfArea(pocketMaster, triArea)
  !   print*,'triArea ',triArea

  ! end subroutine makePocketZip

  ! !
  ! ! ======================================
  ! subroutine pocketZip(s)

  !   use overset
  !   use kdtree2_module
  !   implicit none

  !   ! Input parameters
  !   type(oversetString), intent(inout), target :: s

  !   ! Local variables
  !   integer(kind=intType) :: i, j, k, ii, im1, ip1, nalloc, idx, nFound, N
  !   integer(kind=intType) :: nNodes, nElems, elem1, elem2, imin
  !   logical :: lastNodeZipper, inTri, overlapFound
  !   real(kind=realType), dimension(3) :: v1, v2, norm, c
  !   real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm, costhetaMax

  !   integer(Kind=intType), dimension(:), allocatable :: nodeMap
  !   type(kdtree2_result), dimension(:), allocatable  :: results
  !   real(kind=realType), dimension(:, :), pointer :: xTmp, normTmp
  !   integer(kind=intType), dimension(:, :), pointer :: connTmp
  !   integer(kind=intType), dimension(:), pointer :: indTmp, pNodesTmp
 
  !   ! ----------------------------------
  !   ! pocketStrings 's' should be all peridic
  !   if (.not.s%isPeriodic) then
  !      print*,' Non-periodic pocketStrings ID: ',s%myID
  !      stop
  !   end if

  !   N = s%nNodes
  !   print *,'N ', N

  !   nAlloc = 25
  !   allocate(results(nAlloc))
  !   allocate(nodeMap(s%nNodes))
  !   nodeMap = 1

  !   ! Find min angled ear
  !   costhetaMax = -Large 
  !   do ii=1, N

  !      im1 = ii - 1
  !      ip1 = ii + 1

  !      ! Peroidic string at end...loop around
  !      if (ii == N) then 
  !         ip1 = 1
  !      else if (ii == 1) then
  !         im1 = s%nNodes
  !      end if

  !      ! Determine the angle between the vectors
  !      v1 = s%x(:, ip1) - s%x(:, ii)
  !      v2 = s%x(:, im1) - s%x(:, ii)
  !      v1nrm = norm2(v1)
  !      v2nrm = norm2(v2)
  !      call cross_prod(v2, v1, norm)
  !      norm = norm / norm2(norm)

  !      ! Interior node norm and potential node norm should point
  !      ! in the same direction.
  !      if (dot_product(norm, s%norm(:, ii)) > zero) then
          
  !         ! Dot product of im1 and ip1 nodes should be close
  !         if (dot_product(s%norm(:, ip1), s%norm(:, im1)) > 0.80) then

  !            costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)

  !            ! cos(theta) is the largest for smallest angle
  !            if (costhetaMax <= costheta) then
  !               costhetaMax = costheta
  !               imin = ii
  !            end if
  !         end if
  !      end if

  !   end do !ii 

  !   ! Zip about imin node
  !   ii = imin
    
  !   im1 = ii - 1
  !   ip1 = ii + 1

  !   ! Peroidic string at end...loop around
  !   if (ii == N) then 
  !      ip1 = 1
  !   else if (ii == 1) then
  !      im1 = s%nNodes
  !   end if
  !   print*,'min: im1, ii, ip1 ', im1, ii, ip1

  !   ! Check that this triangle does not contain any other 
  !   ! pocketStrings nodes

  !   c = half*(s%x(:, ip1) + s%x(:, im1))
  !   r2 = (c(1) - s%x(1, ii))**2 +  (c(2) - s%x(2, ii))**2 +  (c(3) - s%x(3, ii))**2
    
  !   r2 = max(r2, (s%x(1, ip1) - s%x(1, im1))**2 + (s%x(2, ip1) - s%x(2, im1))**2 + &
  !        (s%x(3, ip1) - s%x(3, im1))**2)

  !   nFound = 0
  !   outerLoop: do 
       
  !      call kdtree2_r_nearest(s%p%tree, c, r2, nfound, nalloc, results) 
  !      if (nFound < nAlloc) then 
  !         exit outerLoop
  !      end if
       
  !      ! Allocate more space and keep going
  !      deallocate(results)
  !      nAlloc = nAlloc * 2
  !      allocate(results(nAlloc))
  !   end do outerLoop

  !   ! We can now be sure that we have all the points inside our
  !   ! ball. Next we proceed to systematically check them. 
  !   overlapFound = .False.
  !   nodeFoundLoop: do k=1, nFound
  !      ! Note that we do check nodes from our own string,
  !      ! except for the the three nodes we're dealing
  !      ! with. Remember that we are working in our parent's
  !      ! ording here.
  !      idx = results(k)%idx 
       
  !      notPartofTriangle: if (idx /= s%pNodes(im1) .and. &
  !           idx /= s%pNodes(ii) .and. idx /= s%pNodes(ip1)) then 
          
  !         ! Only check if the node normal of the point we're
  !         ! checking is in the same direction as the triangle. 
  !         if (dot_product(s%norm(:, ii), s%p%norm(:, idx)) > zero) then 
             
  !            ! Finally do the actual trianlge test
  !            call pointInTriangle(s%x(:, ip1), s%x(:, ii), s%x(:, im1), &
  !                 s%p%x(:, idx), inTri)
  !            if (inTri) then 
  !               ! As soon as 1 is in the triangle, we know the
  !               ! gap string is no good. 
  !               overlapFound = .True. 
  !               exit nodeFoundLoop
  !            end if
  !         end if
  !      end if notPartofTriangle
  !   end do nodeFoundLoop

  !   if (.not. overlapFound) then 
  !      ! This triangle is good!
       
  !      s%p%nTris = s%p%nTris+ 1
  !      s%p%tris(:, s%p%nTris) = (/s%pNodes(ip1), s%pNodes(ii),s%pNodes(im1)/)
  !      nodeMap(ii) = 0

  !      ! The two shorted string edges have been used for selfZip
  !      ! 
  !      elem1 = s%p%nte(2, s%pNodes(ii))
  !      elem2 = s%p%nte(3, s%pNodes(ii))
  !      s%p%elemUsed(elem1) = 1
  !      s%p%elemUsed(elem2) = 1
  !   end if

  !   ! Modify the pocketStrings to remove the two elements and the node
  !   ! that got eliminated due to pocketZipping. 

  !   ! Save pointers to existing data
  !   nNodes = s%nNodes
  !   nElems = s%nElems
  !   xTmp => s%x
  !   normTmp => s%norm
  !   indTmp => s%ind
  !   connTmp => s%conn
  !   pNodesTmp => s%pNodes
    
  !   ! nodeMap(imin) = 0, and 1 for the rest nodes. Create new nodeMap
  !   ! by taking off node 'imin'.
  !   j = 0
  !   do i=1, s%nNodes
  !      if (nodeMap(i) == 1) then 
  !         j = j + 1
  !         nodeMap(i) = j
  !      end if
  !   end do
    
  !   ! Update the number of nodes/elems in our shorted chain. Every
  !   ! zipper reduces the number of nodes and number of elems by 1
  !   s%nNodes = s%nNodes - 1
  !   s%nElems = s%nElems - 1

  !   allocate(s%x(3, s%nNodes), s%norm(3, s%nNodes), s%ind(s%nNodes), &
  !            s%pNodes(s%nNodes), s%conn(2, s%nElems))

  !   do i=1, nNodes
  !      if (nodeMap(i) /= 0) then 
  !         s%x(:, nodeMap(i)) = xTmp(:, i)
  !         s%norm(:, nodeMap(i)) = normTmp(:, i)
  !         s%ind(nodeMap(i)) = indTmp(i)
  !         s%pNodes(nodeMap(i)) = pNodesTmp(i)

  !         ! Update string's parent's child node data
  !         s%p%cNodes(:, s%pNodes(nodeMap(i))) = (/s%myID, nodeMap(i), s%nNodes/)
  !      end if
  !   end do

  !   ! Since we know the string was in order, we can simply redo the connectivity
  !   do i=1, s%nElems
  !      s%conn(:, i) = (/i, i+1/)
  !   end do
  !   if (s%isPeriodic) then 
  !      s%conn(2, s%nElems) = 1
  !   end if

  !   ! Dellocate the existing memory
  !   deallocate(xTmp, normTmp, indTmp, connTmp, pNodesTmp)

  ! end subroutine pocketZip

  ! subroutine computeTriSurfArea(master, area)

  !   ! Computes area sum of all triangles belonging to object master
  !   use overset
  !   implicit none
 
  !   ! Input parameters
  !   type(oversetString), intent(in) :: master
  !   real(kind=realType), intent(out) :: area

  !   ! Local variables
  !   integer(kind=intType) :: i, n1, n2, n3
  !   real(kind=realType), dimension(3) :: v1, v2, norm

  !   area = 0.0
  !   do i=1, master%nTris
  !      n1 = master%tris(1, i)
  !      n2 = master%tris(2, i)
  !      n3 = master%tris(3, i)

  !      v1 = master%x(:, n2) - master%x(:, n1)
  !      v2 = master%x(:, n3) - master%x(:, n1)
  !      call cross_prod(v1, v2, norm)
  !      area = area + half*norm2(norm)
  !   end do

  ! end subroutine computeTriSurfArea
end module stringOps
