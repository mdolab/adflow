module stringOps

  ! Import oversetString becuase every routine uses this.
  use oversetData, only : oversetString
  contains

  subroutine nullifyString(string)

    use constants
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
         string%tris, &
         string%surfCellID)

  end subroutine nullifyString

  subroutine deallocateString(string)

    use constants
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

    if (associated(string%surfCellID)) &
         deallocate(string%surfCellID)

    call nullifyString(string)

  end subroutine deallocateString

  subroutine setStringPointers(string)

    use constants
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

  subroutine createOrderedStrings(master, strings, nString)

    use constants
    implicit none

    ! Input/Output
    type(oversetString) :: master
    type(oversetString), dimension(:), allocatable :: strings

    ! Working
    integer(kind=intType) :: nElems, nNodes, curElem, nString, iStart,i, firstElem
    type(oversetString), pointer  :: stringsLL, str

    ! The next step is to create ordered strings based on the
    ! connectivity. This is a purely logical operation. We don't know
    ! how many actual strings we will need so we will use a linked
    ! list as we go.
    call createNodeToElem(master)

    ! Allocate some additional arrays we need for doing the chain
    ! searches.
    nElems = master%nElems
    nNodes = master%nNodes
    allocate(master%elemUsed(nElems), master%subStr(2, nElems), &
         master%cNodes(2, nNodes))

    master%cNodes = 0
    master%elemUsed = 0
    curElem = 1
    nString = 0
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
       if (nString == 0) then
          allocate(stringsLL)
          nString = 1
          stringsLL%next => stringsLL
          str => stringsLL
       else
          allocate(str%next)
          str%next%next => stringsLL
          str => str%next
          nString = nString + 1
       end if

       ! Create a substring from master based on the elements we
       ! have in the buffer
       call createSubStringFromElems(master, str, nString)

       ! Scan through until we find the next unused element:
       do while((master%elemUsed(curElem) == 1) .and. (curElem < master%nElems))
          curElem = curElem + 1
       end do
    end do

    ! Put the strings into an regular array which will be easier to
    ! manipulate.
    allocate(strings(nString))
    str => stringsLL
    i = 0
    do while (i < nString)
       i = i + 1
       strings(i) = str ! This is derived type assigment.
       call nullifyString(str)
       str => str%next
    end do

  end subroutine createOrderedStrings

  subroutine performSelfZip(master, strings, nStrings, debugZipper)

    use constants
    use kdtree2_module
    use inputOverset, only : selfZipCutoff
    implicit none

    ! Input/Output
    type(oversetString) :: master
    type(oversetString), dimension(nStrings), target ::  strings
    integer(kind=intType) :: nStrings
    logical, intent(in) :: debugZipper

    ! Workging
    type(oversetString), pointer :: str
    real(kind=realType) ::  cutOff
    integer(kind=intType) :: i, j, nZipped

    ! Now determine if there are any "holes" or periodic strings
    ! without anything inside of it. Ie closed loops. If there isn't
    ! we can self zip. Otherwise, we falg it so that it isn't touched
    ! and automatically pocket zipped at th end.

    do i=1, nStrings
       str => strings(i)
       str%isPocket = .True.
       do j=1, str%nNodes
          if (str%otherID(1, j) /= -1) then
             str%isPocket = .False.
          end if
       end do

       if (.not. str%isPocket) then
          zipperLoop: do j=1, 5
             if (j== 1) then
                cutOff = selfZipCutoff
             else
                cutOff = 90_realType
             end if
             call selfZip(strings(i), cutOff, nZipped)
             if (nZipped == 0) then
                exit zipperLoop
             end if
          end do zipperLoop
       end if
    end do

    ! Now we need to redo the string matching becuase the self-zip
    ! shortened the strings
    call stringMatch(strings, nStrings, debugZipper)

  end subroutine performSelfZip

  subroutine reduceGapString(string)

    ! Generic routine for removing duplicate nodes on the given
    ! string. The string is returned with the nodes and connectivities
    ! adjusted accordingly.

    use constants
    use utils, only : pointReduce, myNOrm2
    implicit none

    ! Input/Ouput
    type(oversetString), intent(inout) :: string

    ! Working:
    real(kind=realType) :: minEdge
    integer(kind=intType) :: nUnqiue, i, n1, n2, nUnique,idx
    integer(kind=intType), dimension(:), allocatable :: link
    real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
    real(kind=realType), dimension(:, :), pointer :: nodeDataPtr
    integer(kind=intType) , dimension(:, :), pointer :: intNodeDataPtr
    integer(kind=intType), dimension(:), allocatable :: normCounter
    real(kind=realType), dimension(:, :), allocatable :: uniqueNorms

    ! We will do a sort of adaptive tolernace here: Get the minium edge
    ! length and base the tolerance on that:

    minEdge = huge(1.0d0)

    do i=1, string%nElems
       n1 = string%conn(1, i)
       n2 = string%conn(2, i)
       minEdge = min(minEdge, mynorm2(string%x(:, n1) - string%x(:, n2)))
    end do

    allocate(link(string%nNodes), uniqueNodes(3, string%nNodes))

    call pointReduce(string%x, string%nNodes, minEdge/1000.0, uniqueNodes, link, nUnique)

    ! Now average the normals for any duplicate nodes. This is to handle any discrepancies
    ! for h-type mesh topologies where the surface block is fully represented by the volume connectivity
    allocate(normCounter(nUnique), uniqueNorms(3, nUnique))
    normCounter(:)=zero
    uniqueNorms(:,:) = zero
    ! sum the norms for the unique node and count how many duplicates there are for a given node
    do i = 1,string%nNodes
       idx = link(i)
       uniqueNorms(:,idx) =  uniqueNorms(:,idx)+ string%nodeData(4:6,i)
       normCounter(idx) = normCounter(idx)+1
    end do

    ! Now divide to get the average and assign back to original data storage
    do i = 1,string%nNodes
       idx = link(i)
       string%nodeData(4:6,i) = uniqueNorms(:,idx)/normCounter(idx)
    end do
    deallocate(normCounter, uniqueNorms)
    ! Averageing is complete

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

    use constants
    use utils, only : terminate
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
                call terminate("makeBoundaryString", "Inconsistent duplicate edge.")
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

    use constants
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

    use constants
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

    use constants
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

    use constants
    use kdtree2_module
    use utils, only : myNorm2, cross_prod
    implicit none

    ! Input/Output
    type(oversetString), intent(inout), target :: s
    integer(Kind=intType), intent(out) :: nZipped
    real(kind=realType), intent(in) :: cutOff

    ! Working
    integer(kind=intType) :: i, j, k,  N, ii, im1, ip1
    logical :: lastNodeZippered,  added
    real(kind=realType), dimension(3) :: v1, v2, norm
    real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm
    integer(Kind=intType), dimension(:), allocatable :: nodeMap
    type(kdtree2_result), dimension(:), allocatable  :: results

    ! Perform self zipping on the supplied string. The string at this
    ! point should be either peroidic or since sinded --- no multiple
    ! loops should be left. Therefore, we can count on the nodes being
    ! in order.

    allocate(results(25))
    allocate(nodeMap(s%nNodes))
    nodeMap = 1

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

    do while (ii <= N)

       ! Peroidic string at end...loop around
       if (s%isPeriodic .and. ii == N) then
          ip1 = 1
       end if

       lastNodeZippered = .False.

       ! Determine the anlge between the vectors
       v1 = s%x(:, ip1) - s%x(:, ii)
       v2 = s%x(:, im1) - s%x(:, ii)
       v1nrm = mynorm2(v1)
       v2nrm = mynorm2(v2)
       call cross_prod(v2, v1, norm)
       norm = norm / mynorm2(norm)

       if (dot_product(norm, s%norm(:, ii)) > zero) then


          ! the dot product of the im1 and ip1 nodes have to be close
          if (dot_product(s%norm(:, ip1), s%norm(:, im1)) > 0.80) then

             costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)

             if (costheta > cosCutoff) then

                call addPotentialTriangle(s, im1, ii, ip1, nodeMap, &
                     results, added)


                if (added) then
                   nZipped = nZipped + 1
                   lastNodeZippered = .True.
                end if
             end if
          end if
       end if

       if (lastNodeZippered) then
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

    call shortenString(s, nodeMap)
    deallocate(results, nodeMap)

  end subroutine selfZip

  subroutine crossZip(str1, N1, N2, str2, N3, N4, debugZipper, failed)

    use constants
    use utils, only : myNorm2, cross_prod

    implicit none

    type(oversetString), intent(inout) :: str1, str2
    integer(kind=intType) :: N1, N2, N3, N4
    logical :: debugZipper, failed
    ! Working
    type(oversetString), pointer :: p
    integer(kind=intType) :: stepsA, stepsB, nStepsA, nStepsB
    integer(kind=intType) :: nTriToAdd, ii, i, j, k, A, B, Ap, Bp
    integer(kind=intType) :: aPrev, bPrev
    real(kind=realType), dimension(3) :: ptA, ptB, ptAp, ptBp
    !real(kind=realType), dimension(3) :: ptAPrev, ptBPRev
    real(kind=realType), dimension(3) :: Aoff, Boff, ApOff, BpOff
    real(kind=realType), dimension(3) :: normA, normB, normAp, normBp
    real(kind=realType), dimension(3) :: perpA, perpB, perpAp, perpBp
    !real(kind=realType), dimension(3) :: normAPrev, normBPrev
    !real(kind=realType), dimension(3) :: perpAPrev, perpBPrev
    real(kind=realType), dimension(3) :: triNorm1, quadNorm1
    real(kind=realType), dimension(3) :: triNorm2, quadNorm2
    logical :: aValid, bValid, advanceA, aPreferred, area1, area2
    logical :: advanceB
    logical :: changeA, changeB
    logical :: aValidPrev, bValidPrev, advanceAPrev, advanceBPrev
    real(kind=realType) ::  sum1, sum2, h, dpa, dpb
    !am real(kind=realType), parameter :: cutOff = 0.95*3
    real(kind=realType), parameter :: cutOff = 0.85*3
    ! First determine the the total number of triangles we will add
    ! total. It is equal to the total number of triangles on each
    ! string. This will form the index on the do loop.
    failed = .False.
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

    ! Initialize these out of bounds incase something goes very wrong.
    APrev = -1
    BPrev = -1

    ! The number of steps we've performed in each edge
    stepsA = 0
    stepsB = 0

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
       triNorm1 = triNorm1 / mynorm2(triNorm1)

       ! Compute the sum of the dot product of the nodal norms with the triNorm
       sum1 = dot_product(triNorm1, normA) + dot_product(triNorm1, normB) + &
            dot_product(triNorm1, normAp)

       call cross_prod(ptB-ptA, ptBp-ptA, triNorm2)
       triNorm2 = triNorm2 / mynorm2(triNorm2)

       sum2 = dot_product(triNorm2, normA) + dot_product(triNorm2, normB) + &
            dot_product(triNorm2, normBp)

       ! Only use this to help pick one if both are still valid:
       if (aValid .and. bValid .and. dot_product(triNorm1, triNorm2) < 0.8) then

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
          ! Flag the nodes as used
          str1%xZipNOdeUsed(A) = 1
          str1%xZipNOdeUsed(Ap) = 1
          str2%xZipNOdeUsed(B) = 1

       else if (bValid .and. .not. aValid) then

          ! We have no choice but to take B+

          call addTri(A, str1, B, str2, Bp, str2)
          str1%xZipNOdeUsed(A) = 1
          str2%xZipNOdeUsed(B) = 1
          str2%xZipNOdeUsed(Bp) = 1

          advanceA = .False.
          advanceB = .True.

       else if (aValid .and. bValid) then

          ! We could take either. Use the preferred triangle.
          if (aPreferred) then

             call addTri(A, str1, B, str2, Ap, str1)

             str1%xZipNOdeUsed(A) = 1
             str1%xZipNOdeUsed(Ap) = 1
             str2%xZipNOdeUsed(B) = 1

             advanceA = .True.
             advanceB = .False.

          else

             call addTri(A, str1, B, str2, Bp, str2)
             str1%xZipNOdeUsed(A) = 1
             str2%xZipNOdeUsed(B) = 1
             str2%xZipNOdeUsed(Bp) = 1

             advanceA = .False.
             advanceB = .True.

          end if

       else

          ! Things are not looking good...but

          if (avalidPRev .and. bvalidPrev) then
             ! We might be able to save it! The last triangle we added
             ! was a choice..both were valid, but we picked one
             ! because it was preferred. Now we know the one we did
             ! pick screwed us for the next triangle...go back and
             ! pick the other one instead!

             ! First 'delete' the triangle by decrementing the tri
             ! counter and edge counters
             p => str1%p
             p%nTris = p%nTris - 1
             p%nEdges = p%nEdges - 3

             ! Now we determine which one was actually added and add
             ! the other one instead

             if (advanceAPrev) then
                ! We need to add the old B triangle instead, which
                ! means the A triangle we had added was bad
                aValidPrev = .False.
                stepsB = stepsB + 1
                stepsA = stepsA - 1

                call addTri(APrev, str1, B, str2, Bp, str2)

                ! Reset the 'A' data by shuffling backwards: The 'A'
                ! data is copied to 'Ap' and the 'A' data is restored from Aprev

                Ap = A
                ptAp= ptA
                normAp = normA
                perpAp = perpA

                A = Aprev
                ptA = str1%x(:, A)
                normA = str1%norm(:, A)
                perpA = str1%perpNorm(:, A)

                ! Increment the 'B' data since we actually used B

                B = Bp
                ptB = ptBp
                normB = normBp
                perpB = perpBp

                ! And get the new data for Bp
                Bp = nextNode(str2, B, .False.)
                ptBp = str2%x(:, Bp)
                normBp = str2%norm(:, Bp)
                perpBp = str2%perpNorm(:, Bp)

                ! We *actually* advanced B so..
                advanceBPrev = .True.
                advanceAPrev = .False.
             else

                ! We need to add the old A triangle, which means the B
                ! triangle we had added was bad
                bValidPrev = .False.
                stepsB = stepsB - 1
                stepsA = stepsA + 1
                call addTri(A, str1, Bprev, str2, Ap, str1)

                ! Reset the 'B' data by shuffling backwards: The 'B'
                ! data is copied to 'Bp' and the 'B' data is restored from Bprev

                Bp = B
                ptBp= ptB
                normBp = normB
                perpBp = perpB

                B = Bprev
                ptB = str2%x(:, B)
                normB = str2%norm(:, B)
                perpB = str2%perpNorm(:, B)

                ! Increment the 'A' data since we actually used A

                A = Ap
                ptA = ptAp
                normA = normAp
                perpA = perpAp

                ! And get the new data for Ap
                Ap = nextNode(str2, A, .True.)
                ptAp = str1%x(:, Ap)
                normAp = str1%norm(:, Ap)
                perpAp = str1%perpNorm(:, Ap)
                ! We *actually* advanced A so..
                advanceAPrev = .True.
                advanceBPrev = .False.

             end if

             ! We *don't* increment ii since this is in essence still
             ! the "last" iteration. We just cycle and try the current
             ! one again.
             if (debugZipper) then
                print *,'Saved cross zip from bad front.'
             end if
             cycle

          end  if

          advanceA = .False.
          advanceB = .False.

          ! Ewww. neither triangle is valid. Do not add the triangle
          ! just return and let the cross zip restart. This should
          ! skip over the bad area and the pocket zip can do the bad
          ! region.
          failed = .True.
          return

       end if

       ! Now we have to shuffle along the string.
       if (advanceA .and. .not.advanceB) then

          stepsA = stepsA + 1

          ! Save a copy of the previous A info
          APrev = A

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

          ! Save a copy of the previous B info incase we need it
          BPrev = B

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
       end if

       ! Save the prevoius valid triangles and what was advanced
       aValidPrev = aValid
       bValidPrev = bValid
       advanceAPrev = advanceA
       advanceBPrev = advanceB

       ! Finally increment the number of triangles we've used so far.
       ii = ii + 1
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

      vecA = vec1 / mynorm2(vec1)
      vecB = vec2 / mynorm2(vec2)

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

      use constants
      use utils, only : myNorm2, cross_prod
      implicit none

      ! Input/Output
      real(kind=realType), intent(in), dimension(3) :: pt1, pt2, pt3
      real(kind=realType) :: triArea

      ! Working
      real(kind=realType), dimension(3) :: norm

      call cross_prod(pt2-pt1, pt3-pt1, norm)
      triArea = half * mynorm2(norm)

    end function triArea

  end subroutine crossZip

  subroutine addTri(A, sA, B, sB, C, sC)

    ! Form a triangle from index 'A' on string 'sA' , index 'B' on
    ! string 'sB' and index 'C' on string 'sC'

    use constants
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

  end subroutine addTri

  subroutine makeCrossZip(p, strings, nStrings, debugZipper)

    use constants
    implicit none

    ! Input/output
    integer(kind=intType), intent(in) :: nStrings
    type(oversetString), intent(inout), target :: p, strings(nStrings)
    type(oversetString), pointer :: s, s1, s2
    logical, intent(in) :: debugZipper
    ! Working
    integer(kind=intType) :: i, iStart, iEnd, jStart, jEnd, iStart_j, iEnd_j
    integer(kind=intType) :: curOtherID, iString, ii, nextI, curIStart, nIElems_j
    integer(kind=intType) :: nIElemsBeg, nJElemsBeg, nElem1, nElem2
    integer(kind=intType) ::iStart_orig, iEnd_orig, jStart_orig, jEnd_orig
    logical :: fullLoop1, fullLoop2, dummy, failed
    ! The purpose of this routine is to determine the ranges on two
    !  paired strings that are continuously paired and suitable for
    !  performing cross zipping.

    ! Allocate arrays to keep track of nodes that have already been
    !  used in cross zipping.
    do i=1, nstrings
       s => strings(i)
       allocate(s%XzipNodeUsed(s%nNodes))
       s%xZipNodeUsed = 0
    end do

    strLoop: do iString=1,nStrings

       ! Skip strings that were pocekts
       if (strings(iString)%isPocket) then
          cycle
       end if

       ! S1 is the curent '1' string we are working with
       s1 => strings(iString)

       ! Find the lowest node number that isn't used:
       curIStart = startNode(s1)
       do while(curIStart > 0)

          if (debugZipper) then
             print *,'------------------------------------------------'
             print *,'Starting string ', s1%myid, 'at index ', curIstart
             print *,'------------------------------------------------'
          end if

          iStart = curIStart
          ! Other ID is the string attached at the current pt.
          curOtherID = s1%otherID(1, iStart)

          if (curOtherID == -1) then
             print *,'*************************************************************************'
             print *,'Error during makeCrossZip: Point ', iStart, 'does not have a matching point'
             print *,'Position: ', s1%x(:, iStart)
             print *,'*************************************************************************'
             stop
          end if

          ! S2 is the current '2' string we are working with
          s2 => strings(curOtherID)
          jStart = s1%otherID(2, iStart)

          ! ---------------- s1 increments -------------
          ! The goal is to increment s1 as far as we can go in the
          ! NEGATIVE direction.
          call traceMatch(s1, iStart, .False., curOtherID, iEnd, fullLoop1)

          if (.not. fullLoop1) then
             ! Now set iStart to iEnd. Basically we start right at the
             ! negative end the chain and traverse in the POSITIVE
             ! direction.
             iStart = iEnd
             call traceMatch(s1, iStart, .True., curOtherID, iEnd, dummy)
          end if

          ! Now, iStart -> iEnd (in the positive order) is the maximum
          ! possible extent that s1 could be connected to s1
          ! over. However, s2 may have something to say about that. We
          ! do the same operation for s2. Note that the orders are reversed.

          ! ---------------- s2 increments -------------
          call traceMatch(s2, jStart, .True., s1%myID, jEnd, fullLoop2)

          ! If the first jnode isnt' actually matched to me, like I am
          ! to him. Therefore skip me, and go to the next one.
          if (jStart == jEnd .and. .not. fullLoop2 .and. &
               s2%otherID(1, jStart) /= s1%myID) then
             s1%xZipNodeUsed(curIStart) = 1
             curIStart = startNode(s1)
             cycle
          end if

          if (.not. fullLoop2) then
             jStart = jEnd
             call traceMatch(s2, jStart, .False., s1%myID, jEnd, dummy)
          end if

          if ((iStart == iEnd .and. .not. fullLoop1)  .or.&
               (jStart == jEnd .and. .not. fullLoop2)) then
             ! Can't go anywhere. Flag this node and the next.

             s1%xZipNodeUsed(curIStart) = 1
             curIStart = startNode(s1)
             cycle
          end if

          if (debugZipper) then
             print *,'Initial Range s1:', istart, iend, fullLoop1
             print *,'Initial Range s2:', jstart, jend, fullLoop2
          end if

          ! Save the original start/endn locations
          iStart_orig = iStart
          iEnd_orig = iEnd
          jStart_orig = jStart
          jEnd_orig = jEnd

          if ((istart == iend .and. fullLoop1) .and. &
               (jstart == jend .and. fullLoop2)) then
             ! s1 fully attached to s2

             call closestSymmetricNode(s1, s2, istart, jstart)
             iEnd = iStart
             jEnd = jStart

          else if((iStart == iEnd .and. fullLoop1) .and. .not. fullLoop2) then

             ! Project jStart and jEnd onto s1
             iStart = s2%otherID(2, jStart)
             iEnd   = s2%otherID(2, jEnd)

          else if((jStart == jEnd .and. fullLoop2) .and. .not. fullLoop1) then

             ! Project iStart and iEnd onto s2
             jStart = s1%otherID(2, iStart)
             jEnd   = s1%otherID(2, iEnd)

          else

             ! part of s1 is attached to part of s2

             nIElemsBeg = elemsForRange(s1, iStart, iEnd, .True.)
             nJElemsBeg = elemsForRange(s2, jStart, jEnd, .False.)

             ! This the "projection" of the 'j' string on the 'i'
             ! string. Basically this is the range the 'j' string
             ! wants to "use up" on the i string.
             iStart_j = s2%otherID(2, jStart)
             iEnd_j   = s2%otherID(2, jEnd)

             ! These could match up with iStart and iEnd or they could
             ! not. That's what we need to determine here.

             ! Need to check if iStart_j is "larger" than istart. We
             ! just increment using nextNode to take care of periodic
             ! boundaries.
             if (iStart_j /= iStart) then
                ! The starting points are different. Increment iStart
                ! until we find
                i = iStart
                do ii=1, nIElemsBeg
                   i = nextNode(s1, i)
                   if (i == iStart_j) then
                      iStart = iStart_j
                      exit
                   end if
                end do
             end if

             if (iEnd_j /= iEnd) then
                ! The starting points are different. Decrement jEnd
                ! until we find
                i = iEnd
                do ii=1, nIElemsBeg
                   i = prevNode(s1, i)
                   if (i == iEnd_j) then
                      iEnd = iEnd_j
                      exit
                   end if
                end do
             end if

             ! Now with the updated range. Project the iRange back to the
             ! the final J range.
             jStart = s1%otherID(2, iStart)
             jEnd   = s1%otherID(2, iEnd)

          end if

          if (debugZipper) then
             print *,'Zipping string: ', s1%myid, ' with ', s2%myid
             print *,'s1 range:', istart, iend
             print *,'s2 range:', jstart, jend
          end if

          ! Before we do the zip, make sure the ranges have not
          ! degenerated to 0 elements
          nElem1 = 1 ! Doesn't matter, just not zero
          nElem2 = 1 ! Doesn't matter, just not zero
          if (.not. fullLoop1) then
             nElem1 = elemsForRange(s1, iStart, iEnd, .True.)
             if (iStart == iEnd) then
                nElem1 = 0
             end if

          end if
          if (.not. fullLoop2) then
             nElem2 = elemsForRange(s2, jStart, jEnd, .False.)
             if (jStart == jEnd) then
                nElem2 = 0
             end if
          end if

          if  (nElem1 > 0 .and. nElem2 > 0) then
             ! Do actual cross zip if we still have elements left on both strings
             call crossZip(s1, iStart, iEnd, s2, jStart, jEnd, debugZipper, failed)

             ! If we succefully cross zippered what we were suppoed to
             ! flag all the nodes from the original region as done.
             if (.not. failed) then
                call flagNodesUsed(s1, iStart_orig, iEnd_orig, .True.)
                call flagNodesUsed(s2, jStart_orig, jEnd_orig, .False.)
             else
                ! UhOh. We got stopped part way through. Flag just the
                ! nodes at the beginning that we didn't use. Leave
                ! those for the pocket. The nodes up to where theh
                ! cross zip stopped were flagged internally in cross zip.
                call flagNodesUsed(s1, iStart_orig, iStart, .True.)
                call flagNodesUsed(s2, jStart_orig, jStart, .False.)
             end if
          else

             ! Flag the full range of elements are consumed even
             ! though we didn't do the cross zip. Leave it for the
             ! pocket zipping.
             call flagNodesUsed(s1, iStart_orig, iEnd_orig, .True.)
             call flagNodesUsed(s2, jStart_orig, jEnd_orig, .False.)
          end if

          ! Find the next starting index:
          curIStart = startNode(s1)

       end do
    end do strLoop

  contains

    function startNode(s)
      ! Determine the lowest index of a non-used xzip node for
      ! string 's'.
      implicit none
      type(oversetString) :: s
      integer(kind=intType) :: startNode, i

      ! This will be the return value if all nodes are used:
      startNode = 0
      nodeLoop: do i=1, s%nNodes
         if (s%xZipNodeUsed(i) == 0) then
            startNode = i
            exit nodeLoop
         end if
      end do nodeLoop
    end function startNode

    function nextNode(s, i)

      implicit none
      type(oversetString), intent(iN) :: s
      integer(kind=intType), intent(in) :: i
      integer(kind=intType) :: nextNode

      ! Normally just increment:
      nextNode = i + 1

      if (i == s%nNodes) then
         if (s%isPeriodic) then
            nextNode = 1
         else
            ! Can't go any further
            nextNode = i
         end if
      end if

      ! If the next node is used. The next node is set the current
      ! one.
      if (s%xZipNodeUsed(nextNode) == 1) then
         nextNode = i
      end if
    end function nextNode

    function simpleNextNode(s, i)

      implicit none
      type(oversetString), intent(iN) :: s
      integer(kind=intType), intent(in) :: i
      integer(kind=intType) :: simpleNextNode

      ! Normally just increment:
      simpleNextNode = i + 1

      if (i == s%nNodes) then
         if (s%isPeriodic) then
            simpleNextNode = 1
         else
            ! Can't go any further
            simpleNextNode = i
         end if
      end if
    end function simpleNextNode

    function prevNode(s, i)

      implicit none
      type(oversetString), intent(iN) :: s
      integer(kind=intType), intent(in) :: i
      integer(kind=intTYpe) :: prevNode
      ! Normally just increment:
      prevNode = i - 1

      if (i == 1) then
         if (s%isPeriodic) then
            prevNode = s%nNodes
         else
            ! Can't go any further
            prevNode = i
         end if
      end if

      ! If the next node is used. The next node is set the current
      ! one.
      if (s%xZipNodeUsed(prevNode) == 1) then
         prevNode = i
      end if
    end function prevNode

    subroutine traceMatch(s, iStart, pos, checkID, iEnd, fullLoop)

      implicit none

      ! Given a starting position 'iStart' on string 's', traverse in
      ! the 'POSitive' or '.not. POSitive' direction checking that the
      ! otherID still matches "checkID". Return the ending position
      ! 'iEnd'.

      ! Input/Output
      type(oversetString) :: s
      integer(kind=intType), intent(in) :: iStart, checkID
      logical, intent(in) :: pos
      integer(kind=intType), intent(out) :: iEnd
      logical, intent(out) :: fullLoop

      ! Working
      integer(kind=intType) :: i, nextI

      i = iStart
      fullLoop = .False.

      traverseLoop: do
         if (pos) then
            nextI = nextNode(s, i)
         else
            nextI = prevNode(s, i)
         end if
         if (nextI == i .or. s%otherID(1, nextI) /= checkID) then
            ! We can't go any further than we already are
            iEnd = i
            exit traverseLoop
         end if

         ! Continue to the next one.
         i = nextI

         if (i == iStart) then
            fullLoop = .True.
            iEnd = i
            exit traverseLoop
         end if
      end do traverseLoop
    end subroutine traceMatch

    subroutine flagNodesUsed(s, N1, N2, pos)

      implicit none

      ! Input/Output
      type(oversetString) :: s
      integer(kind=intType), intent(in) :: N1, N2
      logical, intent(in) :: pos

      ! Working
      integer(kind=intType) :: nSteps, i, nextI

      if (pos) then
         if (N2 > N1) then
            nSteps = N2 - N1
         else if (N2 < N1) then
            nSteps = N2 + s%nNodes - N1
         else ! N1 == N2
            nSteps = s%nElems
         end if
      else
         if (N1 < N2) then
            nSteps = N1 + s%nNodes - N2
         else if (N1 > N2) then
            nSteps = N1 - N2
         else ! N3 == N4
            nSteps =  s%nElems
         end if
      end if

      s%xZipNodeUsed(N1) = 1
      i = N1
      do ii=1, nSteps
         if (pos) then
            nextI = nextNode(s, i)
         else
            nextI = prevNode(s, i)
         end if

         s%xZipNodeUsed(nextI) = 1
         i = nextI
      end do
    end subroutine flagNodesUsed

    function elemsForRange(s, N1, N2, pos)
      ! Determine the number of elements between N1 and N2 for for the
      ! "POSitive" or "not POSIitive" (negative) direction.

      implicit none
      type(oversetString) :: s
      integer(kind=intType), intent(in) :: N1, N2
      logical :: pos
      integer(kind=intType) :: elemsForRange

      if (.not. s%isPeriodic) then
         if (pos) then
            elemsForRange = N2 - N1
         else
            elemsForRange = N1 - N2
         end if
      else ! Periodic
         if (pos) then
            if (N2 == N1) then
               elemsForRange = s%nElems
            else if (N2 > N1) then
               elemsForRange = N2 - N1
            else
               elemsForRange = N2 + s%nNodes - N1
            end if
         else
            if (N1 == N2) then
               elemsForRange = s%nElems
            else if (N1 > N2) then
               elemsForRange = N1 - N2
            else
               elemsForRange = N1 + s%nNodes - N2
            end if
         end if
      end if
    end function elemsForRange

  end subroutine makeCrossZip

  subroutine makePocketZip(p, strings, nStrings, pocketMaster, debugZipper)
    use constants
    use oversetData, only : oversetString, oversetEdge
    use oversetUtilities, only : qsortEdgeType
    use kdtree2_module
    implicit none

    ! Input/output
    integer(kind=intType), intent(in) :: nStrings
    type(oversetString), intent(in) :: p, strings(nStrings)
    type(oversetString) :: pocketMaster
    logical, intent(in) :: debugZipper

    ! Local variables
    integer(kind=intType) :: i, j, nsum1, nsum2, ndiff1, ndiff2, ipedge, icur
    integer(kind=intType) :: n1, n2, npolyEdges
    integer(kind=intType) :: nNodes1, nNodes2, cn1, cn2, str1, str2
    type(oversetEdge), allocatable, dimension(:) :: polyEdges
    type(oversetEdge) :: e1, e2
    type(oversetString), pointer  :: stringsLL, str
    integer(kind=intType) :: npocketEdges, nFullStrings, nNodes
    integer(kind=intType) :: ip, curElem, nElems, iStart, firstElem
    type(oversetString), allocatable, dimension(:), target :: pocketStringsArr

    ! ---------------------------------------------------------------
    ! PocketZip 1:
    ! First sort the edges.
    ! ---------------------------------------------------------------
    call qsortEdgeType(p%Edges, p%nEdges)

    ! Now gather up the left-over edges for pocket zipping.

    ! Over estimate of remaining pocket edges to zip
    allocate(polyEdges(p%nEdges))

    ! Eliminate the edges going through the ordered edges.
    ! The sorted opposite edges are canceled in pairs.
    npolyEdges = 0
    i = 1
    do while (i <= p%nEdges)

       if (i == p%nEdges) then
          ! This must be the last free edge:
          e1 = p%Edges(i)
          npolyEdges = npolyEdges + 1
          polyEdges(npolyEdges) = e1
          i = i + 1
          cycle
       end if

       ! Two edges in sequence
       e1 = p%Edges(i)
       e2 = p%Edges(i+1)

       ! First determine if e1 is at the end of two single ended
       ! chains. In this case the edge *will* not be paired and that's
       ! correct.

       str1    = p%cNodes(1, e1%n1) ! node1's child fullStrings ID
       cn1     = p%cNodes(2, e1%n1) ! node1's child fullStrings node index
       nNodes1 = strings(str1)%nNodes ! node1's child fullStrings nNodes size

       str2    = p%cNodes(1, e1%n2) ! node2's child fullStrings ID
       cn2     = p%cNodes(2, e1%n2) ! node2's child fullStrings node index
       nNodes2 = strings(str2)%nNodes ! node2's child fullStrings nNodes size

       if (str1 /= str2 ) then
          if (.not.strings(str1)%isperiodic .and. &
               .not.strings(str2)%isperiodic .and. &
               (cn1==1 .or. cn1==nNodes1) .and. (cn2==1 .or. cn2==nNodes2)) then
             ! Increment just 1 in 1 to skip over edge e1.
             i = i + 1
             cycle
          end if
       end if

       ! The sum and difference:
       nsum1 = e1%n1 + e1%n2
       nsum2 = e2%n1 + e2%n2

       ndiff1 = e1%n2 - e1%n1
       ndiff2 = e2%n2 - e2%n1

       if (nsum1 == nsum2 .and. ndiff1 + ndiff2 == 0) then
          ! These edges cancel. Great.
          i = i + 2
          cycle
       else
          ! Add just the first edge
          npolyEdges = npolyEdges + 1
          polyEdges(npolyEdges) = e1
          i = i + 1
       end if
    end do

    ! Define pocketMaster string
    call nullifyString(pocketMaster)
    pocketMaster%myID = 88
    pocketMaster%nElems = nPolyEdges
    pocketMaster%nNodes = nPolyEdges*2
    pocketMaster%nEdges = 0
    allocate(pocketMaster%nodeData(10, 2*nPolyEdges), &
         pocketMaster%intNodeData(3, 2*nPolyEdges), &
         pocketMaster%conn(2, nPolyEdges))

    ! Dump the data into the pocketMaster
    do i=1, nPolyEdges
       pocketMaster%nodeData(:, 2*i-1) = p%nodeData(:, polyEdges(i)%n1)
       pocketMaster%intNodeData(:, 2*i-1) = p%intNodeData(:, polyEdges(i)%n1)

       pocketMaster%nodeData(:, 2*i) = p%nodeData(:, polyEdges(i)%n2)
       pocketMaster%intNodeData(:, 2*i) = p%intNodeData(:, polyEdges(i)%n2)
       pocketMaster%conn(:, i) = (/2*i, 2*i-1/)
    end do

    call setStringPointers(pocketMaster)
    call reduceGapString(pocketMaster)
    call createNodeToElem(pocketMaster)

    ! The next step is to create ordered strings based on the
    ! connectivity. This is a purely logical operation. We don't know
    ! how many actual strings we will need so we will use a linked
    ! list as we go.

    ! Allocate some additional arrays we need for doing the chain
    ! searches.
    nElems = pocketMaster%nElems
    nNodes = pocketMaster%nNodes
    allocate(pocketMaster%elemUsed(nElems), pocketMaster%subStr(2, nElems), &
         pocketMaster%cNodes(2, nNodes))

    pocketMaster%cNodes = 0
    pocketMaster%elemUsed = 0
    curElem = 1
    nFullStrings = 0

    do while (curElem < pocketMaster%nElems)

       ! Arbitrarily get the first node for my element:
       iStart = pocketMaster%conn(1, curElem)
       nElems = pocketMaster%nte(1, iStart)

       firstElem = pocketMaster%nte(2, iStart)
       pocketMaster%subStr(1, 1) = firstElem
       call doChain(pocketMaster, iStart, 1)

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
       call createSubStringFromElems(pocketMaster, str, nFullStrings)

       ! Scan through until we find the next unused element:
       do while((pocketMaster%elemUsed(curElem) == 1) .and. &
            (curElem < pocketMaster%nElems))
          curElem = curElem + 1
       end do
    end do

    ! Temporary strings array for plotting and pocketZipping
    allocate(pocketStringsArr(nFullStrings))
    str => stringsLL
    i = 0
    do while(i < nFullStrings)
       i = i + 1
       pocketStringsArr(i) = str ! Derived type assignment
       call nullifyString(str)
       str => str%next
    end do

    ! Allocate space for pocket triangles.
    ! (n-sided polygon -> n-2 triangles)
    allocate(pocketMaster%tris(3, 10*pocketMaster%nElems))
    allocate(pocketMaster%edges(4*pocketMaster%nElems))
    pocketMaster%nTris = 0

    ! Build the pocketMaster tree
    pocketMaster%tree => kdtree2_create(pocketMaster%x, sort=.True.)

    if (debugZipper) then
       open(unit=101, file="strings_pocket.dat", form='formatted')
       write(101,*) 'TITLE = "PocketStrings Data" '

       write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
            "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
       do i=1, nFullStrings
          ! Temporarily allocate otherID
          allocate(pocketStringsArr(i)%otherID(2, pocketStringsArr(i)%nNodes))
          pocketStringsArr(i)%otherID = -1

          call writeOversetString(pocketStringsArr(i), pocketStringsArr, &
               nFullStrings, 101)
       end do
       close(101)
    end if

    ! Loop over pocketStrings and begin pocketZip starting
    ! from smallest convex ear.
    do i=1,nFullStrings
       if (debugZipper) then
          print *,'Pocket Zipping String ', i, ' of ', nFullStrings
       end if
       pocketZiploop: do while (pocketStringsArr(i)%nNodes > 2)
          ! Each pass zips one triangle. Keep zipping
          ! until last triangle is zipped in the pocket polygon.
          call pocketZip(pocketStringsArr(i))
       end do pocketZiploop
    end do

    ! Destroy the strings array
    do i=1, nFullStrings
       call deallocateString(pocketStringsArr(i))
    end do
    deallocate(pocketStringsArr, polyEdges)

  end subroutine makePocketZip

  subroutine pocketZip(s)

    use constants
    use kdtree2_module
    use utils, only : mynorm2, cross_prod

    implicit none

    ! Input parameters
    type(oversetString), intent(inout), target :: s

    ! Local variables
    integer(kind=intType) :: i, j, k, ii, im1, ip1, N
    integer(kind=intType) :: nNodes, nElems, iimin
    real(kind=realType), dimension(3) :: v1, v2, norm, c
    real(kind=realType) :: cosCutoff, cosTheta, r2, v1nrm, v2nrm, costhetaMax
    real(kind=realType) :: dp, dpMax

    integer(Kind=intType), dimension(:), allocatable :: nodeMap, badNode
    type(kdtree2_result), dimension(:), allocatable  :: results
    logical :: added, iiMinSet
    real(kind=realType), parameter:: fact=0.95_realType
    N = s%nNodes
    allocate(results(25), nodeMap(N), badNode(N))
    nodeMap = 1
    badNode = 0 ! Will become 1 if bad
    outerZiploop: do

       ! No choice for the last triangle:
       if (N==3) then
          ii = 1
          im1 = prevNode(ii)
          ip1 = nextNode(ii)
          ! We don't call addPotentialTriangle because we don't have a
          ! choice anymore. Just call the raw addTri command
          call addTri(ip1, s, ii, s, im1, s)
          ! and flag the node as gone
          nodeMap(ii) = 0
          exit outerZipLoop
       end if

       iiMinSet = .False.

       ! First find the largest dot product:
       dpMax = -one
       nodeloop1: do ii=1, N

          if (badNode(ii) == 1) then
             cycle nodeLoop1
          end if

          ip1 = nextNode(ii)
          im1 = prevNode(ii)

          ! Determine the angle between the vectors
          v1 = s%x(:, im1) - s%x(:, ii)
          v2 = s%x(:, ip1) - s%x(:, ii)
          v1nrm = mynorm2(v1)
          v2nrm = mynorm2(v2)
          call cross_prod(v2, v1, norm)
          norm = norm / mynorm2(norm)
          dpMax = max(dpmax, dot_product(norm, s%norm(:, ii)))
       end do nodeloop1

       ! Next find the largest cosTheta that is winthin a factor
       !  of dpMax
       costhetaMax = -Large
       nodeloop2: do ii=1, N

          if (badNode(ii) == 1) then
             cycle nodeLoop2
          end if

          ip1 = nextNode(ii)
          im1 = prevNode(ii)

          ! Determine the angle between the vectors
          v1 = s%x(:, im1) - s%x(:, ii)
          v2 = s%x(:, ip1) - s%x(:, ii)
          v1nrm = mynorm2(v1)
          v2nrm = mynorm2(v2)
          call cross_prod(v2, v1, norm)
          norm = norm / mynorm2(norm)
          dp = dot_product(norm, s%norm(:, ii))
          if (dp > dpMax*fact) then ! We take this
             costheta = dot_product(v1, v2)  / (v1nrm * v2nrm)
             if (cosTheta > cosThetaMax) then
                costhetaMax = costheta
                iiMinSet = .True.
                iimin = ii
             end if
          end if
       end do nodeloop2

       if (iiMinSet) then
          ! Zip about node "iimin" if it was set:
          ii = iimin
          ip1 = nextNode(ii)
          im1 = prevNode(ii)
          call addPotentialTriangle(s, ip1, ii, im1, nodeMap, results, added)
          if (added) then
             ! This triangle was good!
             exit outerZipLoop
          else
             ! Bad node. Need to cycle through rest of pocket nodes.
             ! Remember this bad node in next cycle.
             badNode(ii) = 1
             cycle outerZiploop
          end if
       else
          ! What does this mean? We didn't find any node to zip. Are they all bad?
          print *,'Problem with pocket zipper. Somehow we were not able to find "&
               &"node to add a triangle on. This should not happen. Contact the "&
               &"Developers.'
          stop
       end if
    end do outerZiploop

    ! Modify the pocketStrings to remove the two elements and the node
    ! that got eliminated due to pocketZipping.
    call shortenString(s, nodeMap)
    deallocate(nodeMap, badNode, results)

  contains
    function nextNode(ii)
      implicit none
      integer(kind=intType) :: ii, nextNode
      nextNode = ii + 1
      if (ii == N) then
         nextNode = 1
      end if
    end function nextNode

    function prevNode(ii)
      implicit none
      integer(kind=intType) :: ii, prevNode
      prevNode = ii - 1
      if (ii == 1) then
         prevNode = N
      end if
    end function prevNode
  end subroutine pocketZip

  subroutine computeTriSurfArea(master, area)

    ! Computes area sum of all triangles belonging to object master
    use constants
    use utils, only : mynorm2, cross_prod
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
       area = area + half*mynorm2(norm)
    end do

  end subroutine computeTriSurfArea

  function triOverlap(pt1, pt2, pt3, str, i1, i2)

    use constants
    use utils, only : mynorm2, cross_prod
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
    triNorm = triNorm / mynorm2(triNorm)

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

  subroutine shortenString(s, nodeMap)

    ! This is an auxilary routine that take a string 's', and a node
    ! map of len s%nNodes, with 1 or 0. A 1 means that the node will
    ! be in the shortened string, 0 means that the node should be
    ! deleted.
    use constants
    implicit none

    ! Input/Output
    type(oversetString) :: s
    integer(kind=intType), dimension(:), intent(inout) :: nodeMap

    ! Working
    integer(kind=intType) :: nNodes, nElems, nRemoved, i, j
    real(kind=realType), dimension(:, :), pointer :: nodeDataTmp
    integer(kind=intType), dimension(:, :), pointer :: connTmp, intNodeDataTmp
    integer(kind=intType), dimension(:), pointer :: pNodesTmp


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
    nRemoved = 0
    do i=1, s%nNodes
       if (nodeMap(i) == 1) then
          j = j + 1
          nodeMap(i) = j
       else
          nRemoved = nRemoved + 1
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
    s%nNodes = s%nNodes - nRemoved
    s%nElems = s%nElems - nRemoved

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

    ! Recreate the node to elem
    if (s%nNodes >=3 ) then
       call createNodeToElem(s)
    end if

  end subroutine shortenString

  subroutine addPotentialTriangle(s, im1, ii, ip1, nodeMap, results, added)

    ! Common routine (for pocketZip and selfZip) to potentially add a
    ! triangle resulting from a single string.
    use constants
    use kdtree2_priority_queue_module
    use kdtree2_module
    implicit none

    ! Input/Output
    type(oversetString) :: s
    integer(kind=intType), intent(in) :: im1, ii, ip1
    integer(kind=intType), intent(inout), dimension(:) :: nodeMap
    type(kdtree2_result), dimension(:), allocatable :: results
    logical, intent(out) :: added
    ! Working:
    real(kind=realType) :: r2
    real(kind=realType), dimension(3) :: v1, v2, norm, c
    integer(kind=intType) ::  nFound, nalloc, idx, k, j, i
    logical :: overlapFound, inTri

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
    added = .False.

    c = half*(s%x(:, ip1) + s%x(:, im1))
    r2 = (c(1) - s%x(1, ii))**2 +  (c(2) - s%x(2, ii))**2 +  (c(3) - s%x(3, ii))**2

    r2 = max(r2, (s%x(1, ip1) - s%x(1, im1))**2 + (s%x(2, ip1) - s%x(2, im1))**2 + &
         (s%x(3, ip1) - s%x(3, im1))**2)

    nFound = 0
    outerLoop: do
       nalloc = size(results)
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
                ! triangle is no good.
                overlapFound = .True.
                exit nodeFoundLoop
             end if
          end if
       end if notPartofTriangle
    end do nodeFoundLoop

    if (.not. overlapFound) then

       ! This triangle is good!
       added = .True.

       ! Call the generic addTri Routine. Here all the ndoes from the
       ! triangle come from the same string.
       call addTri(ip1, s, ii, s, im1, s)

       ! Flag this node as gone
       nodeMap(ii) = 0

    end if
  end subroutine addPotentialTriangle

  subroutine closestSymmetricNode(s1, s2, i, j)
    use constants
    use utils, only : mynorm2

    implicit none

    ! Input/Output
    type(oversetString) :: s1, s2
    integer(kind=intType), intent(out) :: i, j
    real(kind=realType) :: minDist, dist
    integer(kind=intType) :: ii

    ! Working:
    minDist = large

    do ii=1, s1%nNodes
       ! "The other index of the matching node on the other string is
       ! me" ie. "I point to you and you point to me"

       if (s2%otherID(2, s1%otherID(2, ii)) == ii) then

          dist = mynorm2(s1%x(:, ii) - s2%x(:, s1%otherID(2, ii)))

          if (dist < minDist) then
             minDist = dist
             i = ii
             j = s1%otherID(2, ii)
          end if
       end if
    end do

  end subroutine closestSymmetricNode

  subroutine stringMatch(strings, nStrings, debugZipper)
    use constants
    use kdtree2_priority_queue_module
    use kdtree2_module
    use utils, only : mynorm2

    implicit none

    ! Input/output
    type(oversetString), dimension(nstrings), target :: strings
    integer(kind=intType), intent(in) :: nStrings
    logical, intent(in) :: debugZipper

    ! Working
    integer(kind=intType) :: i, j, k, idx, oid(4)
    integer(kind=intType) ::  nAlloc, nUnique, nSearch
    type(kdtree2_result), allocatable, dimension(:) :: results
    type(oversetString), pointer ::  str, master
    logical :: checkLeft, checkRight, concave
    logical :: checkLeft2, checkRight2, concave2
    logical :: leftOK, rightOK, isEndNode
    real(kind=realType), dimension(3) :: xj, xjp1, xjm1, normj
    real(kind=realType), dimension(3) :: xk, xkp1, xkm1, normk
    real(kind=realType), dimension(3) :: myPt, otherPt, eNorm
    real(kind=realType) ::  fact, dStar, curDist, minDist, edgeLength
    integer(kind=intTYpe) :: otherID, otherIndex, closestOtherIndex, closestOtherString
    integer(kind=intType) :: id, index
    real(kind=realType) :: timeA,  pt(3),   v(3), cosTheta,  cutOff, dist, maxH, ratio

    if (nStrings == 0) then
       return
    end if

    ! Now make we determine the nearest point on another substring
    ! for each point.
    nAlloc = 50
    allocate(results(nAlloc))
    master => strings(1)%p

    ! Loop over the fullStrings
    do i=1, nStrings
       str => strings(i) ! Easier readability

       ! No need to do anything with the pocket string.
       if (str%isPocket) then
          cycle
       end if

       ! Allocate space for otherID as it is not done yet
       if (associated(str%otherID)) then
          deallocate(str%otherID)
       end if

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
          isEndNode = .False.
          if (.not. (checkLeft .eqv. checkRight)) then
             ! Since we don't need to check one side, this means we're
             ! at the end of the chain. This is important since this
             ! node *MUST* be attached to another node on another
             ! chain at the end
             isEndNode = .True.
          end if

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
                ! Check 4d: If this is an end node, we need to check
                ! if the potential canditate is also a end node
                ! ---------------------------------------------
                if (isEndNode) then
                   if (checkRight2 .eqv. checkLeft2) then
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
                v = v/mynorm2(v)

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

       ! Do a sanity check to fix some extraordinary cases. If a node
       ! hasn't found a neighbouring string but each of the two nodes
       ! either side have, and they found the *same* string, just
       ! accept that.

       do j=3, str%nNodes-2
          if (str%otherID(1, j) == -1) then
             ! Bad node:
             oid(1) = str%otherID(1, j-2)
             oid(2) = str%otherID(1, j-1)
             oid(3) = str%otherID(1, j+1)
             oid(4) = str%otherID(1, j+2)

             if (oid(1) /= -1 .and. &
                  oid(1) == oid(2) .and. &
                  oid(1) == oid(3) .and. &
                  oid(1) == oid(4)) then

                if (debugZipper) then
                   print *,'****************************************************************'
                   print *,'Warning: Fixing a bad association on string ', i, 'at index', j
                   print *,'****************************************************************'
                end if

                ! We have a '-1' surrounded by the same gap string

                ! Set the stringID
                str%otherID(1, j) = oid(1)

                ! Estimate what the other index should be. Since this
                ! is in the middle of the string, the exact index
                ! shouldn't matter.
                str%otherID(2, j) = str%otherID(2, j-1)
             end if
          end if
       end do

    end do
  end subroutine stringMatch

  subroutine writeOversetString(str, strings, n, fileID)

    use constants
    use utils, only : mynorm2
    implicit none

    type(oversetString), intent(in) :: str
    type(oversetString), intent(in), dimension(n) :: strings
    integer(kind=intType), intent(in) :: fileID, n
    integer(kind=intType) :: i, j, id, index
    real(kind=realType), dimension(3) :: myPt, otherPT, vec
    real(kind=realType) :: maxH, dist, ratio

    character(80) :: zoneName


    write (zoneName,"(a,I5.5)") "Zone T=gap_", str%myID
    write (fileID, *) trim(zoneName)

    write (fileID,*) "Nodes = ", str%nNodes, " Elements= ", str%nElems, " ZONETYPE=FELINESEG"
    write(fileID, *) "DATAPACKING=BLOCK"
13  format (E20.12)

    ! Nodes
    do j=1,3
       do i=1, str%nNodes
          write(fileID,13) str%x(j, i)
       end do
    end do

    ! Node normal
    do j=1,3
       do i=1, str%nNodes
          write(fileID,13) str%norm(j, i)
       end do
    end do

    ! Vector between closest points
    do j=1,3
       do i=1, str%nNodes
          myPt = str%x(:, i)
          id = str%otherID(1, i)
          if (id /= -1) then
             index = str%otherID(2, i)
             otherPt = strings(id)%x(:, index)
             vec = otherPt - myPt
          else
             vec = zero
          end if

          write(fileID,13) vec(j)
       end do
    end do

    ! global node ID
    do i=1, str%nNodes
       write(fileID,13) real(str%ind(i))
    end do

    ! gapID
    do i=1, str%nNodes
       write(fileID,13) real(str%myID)
    end do

    ! gap Index
    do i=1, str%nNodes
       write(fileID,13) real(i)
    end do

    if (associated(str%otherID)) then
       ! otherID
       do i=1, str%nNodes
          write(fileID,13) real(str%otherID(1, i))
       end do

       ! other Index
       do i=1, str%nNodes
          write(fileID,13) real(str%otherID(2, i))
       end do
    else
       do i=1, 2*str%nNodes
          write(fileID,13) zero
       end do
    end if


    do i=1, str%nNodes
       myPt = str%x(:, i)
       id = str%otherID(1, i)
       if (id /= -1) then
          index = str%otherID(2, i)
          otherPt = strings(id)%x(:, index)
          dist = mynorm2(myPt - otherPt)
          maxH = max(str%h(i), strings(id)%h(index))
          ratio = dist/maxH
       else
          ratio = zero
       end if

       write(fileID,13) ratio
    end do

15  format(I5, I5)
    do i=1, str%nElems
       write(fileID, 15) str%conn(1, i), str%conn(2, i)
    end do

  end subroutine writeOversetString

  subroutine writeOversetMaster(str,fileID)

    use constants
    use utils, only : mynorm2
    implicit none

    type(oversetString), intent(in) :: str
    integer(kind=intType), intent(in) :: fileID
    integer(kind=intType) :: i, j, id, index
    real(kind=realType), dimension(3) :: myPt, otherPT, vec
    real(kind=realType) :: maxH, dist, ratio

    character(80) :: zoneName


    write (zoneName,"(a,I5.5)") "Zone T=gap_", str%myID
    write (fileID, *) trim(zoneName)

    write (fileID,*) "Nodes = ", str%nNodes, " Elements= ", str%nElems, " ZONETYPE=FELINESEG"
    write(fileID, *) "DATAPACKING=BLOCK"
13  format (E20.12)

    ! Nodes
    do j=1,3
       do i=1, str%nNodes
          write(fileID,13) str%x(j, i)
       end do
    end do

15  format(I5, I5)
    do i=1, str%nElems
       write(fileID, 15) str%conn(1, i), str%conn(2, i)
    end do

  end subroutine writeOversetMaster


  subroutine writeOversetTriangles(string, fileName, startTri, endTri)

    use constants
    implicit none

    type(oversetString), intent(inout) :: string
    integer(kind=intType), intent(in) :: startTri, endTri
    character(*) :: fileName
    integer(kind=intType) :: i, j
    character(80) :: zoneName

    open(unit=101, file=trim(fileName), form='formatted')
    write(101,*) 'TITLE = "Triangles"'
    write(101,*) 'Variables = "X", "Y", "Z"'

    write (zoneName,"(a,I5.5)") "Zone T=triangles_", string%myID
    write (101, *) trim(zoneName)

    write (101,*) "Nodes = ", string%nNodes, " Elements= ", (endTri-startTri+1), " ZONETYPE=FETRIANGLE"
    write (101,*) "DATAPACKING=POINT"
13  format (E20.12)

    ! Write all the coordinates
    do i=1, string%nNodes
       do j=1, 3
          write(101,13, advance='no') string%x(j, i)
       end do
       write(101,"(1x)")
    end do

15  format(I7, I7, I7)
    do i=startTri, endTri
       write(101, 15) string%tris(1, i), string%tris(2, i), string%tris(3, i)
    end do
    close(101)
  end subroutine writeOversetTriangles

  subroutine writeZipperDebug(str)

    ! Save the state of an unsplit string such that it can be debugged
    ! later without running overset interpolation.

    use constants
    implicit none

    type(oversetString) :: str
    integer(kind=intType) :: i, j

    open(unit=101, file="debug.zipper", form='formatted')
    write(101, *) str%nNodes
    write(101, *) str%nElems
    do i=1, str%nNodes
       do j=1, 10
          write (101,*) str%nodeData(j, i)
       end do
    end do

    do i=1, str%nElems
       do j=1, 2
          write (101,*) str%conn(j, i)
       end do
    end do

    do i=1, str%nNodes
       do j=1, 3
          write (101,*) str%intNodeData(j, i)
       end do
    end do
    close(101)
  end subroutine writeZipperDebug

  subroutine loadZipperDebug(fileName, str)

    ! Save the state of an unsplit string such that it can be debugged
    ! later without running overset interpolation.

    use constants
    implicit none

    character(*), intent(in) :: fileName
    type(oversetString) :: str
    integer(kind=intType) :: i, j

    open(unit=101, file=fileName, form='formatted')
    read(101, *) str%nNodes
    read(101, *) str%nElems
    call nullifyString(str)

    allocate(str%nodeData(10, str%nNodes))
    allocate(str%conn(2, str%nElems))
    allocate(str%intNodeData(3, str%nNodes))

    do i=1, str%nNodes
       do j=1, 10
          read (101,*) str%nodeData(j, i)
       end do
    end do

    do i=1, str%nElems
       do j=1, 2
          read (101,*) str%conn(j, i)
       end do
    end do

    do i=1, str%nNodes
       do j=1, 3
          read (101,*) str%intNodeData(j, i)
       end do
    end do
    close(101)

    call setStringPointers(str)

  end subroutine loadZipperDebug

  subroutine pointInTriangle(x1, x2, x3, pt, inTri)

    use constants
    use utils, only : cross_prod
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
    use utils, only : cross_prod
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

    use constants
    use utils, only : cross_prod
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
    logical ::  nodeInFrontOfEdges
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

    use constants
    use utils, only : mynorm2, cross_prod

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
    uNrm = mynorm2(u)
    u = u/uNrm

    call cross_prod(norm, u, v)
    v = v /mynorm2(v)

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

    use constants
    use utils, only : mynorm2, cross_prod
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
    uNrm = mynorm2(u)
    u = u/uNrm

    call cross_prod(norm, u, v)
    v = v /mynorm2(v)

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


end module stringOps
