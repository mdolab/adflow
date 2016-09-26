module zipperMesh
contains

  !
  !      createZipperMesh zips multiple overlapping surface meshes.         
  !      First, it eliminates overlapping quads and then stiches the        
  !      non-overlapping surface mesh boundaries with triangular            
  !      surface grids.                                                     
  !      In an overset framework the overlapping surface grids give         
  !      wrong estimate of airloads. Zipper mesh overcomes this by          
  !      creating a more accurate representation of the surface in the      
  !      mesh overlap regions.                                              
  !      ref:  "Enhancements to the Hybrid Mesh Approach to Surface Loads   
  !             Integration on Overset Structured Grids", William M. Chan   
  !      http://people.nas.nasa.gov/~wchan/publications/AIAA-2009-3990.pdf  
  !

  subroutine createZipperMesh(level, sps, oWallSendList, oWallRecvList, &
       nOwallSend, nOwallRecv, size1, size2, work, nWork)

    use constants
    use communication, only : myID, adflow_comm_world, nProc, recvRequests, &
         sendRequests, commPatternCell_2nd, internalCell_2nd
    use blockPointers
    use overset, only : oversetString, oversetWall
    use wallDistanceData, only : xVolumeVec, IS1
    use stringops
    use inputOverset
    use adtapi
    use utils, only : setPointers, EChk, isWallType
    use adjointvars, only :nNodesLocal
    use oversetUtilities
    use oversetCommUtilities, only : exchangeSurfaceIBlanks, recvOwall, sendOWall
    use oversetPackingRoutines
    use oversetInitialization
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level, sps, nOwallSend, nOwallRecv, nWork
    integer(kind=intType), intent(in) :: size1, size2
    integer(kind=intType), intent(in) :: oWallSendList(2, size1)
    integer(kind=intType), intent(in) :: oWallRecvList(2, size2)
    integer(kind=intType), intent(in) :: work(4, nWork)

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, jj, kk, iStart, iSize
    integer(kind=intType) :: iDom, jDom, nn, mm, n, ierr, iProc, iWork, nStrings
    integer(kind=intType), dimension(:,:), allocatable :: tmpInt2D
    logical, dimension(:), allocatable :: oWallReady
    type(oversetWall), dimension(:), allocatable :: oWalls
    integer(kind=intType), dimension(:), allocatable :: intRecvBuf
    type(oversetString) :: master, pocketMaster
    integer(kind=intType), dimension(:), allocatable :: nodeIndices

    type(oversetString), dimension(:), allocatable, target :: strings
    integer(kind=intType) :: nFullStrings, nUnique

    ! MPI/Communication related
    integer status(MPI_STATUS_SIZE) 
    integer(kind=intType), dimension(:, :), allocatable :: bufSizes
    integer(kind=intType), dimension(:, :), allocatable :: recvInfo
    integer(kind=intType) :: sendCount, recvCount, index

    ! Wall search related
    integer(kind=intType) :: ncells
    type(oversetWall), dimension(:), allocatable, target :: walls
    type(oversetWall),  target :: fullWall

    ! -------------------------------------------------------------------
    ! Step 1: Eliminate any gap overlaps between meshes
    ! -------------------------------------------------------------------

    ! Determine the average area of surfaces on each cluster. This
    ! will be independent of any block splitting distribution. 

    call determineClusterAreas()

    ! Set the boundary condition blank values
    call initBCDataIBlank(level, sps)

    ! Create the surface deltas
    call surfaceDeviation(level, sps)

    ! Alloc data for the OWalls
    nn = max(nProc, 2*nOWallSend, 2*nOwallRecv)
    allocate(tmpInt2D(nDomTotal, 2), bufSizes(nDomTotal, 2), &
         oWallReady(nDomTotal), recvInfo(2, nn))

    tmpInt2D = 0
    ! Need to get the differt sizes for the OWalls since they are now
    ! based on the primal mesh as opposed to do the dual mesh as
    ! previously done
    do nn=1, nDom
       call setPointers(nn, level, sps)
       iDom = cumDomProc(myid) + nn
       call getOWallBufferSizes  (il, jl, kl, tmpInt2D(iDom, 1), tmpInt2D(iDom, 2), .False.)
    end do

    ! Make sure everyone has the sizes
    call mpi_allreduce(tmpInt2D, bufSizes, 2*nDomTotal, adflow_integer, MPI_SUM, &
         adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    deallocate(tmpInt2D)

    ! allocate oWalls for the primal mesh
    allocate(oWalls(nDomTotal))
    do iDom=1, nDomtotal
       if (bufSizes(iDom, 1) == 0) then 
          oWallReady(iDom) = .True.
       else
          oWallReady(iDom) = .False.
       end if
    end do

    ! Initialize the primal walls
    do nn=1, nDom
       call setPointers(nn, level, sps)
       iDom = cumDomProc(myid) + nn
       call initializeOWall(oWalls(iDom), .False., clusters(iDom))
       call packOWall(oWalls(iDom), .False.)
       oWallReady(iDom) = .True. 
    end do

    ! Post all the oWall iSends
    sendCount = 0
    do jj=1, nOWallSend
       iProc = oWallSendList(1, jj)
       iDom = oWallSendList(2, jj)
       call sendOWall(oWalls(iDom), iDom, iProc, 0, sendCount)
    end do

    recvCount = 0
    do jj=1, nOWallRecv
       iProc = oWallRecvList(1, jj)
       iDom = oWallRecvList(2, jj)
       call recvOWall(oWalls(iDom), iDom, iProc, 0, &
            bufSizes(iDom, 1), bufSizes(iDom, 2), recvCount, recvInfo)
    end do

    ! Complete all the recives and sends
    do i=1, recvCount
       call mpi_waitany(recvCount, recvRequests, index, status, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    do i=1,sendCount
       call mpi_waitany(sendCount, sendRequests, index, status, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    ! Unpack any blocks we received if necessary:
    do i=1, recvCount

       ! Global domain index of the recv that finished
       iDom = recvInfo(1, i)
       if (.not. oWalls(iDom)%allocated) then 
          call unpackOWall(oWalls(iDom))
       end if
    end do

    ! Determine the size of the buffer we need locally for the
    ! receives. We do this outside the main iteration loop since it
    ! always the same size.
    ii = 0
    do jj=1, noWallSend
       ! These blocks are by definition local. 
       iDom = oWallSendList(2, jj)
       ii = ii + oWalls(iDom)%maxCells
    end do
    allocate(intRecvBuf(max(1, ii)))

    ! ------------------------ Performing Searches ----------------        
    do iWork=1, nWork
       iDom = work(1, iWork)
       jDom = work(2, iWork)
       call wallSearch(oWalls(iDom), oWalls(jDom))
    end do

    ! ------------------------ Receiving iBlank back ---------------

    sendCount = 0
    do jj=1, noWallRecv
       iProc = oWallRecvList(1, jj)
       iDom = oWallRecvList(2, jj)
       sendCount = sendCount + 1
       call mpi_isend(oWalls(iDom)%iBlank, oWalls(iDom)%maxCells, adflow_integer, &
            iproc, iDom, adflow_comm_world, sendRequests(sendCount), ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    recvCount = 0
    iStart = 1
    do jj=1, noWallSend
       iProc = oWallSendList(1, jj)
       iDom = oWallSendList(2, jj)
       iSize = oWalls(iDom)%maxCells
       recvCount = recvCount + 1       

       call mpi_irecv(intRecvBuf(iStart), iSize, adflow_integer, &
            iProc, iDom, adflow_comm_world, recvRequests(recvCount), ierr)
       call ECHK(ierr, __FILE__, __LINE__)
       iStart = iStart + iSize
    end do

    ! Now wait for the sends and receives to finish
    do i=1, sendCount
       call mpi_waitany(sendCount, sendRequests, index, status, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    do i=1, recvCount
       call mpi_waitany(recvCount, recvRequests, index, status, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
    end do

    ! Process the oWalls we own locally
    do nn=1, nDom
       call setPointers(nn, level, sps)
       iDom = cumDomProc(myid) + nn
       ii = 0
       do mm=1, nBocos
          if (isWallType(BCType(mm))) then 
             do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                   ii = ii +1
                   if (oWalls(iDom)%iBlank(ii) == -2) then 
                      BCData(mm)%iblank(i,j) = -2
                   end if
                end do
             end do
          end if
       end do
    end do

    ! And update based on the data we received from other processors
    ii = 0
    do kk=1, noWallSend

       iDom = oWallSendList(2, kk)
       nn = iDom - cumDomProc(myid)

       ! Set the block pointers for the local block we are dealing
       ! with:
       call setPointers(nn, level, sps)
       do mm=1, nBocos
          if (isWallType(BCType(mm))) then 
             do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                   ii = ii + 1
                   if (intRecvBuf(ii) == -2) then 
                      BCData(mm)%iBlank(i  , j) = -2
                   end if
                end do
             end do
          end if
       end do
    end do

    ! Ditch our owalls
    call deallocateOWalls(OWalls, nDomTotal)
    deallocate(oWalls, intRecvBuf)

    ! We are still left with an issue since we had only worked with
    ! the owned cells, we don't know if if a halo surface iblank was
    ! modified by another processor. The easiest way to deal with
    ! this is do just do a halo-surface-iblank exchange. Essentially
    ! we will put the surface iblankvalues into the volume iblank,
    ! communicate and the extract again. Easy-peasy
    call exchangeSurfaceIBlanks(level, sps, commPatternCell_2nd, internalCell_2nd)

    ! Before we continue, we do a little more
    ! processing. bowTieElimination tries to eliminate cells
    call bowTieAndIsolationElimination(level, sps)

    ! -------------------------------------------------------------------
    ! Step 2: Identify gap boundary strings and split the strings to 
    !         sub-strings.
    ! -------------------------------------------------------------------

    if (debugZipper) then 
       call writeWalls
    end if

    call makeGapBoundaryStrings(level, sps, master)

    if (myid == 0) then 

       if (debugZipper) then 
          call writeZipperDebug(master)
       end if

       call createOrderedStrings(master, strings, nStrings)
       call performSelfZip(master, strings, nStrings)

       ! Write out any self-zipped triangles
       if (debugZipper) then 
          call writeOversetTriangles(master, "selfzipTriangulation.dat")
       end if

       call makeCrossZip(master, strings, nStrings)

       if (debugZipper) then 
          call writeOversetTriangles(master, "fullTriangulation.dat")
       end if

       ! ---------------------------------------------------------------
       ! Sort through zipped triangle edges and the edges which have not
       ! been used twice (orphan edges) will be ultimately gathered to 
       ! form polygon pockets to be zipped.
       call makePocketZip(master, strings, nStrings, pocketMaster)

       if (debugZipper) then 
          call writeOversetTriangles(pocketMaster, "pocketTriangulation.dat")
       end if


       ! (3 nodes per triangle and 3 DOF per ndoe)
       allocate(nodeIndices(3*3*(master%ntris + pocketMaster%ntris)))

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
    else
       ! Other procs don't have any triangles :-(
       allocate(nodeIndices(0))
    end if

    ! Clean up memory on the root proc
    if (myid == 0) then 
       do i=1, nFullStrings
          call deallocateString(strings(i))
       end do
       deallocate(strings)
       call deallocateString(master)
       call deallocateString(pocketMaster)
    end if

    ! ------------------------------------------------

    ! These three vectors are the MPI vectors with only non-zero size on
    ! the root proc. This is where we scatter the nodes, pressure
    ! tractions and viscous tractions into. 
    call VecCreateMPI(ADFLOW_COMM_WORLD, size(nodeIndices), PETSC_DETERMINE, &
         localZipperNodes, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecDuplicate(localZipperNodes, localZipperTp, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecDuplicate(localZipperNodes, localZipperTv, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! This is a generic global nodal vector that can store the nodal
    ! values, pressure tractions or viscous tractions. 
    call VecCreateMPI(ADFLOW_COMM_WORLD, 3*nNodesLocal(1), PETSC_DETERMINE, &
         globalNodalVec, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Now create the general scatter that goes from the
    ! globalNodalVector to the local vectors. 
    call ISCreateGeneral(adflow_comm_world, size(nodeIndices), &
         nodeIndices, PETSC_COPY_VALUES, IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecScatterCreate(globalNodalVec, IS1, localZipperNodes, PETSC_NULL_OBJECT, &
         nodeZipperScatter, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call ISDestroy(IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Free the remaining memory
    deallocate(nodeIndices)

  end subroutine createZipperMesh
  !
  !       determineClusterArea determine the average cell surface area   
  !       for all blocks in a particular cluster. This is used for       
  !       determine blanking preference for overlapping surface cells.   

  subroutine determineClusterAreas

    use constants
    use blockPointers, only : nDom, BCData, nBocos, BCType
    use communication, only : adflow_comm_world, myid
    use overset, onlY : clusterAreas, nClusters, clusters, cumDomProc
    use surfaceFamilies, onlY : totalWallFamilies, wallFamilies
    use surfaceUtils, only : getSurfaceSize, getSurfacePoints, setFamilyInfo
    use utils, only : setPointers, EChk, isWallType
    implicit none

    ! Working
    integer(kind=intType) :: i, j, mm, nn, clusterID, ierr, nPts, nCells
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii
    integer(kind=intType) :: lower_left, lower_right, upper_left, upper_right
    real(kind=realType), dimension(:), allocatable :: localAreas
    integer(kind=intType), dimension(:), allocatable :: localCount, globalCount
    real(kind=realType), dimension(:, :), allocatable :: pts
    real(kind=realType) :: fact , v1(3), v2(3), sss(3), da

    if (allocated(clusterAreas)) then 
       ! We only ever do this once!
       return
    end if

    allocate(clusterAreas(nClusters), localAreas(nClusters), &
         localCount(nClusters), globalCount(nClusters))

    localAreas = zero
    localCount = 0

    call getSurfaceSize(nPts, nCells, wallFamilies, totalWallFamilies)
    call setFamilyInfo(wallfamilies, totalWallFamilies)
    allocate(pts(3, nPts))
    call getSurfacePoints(pts, nPts, 1_intType)

    ii = 0 
    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)

       clusterID = clusters(cumDomProc(myid) + nn)

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,nBocos
          if( isWalltype(BCType(mm))) then            

             ! Store the cell range of the subfaces a bit easier.
             ! As only owned faces must be considered the nodal range
             ! in BCData must be used to obtain this data.

             jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

             ! Compute the dual area at each node. Just store in first dof
             do j=jBeg, jEnd ! This is a face loop
                do i=iBeg, iEnd ! This is a face loop 

                   ! Compute Normal
                   lower_left  = ii + (j-jBeg)*(iEnd-iBeg+2) + i-iBeg + 1
                   lower_right = lower_left + 1
                   upper_left  = lower_right + iend - ibeg + 1
                   upper_right = upper_left + 1

                   v1(:) = pts(:, upper_right) -   pts(:, lower_left)
                   v2(:) = pts(:, upper_left) -  pts(:, lower_right)

                   ! Cross Product
                   call cross_prod(v1, v2, sss)
                   da = fourth*(sss(1)**2 + sss(2)**2 + sss(3)**2)
                   localAreas(clusterID) = localAreas(clusterID) + da
                   localCount(clusterID) = localCount(clusterID) + 1
                end do
             end do

             ! Note how iBeg,iBeg is defined above... it is one MORE
             ! then the starting node (used for looping over faces, not
             ! nodes)
             ii = ii + (jEnd-jBeg+2)*(iEnd-iBeg+2)

          end if
       end do bocos
    end do domains

    ! All reduce sum for the localAreas to get clusterAreas and
    ! localCount to get globalCount

    call mpi_allreduce(localAreas, clusterAreas, nClusters, adflow_real, &
         MPI_SUM, adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call mpi_allreduce(localCount, globalCount, nClusters, adflow_integer, &
         MPI_SUM, adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Final get the average global area
    do i=1, nClusters
       clusterAreas(i) = clusterAreas(i)/max(globalCount(i), 1)
    end do

    deallocate(localAreas, localCount, globalCount)
  end subroutine determineClusterAreas

  subroutine initBCDataiBlank(level, sps)

    use constants
    use blockPointers
    use communication
    use utils, only : setPointers, isWallType
    use oversetUtilities, only : flagForcedReceivers
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level, sps

    ! Local variables
    integer(kind=intType) :: mm, nn, i, j, k, iBeg, iEnd, jBeg, jEnd
    logical :: side(4)

    integer(kind=intType), dimension(:, :), pointer :: ibp, gcp, frx
    integer(kind=intType), dimension(:, :, :), pointer :: forcedRecv
    integer(kind=intType), dimension(:, :), allocatable :: toFlip

    ! This routine initializes the surface cell iblank based on the
    ! volume iblank. It is not a straight copy since we a little
    ! preprocessing 

    ! This is a little trickier than it seems. The reason is that we
    ! will allow a limited number of interpolated cells to be used
    ! directly in the integration provided the meet certain criteria.

    ! Consider the following
    !  +------+--------+-------+
    !  |ib=1  |  ib=1  | ib= 1 |
    !  |      |        |       |
    !  |      |        |       |
    !  +------+========+-------+
    !  |ib=1  || ib=-1 || ib=1 |
    !  |      ||       ||      |
    !  |      ||       ||      |
    !==+======+--------+=======+==
    !  |ib=-1 |  ib=-1 | ib=-1 |
    !  |      |        |       |
    !  |      |        |       |
    !  +------+--------+-------+
    !
    ! The boundary between real/interpolated cells is marked by double
    ! lines. For zipper mesh purposes, it is generally going to be
    ! better to treat the center cell, as a regular force integration
    ! cell (ie surface iblank=1). The criteria for selection of these
    ! cells is:

    ! 1. The cell must not have been a forced receiver (ie at overset outer
    ! bound)
    ! 2. Any pair *opposite* sides of the cell must be compute cells.

    ! This criterial allows one-cell wide 'slits' to be pre-eliminated.


    domainLoop: do nn=1, nDom
       call setPointers(nn, level, sps)
       allocate(forcedRecv(1:ie, 1:je, 1:ke))
       call flagForcedReceivers(forcedRecv)

       bocoLoop: do mm=1, nBocos
          wallType: if (isWallType(BCType(mm))) then
             select case (BCFaceID(mm))
             case (iMin)
                ibp => iblank(2, :, :)
                gcp => globalCell(2, :, :)
                frx => forcedRecv(2, :, :)
             case (iMax)
                ibp => iblank(il, :, :)
                gcp => globalCell(il, :, :)
                frx => forcedRecv(il, :, :)
             case (jMin)
                ibp => iblank(:, 2, :)
                gcp => globalCell(:, 2, :)
                frx => forcedRecv(:, 2, :)
             case (jMax)
                ibp => iblank(:, jl, :)
                gcp => globalCell(:, jl, :)
                frx => forcedRecv(:, jl, :)
             case (kMin)
                ibp => iblank(:, :, 2)
                gcp => globalCell(:, :, 2)
                frx => forcedRecv(:, :, 2)
             case (kMax)
                ibp => iblank(:, :, kl)
                gcp => globalCell(:, :, kl)
                frx => forcedRecv(:, :, kl)
             end select

             ! -------------------------------------------------
             ! Step 1: Set the (haloed) cell iBlanks directly from
             ! the volume iBlanks
             ! -------------------------------------------------
             jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

             ! Just set the cell iblank directly from the cell iblank
             ! above it. Remember the +1 in ibp is for the pointer
             ! offset. These ranges *ALWAYS* give 1 level of halos
             do j=jBeg-1, jEnd+1
                do i=iBeg-1, iEnd+1
                   BCData(mm)%iBlank(i,j) = ibp(i+1, j+1)
                end do
             end do

             ! Make bounds a little easier to read. Owned cells only
             ! from now on.

             ! -------------------------------------------------
             ! Step 2: Slit elimination
             ! -------------------------------------------------

             ! Now we loop back through the cells again. For
             ! interpolated cells with iblank=-1 we see if it satifies
             ! the criteria above. If so we flag it wil "toFlip" =
             ! 1. Note that we can't set a particular iblank directly
             ! since that could cause an "avalance" effect with the
             ! later iterations using the updated iblank from a previous
             ! iteration.
             allocate(toFlip(iBeg:iEnd, jBeg:jEnd))
             toFlip = 0
             do j=jBeg, jEnd
                do i=iBeg, iEnd

                   ! We *might* add it if the interpolated cell is
                   ! touching two real cell on opposite sides.
                   if (BCData(mm)%iBlank(i, j) == -1 .and. validCell(i, j)) then

                      ! Reset the side flag
                      side = .False.

                      if (validCell(i-1, j) .and. BCData(mm)%iBlank(i-1, j) == 1) then
                         side(1) = .True.
                      end if

                      if (validCell(i+1, j) .and. BCData(mm)%iBlank(i+1, j) == 1) then
                         side(2) = .True.
                      end if

                      if (validCell(i, j-1) .and. BCData(mm)%iBlank(i, j-1) ==1 ) then
                         side(3) = .True.
                      end if

                      if (validCell(i, j+1) .and. BCData(mm)%iBlank(i, j+1) == 1) then
                         side(4) = .True.
                      end if

                      if ((side(1) .and. side(2)) .or. (side(3) .and. side(4)))  then
                         toFlip(i,j) = 1
                      end if
                   end if
                end do
             end do

             ! Now just set the cell surface iblank to 1 if we
             ! determined above we need to flip  the cell

             do j=jBeg, jEnd
                do i=iBeg, iEnd
                   if (toFlip(i, j) == 1) then
                      BCData(mm)%iBlank(i, j) = 1
                   end if
                end do
             end do
             deallocate(toFlip)
          end if wallType
       end do bocoLoop
       deallocate(forcedRecv)
    end do domainLoop

  contains
    function validCell(i, j)
      implicit none
      integer(kind=intType), intent(in) :: i, j
      logical :: validCell

      ! for our purposes here, a valid cell is one that:
      ! 1. Is not a boundary halo. ie has globalCell >= 0
      ! 2. It is not a force receiver.

      validCell = .False.
      if (gcp(i+1, j+1) >= 0 .and. frx(i, j) == 0) then
         validCell = .True.
      end if
    end function validCell
  end subroutine initBCDataiBlank

  !
  !      surfaceDeviation computes an approximation of the maximum          
  !      deviation a surface could be as compared to an underlying "exact"  
  !      surface. The purpose is to compute an adaptive "near wall distance"
  !      value that can be used to determine if a point is "close" to a     
  !     * wall. 
  !

  subroutine surfaceDeviation(level, sps)

    use constants
    use blockPointers, only :BCdata, x, nBocos, nDom, BCType, il, jl, kl, BCFaceID
    use utils, only : setPointers, myNorm2, isWallType
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level, sps

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, jj, kk, nn, iBeg, iEnd, jBeg, jEnd, mm
    real(kind=realType) ::  deviation
    real(kind=realType), dimension(:, :, :), pointer :: xx

    ! Loop over blocks
    do nn=1, nDom
       call setPointers(nn, level, sps)

       bocoLoop: do mm=1, nBocos
          wallType: if (isWallType(BCType(mm))) then

             ! Extract pointers for the primal wall mesh

             select case (BCFaceID(mm))
             case (iMin)
                xx => x(1, :, :, :)
             case (iMax)
                xx => x(il, :, :, :)
             case (jMin)
                xx => x(:, 1, :, :)
             case (jMax)
                xx => x(:, jl, :, :)
             case (kMin)
                xx => x(:, :, 1, :)
             case (kMax)
                xx => x(:, :, kl, :)
             end select

             jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd  

             ! The procedure goes in 2 passes. The first pass checks all
             ! the i-direction edges, and the second all the j-direction
             ! edges. For every edge, we estimate the max deviation
             ! along that edge and then the surface will use the maximum
             ! deviation from each of the 4 edges. Ie we scatter the
             ! edge deviation to the two cells next to it. We only do
             ! the real cells here. Boundary halos get -one set (below)
             ! and then actual compute halos are set with an exchange. 

             bcData(mm)%delta = -one

             ! ------------------
             ! Check the i-edges
             ! ------------------
             do j=jBeg, jEnd       ! <------- Node loop
                do i=iBeg+1, iEnd  ! <------- Face Loop

                   ! We will creating a local cubic approximation of the
                   ! local edge. This will use node i-2, i-1, i, and
                   ! i+1. However, due to the pointer offset, these are
                   ! all shifted by 1 to get: i-1, i, i+1, i+2

                   deviation = checkDeviation(xx(i-1, j, :), xx(i, j, :), xx(i+1, j, :), &
                        xx(i+2, j, :))

                   ! Cell to the bottom:
                   if (j-1 >= jBeg+1) then 
                      bcData(mm)%delta(i, j-1) = max(bcData(mm)%delta(i, j-1), deviation)
                   end if

                   ! Cell to the top:
                   if (j+1 <= jEnd) then 
                      bcData(mm)%delta(i, j+1) = max(bcData(mm)%delta(i, j+1), deviation)
                   end if
                end do
             end do

             ! -----------------
             ! Check the j-edges
             ! -----------------
             do j=jBeg+1, jEnd   ! <------- Face loop
                do i=iBeg, iEnd  ! <------- Node Loop

                   ! We will creating a local cubic approximation of the
                   ! local edge. This will use node j-2, j-1, j, and
                   ! j+1. However, due to the pointer offset, these are
                   ! all shifted by 1 to get: j-1, j, j+1, j+2

                   deviation = checkDeviation(xx(i, j-1, :), xx(i, j, :), xx(i, j+1, :), &
                        xx(i, j+2, :))

                   ! Cell to the left:
                   if (i-1 >= iBeg+1) then 
                      bcData(mm)%delta(i-1, j) = max(bcData(mm)%delta(i-1, j), deviation)
                   end if

                   ! Cell to the right:
                   if (i+1 <= iEnd) then 
                      bcData(mm)%delta(i+1, j) = max(bcData(mm)%delta(i+1, j), deviation)
                   end if

                end do
             end do
          end if wallType
       end do bocoLoop
    end do

    ! Exchange so that halos get correct values set as well. THIS IS BROKEN FIX IT!
    !call exchangeSurfaceDelta(level, sps, commPatternCell_1st, internalCell_1st)

    ! Now make one pass back and compute a delta for the nodes. Of
    ! course, this technically makes no sense: The nodes should
    ! exactly
    ! using this as a surrogate for what is near a surface, it make a
    ! little sense. Essentially we go through the nodes, and take the
    ! max deviation from the cells surrpounding it. 

    ! Loop over blocks
    do nn=1, nDom
       call setPointers(nn, level, sps)

       bocoLoop2: do mm=1, nBocos
          wallType2: if (isWallType(BCType(mm))) then

             jBeg = BCdata(mm)%jnBeg; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg; iEnd = BCData(mm)%inEnd  

             do j=jBeg, jEnd
                do i=iBeg, iEnd

                   ! Since we are taking the max and the boundary halos
                   ! have a value of -one it's ok to blindy just take
                   ! the max from each of the 4 cells surrounding each node. 

                   BCData(mm)%deltaNode(i, j) = max(&
                        BCData(mm)%delta(i  , j  ), &
                        BCData(mm)%delta(i+1, j  ), &
                        BCData(mm)%delta(i  , j+1), &
                        BCData(mm)%delta(i+1, j+1))
                end do
             end do
          end if wallType2
       end do bocoLoop2
    end do

  end subroutine surfaceDeviation

  function checkDeviation(P0, P1, P2, P3)

    ! Find the maximum deviation between a local cubic approximation
    ! formed by nodes P0, P1, P2 and P3, with the linear approximation
    ! formed by nodes P1 and P2. 

    ! See this article for the implementation. 
    ! https://en.wikipedia.org/wiki/Centripetal_Catmull-Rom_spline

    use constants
    use utils, only : myNorm2
    implicit none

    ! Input Parameters
    real(kind=realType), intent(in), dimension(3) :: P0, P1, P2, P3

    ! Function value
    real(kind=realType) ::  checkDeviation

    ! Working Parameters
    real(kind=realType) :: t0, t1, t2, t3
    real(kind=realType), dimension(3) :: A1, A2, A3, B1, B2
    real(kind=realType), parameter :: alpha=half
    integer(kind=intType), parameter :: N=20
    integer(kind=intType) :: i
    real(kind=realType) :: t, P(3), Q(3), s

    t0 = zero
    t1 = t0 + mynorm2(P1-P0)**alpha
    t2 = t1 + mynorm2(P2-P1)**alpha
    t3 = t2 + mynorm2(P3-P2)**alpha

    ! Normalize
    t1 = t1/t3
    t2 = t2/t3
    t3 = one

    ! Loop over the number of points to check. We need to go between t2
    ! and t3. No need to check the first and last since the devaition
    ! there is zero by construction.
    checkDeviation = zero

    do i=1, N
       s = (i-one)/(N-one)
       t = (one-s)*t1 + s*t2


       ! Spline pt
       A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
       A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
       A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1

       B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
       B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2

       P = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2

       ! Now project the cubic point onto the line to get point Q

       Q = P1 + dot_product(P-P1, P2-P1)/dot_product(P2-P1, P2-P1) * (P2 - P1)

       ! Just get the distance between the two points.
       checkDeviation = max(checkDeviation, mynorm2(Q-P))

    end do

  end function checkDeviation

  subroutine writeWalls

    use communication
    use overset
    use constants
    use blockPointers
    use utils, only : setPointers
    implicit none

    character(80) :: fileName, zoneName
    integer(kind=intType) :: i, j, nn, iDom, iBeg, iEnd, jBeg, jEnd, mm, iDim
    real(kind=realType), dimension(:, :, :), pointer :: xx

    write (fileName,"(a,I5.5,a)") "wall_", myid, ".dat"

    open(unit=101,file=trim(fileName),form='formatted')
    write(101,*) 'TITLE = "mywalls"'
    write(101,*) 'Variables = "X", "Y", "Z", "CellIBlank"'

    do nn=1,nDom
       iDom = nn + cumDomProc(myid)
       call setPointers(nn, 1, 1)
       if (nBocos > 0) then 
          do mm=1, nBocos
             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
             if (BCType(mm) == EulerWall .or. &
                  BCType(mm) == NSWallAdiabatic .or. &
                  BCType(mm) == NSWallIsothermal) then
                select case (BCFaceID(mm))
                case (iMin)
                   xx => x(1,:,:,:)
                case (iMax)
                   xx => x(il,:,:,:)
                case (jMin)
                   xx => x(:,1,:,:)
                case (jMax)
                   xx => x(:,jl,:,:)
                case (kMin)
                   xx => x(:,:,1,:)
                case (kMax)
                   xx => x(:,:,kl,:)
                end select

                write(zoneName, "(a,I5.5,a,I5.5)") "Zone", iDom, "_Proc_", myid
110             format('ZONE T=',a, " I=", i5, " J=", i5)
                write(101, 110), trim(zoneName), iEnd-iBeg+1, jEnd-jBeg+1
                write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"

13              format (E14.6)
                do iDim=1,3
                   do j=jBeg, jEnd
                      do i=iBeg, iEnd
                         write(101, *) xx(i+1, j+1, iDim)
                      end do
                   end do
                end do

                do j=jBeg+1, jEnd
                   do i=iBeg+1, iEnd
                      write(101, *) BCData(mm)%iBlank(i, j)
                   end do
                end do
             end if
          end do
       else
          ! Write dummy zone
          write(zoneName, "(a,I1,a,I1)") "Zone", 0, "_Proc_", myid
          write(101, 110), trim(zoneName), 1, 1
          write (101,*) "DATAPACKING=POINT"
          write(101, *) zero, zero, zero, one, one
       end if
    end do
    close(101)
  end subroutine writeWalls

end module zipperMesh
