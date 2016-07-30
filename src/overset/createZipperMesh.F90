!
!     **********************************************************************
!     *                                                                    *
!     * createZipperMesh zips multiple overlapping surface meshes.         *
!     * First, it eliminates overlapping quads and then stiches the        *
!     * non-overlapping surface mesh boundaries with triangular            *
!     * surface grids.                                                     *
!     *                                                                    *
!     * In an overset framework the overlapping surface grids give         *
!     * wrong estimate of airloads. Zipper mesh overcomes this by          *
!     * creating a more accurate representation of the surface in the      *
!     * mesh overlap regions.                                              *
!     *                                                                    *
!     * ref:  "Enhancements to the Hybrid Mesh Approach to Surface Loads   *
!     *        Integration on Overset Structured Grids", William M. Chan   *
!     * http://people.nas.nasa.gov/~wchan/publications/AIAA-2009-3990.pdf  *
!     **********************************************************************
!

subroutine createZipperMesh(level, sps, oWallSendList, oWallRecvList, &
     nOwallSend, nOwallRecv, size1, size2, work, nWork)

  use communication
  use blockPointers
  use overset
  use BCTypes
  use inputTimeSpectral
  use wallDistanceData, only : xVolumeVec, IS1
  use stringops
  use inputOverset
  use adtapi
  use adjointvars
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
  logical :: isWallType
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
  call mpi_allreduce(tmpInt2D, bufSizes, 2*nDomTotal, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
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
     call mpi_isend(oWalls(iDom)%iBlank, oWalls(iDom)%maxCells, sumb_integer, &
          iproc, iDom, sumb_comm_world, sendRequests(sendCount), ierr)
     call ECHK(ierr, __FILE__, __LINE__)
  end do

  recvCount = 0
  iStart = 1
  do jj=1, noWallSend
     iProc = oWallSendList(1, jj)
     iDom = oWallSendList(2, jj)
     iSize = oWalls(iDom)%maxCells
     recvCount = recvCount + 1       

     call mpi_irecv(intRecvBuf(iStart), iSize, sumb_integer, &
          iProc, iDom, sumb_comm_world, recvRequests(recvCount), ierr)
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
     call writeWalls(oWalls, size(oWalls))
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
  call VecCreateMPI(SUMB_COMM_WORLD, size(nodeIndices), PETSC_DETERMINE, &
       localZipperNodes, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(localZipperNodes, localZipperTp, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(localZipperNodes, localZipperTv, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! This is a generic global nodal vector that can store the nodal
  ! values, pressure tractions or viscous tractions. 
  call VecCreateMPI(SUMB_COMM_WORLD, 3*nNodesLocal(1), PETSC_DETERMINE, &
       globalNodalVec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now create the general scatter that goes from the
  ! globalNodalVector to the local vectors. 
  call ISCreateGeneral(sumb_comm_world, size(nodeIndices), &
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

