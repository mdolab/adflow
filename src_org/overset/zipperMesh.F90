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

  subroutine createZipperMesh(zipperFamList, nZipFam)

    use constants
    use communication, only : myID, adflow_comm_world, nProc, recvRequests, &
         sendRequests, commPatternCell_2nd, internalCell_2nd
    use blockPointers, only : nDom, BCData, nBocos, BCType, il, jl, kl
    use oversetData, only : oversetString, oversetWall, CSRMatrix, cumDomProc, nDomTotal, &
         clusters, overlapMatrix, &
         oversetPresent, zipperMesh, zipperMeshes
    use wallDistanceData, only : xVolumeVec, IS1, IS2
    use utils, only : setPointers, EChk
    use adjointvars, only :nNodesLocal
    use sorting, only : famInList
    use oversetUtilities, only : getWorkArray, deallocateOSurfs, transposeOverlap
    use oversetCommUtilities, only : exchangeSurfaceIBlanks, recvOSurf, sendOSurf, &
         getOSurfCommPattern
    use oversetPackingRoutines, only : packOSurf, unpackOSurf, getOSurfBufferSizes
    use oversetInitialization, only : initializeOSurf
    use inputOverset, only : debugZipper, useZipperMesh
    use surfaceFamilies, only : BCFamExchange, famNames, BCFamGroups
    use stringOps
    use gapBoundaries
    use wallSearches, only : wallSearch

#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in), dimension(nZipFam) :: zipperFamList
    integer(kind=intType), intent(in) :: nZipFam

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, jj, kk, iStart, iSize, sps, level, iStr
    integer(kind=intType) :: iDom, jDom, nn, mm, n, ierr, iProc, iWork, nodeFam(3)
    integer(kind=intType) :: nNodeTotal, nTriTotal, offset, iBCGroup, nFam
    integer(kind=intType), dimension(:), allocatable :: famList

    integer(kind=intType), dimension(:,:), allocatable :: tmpInt2D, work
    logical, dimension(:), allocatable :: oSurfReady
    type(oversetWall), dimension(:), allocatable :: oSurfs
    integer(kind=intType), dimension(:), allocatable :: intRecvBuf
    type(oversetString), target :: master, pocketMaster
    type(oversetString), pointer :: str
    integer(kind=intType), dimension(:), allocatable :: nodeIndices

    integer(kind=intType) :: nFullStrings, nUnique
    integer(kind=intType) :: nOSurfRecv, nOSurfSend
    integer(kind=intType) , dimension(:,:), allocatable :: oSurfRecvList, oSurfSendList
    ! MPI/Communication related
    integer mpiStatus(MPI_STATUS_SIZE)
    integer(kind=intType), dimension(:, :), allocatable :: bufSizes
    integer(kind=intType), dimension(:, :), allocatable :: recvInfo
    integer(kind=intType) :: sendCount, recvCount, index
    type(CSRMatrix), pointer :: overlap
    type(CSRMatrix) :: overlapTranspose
    type(zipperMesh), pointer :: zipper
    ! Wall search related
    integer(kind=intType) :: ncells
    type(oversetWall), dimension(:), allocatable, target :: walls
    type(oversetWall),  target :: fullWall

    if (.not. oversetPresent .or. (.not. useZipperMesh)) then
       ! Not overset so we don't can't have a zipper.
       return
    end if

    ! Zipper is only implemented for single grid and 1 spectral
    ! instance (ie not time spectral).
    level = 1
    sps = 1

    call initBCDataIblank(level, sps)

    ! We build the zipper meshes *independently* for each BCType.
    BCGroupLoop: do iBCGroup=1, nFamExchange

       ! Set a pointer to the zipper we are working on to make code
       ! easier to read
       zipper => zipperMeshes(iBCGroup)

       ! Deallocate if already exists
       if (zipper%allocated) then
          call VecScatterDestroy(zipper%scatter, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
          call VecDestroy(zipper%localVal, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
          zipper%allocated = .False.
       end if

       ! Note that the zipper%conn could be allocated with size of
       ! zero, but the vecScatter and Vec are not petsc-allocated.
       if (allocated(zipper%conn)) then
          deallocate(zipper%conn, zipper%fam, zipper%indices)
       end if

       ! Before we can proceed with the zipper, we need to generate
       ! intersection of the zipperFamList() with the families on this
       ! BCGroup. This would be so much easier in Python...

       if (allocated(famList)) then
          deallocate(famList)
       end if

       nFam = 0
       do i=1, size(BCFamGroups(iBCGroup)%famList)
          do j=1, size(zipperFamlist)
             if (BCFamGroups(iBCGroup)%famList(i) == zipperFamList(j)) then
                nFam = nFam + 1
             end if
          end do
       end do

       ! If nFam is zero, no need ot do anything for this zipper. Just
       ! allocated zero-sized arrays so we know the size is 0.
       if (nFam == 0) then
          allocate(zipper%conn(3, 0), zipper%fam(0), zipper%indices(0))

          cycle ! to the next BCGroup
       else
          ! Do a second pass and fill up the famList
          allocate(famList(nFam))
          nFam = 0
          do i=1, size(BCFamGroups(iBCGroup)%famList)
             do j=1, size(zipperFamlist)
                if (BCFamGroups(iBCGroup)%famList(i) == zipperFamList(j)) then
                   nFam = nFam + 1
                   famList(nFam) = zipperFamList(j)
                end if
             end do
          end do
       end if

       if (debugZipper .and. myid == 0) then
          write(*,"(a)",advance="no") '-> Creating zipper for families : '
          do i=1, size(famList)
             write(*,"(a,1x)",advance="no") trim(famNames(famList(i)))
          end do
          print "(1x)"
       end if

       overlap => overlapMatrix(level, sps)
       call transposeOverlap(overlap, overlapTranspose)
       ! -------------------------------------------------------------------
       ! Step 1: Eliminate any gap overlaps between meshes
       ! -------------------------------------------------------------------

       ! Determine the average area of surfaces on each cluster. This
       ! will be independent of any block splitting distribution.
       call determineClusterAreas(famList)

       ! Set the boundary condition blank values
       call slitElimination(famList, level, sps)

       ! Create the surface deltas
       call surfaceDeviation(famList, level, sps)

       ! ================ WE NORMALLY GOT THIS FROM OVERSETCOMM =============
       ! Get the OSurf buffer sizes becuase we need that for getOSurfCommPattern
       allocate(bufSizes(nDomTotal, 6), tmpInt2D(nDomTotal, 6))
       tmpInt2D = 0
       do nn=1, nDom
          call setPointers(nn, level, sps)
          iDom = cumDomProc(myid) + nn
          call getOSurfBufferSizes(famList, il, jl, kl, tmpInt2D(iDom, 5), &
               tmpInt2D(iDom, 6), .True.)
       end do

       call mpi_allreduce(tmpInt2D, bufSizes, 6*nDomTotal, adflow_integer, MPI_SUM, &
            adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)

       ! Get the basic surface comm pattern.
       ! For sending, the worse case is sending all my blocks/fringes/walls to
       ! everyone but myself:
       ii = nDom*(nProc-1)
       allocate(oSurfSendList(2, ii))

       ! For receiving, the worse receive is all the blocks/fringes/wall I
       ! don't already have:
       ii = nDomTotal - nDom
       allocate(oSurfRecvList(2, ii))

       call getOSurfCommPattern(overlap, overlapTranspose, &
            oSurfSendList, nOSurfSend, oSurfRecvList, nOSurfRecv, bufSizes(:, 6))
       deallocate(tmpInt2D, bufSizes)
       ! ========================================================================

       ! Alloc data for the OSurfs
       nn = max(nProc, 2*nOSurfSend, 2*nOSurfRecv)
       allocate(tmpInt2D(nDomTotal, 2), bufSizes(nDomTotal, 2), &
            oSurfReady(nDomTotal), recvInfo(2, nn))

       tmpInt2D = 0
       ! Need to get the differt sizes for the OSurfs since they are now
       ! based on the primal mesh as opposed to do the dual mesh as
       ! previously done
       do nn=1, nDom
          call setPointers(nn, level, sps)
          iDom = cumDomProc(myid) + nn
          call getOSurfBufferSizes(famList, il, jl, kl, tmpInt2D(iDom, 1), &
               tmpInt2D(iDom, 2), .False.)
       end do

       ! Make sure everyone has the sizes
       call mpi_allreduce(tmpInt2D, bufSizes, 2*nDomTotal, adflow_integer, MPI_SUM, &
            adflow_comm_world, ierr)
       call ECHK(ierr, __FILE__, __LINE__)
       deallocate(tmpInt2D)

       ! allocate oSurfs for the primal mesh
       allocate(oSurfs(nDomTotal))
       do iDom=1, nDomtotal
          if (bufSizes(iDom, 1) == 0) then
             oSurfReady(iDom) = .True.
          else
             oSurfReady(iDom) = .False.
          end if
       end do

       ! Initialize the primal walls
       do nn=1, nDom
          call setPointers(nn, level, sps)
          iDom = cumDomProc(myid) + nn
          call initializeOSurf(famList, oSurfs(iDom), .False., clusters(iDom))
          call packOSurf(famList, oSurfs(iDom), .False.)
          oSurfReady(iDom) = .True.
       end do

       ! Post all the oSurf iSends
       sendCount = 0
       do jj=1, nOSurfSend
          iProc = oSurfSendList(1, jj)
          iDom = oSurfSendList(2, jj)
          call sendOSurf(oSurfs(iDom), iDom, iProc, 0, sendCount)
       end do

       recvCount = 0
       do jj=1, nOSurfRecv
          iProc = oSurfRecvList(1, jj)
          iDom = oSurfRecvList(2, jj)
          call recvOSurf(oSurfs(iDom), iDom, iProc, 0, &
               bufSizes(iDom, 1), bufSizes(iDom, 2), recvCount, recvInfo)
       end do

       ! Complete all the recives and sends
       do i=1, recvCount
          call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
       end do

       do i=1,sendCount
          call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
       end do

       ! Unpack any blocks we received if necessary:
       do i=1, recvCount

          ! Global domain index of the recv that finished
          iDom = recvInfo(1, i)
          if (.not. oSurfs(iDom)%allocated) then
             call unpackOSurf(oSurfs(iDom))
          end if
       end do

       ! Determine the size of the buffer we need locally for the
       ! receives. We do this outside the main iteration loop since it
       ! always the same size.
       ii = 0
       do jj=1, noSurfSend
          ! These blocks are by definition local.
          iDom = oSurfSendList(2, jj)
          ii = ii + oSurfs(iDom)%maxCells
       end do
       allocate(intRecvBuf(max(1, ii)))

       ! ------------------------ Performing Searches ----------------
       call getWorkArray(overlap, work)

       do iWork=1, size(work,2)
          iDom = work(1, iWork)
          jDom = work(2, iWork)
          call wallSearch(oSurfs(iDom), oSurfs(jDom))
       end do

       ! ------------------------ Receiving iBlank back ---------------

       sendCount = 0
       do jj=1, noSurfRecv
          iProc = oSurfRecvList(1, jj)
          iDom = oSurfRecvList(2, jj)
          sendCount = sendCount + 1
          call mpi_isend(oSurfs(iDom)%iBlank, oSurfs(iDom)%maxCells, adflow_integer, &
               iproc, iDom, adflow_comm_world, sendRequests(sendCount), ierr)
          call ECHK(ierr, __FILE__, __LINE__)
       end do

       recvCount = 0
       iStart = 1
       do jj=1, noSurfSend
          iProc = oSurfSendList(1, jj)
          iDom = oSurfSendList(2, jj)
          iSize = oSurfs(iDom)%maxCells
          recvCount = recvCount + 1

          call mpi_irecv(intRecvBuf(iStart), iSize, adflow_integer, &
               iProc, iDom, adflow_comm_world, recvRequests(recvCount), ierr)
          call ECHK(ierr, __FILE__, __LINE__)
          iStart = iStart + iSize
       end do

       ! Now wait for the sends and receives to finish
       do i=1, sendCount
          call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
       end do

       do i=1, recvCount
          call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
          call ECHK(ierr, __FILE__, __LINE__)
       end do

       ! Process the oSurfs we own locally
       do nn=1, nDom
          call setPointers(nn, level, sps)
          iDom = cumDomProc(myid) + nn
          ii = 0
          do mm=1, nBocos
             if (famInList(BCData(mm)%famID, famList)) then
                do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                   do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                      ii = ii +1
                      if (oSurfs(iDom)%iBlank(ii) == -2) then
                         BCData(mm)%iblank(i,j) = -2
                      end if
                   end do
                end do
             end if
          end do
       end do

       ! And update based on the data we received from other processors
       ii = 0
       do kk=1, noSurfSend

          iDom = oSurfSendList(2, kk)
          nn = iDom - cumDomProc(myid)

          ! Set the block pointers for the local block we are dealing
           !with:
          call setPointers(nn, level, sps)
          do mm=1, nBocos
             if (famInList(BCData(mm)%famID, famList)) then
                do j=BCData(mm)%jnBeg+1, BCData(mm)%jnEnd
                   do i=BCData(mm)%inBeg+1, BCData(mm)%inEnd
                      ii = ii + 1
                      if (intRecvBuf(ii) == -2) then
                         BCData(mm)%iBlank(i  , j) = -2
                      end if
                   enddo
                end do
             end if
          end do
       end do
       ! Release some unnecessary memory
       deallocate(bufSizes, oSurfSendList, oSurfRecvList, oSurfReady, recvInfo, work)

       ! Ditch our oSurfs
       call deallocateOSurfs(OSurfs, nDomTotal)
       deallocate(oSurfs, intRecvBuf)

       ! Before we continue, we do a little more
       ! processing. bowTieElimination tries to eliminate cells
       call bowTieAndIsolationElimination(famList, level, sps)

       ! -------------------------------------------------------------------
       ! Step 2: Identify gap boundary strings and split the strings to
       !         sub-strings.
       ! -------------------------------------------------------------------

       if (debugZipper) then
          call writeWalls(famList)
       end if

       call makeGapBoundaryStrings(famList, level, sps, master)

       rootProc: if (myid == 0) then
          if (debugZipper) then
             call writeZipperDebug(master)
          end if

          ! Run the common core zipper routines
          call zipperCore(master, pocketMaster, debugZipper)

          ! Before we create the final data structures for the zipper
          ! mesh; we will combine the master and pocketMaster
          ! strings, and unique-ify the indices to help keep amount of data
          ! transfer to a minimum. We also determine at this point which
          ! triangle belongs to which family.

          call setStringPointers(master)
          nTriTotal = master%ntris + pocketMaster%ntris
          nNodeTotal = master%nNodes + pocketMaster%nNodes
          allocate(zipper%conn(3, nTriTotal), zipper%fam(nTriTotal), zipper%indices(nNodeTotal))

          ii = 0
          jj = 0
          outerLoop: do iStr=1,2
             ! Select which of the two we are dealing with
             if (iStr ==1) then
                str => master
                offset = 0
             else
                str => pocketMaster
                offset = master%nNodes
             end if

             ! Loop over the number of triangles.
             do i=1, str%nTris

                ! Extract the family from the nodes
                do j=1,3
                   nodeFam(j) = str%family(str%tris(j, i))
                end do

                ! Increment the running counter for all triangles.
                ii = ii + 1

                ! Set the family
                zipper%Fam(ii) = selectNodeFamily(nodeFam)

                ! Set the connectivity.
                zipper%conn(:, ii) = str%tris(:, i) + offset
             end do

             ! Loop over the nodal indices and add
             do j=1, str%nNodes
                ! And set the global indices that the zipper needs. Note
                ! that we are doing zipperIndices as a scalar and that
                ! the indices refer to the global nodes so are already in
                ! 0-based ordering.
                jj = jj + 1
                zipper%indices(jj) = str%ind(j)
             end do
          end do outerLoop

          call deallocateString(master)
          call deallocateString(pocketMaster)
       else
          ! Other procs don't have any triangles *sniffle* :-(
          allocate(zipper%conn(3, 0), zipper%fam(0), zipper%indices(0))
       end if rootProc

       call VecCreateMPI(ADFLOW_COMM_WORLD, size(zipper%indices), PETSC_DETERMINE, &
            zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Now create the general scatter that goes from the
       ! globalNodalVector to the local vectors.
       call ISCreateGeneral(adflow_comm_world, size(zipper%indices), &
            zipper%indices-1, PETSC_COPY_VALUES, IS1, ierr)
       call EChk(ierr,__FILE__,__LINE__)


#if PETSC_VERSION_GE(3,8,0)
       call VecScatterCreate(BCFamExchange(iBCGroup, sps)%nodeValLocal, IS1, &
            zipper%localVal, PETSC_NULL_IS, zipper%scatter, ierr)
#else
       call VecScatterCreate(BCFamExchange(iBCGroup, sps)%nodeValLocal, IS1, &
            zipper%localVal, PETSC_NULL_OBJECT, zipper%scatter, ierr)
#endif
       call EChk(ierr,__FILE__,__LINE__)

       call ISDestroy(IS1, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Flag to keep track of allocated PETSc object.
       zipper%allocated = .True.

    end do BCGroupLoop

  contains
    function selectNodeFamily(nodeFam)
      implicit none
      integer(kind=intType), dimension(3) :: nodeFam
      integer(kind=intType) :: selectNodeFamily

      ! We have a few cases to check:

      if (nodeFam(1) == nodeFam(2) .and. nodeFam(1) == nodeFam(3)) then
         ! Case 1: All of the nodes have the same family. This is what
         ! *should* happen all the time:

         selectNodeFamily = nodeFam(1)

      else if (nodeFam(1) == nodeFam(2)) then

         ! Case 2a: First two nodes are the same. Take the family from 1.

         selectNodeFamily = nodeFam(1)

      else if (nodeFam(2) == nodeFam(3)) then

         ! Case 2b: Last two nodes are the same. Take the family from 2.

         selectNodeFamily = nodeFam(2)

      else if (nodeFam(1) == nodeFam(3)) then

         ! Case 2b: First and last nodes are the same. Take the family from 1.

         selectNodeFamily = nodeFam(1)

      else

         ! All nodes are different. We arbitrarily take the first and
         ! print a warning becuase this should not happen.

         selectNodeFamily = nodeFam(1)
         print *,'Family for triangle could not be uniquely determined. Nodes are from 3 different families!'
      end if

    end function selectNodeFamily
  end subroutine createZipperMesh

  subroutine checkZipper(fileName)

    ! Special routine for checking zipper mesh loaded from debug file
    use constants
    use stringOps
    implicit none

    ! Input/Output
    character(*), intent(in) :: fileName

    ! Working
    type(oversetString) :: master, pocketMaster

    call loadZipperDebug(fileName, master)
    call zipperCore(master, pocketMaster, .True.)

    print *, 'Zipper successfully completed'
  end subroutine checkZipper

  subroutine zipperCore(master, pocketMaster, debugZipper)

    ! Common routine for creating the zipper from a given master
    use constants
    use stringOps
    use kdtree2_module
    implicit none

    ! Input/Output
    type(oversetString), intent(inout) :: master, pocketMaster
    logical, intent(in):: debugZipper

    ! Local
    type(oversetString), dimension(:), allocatable, target :: strings
    type(oversetString), pointer :: str
    integer(kind=intType) :: nStrings, i, j, nTriSelf

    if (debugZipper) then
       open(unit=101, file="master_beforeStrings.dat", form='formatted')
       write(101,*) 'TITLE = "Master Data" '
       write(101,*) 'Variables = "X" "Y" "Z"'
       call writeOversetMaster(master, 101)
       close(101)
    end if

    call createOrderedStrings(master, strings, nStrings)


    master%myID = 99

    ! Allocate space for the maximum number of directed edges. This
    ! is equal to the initial number of edges (nElems) plus 3 times
    ! the number of triangles we will add, which is also nElems. Now,
    ! we will probably not actualy have that many since we save a
    ! triangle and 3 edges for every self zip that is
    ! applied. Therefore we know this will always be enough
    allocate(master%edges(4*master%nElems))

    master%nEdges = 0

    do i=1, nStrings
       str => strings(i)
       do j=1, str%nElems
          master%nEdges = master%nEdges + 1
          master%edges(master%nEdges)%n1 = str%p%conn(1, str%pElems(j)) !<-- first node
          master%edges(master%nEdges)%n2 = str%p%conn(2, str%pElems(j)) !<-- second node
       end do
    end do

    ! Allocate space for the triangles. Again, this can be at most,
    ! nElems, but the total number of elements will most likely be
    ! smaller due to self zipping. If someone puts
    allocate(master%tris(3, master%nElems))
    master%nTris = 0

    ! Build the master tree
    master%tree => kdtree2_create(master%x, sort=.True.)

    ! Perform the string association:
    call stringMatch(strings, nStrings, debugZipper)

    if (debugZipper) then
       open(unit=101, file="strings_beforeSelfZip.dat", form='formatted')
       write(101,*) 'TITLE = "Gap Strings Data" '
       write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
            "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
       do i=1, nStrings
          call writeOversetString(strings(i), strings, nStrings, 101)
       end do
       close(101)
    end if

    call performSelfZip(master, strings, nStrings, debugZipper)

   ! Write out any self-zipped triangles
    nTriSelf = master%nTris

    if (debugZipper) then
       call writeOversetTriangles(master, "selfzipTriangulation.dat", 1, master%nTris)
    end if

    ! Write out the gaps AFTER the self zip
    if (debugZipper) then
       open(unit=101, file="strings_afterSelfZip.dat", form='formatted')
       write(101,*) 'TITLE = "Gap Strings Data" '
       write(101,*) 'Variables = "X" "Y" "Z" "Nx" "Ny" "Nz" "Vx" "Vy" "Vz" "ind" &
            "gapID" "gapIndex" "otherID" "otherIndex" "ratio"'
       do i=1, nStrings
          call writeOversetString(strings(i), strings, nStrings, 101)
       end do
       close(101)
    end if

    ! Now do the cross zip
    call makeCrossZip(master, strings, nStrings, debugZipper)

    ! And write out the triangle from the cross zip
    if (debugZipper) then
       call writeOversetTriangles(master, "crossZipTriangulation.dat", nTriSelf+1, master%nTris)
    end if

    ! ---------------------------------------------------------------
    ! Sort through zipped triangle edges and the edges which have not
    ! been used twice (orphan edges) will be ultimately gathered to
    ! form polygon pockets to be zipped.
    if (debugZipper) then
       print *,'Doing pocket zip'
    end if
    call makePocketZip(master, strings, nStrings, pocketMaster, debugZipper)

    if (debugZipper) then
       call writeOversetTriangles(pocketMaster, "pocketTriangulation.dat", 1, pocketMaster%nTris)
    end if

    ! Clean up the reminder of the sting memory on the root proc
    do i=1, nStrings
       call deallocateString(strings(i))
    end do
    deallocate(strings)

  end subroutine zipperCore

  !
  !       determineClusterArea determine the average cell surface area
  !       for all blocks in a particular cluster. This is used for
  !       determine blanking preference for overlapping surface cells.

  subroutine determineClusterAreas(famList)

    use constants
    use blockPointers, only : nDom, BCData, nBocos, BCType
    use communication, only : adflow_comm_world, myid
    use oversetData, onlY : clusterAreas, nClusters, clusters, cumDomProc
    use utils, only : setPointers, EChk, setBCPointers, cross_prod
    use BCPointers, only : xx
    use sorting, only : famInList
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in), dimension(:) :: famList

    ! Working
    integer(kind=intType) :: i, j, mm, nn, clusterID, ierr, nPts, nCells
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd
    real(kind=realType), dimension(:), allocatable :: localAreas
    integer(kind=intType), dimension(:), allocatable :: localCount, globalCount
    real(kind=realType), dimension(:, :), allocatable :: pts
    real(kind=realType) :: fact , v1(3), v2(3), sss(3), da

    if (allocated(clusterAreas)) then
       ! We only ever do this once!
       deallocate(clusterAreas)
    end if

    allocate(clusterAreas(nClusters), localAreas(nClusters), &
         localCount(nClusters), globalCount(nClusters))

    localAreas = zero
    localCount = 0

    domains: do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)

       clusterID = clusters(cumDomProc(myid) + nn)

       ! Loop over the number of boundary subfaces of this block.
       bocos: do mm=1,nBocos
          famInclude: if (famInList(BCData(mm)%famID, famList)) then
             ! Store the cell range of the subfaces a bit easier.
             ! As only owned faces must be considered the nodal range
             ! in BCData must be used to obtain this data.

             jBeg = BCData(mm)%jnBeg + 1; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg + 1; iEnd = BCData(mm)%inEnd

             call setBCPointers(mm, .True.)

             ! Compute the dual area at each node. Just store in first dof
             do j=jBeg, jEnd ! This is a face loop
                do i=iBeg, iEnd ! This is a face loop

                   v1(:) = xx(i+1, j+1, :) - xx(i,   j,  :)
                   v2(:) = xx(i  , j+1, :) - xx(i+1, j,  :)

                   ! Cross Product
                   call cross_prod(v1, v2, sss)
                   da = half*sqrt(sss(1)**2 + sss(2)**2 + sss(3)**2)
                   localAreas(clusterID) = localAreas(clusterID) + da
                   localCount(clusterID) = localCount(clusterID) + 1
                end do
             end do
          end if famInclude
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
    use utils, only : setPointers
    use oversetUtilities, only : flagForcedRecv
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in) :: level, sps

    ! Local variables
    integer(kind=intType) :: mm, nn, i, j, k, iBeg, iEnd, jBeg, jEnd
    logical :: side(4)

    integer(kind=intType), dimension(:, :), pointer :: ibp, gcp, frx
    integer(kind=intType), dimension(:, :), allocatable :: toFlip

    ! This routine initializes the surface cell iblank based on the
    ! volume iblank.

    domainLoop: do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Setting the surface IBlank array is done for *all* bocos.
       bocoLoop: do mm=1, nBocos
          select case (BCFaceID(mm))
          case (iMin)
             ibp => iblank(2, :, :)
          case (iMax)
             ibp => iblank(il, :, :)
          case (jMin)
             ibp => iblank(:, 2, :)
          case (jMax)
             ibp => iblank(:, jl, :)
          case (kMin)
             ibp => iblank(:, :, 2)
          case (kMax)
             ibp => iblank(:, :, kl)
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
       end do bocoLoop
    end do domainLoop
  end subroutine initBCDataiBlank


  subroutine slitElimination(famList, level, sps)

    use constants
    use blockPointers
    use communication
    use utils, only : setPointers
    use oversetUtilities, only : flagForcedRecv
    use sorting, only : famInList
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in), dimension(:) :: famList
    integer(kind=intType), intent(in) :: level, sps

    ! Local variables
    integer(kind=intType) :: mm, nn, i, j, k, iBeg, iEnd, jBeg, jEnd
    logical :: side(4)

    integer(kind=intType), dimension(:, :), pointer :: ibp, gcp, frx
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

    call flagForcedRecv()
    domainLoop: do nn=1, nDom
       call setPointers(nn, level, sps)

       ! Setting the surface IBlank array is done for *all* bocos.
       bocoLoop: do mm=1, nBocos
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

          ! Only do the slit elimination if we actually care about
          ! this surface for the zipper
          famInclude: if (famInList(BCData(mm)%famID, famList)) then

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
          end if famInclude
       end do bocoLoop
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
      if (gcp(i+1, j+1) >= 0 .and. frx(i+1, j+1) == 0) then
         validCell = .True.
      end if
    end function validCell
  end subroutine slitElimination


  !
  !      surfaceDeviation computes an approximation of the maximum
  !      deviation a surface could be as compared to an underlying "exact"
  !      surface. The purpose is to compute an adaptive "near wall distance"
  !      value that can be used to determine if a point is "close" to a
  !     * wall.
  !

  subroutine surfaceDeviation(famList, level, sps)

    use constants
    use blockPointers, only :BCdata, x, nBocos, nDom, BCType, il, jl, kl, BCFaceID
    use utils, only : setPointers, myNorm2, setBCPointers
    use sorting, only : famInList
    use BCPointers, only : xx
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in), dimension(:) :: famList
    integer(kind=intType), intent(in) :: level, sps

    ! Local Variables
    integer(kind=intType) :: i, j, k, ii, jj, kk, nn, iBeg, iEnd, jBeg, jEnd, mm
    real(kind=realType) ::  deviation

    ! Loop over blocks
    do nn=1, nDom
       call setPointers(nn, level, sps)

       bocoLoop: do mm=1, nBocos
          famInclude: if (famInList(BCData(mm)%famID, famList)) then

             call setBCPointers(mm, .True.)
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
          end if famInclude
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
          famInclude2: if (famInList(BCData(mm)%famID, famList)) then

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
          end if famInclude2
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

  subroutine writeWalls(famList)


    !use oversetData
    use constants
    use blockPointers
    use utils, only : setPointers, setBCPointers
    use BCPointers, only : xx
    use sorting, only : famInList
    use communication, only : myid, adflow_comm_world, nProc
    use utils, only : EChk
    implicit none
    integer(kind=intType), intent(in), dimension(:) :: famList
    character(80) :: fileName, zoneName
    integer(kind=intType) :: i, j, nn, iDom, iBeg, iEnd, jBeg, jEnd, mm, iDim
    integer(kind=intType) :: nNode, nCell, tag, ierr, ii, iProc, nLocalBoco, tmp(5)
    real(kind=realType), dimension(:), allocatable :: xBuffer
    integer(kind=intType), dimension(:), allocatable :: iblankBuffer, bocosPerProc
    integer(kind=intType), dimension(:,:,:), allocatable :: faceInfo
    integer, dimension(mpi_status_size) :: mpiStatus

    ! Write a gathered surface tecplot file.
    if (myid == 0) then
       write (fileName,"(a)") "zipper_wall.dat"

       open(unit=101,file=trim(fileName),form='formatted')
       write(101,*) 'TITLE = "zipper walls"'
       write(101,*) 'Variables = "X", "Y", "Z", "CellIBlank"'
    end if

    ! Before we start, the root processor needs to know how many
    ! receives we can expect from each processor.

    nLocalBoco = 0
    do nn=1,nDom
       call setPointers(nn, 1, 1)
       do mm=1, nBocos
          famInclude: if (famInList(BCData(mm)%famID, famList)) then
             nLocalBoco = nLocalBoco + 1
          end if famInclude
       end do
    end do

    if (myid == 0) then
       allocate(bocosPerProc(0:nProc-1))
    end if

    call mpi_gather(nLocalBoco, 1, adflow_integer, bocosPerProc, 1, &
         adflow_integer, 0, adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)
    ! Now setup the info array on the root proc:
    if (myid == 0) then
       allocate(faceInfo(5, maxval(bocosPerProc), 0:nProc-1))
    end if

    if (myid /= 0) then
       tag = 0
       do nn=1, nDom
          call setPointers(nn, 1, 1)
          do mm=1, nBocos
             jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
             famInclude2: if (famInList(BCData(mm)%famID, famList)) then
                tag = tag + 1
                tmp = (/iBeg, iEnd, jBeg, jEnd, nBkGlobal/)
                call mpi_send(tmp, 5, &
                     adflow_integer, 0, tag, adflow_comm_world, ierr)
                call ECHK(ierr, __FILE__, __LINE__)
             end if famInclude2
          end do
       end do
    else
       ! Receive the size info:
       do iProc=1, nProc-1
          do tag=1, bocosPerProc(iProc)
             call mpi_recv(tmp, 5, adflow_integer, iProc, tag, &
                  adflow_comm_world, mpiStatus, ierr)
             call ECHK(ierr, __FILE__, __LINE__)
             faceInfo(:, tag, iProc) = tmp
          end do
       end do
    end if

    ! Need a barrier between the two sets of the comms just in case.
    call MPI_Barrier(adflow_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    tag = 0
    do nn=1,nDom
       call setPointers(nn, 1, 1)
       do mm=1, nBocos
          jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
          iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
          famInclude3: if (famInList(BCData(mm)%famID, famList)) then
             call setBCPointers(mm, .True.)

             nNode = (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
             nCell = (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
             allocate(iBlankBuffer(nCell), xBuffer(3*nNode))
             ii = 0
             do j=jBeg+1, jEnd
                do i=iBeg+1, iEnd
                   ii = ii + 1
                   iBlankBuffer(ii) = BCData(mm)%iBlank(i, j)
                end do
             end do

             ii = 0
             do iDim=1,3
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      xBuffer(ii) = xx(i+1, j+1, iDim)
                   end do
                end do
             end do

             if (myid == 0) then
                ! We can write it directly:
                call writeZone(iBeg, iEnd, jBeg, jEnd, nBkGlobal, xBuffer, iBlankBuffer)
             else

                tag = tag + 1
                call mpi_send(iBlankBuffer, nCell, adflow_integer, 0, tag, &
                     adflow_comm_world, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                tag = tag + 1
                call mpi_send(xBuffer, nNode*3, adflow_real, 0, tag, &
                     adflow_comm_world, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

             end if
             deallocate(iBlankBuffer, xBuffer)
          end if famInclude3
       end do
    end do

    ! Complete the receives and writes on the root proc:
    if (myid == 0) then
       ! Receive the nodes and iblank info
       do iProc=1, nProc-1
          do tag=1, bocosPerProc(iProc)

             iBeg = faceInfo(1, tag, iProc)
             iEnd = faceInfo(2, tag, iProc)

             jBeg = faceInfo(3, tag, iProc)
             jEnd = faceInfo(4, tag, iProc)

             nBkGlobal =faceInfo(5, tag, iProc)

             nNode = (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
             nCell = (iEnd - iBeg + 1)*(jEnd - jBeg + 1)

             allocate(iBlankBuffer(nCell), xBuffer(3*nNode))

             call mpi_recv(iBlankBuffer, nCell, adflow_integer, iProc, (2*tag-1), &
                  adflow_comm_world, mpiStatus, ierr)
             call ECHK(ierr, __FILE__, __LINE__)

             call mpi_recv(xBuffer, nNode*3, adflow_real, iProc, 2*tag, &
                  adflow_comm_world, mpiStatus, ierr)
             call ECHK(ierr, __FILE__, __LINE__)

             call writeZone(iBeg, iEnd, jBeg, jEnd, nBkGlobal, xBuffer, iBlankBuffer)

             deallocate(iBlankBuffer, xBuffer)

          end do
       end do
       deallocate(faceInfo, bocosPerProc)
       close(101)
    end if

    contains

      subroutine writeZone(iBeg, iEnd, jBeg, jEnd, nBkGlobal, xx, iblank)

        implicit none
        ! Input
        integer(kind=intType), intent(in) :: iBeg, iEnd, jBeg, jEnd, nBkGlobal
        real(kind=realType), intent(in), dimension(:) :: xx
        integer(kind=intType), intent(in), dimension(:) :: iblank

        character(80) :: zoneName
        integer(kind=intType) :: iDim

        write(zoneName, "(a,I5.5)") "Zone_", nBkGlobal
110     format('ZONE T=',a, " I=", i5, " J=", i5)
        write(101, 110) trim(zoneName), iEnd-iBeg+1, jEnd-jBeg+1
        write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"
13      format (E20.12)

        ! The 3 is for the three coordinate directions
        nNode = (iEnd - iBeg + 1)*(jEnd - jBeg + 1)
        nCell = (iEnd - iBeg)*(jEnd - jBeg)

        do i=1, 3*nNode
           write(101, 13) xx(i)
        end do

        do i=1, nCell
           write(101, *) iBlank(i)
        end do
      end subroutine writeZone
    end subroutine writeWalls

  subroutine bowTieAndIsolationElimination(famList, level, sps)

    use constants
    use blockPointers
    use communication
    use utils, only : setPointers
    use oversetCommUtilities, only : exchangeSurfaceIBlanks
    use sorting, only : famInList
    implicit none

    ! Input Parameters
    integer(kind=intType), intent(in), dimension(:) :: famList
    integer(kind=intType), intent(in) :: level, sps

    ! Local variables
    integer(kind=intType) :: mm, nn, i, j, k, e, iBeg, iEnd, jBeg, jEnd
    logical :: side(4)

    integer(kind=intType), dimension(:, :), pointer :: ibp, gcp
    integer(kind=intType), dimension(:, :), allocatable :: toFlip, nE, nC

    ! This routine initializes the surface cell iblank based on the
    ! volume iblank. It is not a straight copy since we a little
    ! preprocessing to eliminate a few particularly nasty cases.
    ! Three analysis are performed:
    ! 1. Bow-tie elimination
    ! 2. Single cell elmination

    bowTieLoop: do E=0, 2
       domainLoop1: do nn=1, nDom
          call setPointers(nn, level, sps)

          bocoLoop1: do mm=1, nBocos
             famInclude1: if (famInList(BCData(mm)%famID, famList)) then

                select case (BCFaceID(mm))
                case (iMin)
                   ibp => iblank(2, :, :)
                   gcp => globalCell(2, :, :)
                case (iMax)
                   ibp => iblank(il, :, :)
                   gcp => globalCell(il, :, :)
                case (jMin)
                   ibp => iblank(:, 2, :)
                   gcp => globalCell(:, 2, :)
                case (jMax)
                   ibp => iblank(:, jl, :)
                   gcp => globalCell(:, jl, :)
                case (kMin)
                   ibp => iblank(:, :, 2)
                   gcp => globalCell(:, :, 2)
                case (kMax)
                   ibp => iblank(:, :, kl)
                   gcp => globalCell(:, :, kl)
                end select

                ! -------------------------------------------------
                ! Step 2: Bow-tie elimination: Elimiate cells
                ! that touch only at a corner.
                ! -------------------------------------------------

                ! Make bounds a little easier to read. Owned cells only
                ! from now on.
                jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
                iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

                ! Allocate two tmporary auxilary arrays 'eN'->
                ! edgeNeighbours and 'cN'-> cornerNeighbous. For every
                ! comupte cell determine the number of compute
                ! neighbours connected along edges and at corners
                allocate(nE(iBeg:iEnd, jBeg:jEnd), nC(iBeg:iEnd, jBeg:jEnd))!, &

                call findBowTies()

                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      if (BCData(mm)%iBlank(i, j) > 0 .and. nC(i,j) >=1 .and. nE(i,j) <= E) then
                         BCData(mm)%iBlank(i, j) = 0
                      end if
                   end do
                end do

                deallocate(nC, nE)
             end if famInclude1
          end do bocoLoop1
       end do domainLoop1

       ! Since we potentially changed iBlanks, we need to updated by
       ! performing an exchange.
       call exchangeSurfaceIBlanks(famList, level, sps, commPatternCell_2nd, internalCell_2nd)

    end do bowTieLoop

    domainLoop2: do nn=1, nDom
       call setPointers(nn, level, sps)

       bocoLoop2: do mm=1, nBocos
          famInclude2: if (famInList(BCData(mm)%famID, famList)) then

             select case (BCFaceID(mm))
             case (iMin)
                gcp => globalCell(2, :, :)
             case (iMax)
                gcp => globalCell(il, :, :)
             case (jMin)
                gcp => globalCell(:, 2, :)
             case (jMax)
                gcp => globalCell(:, jl, :)
             case (kMin)
                gcp => globalCell(:, :, 2)
             case (kMax)
                gcp => globalCell(:, :, kl)
             end select

             ! Make bounds a little easier to read. Owned cells only
             ! from now on.
             jBeg = BCData(mm)%jnBeg+1 ; jEnd = BCData(mm)%jnEnd
             iBeg = BCData(mm)%inBeg+1 ; iEnd = BCData(mm)%inEnd

             ! -------------------------------------------------
             ! Step 3: Single-cell elimination: Elimiate cells
             ! that do not touch any other cells.
             ! -------------------------------------------------

             allocate(nE(iBeg:iEnd, jBeg:jEnd), nC(iBeg:iEnd, jBeg:jEnd))

             call setNeighbourCounts()

             ! This is easy, if a compute cell is stil around with no
             ! neighbours, kill it
             do j=jBeg, jEnd
                do i=iBeg, iEnd
                   if (BCData(mm)%iBlank(i, j) == 1 .and. nE(i,j) == 0 .and. nC(i,j) == 0) then
                      BCData(mm)%iBlank(i, j) = 0
                   end if
                end do
             end do

             deallocate(nE, nC)

          end if famInclude2
       end do bocoLoop2
    end do domainLoop2

    ! Again, since we potentially changed iBlanks, we need to updated by
    ! performing an exchange.
    call exchangeSurfaceIBlanks(famList, level, sps, commPatternCell_2nd, internalCell_2nd)
  contains
    subroutine findBowTies

      implicit none
      ! For every compute determine the number of compute neighbours
      ! connected along edges (nE) and the number of bow-ties (in nC)

      integer(kind=intType) :: i, j, e(4)

      nE = 0
      nC = 0

      do j=jBeg, jEnd
         do i=iBeg, iEnd
            if (BCData(mm)%iBlank(i, j) >= 0 ) then
               e = 0

               !    |  e3  |
               !  --c4-----c3-
               !    |      |
               ! e4 |   x  | e2
               !    |      |
               ! -- c1-----c2-
               !    |  e1  |

               ! Set the status of each of the 4 edges:

               if (BCData(mm)%iBlank(i, j-1) == 1) &
                    e(1) = 1

               if (BCData(mm)%iBlank(i+1, j) == 1) &
                    e(2) = 1

               if (BCData(mm)%iBlank(i, j+1) == 1) &
                    e(3) = 1

               if (BCData(mm)%iBlank(i-1, j) == 1) &
                    e(4) = 1

               ! Check the 4 corner neighbours for bow-tie status
               if (BCData(mm)%iBlank(i-1, j-1) == 1 .and. e(4) == 0 .and. e(1) == 0) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i+1, j-1) == 1 .and. e(1) == 0 .and. e(2) == 0) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i+1, j+1) == 1 .and. e(2) == 0 .and. e(3) == 0) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i-1, j+1) == 1 .and. e(3) == 0 .and. e(4) == 0) &
                    nC(i, j) = nC(i, j) + 1

               nE(i, j) = sum(e)
            end if
         end do
      end do
    end subroutine findBowTies

    subroutine setNeighbourCounts

      implicit none
      ! For every comute determine the number of compute neighbours
      ! connected along edges and at corners

      integer(kind=intType) :: i, j

      nE = 0
      nC = 0

      do j=jBeg, jEnd
         do i=iBeg, iEnd
            if (BCData(mm)%iBlank(i, j) == 1) then

               !    |  e3  |
               !  --c4-----c3-
               !    |      |
               ! e4 |   x  | e2
               !    |      |
               ! -- c1-----c2-
               !    |  e1  |

               ! Set the status of each of the 4 edges:

               if (BCData(mm)%iBlank(i, j-1) == 1) &
                    nE(i, j) = nE(i, j) + 1

               if (BCData(mm)%iBlank(i+1, j) == 1) &
                    nE(i, j) = nE(i, j) + 1

               if (BCData(mm)%iBlank(i, j+1) == 1) &
                    nE(i, j) = nE(i, j) + 1

               if (BCData(mm)%iBlank(i-1, j) == 1) &
                    nE(i, j) = nE(i, j) + 1

               ! Check the 4 corner neighbours for compute neighbour
               if (BCData(mm)%iBlank(i-1, j-1) == 1) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i+1, j-1) == 1) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i+1, j+1) == 1) &
                    nC(i, j) = nC(i, j) + 1

               if (BCData(mm)%iBlank(i-1, j+1) == 1) &
                    nC(i, j) = nC(i, j) + 1
            end if
         end do
      end do
    end subroutine setNeighbourCounts
  end subroutine bowTieAndIsolationElimination
end module zipperMesh
