subroutine setSurfaceFamilyInfo

  use constants
  use su_cgns
  use blockPointers, onlY : nDom, flowDoms, nBocos, cgnsSubFace, BCType
  use cgnsGrid, onlY : cgnsDoms
  use communication, only : myid, sumb_comm_world
  use inputTimeSpectral, only : nTimeIntervalsSpectral
  use surfaceFamilies, only : wallExchange, familyExchanges, famNames, &
       famGroups, wallFamilies, famIsWall, totalFamilies, totalWallFamilies
  use utils, only : setPointers, EChk, pointReduce, terminate, convertToLowerCase
  use sorting, only : qsortStrings, bsearchStrings
  implicit none
  !
  !      Local variables.
  !
  integer :: ierr

  integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax, i, j, k, nFam, famID, cgb, iFam
  integer(kind=intType) :: sps, currentFamily, curFam(1)
  character(maxCGNSNameLen), dimension(25) :: defaultFamName 
  character(maxCGNSNameLen) :: curStr, family
  character(maxCGNSNameLen), dimension(:), allocatable :: fulLFamList, uniqueFamList
  integer(kind=intType), dimension(:), allocatable :: localFlag

  ! Process out the family information. The goal here is to
  ! assign a unique integer to each family in each boundary
  ! condition. The CGNS grid has all the information we need.

  ! Firstly make sure that there is actual family specified for
  ! each BC. If there isn't, we will provide one for you. 
  defaultFamName(BCAxisymmetricWedge) = 'axi'
  defaultFamName(BCDegenerateLine) = 'degenerate'
  defaultFamName(BCDegeneratePoint) ='degenerate'
  defaultFamName(BCDirichlet) = 'dirichlet'
  defaultFamName(BCExtrapolate) = 'extrap'
  defaultFamName(BCFarfield) = 'far'
  defaultFamName(BCGeneral) = 'general'
  defaultFamName(BCInflow) = 'inflow'
  defaultFamName(BCInflowSubsonic) = 'inflow'
  defaultFamName(BCInflowSupersonic) = 'inflow'
  defaultFamName(BCNeumann) = 'neumann'
  defaultFamName(BCOutflow) = 'outflow'
  defaultFamName(BCOutflowSubsonic) = 'outflow'
  defaultFamName(BCOutflowSupersonic)  ='outflow'
  defaultFamName(BCSymmetryPlane) = 'sym'
  defaultFamName(BCSymmetryPolar) = 'sympolar'
  defaultFamName(BCTunnelInflow) = 'inflow'
  defaultFamName(BCTunnelOutflow) = 'outflow'
  defaultFamName(BCWall) = 'wall'
  defaultFamName(BCWallInviscid) = 'wall'
  defaultFamName(BCWallViscous) = 'wall'
  defaultFamName(BCWallViscousHeatFlux) = 'wall'
  defaultFamName(BCWallViscousIsothermal) = 'wall'
  defaultFamName(UserDefined) = 'userDefined'

101 format("CGNS Block ",I4,", boundary condition ",I4, ", of type ",a, &
       " does not have a family. Based on the boundary condition type," &
       " a name of: '", a, "' will be used.")

  nFam = 0
  do i=1, size(cgnsDoms)
     do j=1, size(cgnsDoms(i)%bocoInfo)
        if (cgnsDoms(i)%bocoInfo(j)%actualFace) then 
           if (trim(cgnsDoms(i)%bocoInfo(j)%wallBCName) == "") then 
              if (myid == 0) then 
                 ! Tell the user we are adding an automatic family name
                 write(*, 101), i, j, trim(BCTypeName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS)), &
                      trim(defaultFamName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS))
              end if
              cgnsDoms(i)%bocoInfo(j)%wallBCName = trim(defaultFamName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS))
           end if
           nFam = nFam + 1
        end if
     end do
  end do

  ! Allocate space for the full family list
  allocate(fullFamList(nFam))
  nFam = 0
  do i=1, size(cgnsDoms)
     do j=1, size(cgnsDoms(i)%bocoInfo)
        if (cgnsDoms(i)%bocoInfo(j)%actualFace) then 
           nFam = nFam + 1
           fullFamList(nfam) = cgnsDoms(i)%bocoInfo(j)%wallBCName
           call convertToLowerCase(fullFamList(nFam))
        end if
     end do
  end do

  ! Now sort the family names:
  call qsortStrings(fullFamList, nFam)

  ! Next we need to generate a unique set of names. 
  allocate(uniqueFamList(nFam))

  curStr = fullFamList(1)
  uniqueFamList(1) = curStr
  j = 1
  i = 1
  do while(i < nFam)

     i = i + 1
     if (fullFamList(i) == curStr) then 
        ! Same str, do nothing. 
     else
        j = j + 1
        curStr = fullFamList(i)
        uniqueFamList(j) = curStr
     end if
  end do


  totalFamilies = j
  ! Now copy the uniqueFamList back to "fullFamList" and allocate
  ! exactly the right size. 
  deallocate(fullFamList)
  allocate(fullFamList(totalFamilies))
  fulLFamList(1:totalFamilies) = uniqueFamList(1:totalFamilies)
  deallocate(uniqueFamList)

  ! Now each block boundary condition can uniquely determine it's
  ! famID. We do all BC on all blocks and levels. 
  nLevels = ubound(flowDoms,2)
  do nn=1, nDom
     call setPointers(nn, 1_intType, 1_intType)
     do mm=1, nBocos

        cgb = flowDoms(nn, 1, 1)%cgnsBlockID
        family = cgnsDoms(cgb)%bocoInfo(cgnsSubface(mm))%wallBCName
        call convertToLowerCase(family)

        famID = bsearchStrings(family, fullFamList, totalFamilies)
        if (famID == 0) then 
           ! Somehow we never found the family...
           call terminate("setSurfaceFamilyInfo", &
                "An error occuring in assigning families")
        end if

        ! Now set the data on each of the level/sps instances
        do sps=1, nTimeIntervalsSpectral
           do level=1,nlevels

              flowDoms(nn, level, sps)%bcData(mm)%famID = famID
              flowDoms(nn, level, sps)%bcData(mm)%family = family

           end do
        end do
     end do
  end do

  allocate(famNames(totalFamilies), famIsWall(totalFamilies), localFlag(totalFamilies))
  localFlag = 0
  famNames = fullFamList

  ! Determine which of the unique families are walls. This is
  ! slightly inefficient but not terribly so.
  do iFam=1, totalFamilies

     do nn=1,nDom
        call setPointers(nn, 1_intType, 1_intType)
        do mm=1, nBocos
           if (flowDoms(nn, 1, 1)%bcData(mm)%famID == iFam .and. &
                (BCType(mm) == EulerWall .or. BCtype(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSwallIsoThermal)) then 
              localFlag(iFam) = 1
           end if
        end do
     end do
  end do

  call mpi_allreduce(localFlag, famIsWall, totalFamilies, sumb_integer, MPI_SUM, sumb_comm_world, ierr)

  ! Save the wall family list. 
  totalWallFamilies = 0
  do i=1,totalFamilies
     if (famIsWall(i) > 0) then 
        totalWallFamilies = totalWallFamilies + 1
     end if
  end do

  allocate(wallFamilies(totalWallFamilies))
  k = 0
  do i=1,totalFamilies
     if (famIsWall(i) > 0) then 
        k = k + 1
        wallFamilies(k) = i
     end if
  end do

  ! Determine the local-to-global mapping for the traction
  ! computation. Before we do so, set the famGroup list to
  ! include all boundary conditions

  allocate(famGroups(totalFamilies))

  ! Finally create the scatter context for each individual family
  ! as well as a special one for all wall families

  allocate(familyExchanges(totalFamilies, ntimeIntervalsSpectral))

  do sps=1, nTimeIntervalsSpectral
     do i=1, totalFamilies
        curFam(1) = i
        nFam = 1
        call createNodeScatterForFamilies(curFam, nFam, familyExchanges(i, sps), sps)
     end do
  end do

  allocate(wallExchange(nTimeIntervalsSpectral))
  do sps=1, nTimeIntervalsSpectral
     call createNodeScatterForFamilies(wallFamilies(sps), totalWallFamilies, wallExchange(sps), sps)
  end do

end subroutine setSurfaceFamilyInfo

subroutine createNodeScatterForFamilies(famList, nFam, famExchange, sps)

  ! The purpose of this routine is to create the appropriate data
  ! structures that allow for the averaging of cell based surface
  ! quantities to node-based quantities. The primary reason for this
  ! is that the viscous stress tensor is not available at halo cells
  ! and therefore it is not possible to create consistent node-based
  ! values locallly. What the scatter does is allows us to sum the
  ! nodal values across processors, average them and finally update
  ! the node based values to be consistent. This operation is
  ! necessary for several operations:

  ! 1. Integration of forces over zipper triangles requires force/area
  ! at nodes.
  ! 2. Lift distributions/slices also requires node-based tractions
  ! 3. Node-based output for tecplot files.
  use constants
  use communication, only : sumb_comm_world, myid, nProc
  use surfaceFamilies, only : famGroups, familyExchange, &
       IS1, IS2, PETSC_COPY_VALUES, PETSC_DETERMINE
  use utils, only : pointReduce, eChk
  implicit none

  ! Input Parameters
  integer(kind=intType) , dimension(nFam), intent(in) :: famList
  integer(kind=intType) , intent(in) :: nFam, sps
  type(familyExchange), intent(inout) :: famExchange

  ! Working param
  integer(kind=intType) :: i,j, ierr, nNodesLocal, nNodesTotal, nCellsLocal
  integer(kind=intType) :: nUnique, iSize, iStart, iEnd, iProc
  real(kind=realType), dimension(:, :), allocatable :: localNodes, allNodes
  real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
  integer(kind=intType), dimension(:), allocatable :: link, localIndices, startIndices, endIndices
  integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
  real(kind=realType) :: tol
  integer(kind=intType) :: status(MPI_STATUS_SIZE)

  ! Note that for overset, the reduction must be done by cluster. 

  ! Get the nodes for the surfaces. Note that the nodeScatter works
  ! over *all* possible boundary condition surfaces. 
  if (allocated(famGroups)) then 
     deallocate(famGroups)
  end if
  allocate(famGroups(nFam), famExchange%famGroups(nFam))
  famGroups = famList
  famExchange%famGroups = famList
  famExchange%nFam = nFam
  famExchange%sps = sps

  call getSurfaceSize(nNodesLocal, nCellsLocal, famList, nFam)
  famExchange%nNodes = nNodesLocal
  ! Allocate space for the nodes/connectivities
  allocate(famExchange%conn(4, nCellsLocal), famExchange%elemFam(nCellsLocal), &
       famExchange%fc(nNodesLocal), famExchange%nodalValues(nNodesLocal, 3))

  call getSurfaceConnectivity(famExchange%conn, nCellsLocal)
  call getSurfaceFamily(famExchange%elemFam, nCellsLocal)


  allocate(localNodes(3, nNodesLocal), nNodesProc(nProc), cumNodesProc(0:nProc))
  call getSurfacePoints(localNodes, nNodesLocal, sps)

  do i=1, nNodesLocal
     famExchange%nodalValues(i, 1:3) = localNodes(1:3, i)
  end do

  ! Determine the total number of nodes on each proc
  call mpi_allgather(nNodesLocal, 1, sumb_integer, nNodesProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Determine cumulative version
  cumNodesProc(0) = 0_intType
  nNodesTotal = 0
  do i=1, nProc
     nNodesTotal = nNodesTotal + nNodesProc(i)
     cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)
  end do

  ! Send all the nodes to everyone
  allocate(allNodes(3, nNodesTotal))
  call mpi_allgatherv(localNodes, nNodesLocal*3, sumb_real, allNodes, &
       nNodesProc*3, cumNodesProc*3, sumb_real, sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Local nodes is no longer necessary
  deallocate(localNodes)

  ! Now point reduce
  allocate(uniqueNodes(3, nNodestotal), link(nNodestotal))
  tol = 1e-12

  call pointReduce(allNodes, nNodesTotal, tol, uniqueNodes, link, nUnique)

  ! We can immediately discard everything but link since we are only
  ! doing logical operations here:
  deallocate(uniqueNodes, allNodes)

  ! Now back out the global indices for our local points
  allocate(localIndices(nNodesLocal))
  do i=1, nNodesLocal
     ! The -1 is to convert to 0-based ordering for petsc
     localIndices(i) = link(cumNodesProc(myid) + i)-1
  end do

  ! Create the basic (scalar) local vector
  call VecCreateMPI(SUMB_COMM_WORLD, nNodesLocal, PETSC_DETERMINE, &
       famExchange%nodeValLocal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPI(SUMB_COMM_WORLD, nCellsLocal, PETSC_DETERMINE, &
       famExchange%localWeight, ierr)
  call EChk(ierr,__FILE__,__LINE__)


  ! Create the basic global vector. This is slightly tricker than it
  ! sounds. We could just make it uniform, but then there would be
  ! more communicaiton than necessary. Instead what we do is determine
  ! the min and max range of local indices on the proc and the one
  ! before it. A little diagram will help
  !
  ! Proc 0 +---------------------+
  ! Proc 1                    +-------------+
  ! Proc 2                                +----------------+
  !
  ! Proc zero has a many global nodes as local since they are by
  ! definition all unqiue. Proc 1 then will start at 1 more than the
  ! proc 0 and continue to it's maximum value. Proc 2 starts at the
  ! end of proc 1 etc. This way the vast majority of the global nodes
  ! are owned locally. 

  ! In order to determine the owning range for each processor, it is
  ! much trickier than it sounds. We do a linear cascasde through the
  ! procs sending the upper range from proc 0 to proc 1, then proc1 to
  ! proc 2 and so on.

  ! Proc zero owns all of it's nodes. 
  if (myid == 0) then
     iStart = 0
     if (nNodesLocal == 0) then 
        iEnd = 0
     else
        iEnd = maxval(localIndices) + 1
     end if
  end if

  do iProc=0, nProc-2
     if (myid == iProc) then 
        ! I need to send my iEnd to proc+1
        call mpi_send(iEnd, 1, sumb_integer, iProc+1, iProc, sumb_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)
     else if(myid == iProc+1) then 

        ! Receive the value from the proc below me:
        call mpi_recv(iEnd, 1, sumb_integer, iProc, iProc, sumb_comm_world, status, ierr)
        call EChk(ierr,__FILE__,__LINE__)
        
        ! On this proc, the start index is the 
        iStart = iEnd
        if (nNodesLOCAl == 0) then 
           iEnd = iStart
        else
           iEnd = max(iStart, maxval(localIndices)+1)
        end if
     end if
  end do

  iSize = iEnd-iStart
  !print *,'isize:', myid, iEnd, iStart, iSize, nnodeslocal, nunique
  ! Create the actual global vec. Note we also include nUnique to make
  ! sure we have all the local sizes correct.
  call VecCreateMPI(SUMB_COMM_WORLD, iSize, nUnique, &
       famExchange%nodeValGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDuplicate(famExchange%nodeValGlobal, famExchange%sumGlobal, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now create the scatter that goes from the local vec to the global
  ! vec.

  ! Indices for the local vector is just a stride, starting at the
  ! offset
  call ISCreateStride(SUMB_COMM_WORLD, nNodesLocal, cumNodesProc(myid), &
       1, IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Indices for the global vector are the "localIndices" we previously
  ! computed. 
  call ISCreateGeneral(sumb_comm_world, nNodesLocal, localIndices, &
       PETSC_COPY_VALUES, IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterCreate(famExchange%nodeValLocal, IS1, famExchange%nodeValGlobal, IS2, &
       famExchange%scatter, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! And dont' forget to destroy the index sets
  call ISDestroy(IS1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call ISDestroy(IS2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  deallocate(localIndices)
  famExchange%allocated = .True.
end subroutine createNodeScatterForFamilies


! The only two groups we need to deal with in fortran directly is the
! fullFamilyList and all walls. We write special routines for them. 
subroutine setFullFamilyList()

  use surfaceFamilies
  implicit none

  integer(kind=intType) :: i
  if (allocated(famGroups)) then 
     deallocate(famGroups)
  end if
  allocate(famGroups(totalFamilies))
  do i=1, totalFamilies
     famGroups(i) = i
  end do

end subroutine setFullFamilyList

subroutine setWallFamilyList()

  use surfaceFamilies
  implicit none
  integer(kind=intType) :: i

  if (allocated(famGroups)) then 
     deallocate(famGroups)
  end if
  allocate(famGroups(totalWallFamilies))
  do i=1, totalWallFamilies
     famGroups(i) = wallFamilies(i)
  end do

end subroutine setWallFamilyList

subroutine setFamilyInfo(famList, n)

  use surfaceFamilies
  implicit none
  integer(kind=intType), intent(in) :: famList(n), n

  if (allocated(famGroups)) then 
     deallocate(famGroups)
  end if
  allocate(famGroups(n))
  famGroups = famList

end subroutine setFamilyInfo

subroutine mapVector(vec1, n1, famList1, nf1, vec2, n2, famList2, nf2)

  ! Map one vector, vec1 of size (3,n1) defined on family list 'famList1' onto
  ! vector, vec2, of size (3, n2) defined on family list 'famList2'

  ! This operation is actually pretty fast since it just requires a
  ! single copy of surface-based data. 
  use constants
  use blockPointers
  use sorting, only : bsearchIntegers
  implicit none

  ! Input/Output
  integer(kind=intType) :: n1, n2, nf1, nf2
  integer(kind=intType), intent(in) :: famList1(nf1), famList2(nf2)
  real(kind=realType), intent(in) :: vec1(3, n1)
  real(kind=realType), intent(inout) :: vec2(3, n2)

  ! Working
  integer(kind=intType) :: k, ii, jj, nn, mm, iSize, iBeg, iEnd, jBeg, jEnd, famID
  logical :: fam1Included, fam2Included

  ii = 0 
  jj = 0
  domains: do nn=1,nDom
     ! Don't set pointers for speed
     
     ! Loop over the number of boundary subfaces of this block.
     bocos: do mm=1,flowDoms(nn, 1, 1)%nBocos
        famId = flowDoms(nn, 1, 1)%BCdata(mm)%famID

        fam1Included = bsearchIntegers(famID, famList1, nf1) > 0
        fam2Included = bsearchIntegers(famID, famList2, nf2) > 0

        jBeg = flowDoms(nn, 1, 1)%bcData(mm)%jnBeg
        jEnd = flowDoms(nn, 1, 1)%bcData(mm)%jnEnd
        
        iBeg = flowDoms(nn, 1, 1)%bcData(mm)%inBeg
        iEnd = flowDoms(nn, 1, 1)%bcData(mm)%inEnd
        iSize = (iEnd-iBeg+1)*(jEnd-jBeg+1)

        if (fam1Included .and. fam2Included) then 
           ! The two lists overlap so copy:
           do k=1, iSize
              vec2(:, k+jj) = vec1(:, k+ii)
           end do
        end if

        ! Finally increment the counters if the face had been inclded
        if (fam1Included) then 
           ii = ii + iSize
        end if

        if (fam2Included) then 
           jj =jj + iSize
        end if

     end do bocos
  end do domains
end subroutine mapVector
