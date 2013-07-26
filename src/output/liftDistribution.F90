subroutine addLiftDistribution(nSegments, dir_vec, distName)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will add the description of a lift distribution   *
  !      *                                                                *
  !      ******************************************************************
  
  use communication
  use liftDistributionData

  implicit none

  ! Input parameters
  character*(*), intent(in) :: distName
  integer(kind=intType), intent(in) :: nSegments
  real(kind=realType), dimension(3) :: dir_vec

  nLiftDists = nLiftDists + 1
  liftDists(nLIftDists)%nSegments = nSegments
  liftDists(nLiftDists)%dir = dir_vec
  liftDists(nLiftDists)%distName = distName
  liftDists(nLIftDists)%dir_ind = maxloc(dir_vec,1)
  
end subroutine addLiftDistribution

subroutine writeLiftDistributions(sps)

  use constants
  use communication
  use liftDistributionData
  use outputMod
  use su_cgns
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: sps

  real(kind=realType), dimension(3) :: xmin, xmax
  real(kind=realType), parameter :: tol=1e-8
  type(liftDist), pointer :: d
  integer(kind=intType) :: i, iDist, sizes(6)
  integer(kind=intType) :: cgnsInd, cgnsSol, cgnsBase, cgnsZone, ierr, coordID, fieldID, solID
  real(kind=realType), dimension(:,:), allocatable :: CoorX, CoorY, CoorZ
  real(kind=realType), dimension(:, :, :), allocatable :: values

  call initializeLiftDistributionData
  call liftDistGatherForcesAndNodes(sps)

  if (myid == 0) then 

     ! Get the bounding box for the entire geometry
     do i=1,3
        xmin(i) = minval(uniqueNodes(i,:))
        xmax(i) = maxval(uniqueNodes(i,:))
     end do
     
     do iDist=1,nLiftDists
        d => liftDists(iDist)
        d%delta = (xMax(d%dir_ind) - xMin(d%dir_ind))/dble((d%nSegments - 1))
        allocate(d%slicePts(3, d%nSegments))
        allocate(d%slices(d%nSegments))
        d%slicePts = zero
        do i=1,d%nSegments
           if (i == 1) then
              d%slicePts(d%dir_ind, i) = d%slicePts(d%dir_ind, i) + (i-1)*d%delta + tol
           else if (i == d%nSegments) then 
              d%slicePts(d%dir_ind, i) = d%slicePts(d%dir_ind, i) + (i-1)*d%delta -tol
           else
              d%slicePts(d%dir_ind, i) = d%slicePts(d%dir_ind, i) + (i-1)*d%delta 
           end if

           call createSlice(d%slices(i), d%slicePts(:, i), d%dir)
           call integrateSlice(d%slices(i))
        end do

        cgnsInd  = fileIDs(sps)
        cgnsBase = cgnsLiftDistBases(sps)
        sizes(1) = d%nSegments
        sizes(2) = 2
        sizes(3) = d%nSegments-1
        sizes(4) = 2-1
        sizes(5) = 0
        sizes(6) = 0

        call cg_zone_write_f(cgnsInd, cgnsBase, d%distName, sizes, &
             Structured, cgnsZone, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_sol_write_f(cgnsInd, cgnsBase, cgnsZone, &
             "data", Vertex, solID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        allocate(CoorX(d%nSegments, 2), CoorY(d%nSegments, 2), CoorZ(d%nSegments, 2))
          
        coorX = zero
        coorY = zero
        coorZ = zero

        if (d%dir_ind == 1) then! X slices
           coorX(:, 1) = d%slicePts(1, :)
           coorX(:, 2) = d%slicePts(1, :)
           coorY(:, 2) = one
        else if (d%dir_ind == 2) then ! Y slices
           coorY(:, 1) = d%slicePts(2, :)
           coorY(:, 2) = d%slicePts(2, :)
           coorZ(:, 2) = one
        else if (d%dir_ind == 3) then ! Z slices
           coorZ(:, 1) = d%slicePts(3, :)
           coorZ(:, 2) = d%slicePts(3, :)
           coorX(:, 2) = one
        end if
        
        ! Write each set of grid coords
        call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
             'CoordinateX', coorX, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
             'CoordinateY', coorY, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f
        
        call cg_coord_write_f(cgnsInd, cgnsBase, cgnsZone, realDouble, &
             'CoordinateZ', coorZ, coordID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! We currently have 4 values
        allocate(values(d%nSegments, 2, 4))
        
        do i=1,d%nSegments
           values(i, :, 1) = d%slices(i)%pL
           values(i, :, 2) = d%slices(i)%pD
           values(i, :, 3) = d%slices(i)%Clp
           values(i, :, 4) = d%slices(i)%Cdp
        end do

        ! Write the 4 of them:
        call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
             "Lift", values(:, :, 1), fieldID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
             "Drag", values(:, :, 2), fieldID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
             "Cl", values(:, :, 3), fieldID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        call cg_field_write_f(cgnsInd, cgnsBase, cgnsZone, solID, realDouble, &
             "Cd", values(:, :, 4), fieldID, ierr)
        if (ierr .eq. CG_ERROR) call cg_error_exit_f

        ! Destroy the lift slices:
        do i=1,d%nSegments
           call destroySlice(d%slices(i))
        end do
        
        ! Deallocate slice list and point list
        deallocate(d%slices, d%slicePts)

        ! Destroy temp variables
        deallocate(CoorX, CoorY, CoorZ, values)
        
     end do
  end if

end subroutine writeLiftDistributions

subroutine addParaSlice(sliceName, pt, direction)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will add a parametric slice to the list of user   *
  !      * supplied slices.                                               *
  !      *                                                                *
  !      ******************************************************************
  use communication
  use liftDistributionData

  implicit none

  ! Input parameters
  character*(*), intent(in) :: sliceName
  real(kind=realType), dimension(3), intent(in) :: pt, direction

  call initializeLiftDistributionData

  if (myid == 0) then
     nParaSlices = nParaSlices + 1
     call createSlice(paraSlices(nParaSlices), pt, direction)
     paraSlices(nParaSlices)%sliceName = sliceName
  end if
end subroutine addParaSlice

subroutine addAbsSlice(sliceName, pt, direction)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will add an absolute slice to the list of user    *
  !      * supplied slices.                                               *
  !      *                                                                *
  !      ******************************************************************
  use communication
  use liftDistributionData
  
  implicit none

  ! Input parameters
  character*(*), intent(in) :: sliceName
  real(kind=realType), dimension(3), intent(in) :: pt, direction

  call initializeLiftDistributionData

  if (myid == 0) then
     nAbsSlices = nAbsSlices + 1
     call createSlice(absSlices(nAbsSlices), pt, direction)
     absSlices(nAbsSlices)%sliceName = sliceName
  end if
  
end subroutine addAbsSlice

subroutine initializeLiftDistributionData
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine initializes and allocates the data requried    *
  !      * for slice options.                                             *
  !      *                                                                *
  !      * The purpose of this routine is t ost step is we need to        *
  !      * produce a finite-element type mesh for the entire wall         *
  !      * surface. There is a specific reason for this: The slicing      *
  !      * algorithm requies the function values to be stored at the      *
  !      * nodes. Since SUmb is cell centered, to produce the nodal       *
  !      * values, averaging must be done. However, an issue aries at     *
  !      * the the block boundaries; we don't have halo information       *
  !      * for the visous component of the force. Without halo            *
  !      * information, the resuling slice may have a discontinuity       *
  !      * at the block interfaces. To fully circumvent this              *
  !      * problem, we take an alternate approach. After                  *
  !      * communicating the nodes of the (wall) patches and their        *
  !      * connecitivity we find the duplicate nodes and then can         *
  !      * produce a compact finite-element mesh repsenting the           *
  !      * suface. Then the cell centered pressure and shear forces       *
  !      * defined on the faces, can be scattered to the nodes and        *
  !      * divided by the dual area to produce the traction forces        *
  !      * at the nodes. The issue of no halo information is              *
  !      * eliminated since the finite element mesh contains all the      *
  !      * nodes/faces. Note that determining the mapping from the        *
  !      * local patch nodes to the compact gloabl list only needs        *
  !      * to be computed once since the suface connectivity won't        *
  !      * change.                                                        *
  !      *                                                                *
  !      ******************************************************************
  !
  use communication
  use blockPointers
  use inputPhysics
  use outputMod
  use liftDistributionData

  implicit none

  ! Working Variables
  integer(kind=intType) :: ierr, i, j
  real(kind=realType), parameter :: tol=1e-8

  interface
     subroutine pointReduce(pts, N, tol, uniquePts, link, nUnique)
       use precision
       implicit none

       real(kind=realType), dimension(:, :) :: pts
       integer(kind=intType), intent(in) :: N
       real(kind=realType), intent(in) :: tol
       real(kind=realType), dimension(:, :) :: uniquePts
       integer(kind=intType), dimension(:) :: link
       integer(kind=intType) :: nUnique
     end subroutine pointReduce

  end interface

  if (liftDistInitialized) then
     return

  else
     ! Step 1a. Number of local nodes and quads on each proc:
     call getForceSize(nNodesLocal, nCellsLocal)

     ! Step 1b. Communicate to the root processor
     allocate(nNodesProc(nProc), nCellsProc(nProc), cumNodesProc(0:nProc), cumCellsProc(0:nProc))

     call mpi_allgather(nNodesLocal, 1, sumb_integer, nNodesProc, 1, sumb_integer, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call mpi_allgather(nCellsLocal, 1, sumb_integer, nCellsProc, 1, sumb_integer, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Step 1c: Fill Coordinates and Connecitivty for each proc. We will
     ! reuse some of the force code for simpliciity, even if it requires
     ! a bit more work 
     
     allocate(localCells(4, nCellsLocal))
     allocate(localNodes(3, nNodesLocal))

     call getForcePoints(localNodes, nNodesLocal, 1) ! Just use first sps
     ! instance here
     call getForceConnectivity(localCells, nCellsLocal)

     ! Step 1c. Communicate the coordinates to the root proc
     cumNodesProc(0) = 0_intType
     cumCellsProc(0) = 0_intType
     nNodesTotal = 0
     nCellsTotal = 0
     do i=1, nProc
        nNodesTotal = nNodesTotal + nNodesProc(i)
        cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)

        nCellsTotal = nCellsTotal + nCellsProc(i)
        cumCellsProc(i) = cumCellsProc(i-1) + nCellsProc(i)
     end do

     ! Modify the local connectivity with the offset of the global nodes:
     ! ie. first proc is offset by 0, second proc my the nodes on the
     ! first, third by the nodes on the first two, etc. ALso, the
     ! connectivity was zero based, we will convert to 1 ordering here
     do i=1,nCellsLocal
        localCells(:, i) = localCells(:, i) + cumNodesProc(myid) + 1
     end do

     if (myid == 0) then
        ! Allocate space for all nodes and cells only on root proc
        allocate(allNodes(3, nNodesTotal))
        allocate(allCells(4, nCellsTotal))
     end if

     ! Gather the coordinates the connectivity
     call mpi_gatherv(&
          localNodes, nNodesLocal*3, sumb_real, &
          allNodes, nNodesProc*3, cumNodesProc*3, sumb_real, 0, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call mpi_gatherv(&
          localCells, nCellsLocal*4, sumb_integer, &
          allCells, nCellsProc*4, cumCellsProc*4, sumb_integer, 0, &
          sumb_comm_world, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Step 1d. Determine the unique coordiantes and the new cell mapping:
     if (myid == 0) then
        ! Maximum space for the unique coordinates and link array
        allocate(uniqueNodes(3, nNodesTotal))
        allocate(link(nNodesTotal))
        allocate(allForcesP(3, nNodesTotal), allForcesV(3, nNodesTotal))

        call pointReduce(allNodes, nNodesTotal, tol, uniqueNodes, link, nUnique)

        ! Loop over the allCells and relpace the node index which
        ! currently refers to the indexing in 'allNodes' and replace with
        ! the ordering in 'uniqueNodes'
        do i=1,nCellsTotal
           do j=1,4
              allCells(j, i) = link(allCells(j, i))
           end do
        end do
     end if

     ! Allocate a few additional arrays that will be required for the forces:
     allocate(localForcesP(3, nNodesLocal), localForcesV(3, nNodesLocal))

     if (myid == 0) then
        allocate(uniqueTractionsP(3, nUnique), uniqueTractionsV(3, nUnique), &
             dualAreas(nUnique))
        allocate(fc(nUnique))
     end if

     ! Data for the marching squares method:
     ! Which edges are cut by the contour
     msCon1( 1 , : ) = (/ 0 , 0 , 0 , 0 , 0 /)
     msCon1( 2 , : ) = (/ 1 , 4 , 0 , 0 , 0 /)
     msCon1( 3 , : ) = (/ 1 , 2 , 0 , 0 , 0 /)
     msCon1( 4 , : ) = (/ 4 , 2 , 0 , 0 , 0 /)
     msCon1( 5 , : ) = (/ 2 , 3 , 0 , 0 , 0 /)
     msCon1( 6 , : ) = (/ 2 , 3 , 0 , 0 , 0 /) ! Should be 2, 3, 1, 4
     msCon1( 7 , : ) = (/ 1 , 3 , 0 , 0 , 0 /)
     msCon1( 8 , : ) = (/ 4 , 3 , 0 , 0 , 0 /)
     msCon1( 9 , : ) = (/ 4 , 3 , 0 , 0 , 0 /)
     msCon1( 10, : ) = (/ 1 , 3 , 0 , 0 , 0 /)
     msCon1( 11 , : ) = (/ 2 , 3 , 0 , 0 , 0 /) ! Should be 2, 3, 1, 4
     msCon1( 12 , : ) = (/ 2 , 3 , 0 , 0 , 0 /)
     msCon1( 13 , : ) = (/ 4 , 2 , 0 , 0 , 0 /)
     msCon1( 14 , : ) = (/ 1 , 2 , 0 , 0 , 0 /)
     msCon1( 15 , : ) = (/ 1 , 4 , 0 , 0 , 0 /)
     msCon1( 16 , : ) = (/ 0 , 0 , 0 , 0 , 0 /)

     msCon2(1, :) = (/1, 2/)
     msCon2(2, :) = (/2, 3/)
     msCon2(3, :) = (/3, 4/)
     msCon2(4, :) = (/4, 1/)

     liftDistInitialized = .True. 
  end if

end subroutine initializeLiftDistributionData

subroutine destroyLiftDistributionData
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine destroys the data allocated in                 *
  !      * initializeLiftDistribution.                                    *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  use communication
  implicit none

  ! Destroy the data allocated in initializeLiftDistribution

  if (liftDistInitialized) then

     if (myid == 0) then
        ! Free all data allocated only on root proc:
        deallocate(allNodes, allCells, link, uniqueNodes)
        deallocate(allForcesP, allForcesV)
        deallocate(uniqueTractionsP, uniqueTractionsV, dualAreas)
        deallocate(fc)
     end if

     ! Free data for localCells and Nodes
     deallocate(localCells, localNodes)

     ! Free the data for sizing
     deallocate(nNodesProc, nCellsProc, cumNodesProc, cumCellsProc)

     ! Deallocate local force data
     deallocate(localForcesP, localForcesV)

     liftDistInitialized = .False.
  end if

end subroutine destroyLiftDistributionData

subroutine liftDistGatherForcesAndNodes(sps)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine takes the nodes and forces from each processor *
  !      * and communicates to the root proc and assembles the global FE  *
  !      * mesh which is then ready for slicing operations.               *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  use communication
  use inputPhysics
  implicit none

  ! Input Param
  integer(kind=intType), intent(in) :: sps

  !orking param
  integer(kind=intType) :: i, j, ierr
  real(kind=realType), dimension(3) :: n1, n2, n3, n4, v1, v2, sss
  real(kind=realType) :: qa
  logical :: forcesTypeSave

  ! First get and communicate the coordinates

  call getForcePoints(localNodes, nNodesLocal, sps)

  call mpi_gatherv(&
       localNodes, nNodesLocal*3, sumb_real, &
       allNodes, nNodesProc*3, cumNodesProc*3, sumb_real, 0, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now get the forces. Make sure we get forces, and not tractions.
  forcesTypeSave = forcesAsTractions
  forcesAsTractions = .False.
  call getForces(localForcesP, localForcesV, nNodesLocal, sps)
  forcesAsTractions = forcesTypeSave

  ! Gather pressure and viscous forces to root proc:
  call mpi_gatherv(&
       localForcesP, nNodesLocal*3, sumb_real, &
       allForcesP, nNodesProc*3, cumNodesProc*3, sumb_real, 0, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call mpi_gatherv(&
       localForcesV, nNodesLocal*3, sumb_real, &
       allForcesV, nNodesProc*3, cumNodesProc*3, sumb_real, 0, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Process the unique set of forces and tractions using the finite elemnent data:
  if (myid == 0) then

     ! First get the unique set of forces by summing into
     ! uniquetractions. At this point, they are forces and not
     ! tractions. Also update the unique node list using allNodes
     uniqueTractionsP = zero
     uniqueTractionsV = zero
     
     do i=1,nNodesTotal
        uniqueTractionsP(:, link(i)) = uniqueTractionsP(:, link(i)) + allForcesP(:, i)
        uniqueTractionsV(:, link(i)) = uniqueTractionsV(:, link(i)) + allForcesV(:, i)
        uniqueNodes(:, link(i)) = allNodes(:, i)
     end do

     ! Next get the dual areas by computing area of each quad:
     dualAreas = zero
     do i=1,nCellsTotal
        n1 = uniqueNodes(:, allCells(1, i))
        n2 = uniqueNodes(:, allCells(2, i))
        n3 = uniqueNodes(:, allCells(3, i))
        n4 = uniqueNodes(:, allCells(4, i))

        ! Recall these are CCW ordering!
        v1 = n3 - n1
        v2 = n4 - n2

        ! Cross product
        sss(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
        sss(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
        sss(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))

        ! Quarter of area
        qa = fourth*sqrt(sss(1)*sss(1) + sss(2)*sss(2) + sss(3)*sss(3))

        ! Scatter to the nodes
        do j=1,4
           dualAreas(allCells(j, i)) = dualAreas(allCells(j, i)) + qa
        end do
     end do

     ! Now compute the tractions by just dividing by the area:
     do i=1,nUnique
        uniqueTractionsP(:, i) = uniquetractionsP(:, i) / dualAreas(i)
        uniqueTractionsV(:, i) = uniquetractionsV(:, i) / dualAreas(i)
     end do
  end if

end subroutine liftDistGatherForcesAndNodes

subroutine createSlice(slc, pt, dir)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine creates a slice on the FE mesh made by the     *
  !      * plane defined pt and dir. (a pt and direction uniquely         *
  !      * defines a plane.) The data is stored in the variable slc       *
  !      * which is a derived slice type.                                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  implicit none

  ! Input param
  type(slice), intent(inout) :: slc
  real(kind=realType), dimension(3), intent(in) :: pt, dir

  ! Working param
  integer(kind=intType) :: i, j, nMax
  integer(kind=intType) :: patchIndices(4), indexSquare, jj, kk, icon, iCoor, num1, num2
  real(kind=realType) :: f(4), d, ovrdnom
  logical :: logic1
  interface
     subroutine reallocateReal2(realArray, newSize1, newSize2, &
          oldSize1, oldSize2, alwaysFreeMem)
       use precision
       implicit none
       real(kind=realType), dimension(:,:), pointer :: realArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
            oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
     end subroutine reallocateReal2
     
     subroutine reallocateInteger2(intArray, newSize1, newSize2, &
          oldSize1, oldSize2, alwaysFreeMem)
       use precision
       implicit none
       integer(kind=intType), dimension(:,:), pointer :: intArray
       integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                            oldSize1, oldSize2
       logical, intent(in) :: alwaysFreeMem
     end subroutine reallocateInteger2
  end interface

  ! First step is to compute the 'function' value that will be used
  ! for the contour.

  ! Equation of plane: ax + by + cz + d = 0
  d = -pt(1)*dir(1) - pt(2)*dir(2) - pt(3)*dir(3)
  ovrdnom = one/sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
  do i=1,nUnique
     ! Now compute the signed distance
     fc(i) = (dir(1)*uniqueNodes(1, i) + dir(2)*uniqueNodes(2, i) + &
          dir(3)*uniqueNodes(3, i) + d)*ovrdnom
  end do

  ! Estimate size of slice by the sqrt of the number of nodes in the
  ! mesh
  nMax = int(sqrt(dble(nNodesTotal)))
  allocate(slc%w(2,nMax), slc%ind(2, nMax))

  iCoor = 0
  ! Loop over the cells
  do i=1,nCellsTotal

     ! Extract the indices and function values at each corner
     do jj=1,4
        patchIndices(jj) = allCells(jj, i)
        f(jj) = fc(patchIndices(jj))
     end do

     ! Based on the values at each corner, determine which
     ! type contour we have
     indexSquare = 1

     if (f(1) .lt. zero) indexsquare = indexsquare + 1
     if (f(2) .lt. zero) indexsquare = indexsquare + 2
     if (f(3) .lt. zero) indexsquare = indexsquare + 4
     if (f(4) .lt. zero) indexsquare = indexsquare + 8

     logic1 = .true.

     kk = 1
     do while (logic1)
        ! This is the edge
        icon = mscon1(indexSquare, kk)

        if (icon == 0) then
           logic1=.false.
        else

           ! num1, num2 are node indices
           num1 = mscon2(icon,1) 
           num2 = mscon2(icon,2)

           iCoor = iCoor + 1
           if (iCoor > nMax) then
              ! Need to reallocate the arrays. Make it double the size
              call reallocateReal2(slc%w, 2, 2*nMax, 2, nMax, .true.)
              call reallocateInteger2(slc%ind, 2, 2*nMax, 2, nMax, .true.)
              nMax = nMax * 2
           end if
           
           ! Weight factors
           slc%w(2, iCoor) = (zero - f(num1))/(f(num2) - f(num1))
           slc%w(1, iCoor) = one - slc%w(2, icoor)

           ! Store the weight factors
           slc%ind(:, iCoor) = (/patchIndices(num1), patchIndices(num2)/)

           kk = kk + 1
        end if
     end do
  end do ! Cell Loop
  
  slc%nNodes = iCoor

end subroutine createSlice

subroutine destroySlice(slc)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine destroys a slice created by the createSlice    *
  !      * routine                                                        *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  implicit none

  ! Input param
  type(slice), intent(inout) :: slc

  ! Deallocate weights and indices if they are already allocated
  if (associated(slc%w)) then
     deallocate(slc%w)
  end if

  if (associated(slc%ind)) then
     deallocate(slc%ind) 
  end if

end subroutine destroySlice

subroutine integrateSlice(slc) 
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine integrates the forces on slice slc and computes*
  !      * the integrated quantities for lift, drag, cl and cd with       *
  !      * contributions from both the pressure and visoucs forces.       *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  use inputPhysics   
  use flowVarRefState
  implicit none

  ! Input Variables
  type(slice) :: slc

  ! Working variables
  integer(kind=intType) :: i, j, i1, i2
  real(kind=realType), dimension(3) :: x1, x2, pT1, pT2, vT1, vT2, pF, vF
  real(kind=realType) :: len, dmax, dist, fact, scaleDim
  integer(kind=intType) :: bestPair(2)
  ! Loop over elements and integrate:

  pF = zero
  vF = zero
  do i=1,slc%nNodes/2

     ! extract nodes:
     i1 = 2*i-1
     i2 = 2*i
     x1 = slc%w(1,i1)*uniqueNodes(:, slc%ind(1,i1)) + slc%w(2,i1)*uniqueNodes(:, slc%ind(2, i1))
     x2 = slc%w(1,i2)*uniqueNodes(:, slc%ind(1,i2)) + slc%w(2,i2)*uniqueNodes(:, slc%ind(2, i2))

     ! extract pressure tractions
     pT1 = slc%w(1,i1)*uniqueTractionsP(:, slc%ind(1,i1)) + slc%w(2,i1)*uniqueTractionsP(:, slc%ind(2, i1))
     pT2 = slc%w(1,i2)*uniqueTractionsP(:, slc%ind(1,i2)) + slc%w(2,i2)*uniqueTractionsP(:, slc%ind(2, i2))

     ! extract viscous tractions
     vT1 = slc%w(1,i1)*uniqueTractionsV(:, slc%ind(1,i1)) + slc%w(2,i1)*uniqueTractionsV(:, slc%ind(2, i1))
     vT2 = slc%w(1,i2)*uniqueTractionsV(:, slc%ind(1,i2)) + slc%w(2,i2)*uniqueTractionsV(:, slc%ind(2, i2))

     ! Length of this segment
     len = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)
     
     ! Integrate the pressure and viscous forces separately
     pF = pF + half*(pT1 + pT2)*len
     vF = vF + half*(vT1 + vT2)*len
  end do

  ! Set the values we can in the slice
  slc%pL = liftDirection(1)*pF(1) + liftDirection(2)*pF(2) + liftDirection(3)*pF(3)
  slc%pD = dragDirection(1)*pF(1) + dragDirection(2)*pF(2) + dragDirection(3)*pF(3)
  slc%vL = liftDirection(1)*vF(1) + liftDirection(2)*vF(2) + liftDirection(3)*vF(3)
  slc%vD = dragDirection(1)*vF(1) + dragDirection(2)*vF(2) + dragDirection(3)*vF(3)

  ! Compute the chord as the max lenght between two nodes...this is
  ! n^2, so should be changed in the future

  dmax = zero
  do i=1,slc%nNodes
     ! extract node:
     x1 = slc%w(1,i)*uniqueNodes(:, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(:, slc%ind(2, i))

     do j=i+1,slc%Nnodes
        ! extract node:
        x2 = slc%w(1,j)*uniqueNodes(:, slc%ind(1,j)) + slc%w(2,j)*uniqueNodes(:, slc%ind(2, j))
        
        dist = sqrt((x1(1)-x2(1))**2 +(x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)
        
        if (dist > dmax) then
           dmax = dist
           bestPair = (/i, j/)
        end if
     end  do
  end do

  ! Set chord
  slc%chord = dmax
  
  ! Compute factor to get coefficient
  scaleDim = pRef/pInf
  fact = two/(gammaInf*pInf*MachCoef*MachCoef &
       *surfaceRef*LRef*LRef*scaleDim)

  ! Take dmax as chord:
  slc%CLp = slc%pL / slc%chord * fact
  slc%CDp = slc%pD / slc%chord * fact
  slc%CLv = slc%vL / slc%chord * fact
  slc%CDv = slc%vD / slc%chord * fact

end subroutine integrateSlice
  
! subroutine writeSliceTecplot(slc, fileName)

!   use liftDistributionData
!   implicit none

!   character*(*), intent(in) :: fileName
!   type(slice), intent(in) :: slc
!   integer(kind=intType) :: i 
!   open (unit=7, file=fileName)
!   write (7,*) "Title = ""Lift Dist Test"""
!   write (7,*) "Variables = ""X"" ""Y"" ""Z"""! ""Tx"" ""Ty"" ""Tz"""
!   write (7,*) "Zone T=""FE Mesh"""
!   write (7,*) "Nodes = ", slc%nNodes, " Elements= ", slc%nNodes/2, " ZONETYPE=FELINESEG"
!   write (7,*) "DATAPACKING=POINT"
! 13 format (E16.10, E16.10, E16.10)
!   do i=1,slc%nNodes
!      write(7,13) &
!           slc%w(1,i)*uniqueNodes(1, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(1, slc%ind(2, i)), &
!           slc%w(1,i)*uniqueNodes(2, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(2, slc%ind(2, i)), &
!           slc%w(1,i)*uniqueNodes(3, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(3, slc%ind(2, i))!, &
!           ! slc%w(1,i)*uniqueTractionsP(1, slc%ind(1,i)) + slc%w(2,i)*uniqueTractionsP(1, slc%ind(2, i)), &
!           ! slc%w(1,i)*uniqueTractionsP(2, slc%ind(1,i)) + slc%w(2,i)*uniqueTractionsP(2, slc%ind(2, i)), &
!           ! slc%w(1,i)*uniqueTractionsP(3, slc%ind(1,i)) + slc%w(2,i)*uniqueTractionsP(3, slc%ind(2, i))
!   end do
          
! 14 format(I5, I5)
!   do i=1, slc%nNodes/2
!      write(7, 14) 2*i-1, 2*i
!   end do
!   close(7)

! end subroutine writeSliceTecplot
