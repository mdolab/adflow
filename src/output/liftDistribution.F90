subroutine addLiftDistribution(nSegments, dir_vec, dir_ind, distName)
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
  integer(kind=intType), intent(in) :: dir_ind

  nLiftDists = nLiftDists + 1
  liftDists(nLIftDists)%nSegments = nSegments
  liftDists(nLiftDists)%dir = dir_vec
  liftDists(nLiftDists)%distName = distName
  liftDists(nLIftDists)%dir_ind = dir_ind

end subroutine addLiftDistribution

subroutine writeSlicesFile(fileName)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will write the user defined slics to an           *
  !      * to the (ascii) tecplot file fileName. ASCII files are          *
  !      * used for siplicity since very little informatin is actually    *
  !      * written.                                                       *
  !      *                                                                *
  !      ******************************************************************
  use communication
  use liftDistributionData
  use outputMod
  use inputTimeSpectral
  use inputPhysics
  use inputIteration
  implicit none

  ! Input Params
  character*(*), intent(in) :: fileName

  ! Working parameters
  integer(kind=intType) :: file, i, sps, nSolVar
  character(len=maxStringLen) :: fname, sliceName
  character(len=7) :: intString
  real(kind=realType), dimension(3) :: pt, dir
  character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
  
  ! Only write if we actually have lift distributions
  testwriteSlices: if(nParaSlices + nAbsSlices > 0) then

     if(myID == 0 .and. printIterations) then
        print "(a)", "#"
        print "(a)", "# Writing slices file(s) ..."
     endif


     do sps=1,nTimeIntervalsSpectral

        ! If it is time spectral we need to agument the filename
        if (equationMode == timeSpectral) then
           write(intString,"(i7)") sps
           intString = adjustl(intString)
           fname = trim(fileName)//"Spectral"//trim(intString)
        else
           fname = fileName
        end if

        file = 11
        ! Open file on root proc:
        if (myid == 0) then 
           open(unit=file, file=trim(fname))

           ! Write Header Information 
           write (file,*) "Title = ""SUmb Slice Data"""
           write (file,"(a)", advance="no") "Variables = "
           write(file,"(a)",advance="no") " ""CoordinateX"" "
           write(file,"(a)",advance="no") " ""CoordinateY"" "
           write(file,"(a)",advance="no") " ""CoordinateZ"" "

           ! Number of additional variables on slices
           call numberOfSurfSolVariables(nSolVar)
           allocate(solNames(nSolVar))
           call surfSolNames(solNames)

           ! Write the rest of the variables
           do i=1,nSolVar
              write(file,"(a,a,a)",advance="no") """",trim(solNames(i)),""" "
           end do
           write(file,"(1x)")
           deallocate(solNames)

           do i=1,nParaSlices
              call integrateSlice(paraSlices(i))
              call writeSlice(paraSlices(i), file, ifzv+nSolVar)
           end do

           do i=1,nAbsSlices
              ! 'Destroy' the slice...just dealloc the data
              call destroySlice(absSlices(i))

              ! Make new one in the same location
              call createSlice(absSlices(i), absSlices(i)%pt, absSlices(i)%dir)

              ! Integrate and write
              call integrateSlice(absSlices(i))
              call writeSlice(absSlices(i), file, ifzv+nSolVar)
           end do

           ! Close file on root proc
           close(file)
        end if
     end do

     if(myID == 0 .and. printIterations) then
        print "(a)", "# Slices file(s) written"
        print "(a)", "#"
     endif
  end if testwriteSlices
end subroutine writeSlicesFile

subroutine writeLiftDistributionFile(fileName)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will write the added lift distributions           *
  !      * to the (ascii) tecplot file fileName. ASCII files are          *
  !      * used for siplicity since very little informatin is actually    *
  !      * written.                                                       *
  !      *                                                                *
  !      ******************************************************************
  use communication
  use liftDistributionData
  use outputMod
  use inputPhysics
  use inputTimeSpectral
  use inputIteration
  implicit none

  ! Input Params
  character*(*), intent(in) :: fileName

  ! Working parameters
  integer(kind=intType) :: file, sps
  character(len=maxStringLen) :: fname
  character(len=7) :: intString

  ! Only write if we actually have lift distributions
  testwriteLiftDists: if(nLiftDists > 0) then

     if(myID == 0 .and. printIterations) then
        print "(a)", "#"
        print "(a)", "# Writing lift distribution file(s) ..."
     endif

     do sps=1,nTimeIntervalsSpectral

        ! If it is time spectral we need to agument the filename
        if (equationMode == timeSpectral) then
           write(intString,"(i7)") sps
           intString = adjustl(intString)
           fname = trim(fileName)//"Spectral"//trim(intString)
        else
           fname = fileName
        end if


        file = 11
        ! Open file on root proc:
        if (myid == 0) then 
           open(unit=file, file=trim(fname))
        end if

        call writeLiftDistributions(sps, file)

        ! Close file on root proc
        if (myid == 0) then
           close(file)
        end if
     end do

     if(myID == 0 .and. printIterations) then
        print "(a)", "# Lift distribution file(s) written"
        print "(a)", "#"
     endif

  end if testwriteLiftDists
end subroutine writeLiftDistributionFile

subroutine writeLiftDistributions(sps, fileID)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine writes the liftdistribution for the specified  *
  !      * spectral instance. It is assumed that the required file handles*
  !      * are already open and can be written to                         *
  !      *                                                                *
  !      ******************************************************************
  use constants
  use communication
  use liftDistributionData
  use outputMod
  use su_cgns
  use cgnsNames
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: sps, fileID

  real(kind=realType), dimension(3) :: xmin, xmax
  real(kind=realType), parameter :: tol=1e-8
  type(liftDist), pointer :: d
  integer(kind=intType) :: i, j, iVar, iDist, sizes(6)
  real(kind=realType), dimension(:,:), allocatable :: values
  character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistNames
  real(kind=realType) :: dmin, dmax, sumL, sumD, span, delta
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

        ! These are the variable names for the lift distribution:
        allocate(liftDistNames(nLiftDistVar))
        ! Set the names here
        liftDistNames(1) = "eta"
        liftDistNames(2) = cgnsCoorX
        liftDistNames(3) = cgnsCoorY
        liftDistNames(4) = cgnsCoorZ
        liftDistNames(5) = "Lift"
        liftDistNames(6) = "Drag"
        liftDistNames(7) = "Normalized Lift"
        liftDistNames(8) = "Normalized Drag"
        liftDistNames(9) = "CL"
        liftDistNames(10) = "CD"
        liftDistNames(11) = "CLp"
        liftDistNames(12) = "CDp"
        liftDistNames(13) = "CLv"
        liftDistNames(14) = "CDv"
        liftDistNames(15) = "Ellptical"
        liftDistNames(16) = "thickness"
        liftDistNames(17) = "twist"
        liftDistNames(18) = "chord"

        ! Only write header info for first distribution only
        if (iDist == 1) then
           write (fileID,*) "Title= ""SUmb Lift Distribution Data"""
           write (fileID,"(a)", advance="no") "Variables = "
           do i=1,nLIftDistVar
              write(fileID,"(a,a,a)",advance="no") """",trim(liftDistNames(i)),""" "
           end do
           write(fileID,"(1x)")
        end if

        write (fileID,"(a,a,a)") "Zone T= """,trim(d%distName),""""
        write (fileID,*) "I= ",d%nSegments
        write (fileID,*) "DATAPACKING=BLOCK"
15      format (E14.6)

        allocate(values(d%nSegments, nLiftDistVar))
        values = zero

        ! Scaled Eta values
        dmin = minVal(d%slicePts(d%dir_ind, :))
        dmax = maxval(d%slicePts(d%dir_ind, :))
        values(:, 1) = d%slicePts(d%dir_ind, :)/dmax

        ! Coordinate Varaibles
        if (d%dir_ind == 1) then! X slices
           values(:, 2) = d%slicePts(1, :)
        else if (d%dir_ind == 2) then ! Y slices
           values(:, 3) = d%slicePts(2, :)
        else if (d%dir_ind == 3) then ! Z slices
           values(:, 4) = d%slicePts(3, :)
        end if

        ! Other variables
        do i=1,d%nSegments
           ! Total lift and drag 
           values(i, 5)  = d%slices(i)%pL + d%slices(i)%vL
           values(i, 6)  = d%slices(i)%pD + d%slices(i)%vD

           ! Total CL and CD
           values(i, 9)  = d%slices(i)%CLp + d%slices(i)%CLv
           values(i, 10) = d%slices(i)%CDp + d%slices(i)%CDv

           ! Pressure lift and drag coefficients
           values(i, 11) = d%slices(i)%CLp
           values(i, 12) = d%slices(i)%CDp

           ! Viscous lift and drag coefficients
           values(i, 13) = d%slices(i)%CLv
           values(i, 14) = d%slices(i)%CDv

           ! t/c, twist, chord
           values(i, 16) = d%slices(i)%thickness
           values(i, 17) = d%slices(i)%twist
           values(i, 18) = d%slices(i)%chord
        end do

        ! Sum up the lift and drag values from the slices:
        sumL = zero
        sumD = zero
        do i=1,d%nSegments-1
           sumL = sumL + half*(values(i, 5) + values(i+1, 5))
           sumD = sumD + half*(values(i, 6) + values(i+1, 6))
        end do


        ! Now compute the normalized lift, drag and elliptical since
        ! we know the sum. Note delta is non-dimensional!
        delta = values(2, 1) - values(1, 1)

        ! This is the "nondimensional" span...it basically takes into account if you have
        ! a wing not at the sym plane
        span = maxval(values(:, 1)) - minval(values(:, 1)) 
        dmin = minval(values(:, 1))
        sumL = sumL * delta
        sumD = sumD * delta

        do i=1,d%nSegments 

           ! Normalized Lift
           values(i, 7) = values(i, 5) / sumL 
           values(i, 8) = values(i, 6) / sumD

           ! elliptical value
           values(i, 15) = four/pi/span*sqrt(one-(values(i, 1)-dmin)**2/span**2)
        end do

        ! Write all variables in block format
        do j=1,nLiftDistVar
           do i=1,d%nSegments
              write(fileID,15) values(i, j)
           end do
        end do

        ! Destroy the lift slices:
        do i=1,d%nSegments
           call destroySlice(d%slices(i))
        end do

        ! Deallocate slice list and point list
        deallocate(d%slices, d%slicePts)

        ! Destroy temp variables
        deallocate(liftDistNames, values)
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
  integer(kind=intType) :: nSurfVariables, nSliceVariables
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


     ! Now determine of cell centered values that will be
     ! communicated. It is equal to 6 plus the number of surface
     ! variables. The 6 are the three components of pressure and
     ! viscous forces. 
     call numberOfSurfSolVariables(nSurfVariables)

     nSliceVariables = nSurfVariables + 6

     ! Allocate a few additional arrays that will be required for the data transfer:
     allocate(localData(nSliceVariables, nCellsLocal))

     ! Step 1d. Determine the unique coordiantes and the new cell mapping:
     if (myid == 0) then
        ! Maximum space for the unique coordinates and link array
        allocate(uniqueNodes(3, nNodesTotal))
        allocate(link(nNodesTotal))

        call pointReduce(allNodes, nNodesTotal, tol, uniqueNodes, link, nUnique)

        ! Loop over the allCells and relpace the node index which
        ! currently refers to the indexing in 'allNodes' and replace with
        ! the ordering in 'uniqueNodes'
        do i=1,nCellsTotal
           do j=1,4
              allCells(j, i) = link(allCells(j, i))
           end do
        end do

        allocate(globalData(nSliceVariables, nCellsTotal))
        allocate(uniqueData(nSliceVariables, nUnique))
        allocate(fc(nUnique), dualAreas(nUnique))
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
        deallocate(globalData, uniqueNodes, uniqueData, fc, dualAreas)
     end if

     ! Free data for localCells and Nodes
     deallocate(localCells, localNodes)

     ! Free the data for sizing
     deallocate(nNodesProc, nCellsProc, cumNodesProc, cumCellsProc)

     ! Deallocate local force data
     deallocate(localData)

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
  integer(kind=intType) :: i, j, jj, ierr, nFields
  real(kind=realType), dimension(3) :: n1, n2, n3, n4, v1, v2, sss, qfp, qfv
  real(kind=realType) :: qa
  logical :: forcesTypeSave

  ! First get and communicate the coordinates

  call getForcePoints(localNodes, nNodesLocal, sps)

  call mpi_gatherv(&
       localNodes, nNodesLocal*3, sumb_real, &
       allNodes, nNodesProc*3, cumNodesProc*3, sumb_real, 0, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Now get the forces and any other surface data that we neded to
  ! communicate. Since all this data is cell centered we do a single
  ! reduction. 
  call computeSliceSurfaceData(sps, nFields)

  ! Gather pressure and viscous forces to root proc:
  call mpi_gatherv(&
       localData, nCellsLocal*nFields, sumb_real, &
       globalData, nCellsProc*nFields, cumCellsProc*nFields, sumb_real, 0, &
       sumb_comm_world, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Process the unique set of forces and tractions using the finite elemnent data:
  if (myid == 0) then
     
     ! Assign the unique nodes:
     do i=1,nNodesTotal
        uniqueNodes(:, link(i)) = allNodes(:, i)
     end do

     ! Zero the entire uniqueData array
     uniqueData = zero

     ! For the forces we a quarter of each cell and sum onto the unique nodes:
     do i=1,nCellsTotal
        qfp = globalData(ifxp:ifzp, i)*fourth
        qfv = globalData(ifxv:ifzv, i)*fourth
        do jj=1,4
           uniqueData(ifxp:ifzp, allCells(jj, i)) = uniqueData(ifxp:ifzp, allCells(jj, i)) + qfp
           uniqueData(ifxv:ifzv, allCells(jj, i)) = uniqueData(ifxv:ifzv, allCells(jj, i)) + qfv
        end do
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
        do j=ifxp,ifzv
           uniqueData(j, i) = uniqueData(j, i) / dualAreas(i)
        end do
     end do

     ! Finally we have the rest of data to average to the nodes. Add
     ! each cell value to the nodes and keep track of the number so we
     ! have an actual average
     dualAreas = zero
     do i=1,nCellsTotal
        do j=ifzv+1,nFields
           do jj=1,4
              uniqueData(j, allCells(jj, i)) = uniqueData(j, allCells(jj, i)) + globalData(j, i)
           end do
        end do
        do jj=1,4
           dualAreas(allCells(jj, i)) = dualAreas(allCells(jj, i)) + one
        end do
     end do

     ! Finally divide by the node count which was stored in dualAreas
     do i=1,nUnique
        do j=ifzv+1,nFields
           uniqueData(j, i) = uniqueData(j, i) / dualAreas(i)
        end do
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

  ! Set the info for the slic:
  slc%pt = pt
  slc%dir = dir

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

  ! Estimate size of slice by the 5 times sqrt of the number of nodes in the
  ! mesh
  nMax = int(sqrt(dble(nNodesTotal)))*5
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
  real(kind=realType) :: len, dmax, dmin, dist, fact, scaleDim, M(3,3)
  real(kind=realType) :: r(3), r_new(3), hyp, te(3), le(3), theta
  integer(kind=intType) :: bestPair(2), dir_ind
  real(kind=realtype), dimension(:,:), allocatable :: tempCoords
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
     pT1 = slc%w(1,i1)*uniqueData(ifxp:ifzp, slc%ind(1,i1)) + slc%w(2,i1)*uniqueData(ifxp:ifzp, slc%ind(2, i1))
     pT2 = slc%w(1,i2)*uniqueData(ifxp:ifzp, slc%ind(1,i2)) + slc%w(2,i2)*uniqueData(ifxp:ifzp, slc%ind(2, i2))

     ! extract viscous tractions
     vT1 = slc%w(1,i1)*uniqueData(ifxv:ifzv, slc%ind(1,i1)) + slc%w(2,i1)*uniqueData(ifxv:ifzv, slc%ind(2, i1))
     vT2 = slc%w(1,i2)*uniqueData(ifxv:ifzv, slc%ind(1,i2)) + slc%w(2,i2)*uniqueData(ifxv:ifzv, slc%ind(2, i2))

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

  ! Compute the chord as the max length between any two nodes...this
  ! is n^2, so should be changed in the future

  dmax = zero
  bestPair = (/1, 1/)
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

  ! Set chord, protected from zero
  slc%chord = max(dmax, 1e-12)

  ! Compute factor to get coefficient
  scaleDim = pRef/pInf
  fact = two/(gammaInf*pInf*MachCoef*MachCoef*scaleDim)

  ! Take dmax as chord and compute coefficients
  slc%CLp = slc%pL / slc%chord * fact
  slc%CDp = slc%pD / slc%chord * fact
  slc%CLv = slc%vL / slc%chord * fact
  slc%CDv = slc%vD / slc%chord * fact

  ! Default values
  slc%twist = zero
  slc%thickness = zero

  if (slc%nNodes == 0) then
     return
  end if

  ! Lastly we need the twist and the twist and the thickness
  i1 = bestPair(1)
  i2 = bestPair(2)
  x1 = slc%w(1,i1)*uniqueNodes(:, slc%ind(1,i1)) + slc%w(2,i1)*uniqueNodes(:, slc%ind(2, i1))
  x2 = slc%w(1,i2)*uniqueNodes(:, slc%ind(1,i2)) + slc%w(2,i2)*uniqueNodes(:, slc%ind(2, i2))

  if (x1(1) > x2(1)) then
     te = x1
     le = x2
  else
     te = x2
     le = x1
  end if

  ! Length of hyptoneuse is the same
  hyp = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

  if (dir_ind == 1) then
     ! Xslice...we don't how what to do here..could be y or z. Don't
     ! do anything.
     slc%twist = zero
  else if (dir_ind == 2) then
     ! Yslice
     theta = asin((le(3)-te(3))/hyp)
     slc%twist = theta*180.0/pi
  else
     ! Zslice
     theta = asin((le(2)-te(2))/hyp)
     slc%twist = theta*180.0/pi
  end if

  ! Finally we need to get the thickness. For this, compute temporary
  ! section nodes and rotate them by the twist values we just computed
  ! and take the max and min

  ! Back out what is the main index of the slice, x, y or z based on
  ! the direction. Not the best approach, but that's ok
  dir_ind = maxloc(abs(slc%dir),1)

  if (dir_ind == 1) then
     M(1,1) = one; M(1,2) = zero; M(1, 3) = zero;
     M(2,1) = zero; M(2,2) = one; M(2, 3) = zero;
     M(3,1) = zero; M(3,2) = zero; M(3,3) = one;
  else if(dir_ind == 2) then
     ! Y-rotation matrix
     M(1,1) = cos(-theta); M(1,2) = zero; M(1, 3) = sin(-theta);
     M(2,1) = zero; M(2,2) = one; M(2, 3) = zero;
     M(3,1) = -sin(-theta); M(3,2) = zero; M(3,3) = cos(-theta);
  else
     ! Z rotation Matrix
     M(1,1) = cos(theta); M(1,2) = -sin(theta); M(1,3) = zero;
     M(2,1) = sin(theta); M(2,2) = cos(theta); M(2,3) = zero;
     M(3,1) = zero; M(3,2) = zero; M(3,3) = one;
  end if

  allocate(tempCoords(3, slc%nNodes))
  do i=1,slc%nNodes
     ! extract node:
     r = slc%w(1,i)*uniqueNodes(:, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(:, slc%ind(2, i)) - te
     r_new = matmul(M, r)
     tempCoords(:, i) = r_new + te
  end do

  ! Now get the max and the min and divide by the chord for t/c
  if (dir_ind == 1) then
     slc%thickness = 0 ! Again, don't know what to do here
  else if(dir_ind == 2) then
     dmax = maxval(tempCoords(3, :))
     dmin = minval(tempCoords(3, :))
     slc%thickness = (dmax-dmin)/hyp
  else if(dir_ind == 3) then
     dmax = maxval(tempCoords(2, :))
     dmin = minval(tempCoords(2, :))
     slc%thickness = (dmax-dmin)/hyp
  end if
  deallocate(tempCoords)

end subroutine integrateSlice

subroutine writeSlice(slc, fileID, nFields)
  ! Write the data in slice 'slc' to openfile ID fileID

  use liftDistributionData
  implicit none

  ! Input Parameters
  type(slice), intent(in) :: slc
  integer(kind=intType), intent(in) :: fileID, nFields

  ! Working Variables
  integer(kind=intType) :: i, j

  write (fileID,"(a,a,a)") "Zone T= """,trim(slc%sliceName),""""

  ! IF we have nodes actually write:
  if (slc%nNodes > 0) then
     write (fileID,*) "Nodes = ", slc%nNodes, " Elements= ", slc%nNodes/2, " ZONETYPE=FELINESEG"
     write (fileID,*) "DATAPACKING=POINT"
13   format (E14.6)

     do i=1,slc%nNodes
        ! Write the coordinates
        do j=1,3
           write(fileID,13, advance='no') &
                slc%w(1,i)*uniqueNodes(j, slc%ind(1,i)) + slc%w(2,i)*uniqueNodes(j, slc%ind(2, i))
        end do

        ! Write field data
        do j=ifzv+1,nFields
           write(fileID,13, advance='no') &
                slc%w(1,i)*uniqueData(j, slc%ind(1,i)) + slc%w(2,i)*uniqueData(j, slc%ind(2, i))
        end do
        write(fileID,"(1x)")
     end do

15   format(I5, I5)
     do i=1, slc%nNodes/2
        write(fileID, 15) 2*i-1, 2*i
     end do
  else ! Write dummy data so the number of zones are the same

     write (fileID,*) "Nodes = ", 2, " Elements= ", 1, " ZONETYPE=FELINESEG"
     write (fileID,*) "DATAPACKING=POINT"
     do i=1,2
        do j=1,3
           write(fileID,13, advance='no') zero
        end do

        do j=ifzv+1,nFields
           write(fileID,13, advance='no') zero
        end do
        write(fileID,"(1x)")
     end do
     write(fileID, 15) 1, 2
  end if
end subroutine writeSlice
