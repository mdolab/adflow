subroutine addParaSlice(sliceName, pt, direction, famList, n)
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
  use surfaceFamilies
  use inputTimeSpectral 
  implicit none

  ! Input parameters
  character*(*), intent(in) :: sliceName
  real(kind=realType), dimension(3), intent(in) :: pt, direction
  integer(kind=intType), intent(in) :: famList(n), n

  ! Working
  integer(kind=intType) :: sps
  interface
     subroutine createSlice(exch, slc, pt, dir, sliceName, famList, nFam)
       use surfaceFamilies
       use liftdistributiondata
       type(familyExchange), target :: exch
       type(slice), intent(inout) :: slc
       real(kind=realType), dimension(3), intent(in) :: pt, dir
       character*(*), intent(in) :: sliceName
       integer(kind=intType), intent(in) :: famList(nFam), nFam
     end subroutine createSlice
  end interface

  if (.not. allocated(paraSlices)) then 
     allocate(paraSlices(nSliceMax, nTimeIntervalsSpectral))
  end if

  ! We have to add a slice for each spectral instance.
  do sps=1, nTimeIntervalsSpectral
     nParaSlices = nParaSlices + 1

     if (nParaSlices > nSliceMax) then 
        print *,'Error: Exceeded the maximum number of slices. Increase nSliceMax'
        stop
     end if

     call createSlice(wallExchange(sps), paraSlices(nParaSlices, sps), pt, &
          direction, sliceName, famList, n)
  end do

end subroutine addParaSlice

subroutine addAbsSlice(sliceName, pt, direction, famList, n)
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
  use surfaceFamilies
  use inputTimeSpectral
  use surfaceFamilies
  implicit none

  ! Input parameters
  character*(*), intent(in) :: sliceName
  real(kind=realType), dimension(3), intent(in) :: pt, direction
  integer(kind=intType), intent(in) :: famList(n), n

  ! Working
  integer(kind=intType) :: sps
  interface
     subroutine createSlice(exch, slc, pt, dir, sliceName, famList, nFam)
       use surfaceFamilies
       use liftdistributiondata
       type(familyExchange), target :: exch
       type(slice), intent(inout) :: slc
       real(kind=realType), dimension(3), intent(in) :: pt, dir
       character*(*), intent(in) :: sliceName
       integer(kind=intType), intent(in) :: famList(nFam), nFam
     end subroutine createSlice
  end interface

  if (.not. allocated(absSlices)) then 
     allocate(absSlices(nSliceMax, nTimeIntervalsSpectral))
  end if

  do sps=1, nTimeIntervalsSpectral
     nAbsSlices = nAbsSlices + 1

     if (nAbsSlices > nSliceMax) then 
        print *,'Error: Exceeded the maximum number of slices. Increase nSliceMax'
        stop
     end if

     call createSlice(wallExchange(sps), absSlices(nAbsSlices, sps), pt, &
          direction, sliceName, famList, n)
  end do

end subroutine addAbsSlice

subroutine addLiftDistribution(nSegments, dir_vec, dir_ind, distName, famList, n)
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
  use surfaceFamilies
  implicit none

  ! Input parameters
  character*(*), intent(in) :: distName
  integer(kind=intType), intent(in) :: nSegments
  real(kind=realType), dimension(3) :: dir_vec
  integer(kind=intType), intent(in) :: dir_ind
  integer(kind=intType), intent(in) :: famList(n), n

  nLiftDists = nLiftDists + 1
  if (nLiftDists > nLiftDistMax) then 
     print *,'Error: Exceeded the maximum number of lift distributions. &
          &Increase nLiftDistMax'
     stop
  end if

  liftDists(nLIftDists)%nSegments = nSegments
  liftDists(nLiftDists)%dir = dir_vec
  liftDists(nLiftDists)%distName = distName
  liftDists(nLIftDists)%dir_ind = dir_ind

  allocate(liftDists(nLiftDists)%famList(n))
  liftDists(nLiftDists)%famList(:) = famList

end subroutine addLiftDistribution

subroutine writeTecplot(sliceFile, writeSlices, liftFile, writeLift, surfFile, writeSurf)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This is the master routine for writing tecplot data from sumb. *
  !      *                                                                *
  !      * This routine will write the slice, lift and surface files      *
  !      * depending on the flags writeSlics, writeLift and writeSurface. *
  !      * The reason for the combined routine is that we can safely only *
  !      * perform the nodal averaging once which is required for all     *
  !      * three output files.                                            *
  !      ******************************************************************

  use inputTimeSpectral
  use surfaceFamilies
  use communication
  implicit none

  ! Input Params
  character*(*), intent(in) :: sliceFile, liftFile, surfFile
  logical, intent(in) :: writeSlices, writeLift, writeSurf

  ! Working
  integer(kind=intType) :: sps
  real(kind=realType) :: time1
  do sps=1, nTimeIntervalsSpectral
     call computeSurfaceOutputNodalData(wallExchange(sps), .True.)
  end do

  if (writeSlices) then 
     call writeSlicesFile(sliceFile, .False.)
  end if

  if (writeLift) then 
     call writeLiftDistributionFile(liftFile, .False.)
  end if

  if (writeSurf) then 
     call writeTecplotSurfaceFile(surfFile, .False.)
  end if

end subroutine writeTecplot

subroutine writeSlicesFile(fileName, updateSurfaceData)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine is intended to be called from python.          *
  !      *                                                                *
  !      * This routine will write the user defined slics to an           *
  !      * to the (ascii) tecplot file fileName. ASCII files are          *
  !      * used for simplicity since very little information is actually  *
  !      * written.                                                       *
  !      *                                                                *
  !      ******************************************************************
  use communication
  use liftDistributionData
  use outputMod
  use inputTimeSpectral
  use inputPhysics
  use inputIteration
  use inputIO
  use surfaceFamilies
  implicit none

  ! Input Params
  character*(*), intent(in) :: fileName
  logical, intent(in) :: updateSurfaceData

  ! Working parameters
  integer(kind=intType) :: file, i, sps, nSolVar, ierr
  character(len=maxStringLen) :: fname
  character(len=7) :: intString
  character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
  integer(kind=intType), allocatable, dimension(:) :: famList
  type(slice) :: globalSlice
  interface
     subroutine createSlice(exch, slc, pt, dir, sliceName, famList, nFam)
       use surfaceFamilies
       use liftdistributiondata
       type(familyExchange), target :: exch
       type(slice), intent(inout) :: slc
       real(kind=realType), dimension(3), intent(in) :: pt, dir
       character*(*), intent(in) :: sliceName
       integer(kind=intType), intent(in) :: famList(nFam), nFam
     end subroutine createSlice
  end interface

  ! Only write if we actually have lift distributions
  testwriteSlices: if(nParaSlices + nAbsSlices > 0) then

     if(myID == 0 .and. printIterations) then
        print "(a)", "#"
        print "(a)", "# Writing slices file(s) ..."
     endif
     call setWallFamilyList(sps)
     do sps=1,nTimeIntervalsSpectral

        ! Gather the forces and nodes
        if (updateSurfaceData) then 
           call computeSurfaceOutputNodalData(wallExchange(sps), .True.)
        end if

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
           write(file,"(a)",advance="no") " ""XoC"" "
           write(file,"(a)",advance="no") " ""YoC"" "
           write(file,"(a)",advance="no") " ""ZoC"" "

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
        end if
        call mpi_bcast(nSolVar, 1, sumb_integer, 0, sumb_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)

        ! Integration is performed in parallel
        do i=1,nParaSlices
           call integrateSlice(paraSlices(i, sps), globalSlice, nSolVar, .True.)
           if (myid == 0) then 
              call writeSlice(globalSlice, file, nSolVar)
           end if
           call destroySlice(globalSlice)
        end do

        do i=1,nAbsSlices
           ! 'Destroy' the slice...just dealloc the allocated data. 
           ! before we do, save the family list
           allocate(famList(size(absSlices(i, sps)%famList)))
           famList = absSlices(i, sps)%famList
           call destroySlice(absSlices(i, sps))

           ! Make new one in the same location
           call createSlice(wallExchange(sps), absSlices(i, sps), absSlices(i, sps)%pt, &
                absSlices(i, sps)%dir, absSlices(i, sps)%sliceName, famList, size(famList))

           call integrateSlice(absSlices(i, sps), globalSlice, nSolVar, .True.)
           if (myid == 0) then 
              call writeSlice(globalSlice, file, nSolVar)
           end if
           call destroySlice(globalSlice)
           deallocate(famList)
        end do

        !Close file on root proc
        if (myid == 0) then 
           close(file)
        end if
     end do

     if(myID == 0 .and. printIterations) then
        print "(a)", "# Slices file(s) written"
        print "(a)", "#"
     endif
  end if testwriteSlices
end subroutine writeSlicesFile

subroutine writeLiftDistributionFile(fileName, updateSurfaceData)
  !
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
  use surfaceFamilies
  implicit none

  ! Input Params
  character*(*), intent(in) :: fileName
  logical :: updateSurfaceData

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

        if (updateSurfaceData) then 
           call computeSurfaceOutputNodalData(wallExchange(sps), .False.)
        end if

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
  use surfaceFamilies
  implicit none

  ! Input parameters
  integer(kind=intType), intent(in) :: sps, fileID

  real(kind=realType), dimension(3) :: xmin, xmax, xmin_local, xmax_local
  real(kind=realType), parameter :: tol=1e-8
  type(liftDist), pointer :: d
  integer(kind=intType) :: i, j, ii, jj, iDist, ierr, bsearchIntegers
  real(kind=realType), dimension(:,:), allocatable :: values
  character(len=maxCGNSNameLen), dimension(:), allocatable :: liftDistNames
  real(kind=realType) :: dmin, dmax, sumL, sumD, span, delta, xCur(3)
  type(slice) :: localSlice, globalSlice
  type(familyExchange), pointer :: exch
  interface
     subroutine createSlice(exch, slc, pt, dir, sliceName, famList, nFam)
       use surfaceFamilies
       use liftdistributiondata
       type(familyExchange), target :: exch
       type(slice), intent(inout) :: slc
       real(kind=realType), dimension(3), intent(in) :: pt, dir
       character*(*), intent(in) :: sliceName
       integer(kind=intType), intent(in) :: famList(nFam), nFam
     end subroutine createSlice
  end interface

  do iDist=1,nLiftDists
     d => liftDists(iDist)
     xmin_local = huge(real(zero))
     xmax_local = -huge(real(zero))
     exch => wallExchange(sps)
     ! Get the bounding box for the entire geometry we have been slicing.
     elemLoop: do i=1, size(exch%conn, 2)
        if (bSearchIntegers(exch%elemFam(i), d%FamList, size(d%famList)) > 0) then 

           ! Extract each of the 4 nodes on this quad:
           do jj=1,4
              xCur = exch%nodalValues(exch%conn(jj, i), 1:3)
              ! Check the max/min on each index
              do ii=1,3
                 xmin_local(ii) = min(xmin_local(ii) , xCur(ii))
                 xmax_local(ii) = max(xmax_local(ii) , xCur(ii))
              end do
           end do
        end if
     end do elemLoop

     ! Globalize all min/max values. 
     call mpi_allreduce(xmin_local, xmin, 3, sumb_real, MPI_MIN, &
          sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call mpi_allreduce(xmax_local, xmax, 3, sumb_real, MPI_MAX, &
          sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     d%delta = (xMax(d%dir_ind) - xMin(d%dir_ind))/dble((d%nSegments - 1))
     allocate(d%slicePts(3, d%nSegments))

     ! Zero out all segments
     d%slicePts = zero

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
     if (myid == 0) then 
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
     end if

     allocate(values(d%nSegments, nLiftDistVar))
     values = zero

     do i=1,d%nSegments
        if (i==1) then
           d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + tol
        else if (i == d%nSegments) then 
           d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + (i-1)*d%delta - tol
        else
           d%slicePts(d%dir_ind, i) = xMin(d%dir_ind) + (i-1)*d%delta 
        end if
     end do

     ! Scaled Eta values
     dmin = minVal(d%slicePts(d%dir_ind, :))
     dmax = maxval(d%slicePts(d%dir_ind, :))

     values(:, 1) = (d%slicePts(d%dir_ind, :) - dmin)/(dmax-dmin)
     ! Coordinate Varaibles
     if (d%dir_ind == 1) then! X slices
        values(:, 2) = d%slicePts(1, :)
     else if (d%dir_ind == 2) then ! Y slices
        values(:, 3) = d%slicePts(2, :)
     else if (d%dir_ind == 3) then ! Z slices
        values(:, 4) = d%slicePts(3, :)
     end if

     do i=1, d%nSegments
        ! Make new one in the same location
        call createSlice(exch, localSlice, d%slicePts(:, i), &
             d%dir, "does_not_matter", d%famList, size(d%famList))
        call integrateSlice(localSlice, globalSlice, 0, .False.)

        ! Total lift and drag 
        values(i, 5)  = globalSlice%pL + globalSlice%vL
        values(i, 6)  = globalSlice%pD + globalSlice%vD

        ! Total CL and CD
        values(i, 9)  = globalSlice%CLp + globalSlice%CLv
        values(i, 10) = globalSlice%CDp + globalSlice%CDv

        ! Pressure lift and drag coefficients
        values(i, 11) = globalSlice%CLp
        values(i, 12) = globalSlice%CDp

        ! Viscous lift and drag coefficients
        values(i, 13) = globalSlice%CLv
        values(i, 14) = globalSlice%CDv

        ! t/c, twist, chord
        values(i, 16) = globalSlice%thickness
        values(i, 17) = globalSlice%twist
        values(i, 18) = globalSlice%chord

        call destroySlice(localSlice)
        call destroySlice(globalSlice)

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
     if (myid == 0) then 
        do j=1,nLiftDistVar
           do i=1,d%nSegments
              write(fileID,15) values(i, j)
           end do
        end do
     end if

     ! Deallocate slice list and point list
     deallocate(d%slicePts)

     ! Destroy temp variables
     deallocate(liftDistNames, values)
  end do

end subroutine writeLiftDistributions

subroutine writeTecplotSurfaceFile(fileName, updateSurfaceData)
  use communication
  use liftDistributionData
  use outputMod
  use inputTimeSpectral
  use inputPhysics
  use inputIteration
  use inputIO
  use surfaceFamilies
  implicit none

  ! Input Params
  character*(*), intent(in) :: fileName
  logical, intent(in) :: updateSurfaceData

  ! Working parameters
  integer(kind=intType) :: i, j, iProc, sps, nSolVar, ierr, fileID, iSize, iFam
  character(len=maxStringLen) :: fname
  character(len=7) :: intString
  integer(Kind=intType) :: nNodes, nCells, bSearchIntegers
  integer(kind=intType), dimension(:), allocatable :: nodeSizes, nodeDisps
  integer(kind=intType), dimension(:), allocatable :: cellSizes, cellDisps
  character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
  real(kind=realType), dimension(:, :), allocatable :: vars
  integer(kind=intType), dimension(:, :), allocatable :: conn
  integer(kind=intType), dimension(:), allocatable :: mask, elemFam

  if(myID == 0 .and. printIterations) then
     print "(a)", "#"
     print "(a)", "# Writing tecplot surface file(s) ..."
  endif

  spectralLoop: do sps=1,nTimeIntervalsSpectral

     ! Gather the forces and nodes
     if (updateSurfaceData) then 
        call computeSurfaceOutputNodalData(wallExchange(sps), .True.)
     end if

     ! If it is time spectral we need to agument the filename
     if (equationMode == timeSpectral) then
        write(intString,"(i7)") sps
        intString = adjustl(intString)
        fname = trim(fileName)//"Spectral"//trim(intString)
     else
        fname = fileName
     end if

     fileID = 11
     ! Open file on root proc:
     if (myid == 0) then 
        open(unit=fileID, file=trim(fname))

        ! Write Header Information 
        write (fileID,*) "Title = ""SUmb Surface Data"""
        write (fileID,"(a)", advance="no") "Variables = "
        write(fileID,"(a)",advance="no") " ""CoordinateX"" "
        write(fileID,"(a)",advance="no") " ""CoordinateY"" "
        write(fileID,"(a)",advance="no") " ""CoordinateZ"" "

        ! Number of surface variables
        call numberOfSurfSolVariables(nSolVar)
        allocate(solNames(nSolVar))
        call surfSolNames(solNames)

        ! Write the rest of the variables
        do i=1,nSolVar
           write(fileID,"(a,a,a)",advance="no") """",trim(solNames(i)),""" "
        end do

        write(fileID,"(1x)")
        deallocate(solNames)
     end if

     call mpi_bcast(nSolVar, 1, sumb_integer, 0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! ================================================================
     !                        Wall Writing
     ! ================================================================

     ! Firstly we write all the walls. We have already computed the
     ! data for those. For now, just dump all the nodal data on each
     ! zone.  We should be smart and do the fancy variable sharing
     ! stuff, but that can be done later. We also have to gather up
     ! the type of each element corresponding to conn.

     ! Gather up the number of nodes to be set to the root proc:
     allocate(nodeSizes(nProc), nodeDisps(0:nProc))
     nodeSizes = 0
     nodeDisps = 0
     call mpi_allgather(wallExchange(sps)%nNodes, 1, sumb_integer, nodeSizes, 1, sumb_integer, &
          sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     nodeDisps(0) = 0
     do iProc=1, nProc
        nodeDisps(iProc) = nodeDisps(iProc-1) + nodeSizes(iProc)
     end do

     iSize = 3 + 6 + nSolVar
     if (myid == 0) then 
        nNodes = sum(nodeSizes)
        allocate(vars(nNodes, iSIze))
     end if

     ! Gather values to the root proc.
     do i=1, iSize
        call mpi_gatherv(wallExchange(sps)%nodalValues(:, i), wallExchange(sps)%nNodes, &
             sumb_real, vars(:, i), nodeSizes, nodeDisps, sumb_real, 0, sumb_comm_world, ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end do
     ! Now gather up the connectivity
     allocate(cellDisps(0:nProc), cellSizes(nProc))

     call mpi_gather(size(wallExchange(sps)%conn, 2), 1, sumb_integer, &
          cellSizes, 1, sumb_integer, 0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     if (myid == 0) then 
        cellDisps(0) = 0
        do iProc=1, nProc
           cellDisps(iProc) = cellDisps(iProc-1) + cellSizes(iProc)
        end do
        nCells = sum(cellSizes)
        allocate(conn(4, nCells))
        allocate(elemFam(nCells))
     end if

     ! We offset the conn array by nodeDisps(iProc) which
     ! automagically adjusts the connectivity to account for the
     ! number of nodes from different processors

     call mpi_gatherv(wallExchange(sps)%conn+nodeDisps(myid), &
          4*size(wallExchange(sps)%conn, 2), sumb_integer, conn, &
          cellSizes*4, cellDisps*4, sumb_integer, 0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call mpi_gatherv(wallExchange(sps)%elemFam, &
          size(wallExchange(sps)%elemFam), sumb_integer, elemFam, &
          cellSizes, cellDisps, sumb_integer, 0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)


     if (myid == 0) then 
        allocate(mask(nCells))
        do iFam=1, totalWallFamilies

           ! Create a temporary mask
           mask = 0
           do i=1, nCells
              ! Check if this elem is to be included
              if (elemFam(i) == wallFamilies(iFam)) then 
                 mask(i) = 1
              end if
           end do

           write (fileID,"(a,a,a)") "Zone T= """,famNames(wallFamilies(iFam)),""""
           write (fileID,*) "Nodes = ", nNodes, " Elements= ",  sum(mask), " ZONETYPE=FEQUADRILATERAL"
           write (fileID,*) "DATAPACKING=BLOCK"

           ! Write the points (indices 1:3)
           do j=1, 3
              do i=1, nNOdes
                 write(fileID,15) vars(i, j)
              end do
           end do

           ! And the rest of the variables
           do j=1,nSolVar
              do i=1, nNodes
                 write(fileID,15) vars(i, j+9)
              end do
           end do

           ! And now the connectivity
           do i=1, nCells
              ! Check if this elem is to be included
              if (mask(i) == 1) then 
                 write(fileID, 16) conn(1, i), conn(2, i), conn(3,i), conn(4, i)
              end if
           end do
        end do
        deallocate(mask, elemFam, vars , conn)
     end if

     deallocate(cellSizes, cellDisps, nodeSizes, nodeDisps)

     ! ================================================================
     !                     All other surfaces
     ! ================================================================

     do iFam=1, totalFamilies

        notAWall: if (.not. (bsearchIntegers(iFam, wallFamilies, totalWallFamilies)>0)) then 
           call computeSurfaceOutputNodalData(familyExchanges(iFam, sps), .False.)

           ! Gather up the number of nodes to be set to the root proc:
           allocate(nodeSizes(nProc), nodeDisps(0:nProc))
           nodeSizes = 0
           nodeDisps = 0
           call mpi_allgather(familyExchanges(iFam, sps)%nNodes, 1, sumb_integer, nodeSizes, 1, sumb_integer, &
                sumb_comm_world, ierr)
           call EChk(ierr,__FILE__,__LINE__)
           nodeDisps(0) = 0
           do iProc=1, nProc
              nodeDisps(iProc) = nodeDisps(iProc-1) + nodeSizes(iProc)
           end do

           iSize = 3 + 6 + nSolVar
           if (myid == 0) then 
              nNodes = sum(nodeSizes)
              allocate(vars(nNodes, iSIze))
           end if

           ! Gather values to the root proc.
           do i=1, iSize
              call mpi_gatherv(familyExchanges(iFam, sps)%nodalValues(:, i), familyExchanges(iFam, sps)%nNodes, &
                   sumb_real, vars(:, i), nodeSizes, nodeDisps, sumb_real, 0, sumb_comm_world, ierr)
              call EChk(ierr,__FILE__,__LINE__)
           end do
           ! Now gather up the connectivity
           allocate(cellDisps(0:nProc), cellSizes(nProc))

           call mpi_gather(size(familyExchanges(iFam, sps)%conn, 2), 1, sumb_integer, &
                cellSizes, 1, sumb_integer, 0, sumb_comm_world, ierr)
           call EChk(ierr,__FILE__,__LINE__)

           if (myid == 0) then 
              cellDisps(0) = 0
              do iProc=1, nProc
                 cellDisps(iProc) = cellDisps(iProc-1) + cellSizes(iProc)
              end do
              nCells = sum(cellSizes)
              allocate(conn(4, nCells))
           end if

           ! We offset the conn array by nodeDisps(iProc) which
           ! automagically adjusts the connectivity to account for the
           ! number of nodes from different processors

           call mpi_gatherv(familyExchanges(iFam, sps)%conn+nodeDisps(myid), &
                4*size(familyExchanges(iFam, sps)%conn, 2), sumb_integer, conn, &
                cellSizes*4, cellDisps*4, sumb_integer, 0, sumb_comm_world, ierr)
           call EChk(ierr,__FILE__,__LINE__)

           ! Not quite finished yet since we will have gathered nodes from
           ! multiple procs we have to adjust the connectivity

           if (myid == 0) then 
              write (fileID,"(a,a,a)") "Zone T= """,trim(famNames(iFam)),""""
              write (fileID,*) "Nodes = ", nNodes, " Elements= ",  nCells, " ZONETYPE=FEQUADRILATERAL"
              write (fileID,*) "DATAPACKING=BLOCK"

              ! Write the points (indices 1:3)
15            format (E14.6)

              do j=1, 3
                 do i=1, nNOdes
                    write(fileID,15) vars(i, j)
                 end do
              end do

              ! And the rest of the variables
              do j=1,nSolVar
                 do i=1, nNodes
                    write(fileID,15) vars(i, j+9)
                 end do
              end do

              ! And now the connectivity
16            format ((I6) (I6) (I6) (I6))
              do i=1, nCells
                 write(fileID, 16) conn(1, i), conn(2, i), conn(3,i), conn(4, i)
              end do
              deallocate(conn, vars)
           end if
           deallocate(cellSizes, cellDisps, nodeSizes, nodeDisps)
        end if notAWall
     end do

     if (myid == 0) then 
        close(fileID)
     end if
  end do spectralLoop

  if(myID == 0 .and. printIterations) then
     print "(a)", "# Tecplot surface file(s) written"
     print "(a)", "#"
  endif

end subroutine writeTecplotSurfaceFile

subroutine initializeLiftDistributionData

  use communication
  use blockPointers
  use inputPhysics
  use outputMod
  use liftDistributionData
  use inputTimeSpectral
  use surfacefamilies
  implicit none

  ! Working Variables
  integer(kind=intType) :: nPts, nCells, sps

  if (liftDistInitialized) then
     return
  else

     ! Data for the marching squares method: Which edges are cut by
     ! the contour. We haven't dealt with the case of an ambiguous
     ! contour which is case 6 and 11. Since most of the slices we are
     ! doing are planes this won't matter.
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

subroutine computeSurfaceOutputNodalData(exch, includeTractions)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This purpose of this subroutine is to compute all nodal values *
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  use communication
  use inputPhysics
  use blockPointers
  use BCTypes
  use surfaceFamilies, only: famGroups
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#else
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#endif

  ! Input Param  
  type(familyExchange) :: exch
  logical :: includeTractions

  ! Working params
  integer(kind=intType) :: i, j, ii, jj, kk, nn, mm, iSol, ierr, nPts, nCells
  integer(kind=intType) :: nFields, nSolVar, iBeg, iEnd, jBeg, jEnd, ind(4), ni, nj
  integer(kind=intType) :: bSearchIntegers
  integer(kind=intType), dimension(3,2) :: cellRangeCGNS
  character(len=maxCGNSNameLen), dimension(:), allocatable :: solNames
  real(kind=realType), dimension(:), allocatable :: buffer
  real(kind=realType), dimension(:), pointer :: weightPtr, localPtr
  real(kind=realType), dimension(:, :), allocatable :: tmp
  logical :: viscousSubFace


  ! Determine the number of surface variables we have
  call numberOfSurfSolVariables(nSolVar)
  allocate(solNames(nSolVar))
  call surfSolNames(solNames)

  ! Reallocate the nodal values if necessary and zero
  if (allocated(exch%nodalValues)) then 
     deallocate(exch%nodalValues)
  end if
  allocate(exch%nodalValues(exch%nNodes, nSolVar+6+3))
  exch%nodalValues = zero

  ! Set the main family group to the family this exchange. 
  call setFamilyInfo(exch%famGroups, exch%nFam)

  ! The tractions have a sepcial routine so call that first before we
  ! mess with the family information. 

  if (includeTractions) then 

     ! And put the traction where it needs to go if necessary.  Note that
     ! this computation will only work corectly when the wallExchange
     ! which include all wallgroups is used. 

     call computeNodalTractions(exch%sps)
     ii = 0
     do nn=1, nDom
        call setPointers(nn, 1_intType, exch%sps)
        do mm=1, nBocos
           iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
           jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd

           if (bsearchIntegers(BCdata(mm)%famID, &
                famGroups, shape(famGroups)) > 0) then 
              do j=jBeg, jEnd
                 do i=iBeg, iEnd
                    ii = ii + 1
                    exch%nodalValues(ii, 4) = BCData(mm)%Tp(i, j, 1)
                    exch%nodalValues(ii, 5) = BCData(mm)%Tp(i, j, 2)
                    exch%nodalValues(ii, 6) = BCData(mm)%Tp(i, j, 3)

                    exch%nodalValues(ii, 7) = BCData(mm)%Tv(i, j, 1)
                    exch%nodalValues(ii, 8) = BCData(mm)%Tv(i, j, 2)
                    exch%nodalValues(ii, 9) = BCData(mm)%Tv(i, j, 3)
                 end do
              end do
           end if
        end do
     end do
  end if

  ! Get the current set of surface points for the family we just set. 
  allocate(tmp(3, exch%nNodes))
  call getSurfacePoints(tmp, exch%nNodes, exch%sps)

  do i=1, exch%nNodes
     exch%nodalValues(i, 1:3) = tmp(1:3, i)
  end do
  deallocate(tmp)
  ! For the remainder of the variables, use arithematic averaging. 
  call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = zero
  ! ii is the index into the pointer array
  ii = 0
  do nn=1, nDom
     call setPointers(nn, 1_intType, exch%sps)
     do mm=1, nBocos
        iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
        jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
        ni = iEnd - iBeg + 1
        nj = jEnd - jBeg + 1
        if (bsearchIntegers(BCdata(mm)%famID, &
             famGroups, shape(famGroups)) > 0) then 
           do j=0,nj-2
              do i=0,ni-2

                 ! Scatter 1 to each node. 
                 ind(1) = ii + (j  )*ni + i + 1
                 ind(2) = ii + (j  )*ni + i + 2 
                 ind(3) = ii + (j+1)*ni + i + 2 
                 ind(4) = ii + (j+1)*ni + i + 1
                 do jj=1,4
                    localPtr(ind(jj)) = localPtr(ind(jj)) + one
                 end do
              end do
           end do
           ii = ii + ni*nj
        end if
     end do
  end do

  call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Globalize the area
  call vecSet(exch%sumGlobal, zero, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
       exch%sumGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Now compute the inverse of the weighting so that we can multiply
  ! instead of dividing.

  call vecGetArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  localPtr = one/localPtr

  call vecRestoreArrayF90(exch%sumGlobal, localPtr, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  varLoop: do iSol=1, nSolVar

     ! Extract the poitner to the local array
     call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! ii is the continous running element pointer
     ii = 0 
     localPtr = zero
     ! Do each variable separately. 
     domainLoop: do nn=1, nDom
        call setPointers(nn, 1, exch%sps)

        bocoLoop: do mm=1, nBocos
           if (bsearchIntegers(BCdata(mm)%famID, &
                famGroups, shape(famGroups)) > 0) then 

              ! storeSurfSolInBuffer needs to know if the subface is
              ! viscous or not. 
              viscousSubFace = .False.
              if (BCType(mm) == NSWallAdiabatic .or. &
                   BCType(mm) == NSWallIsoThermal) then 
                 viscousSubFace = .true. 
              end if

              ! Determine the cell range *in the original cgns
              ! ordering*. The reason for this is the storeSurfSolInBufer
              ! is normally used for writing CGNS and thus it is that
              ! ording that is important. However, that isn't too hard to
              ! deal with. 

              jBeg = BCData(mm)%jnBeg + 1
              jEnd = BCData(mm)%jnEnd
              iBeg = BCData(mm)%inBeg + 1
              iEnd = BCData(mm)%inEnd

              ! Dummy value of 1 for the face values not set.
              cellRangeCGNS = 1
              select case(BCFaceID(mm))
              case (iMin, iMax)
                 ! I range meaningless
                 cellRangeCGNS(2, 1) = iBeg + jBegOr - 1 
                 cellRangeCGNS(2, 2) = iEnd + jBegOr - 1 

                 cellRangeCGNS(3, 1) = jBeg + kBegOr - 1
                 cellRangeCGNS(3, 2) = jEnd + kBegOr - 1 

              case (jMin, jMax)
                 ! J range meaningless
                 cellRangeCGNS(1, 1) = iBeg + iBegOr - 1 
                 cellRangeCGNS(1, 2) = iEnd + iBegOr - 1 

                 cellRangeCGNS(3, 1) = jBeg + kBegOr - 1
                 cellRangeCGNS(3, 2) = jEnd + kBegOr - 1 

              case (kMin, kMax)
                 ! J range meaningless
                 cellRangeCGNS(1, 1) = iBeg + iBegOr - 1 
                 cellRangeCGNS(1, 2) = iEnd + iBegOr - 1 

                 cellRangeCGNS(2, 1) = jBeg + jBegOr - 1
                 cellRangeCGNS(2, 2) = jEnd + jBegOr - 1 
              end select

              ! Allocate enough space for the 1D buffer
              allocate(buffer((iEnd-iBeg+1)*(jEnd-jBeg+1)))
              kk = 0
              call storeSurfsolInBuffer(exch%sps, buffer, kk, nn, BCfaceID(mm), &
                   cellRangeCGNS, solNames(iSol), viscousSubFace, .False.)

              ! Now since the storeSurfSol just put things in a flat
              ! array and are face based, here we take the 1D face data
              ! and scatter to the nodes. 

              iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
              jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
              ni = iEnd - iBeg + 1
              nj = jEnd - jBeg + 1
              jj = 0
              do j=0,nj-2
                 do i=0,ni-2
                    jj = jj + 1
                    ! Scatter value to each node
                    ind(1) = ii + (j  )*ni + i + 1
                    ind(2) = ii + (j  )*ni + i + 2 
                    ind(3) = ii + (j+1)*ni + i + 2 
                    ind(4) = ii + (j+1)*ni + i + 1
                    do kk=1,4
                       localPtr(ind(kk)) = localPtr(ind(kk)) + buffer(jj)
                    end do
                 end do
              end do

              ii = ii + ni*nj
              deallocate(buffer)
           end if
        end do bocoLoop
     end do domainLoop

     ! Return our pointers
     call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Globalize the function value
     call vecSet(exch%nodeValGlobal, zero, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterBegin(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, exch%nodeValLocal, &
          exch%nodeValGlobal, ADD_VALUES, SCATTER_FORWARD, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Now divide by the weighting sum. We can do this with a
     ! vecpointwisemult since we already divided by the weight. 
     call vecPointwiseMult(exch%nodeValGlobal, exch%nodeValGlobal, &
          exch%sumGlobal, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Temporarly place the nodalValues into the array since that is
     ! where we want the data placed.
     call VecPlaceArray(exch%nodeValLocal, exch%nodalValues(:, iSol+9), ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Push back to the local values
     call VecScatterBegin(exch%scatter, exch%nodeValGlobal, &
          exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecScatterEnd(exch%scatter, exch%nodeValGlobal, &
          exch%nodeValLocal, INSERT_VALUES, SCATTER_REVERSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Restore the array we replaed.
     call VecResetArray(exch%nodeValLocal, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  end do varLoop
  deallocate(solNames)
end subroutine computeSurfaceOutputNodalData


subroutine createSlice(exch, slc, pt, dir, sliceName, famList, nFam)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine creates a slice on a plane defined by pt and   *
  !      * and dir. It only uses the families specified in the famList.   *
  !      * sps define which specral instance to use. 
  !      *                                                                *
  !      ******************************************************************
  !
  use liftDistributionData
  implicit none

  ! Input param
  type(familyExchange), target :: exch
  type(slice), intent(inout) :: slc
  real(kind=realType), dimension(3), intent(in) :: pt, dir
  character*(*), intent(in) :: sliceName
  integer(kind=intType), intent(in) :: famList(nFam), nFam

  ! Working param
  integer(kind=intType) :: i, nMax, nUnique, oldInd, newInd, bSearchIntegers
  integer(kind=intType) :: patchIndices(4), indexSquare, jj, kk, icon, iCoor, num1, num2
  real(kind=realType) :: f(4), d, ovrdnom, tol
  logical :: logic1, foundFam
  real(kind=realType), dimension(:, :), pointer :: tmpWeight, dummy, tmpNodes
  integer(kind=intType), dimension(:, :), pointer :: tmpInd
  integer(kind=intType), dimension(:), allocatable :: link

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

  ! Allocate the family list this slice is to use:
  allocate(slc%famList(nFam))
  slc%famList = famList
  slc%sliceName = sliceName
  ! Set the info for the slice:
  slc%pt = pt
  slc%dir = dir
  slc%nNodes = 0
  slc%exch => exch

  ! First step is to compute the 'function' value that will be used
  ! for the contour.

  ! Equation of plane: ax + by + cz + d = 0
  d = -pt(1)*dir(1) - pt(2)*dir(2) - pt(3)*dir(3)
  ovrdnom = one/sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)

  ! Compute the distance function on all possible surfaces on this
  ! processor.
  do i=1, exch%nNodes
     ! Now compute the signed distance
     exch%fc(i) = (dir(1)*exch%nodalValues(i, 1) + dir(2)*exch%nodalValues(i, 2) + &
          dir(3)*exch%nodalValues(i, 3) + d)*ovrdnom
  end do

  ! Estimate size of slice by the 5 times sqrt of the number of nodes in the
  ! mesh. Exact size doesn't matter as we realloc if necessary. 
  nMax = int(sqrt(dble(exch%nNodes)))
  allocate(tmpWeight(2,nMax), tmpInd(2, nMax), tmpNOdes(3, nMax))

  iCoor = 0
  oldInd = 1

  ! Loop over all elements
  elemLoop: do i=1, size(exch%conn, 2)

     ! Determine if this the family ID of this patch should be
     ! included. We will be a little cheeky here. Since most of the
     ! elements next to each other have the same family, 99% of the
     ! time we can skip the search and jsut look at the index of last
     ! search to see if that entry in the family list matches our
     ! current cell. This reduces the cost to O(1) since if we know
     ! our cell famID is in the list we're done. If the old index
     ! doesn't work, we can still search and find it or it may not
     ! actually be there. Of course for cells that are not in it, you
     ! always have to search. 

     foundFam = slc%famList(oldInd) == exch%elemFam(i)
     if (.not. foundFam) then 
        newInd = bsearchIntegers(exch%elemFam(i), slc%famList, size(slc%famList))
        if (newInd > 0) then
           foundFam = .True.
           oldInd = newInd
        end if
     end if

     includeElem: if (foundFam) then
        ! Extract the indices and function values at each corner
        do jj=1,4
           patchIndices(jj) = exch%conn(jj, i)
           f(jj) = exch%fc(patchIndices(jj))
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
                 call reallocateReal2(tmpWeight, 2, 2*nMax, 2, nMax, .true.)
                 call reallocateReal2(tmpNodes, 3, 2*nMax, 3, nMax, .true.)
                 call reallocateInteger2(tmpInd, 2, 2*nMax, 2, nMax, .true.)
                 nMax = nMax * 2
              end if

              ! Weight factors
              tmpWeight(2, iCoor) = (zero - f(num1))/(f(num2) - f(num1))
              tmpWeight(1, iCoor) = one - tmpWeight(2, icoor)

              ! Store the weight factors
              tmpInd(:, iCoor) = (/patchIndices(num1), patchIndices(num2)/)

              ! Store the physical nodes so we know how to reduce
              tmpNodes(:, iCoor) = &
                   tmpWeight(1, iCoor)*exch%nodalValues(tmpInd(1, iCoor), 1:3) + &
                   tmpWeight(2, iCoor)*exch%nodalValues(tmpInd(2, iCoor), 1:3)

              kk = kk + 1
           end if
        end do
     end if includeElem
  end do  ElemLoop

  ! To save space, we can compact out the doubly defined nodes that
  ! were created during the slicing process. Then we can allocate the
  ! final weight array and index array to be the exact correct
  ! length

  allocate(dummy(3, iCoor), link(iCoor))
  tol=1e-12
  call pointReduce(tmpNodes, iCoor, tol, dummy, link, nUnique)
  allocate(slc%w(2, nUnique), slc%ind(2, nUnique), slc%conn(2, iCoor/2))
  slc%nNodes = nUnique

  ! Modify the data accordingly
  do i=1, iCoor
     slc%w(:, link(i)) = tmpWeight(:, i)
     slc%ind(:, link(i)) = tmpInd(:, i)
  end do

  ! The connectivity is actually link reshaped since the original conn
  ! would have been (1,2), (3,4), (5,6) etc.
  do i=1, iCoor/2
     slc%conn(1, i) = link(2*i-1)
     slc%conn(2, i) = link(2*i)
  end do

  deallocate(tmpNodes, tmpWeight, tmpInd, dummy, link)
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
  if (allocated(slc%w)) then
     deallocate(slc%w)
  end if

  if (allocated(slc%ind)) then
     deallocate(slc%ind) 
  end if

  if (allocated(slc%conn)) then
     deallocate(slc%conn)
  end if

  if (allocated(slc%famList)) then 
     deallocate(slc%famList)
  end if

  if (allocated(slc%vars)) then 
     deallocate(slc%vars)
  end if

end subroutine destroySlice

subroutine integrateSlice(lSlc, gSlc, nFields, doConnectivity)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * This subroutine integrates the forces on slice slc and computes*
  !      * the integrated quantities for lift, drag, cl and cd with       *
  !      * contributions from both the pressure and visoucs forces.       *
  !      * It optionally interpolates solution variables as well.         *
  !      ******************************************************************
  !
  use liftDistributionData
  use inputPhysics   
  use flowVarRefState
  use communication
  implicit none

  ! Input Variables
  type(slice) :: lSlc, gSlc
  integer(kind=intType), intent(in) :: nFields
  logical, intent(in) :: doConnectivity

  ! Working variables
  integer(kind=intType) :: i, j, i1, i2
  real(kind=realType), dimension(3) :: x1, x2, pT1, pT2, vT1, vT2, pF, vF
  real(kind=realType) :: len, dmax, dmin, dist, fact, M(3,3), tmp(4)
  real(kind=realType) :: r(3), r_new(3), hyp, te(3), le(3), theta, w1, w2
  integer(kind=intType) :: bestPair(2), dir_ind, iProc, ierr, iSize
  real(kind=realtype), dimension(:,:), allocatable :: tempCoords
  real(kind=realtype), dimension(:,:), allocatable :: localVals
  integer(kind=intType), dimension(:), allocatable :: sliceNodeSizes, sliceCellSizes
  integer(kind=intType), dimension(:), allocatable :: nodeDisps, cellDisps

  ! Copy the info related to slice
  gSlc%sliceName = trim(lSlc%sliceName)
  gSlc%pt = lSlc%pt
  gSlc%dir = lSlc%dir

  pF = zero
  vF = zero
  iSize = 3 + 6 + nFields

  allocate(localVals(iSize,lSlc%nNodes))

  ! Compute all the local variables we need.

  ! Interpolate the required values
  do i=1, lSlc%nNodes
     i1 = lSlc%ind(1, i)
     i2 = lSlc%ind(2, i)
     w1 = lSlc%w(1, i)
     w2 = lSlc%w(2, i)
     localVals(1:iSize, i) = w1*lslc%exch%nodalValues(i1, 1:iSize) + w2*lslc%exch%nodalValues(i2, 1:iSize)
  end do

  do i=1, size(lSlc%conn, 2)

     ! extract nodes:
     i1 = lslc%conn(1, i)
     i2 = lslc%conn(2, i)

     x1 = localVals(1:3, i1)
     x2 = localVals(1:3, i2)

     ! extract pressure tractions
     pT1 = localVals(4:6, i1)
     pT2 = localVals(4:6, i2)

     ! extract viscous tractions
     vT1 = localVals(7:9, i1)
     vT2 = localVals(7:9, i2)

     ! Length of this segment
     len = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

     ! Integrate the pressure and viscous forces separately
     pF = pF + half*(pT1 + pT2)*len
     vF = vF + half*(vT1 + vT2)*len

  end do

  ! That is as far as we can go in parallel. We now have to gather up
  ! pL, pD, vL vD as well as the nodes to the root proc. 

  ! Gather up the number of nodes to be set to the root proc:
  allocate(sliceNodeSizes(nProc), nodeDisps(0:nProc))
  sliceNodeSizes = 0
  nodeDisps = 0
  call mpi_allgather(lSlc%nNodes,1, sumb_integer, sliceNodeSizes, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  nodeDisps(0) = 0
  do iProc=1, nProc
     nodeDisps(iProc) = nodeDisps(iProc-1) + sliceNodeSizes(iProc)*iSize
  end do

  if (myid == 0) then 
     gSlc%nNodes = sum(sliceNodeSizes)
     allocate(gSlc%vars(iSize, gSlc%nNodes))

  end if

  call mpi_gatherv(localVals, iSize*lSlc%nNodes, sumb_real, gSlc%vars, sliceNodeSizes*iSize, &
       nodeDisps, sumb_real, 0, sumb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! We may also need to gather the connectivity if the slice will have
  ! to be written to a file. 
  if (doConnectivity) then 
     i = size(lslc%conn, 2)
     allocate(cellDisps(0:nProc), sliceCellSizes(nProc))

     call mpi_gather(i, 1, sumb_integer, sliceCellSizes, 1, sumb_integer, &
          0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)


     if (myid == 0) then 
        cellDisps(0) = 0
        do iProc=1, nProc
           cellDisps(iProc) = cellDisps(iProc-1) + sliceCellSizes(iProc)*2
        end do
        allocate(gSlc%conn(2, sum(sliceCellSizes)))
     end if

     ! We offset the conn array by nodeDisps(iProc) which
     ! automagically adjust the connectivity to account for the
     ! number of nodes from different processors

     call mpi_gatherv(lSlc%conn+nodeDisps(myid)/iSize, 2*size(lSlc%conn, 2), sumb_integer, gSlc%conn, &
          sliceCellSizes*2, cellDisps, sumb_integer, 0, sumb_comm_world, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Not quite finished yet since we will have gathered nodes from
     ! multiple procs we have to adjust the connectivity

     deallocate(sliceCellSizes, cellDisps)
  end if

  deallocate(sliceNodeSizes, nodeDisps)
  ! Set the local values we can in the slice
  lSlc%pL = liftDirection(1)*pF(1) + liftDirection(2)*pF(2) + liftDirection(3)*pF(3)
  lSlc%pD = dragDirection(1)*pF(1) + dragDirection(2)*pF(2) + dragDirection(3)*pF(3)
  lSlc%vL = liftDirection(1)*vF(1) + liftDirection(2)*vF(2) + liftDirection(3)*vF(3)
  lSlc%vD = dragDirection(1)*vF(1) + dragDirection(2)*vF(2) + dragDirection(3)*vF(3)

  ! Reduce the lift/drag values
  call mpi_reduce((/lSlc%pL, lSlc%pD, lSlc%vL, lSlc%vD/), tmp, 4, sumb_real, MPI_SUM, &
       0, sumb_comm_world, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (myid == 0) then 
     gSlc%pL = tmp(1)
     gSlc%pD = tmp(2)
     gSlc%vL = tmp(3)
     gSlc%vD = tmp(4)

     ! Compute the chord as the max length between any two nodes...this
     ! is n^2, so should be changed in the future

     dmax = zero
     bestPair = (/1, 1/)
     do i=1,size(gSlc%vars, 2)
        ! extract node:
        x1 = gSlc%vars(1:3, i)

        do j=i+1,size(gSlc%vars, 2)
           ! extract node:
           x2 = gSlc%vars(1:3, j)

           dist = sqrt((x1(1)-x2(1))**2 +(x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

           if (dist > dmax) then
              dmax = dist
              bestPair = (/i, j/)
           end if
        end  do
     end do

     ! Set chord, protected from zero
     gSlc%chord = max(dmax, 1e-12)

     ! Compute factor to get coefficient
     fact = two/(gammaInf*pInf*MachCoef*MachCoef*pRef)

     ! Take dmax as chord and compute coefficients
     gSlc%CLp = gSlc%pL / gSlc%chord * fact
     gSlc%CDp = gSlc%pD / gSlc%chord * fact
     gSlc%CLv = gSlc%vL / gSlc%chord * fact
     gSlc%CDv = gSlc%vD / gSlc%chord * fact

     ! Default values
     gSlc%twist = zero
     gSlc%thickness = zero

     if (gSlc%nNodes == 0) then
        return
     end if

     ! Lastly we need the twist and the twist and the thickness
     i1 = bestPair(1)
     i2 = bestPair(2)
     x1 = gSlc%vars(1:3, i1)
     x2 = gSlc%vars(1:3, i2)

     if (x1(1) > x2(1)) then
        te = x1
        le = x2
     else
        te = x2
        le = x1
     end if

     ! Save the leading and trailing edges so we can do scaled output
     ! later
     gSlc%le = le
     gSlc%te = te

     ! Finally we need to get the thickness. For this, compute temporary
     ! section nodes and rotate them by the twist values we just computed
     ! and take the max and min

     ! Back out what is the main index of the slice, x, y or z based on
     ! the direction. Not the best approach, but that's ok
     dir_ind = maxloc(abs(gSlc%dir),1)

     ! Length of hyptoneuse is the same
     hyp = sqrt((x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 + (x1(3)-x2(3))**2)

     if (dir_ind == 1) then
        ! Xslice...we don't how what to do here..could be y or z. Don't
        ! do anything.
        gSlc%twist = zero
     else if (dir_ind == 2) then
        ! Yslice
        theta = asin((le(3)-te(3))/hyp)
        gSlc%twist = theta*180.0/pi
     else
        ! Zslice
        theta = asin((le(2)-te(2))/hyp)
        gSlc%twist = theta*180.0/pi
     end if

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

     allocate(tempCoords(3, size(gSlc%vars, 2)))
     do i=1, size(gSlc%vars, 2)
        ! extract node:
        r = gSlc%vars(1:3, i) - te
        r_new = matmul(M, r)
        tempCoords(:, i) = r_new + te
     end do

     ! Now get the max and the min and divide by the chord for t/c
     if (dir_ind == 1) then
        gSlc%thickness = 0 ! Again, don't know what to do here
     else if(dir_ind == 2) then
        dmax = maxval(tempCoords(3, :))
        dmin = minval(tempCoords(3, :))
        gSlc%thickness = (dmax-dmin)/hyp
     else if(dir_ind == 3) then
        dmax = maxval(tempCoords(2, :))
        dmin = minval(tempCoords(2, :))
        gSlc%thickness = (dmax-dmin)/hyp
     end if
     deallocate(tempCoords)
  end if
end subroutine integrateSlice

subroutine writeSlice(slc, fileID, nFields)
  ! Write the data in slice 'slc' to openfile ID fileID

  use liftDistributionData
  use inputIO
  implicit none

  ! Input Parameters
  type(slice), intent(in) :: slc
  integer(kind=intType), intent(in) :: fileID, nFields

  ! Working Variables
  integer(kind=intType) :: i, j
  real(kind=realType) :: tmp, tx, ty, tz
  write (fileID,"(a,a,a)") "Zone T= """,trim(slc%sliceName),""""

  ! IF we have nodes actually write:
  if (slc%nNodes > 0) then
     write (fileID,*) "Nodes = ", slc%nNodes, " Elements= ",  size(slc%conn, 2), " ZONETYPE=FELINESEG"
     write (fileID,*) "DATAPACKING=POINT"
13   format (E14.6)

     do i=1,slc%nNodes
        ! Write the coordinates
        do j=1,3
           write(fileID,13, advance='no') slc%vars(j, i)
        end do

        ! Write the scaled coordiantes with the LE at (0,0,0)
        do j=1,3
           tmp = slc%vars(j, i)
           write(fileID,13, advance='no') (tmp - slc%le(j))/slc%chord
        end do

        ! Write field data. Starts at 9 (after 3 coordindates and the 6 tractions)
        do j=1,nFields
           write(fileID,13, advance='no') slc%vars(9+j, i)
        end do

        write(fileID,"(1x)")
     end do

15   format(I5, I5)
     do i=1, size(slc%conn, 2)
        write(fileID, 15) slc%conn(1, i), slc%conn(2, i)
     end do
  else ! Write dummy data so the number of zones are the same

     write (fileID,*) "Nodes = ", 2, " Elements= ", 1, " ZONETYPE=FELINESEG"
     write (fileID,*) "DATAPACKING=POINT"
     do i=1,2
        do j=1,6
           write(fileID,13, advance='no') zero
        end do

        do j=7,nFields
           write(fileID,13, advance='no') zero
        end do

        write(fileID,"(1x)")
     end do
     write(fileID, 15) 1, 2
  end if
end subroutine writeSlice



! subroutine computeNodalAveragingVector(exch, avgType, sps, famList, nFamList)

!   use liftDistributionData
!   use communication
!   use inputPhysics

!   implicit none

! #define PETSC_AVOID_MPIF_H
! #include "include/petscversion.h"
! #if PETSC_VERSION_MINOR > 5
! #include "petsc/finclude/petscsys.h"
! #include "petsc/finclude/petscvec.h"
! #include "petsc/finclude/petscvec.h90"
! #else
! #include "include/finclude/petscsys.h"
! #include "include/finclude/petscvec.h"
! #include "include/finclude/petscvec.h90"
! #endif
!   ! Input
!   type(familyExchange), intent(inout) :: exch
!   integer(kind=intType), intent(in) :: avgType, sps, famList(nFamlist), nFamList

!   ! Working
!   integer(kind=intType) :: i, j, ierr, bsearchintegers
!   real(kind=realType), dimension(:), pointer :: weightPtr, localPtr
!   real(kind=realType), dimension(3) ::sss, v1, v2, n1, n2, n3, n4
!   real(kind=realType) :: cellArea, weight

!   ! There are three kinds of averaging we could do:
!   ! 1. Arithematic
!   ! 2. Area weighted
!   ! 3. Inverse area weighted

!   ! Now compute the nodal area for each quad scattering them into the
!   ! localPtr as we go
!   call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   call vecGetArrayF90(exch%localWeight, weightPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   localPtr = zero

!   do i=1, size(conn, 2)

!      if (bsearchIntegers(elemFam(i), famList, nfamList)>  0) then 

!         n1 = uniqueNodes(:, conn(1, i), sps)
!         n2 = uniqueNodes(:, conn(2, i), sps)
!         n3 = uniqueNodes(:, conn(3, i), sps)
!         n4 = uniqueNodes(:, conn(4, i), sps)

!         ! Recall these are CCW ordering!
!         v1 = n3 - n1
!         v2 = n4 - n2

!         ! Cross product
!         sss(1) = half*(v1(2)*v2(3) - v1(3)*v2(2))
!         sss(2) = half*(v1(3)*v2(1) - v1(1)*v2(3))
!         sss(3) = half*(v1(1)*v2(2) - v1(2)*v2(1))

!         cellArea = sqrt(sss(1)*sss(1) + sss(2)*sss(2) + sss(3)*sss(3))
!         if (avgType == 1) then
!            weight = one
!         else if(avgType == 2) then 
!            weight = cellArea*fourth
!         else if(avgType == 3) then 
!            weight = fourth/(max(cellArea, 1e-30))
!         end if

!         ! Store the weight value for future use.
!         weightPtr(i) = weight

!         ! Scatter to the nodes
!         do j=1,4
!            localPtr(conn(j, i)) = localPtr(conn(j, i)) + weight
!         end do
!      end if
!   end do

!   call vecRestoreArrayF90(nodeValLocal, localPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   call vecRestoreArrayF90(localWeight, weightPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Globalize the weights
!   call vecSet(sumGlobal, zero, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   call VecScatterBegin(tracScatter, nodeValLocal, sumGlobal, ADD_VALUES, &
!        SCATTER_FORWARD, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   call VecScatterEnd(tracScatter, nodeValLocal, sumGlobal, ADD_VALUES, &
!        SCATTER_FORWARD, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   ! Now compute the inverse of the weighting so that we can multiply
!   ! instead of dividing.

!   call vecGetArrayF90(sumGlobal, localPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

!   localPtr = one/localPtr

!   call vecRestoreArrayF90(sumGlobal, localPtr, ierr)
!   call EChk(ierr,__FILE__,__LINE__)

! end subroutine computeNodalAveragingVector
